# fits_metadata_extractor/search.py

import logging
import pandas as pd
import json
from shapely.geometry import Point, Polygon
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np

def is_point_in_convex_polygon(ra, dec, polygon_coords, polygon_frame='icrs'):
    """
    Determines if a point (ra, dec) in 'icrs' frame is inside a convex polygon defined by polygon_coords,
    which are in 'polygon_frame'.
    """
    try:
        if not polygon_coords or len(polygon_coords) < 3:
            logging.error("Polygon must have at least 3 vertices.")
            return False

        # Convert polygon_coords to 'icrs' frame if needed
        if polygon_frame.lower() == 'galactic':
            # Convert polygon_coords from Galactic to ICRS
            l_vals, b_vals = zip(*polygon_coords)
            sky_coords = SkyCoord(l=l_vals * u.deg, b=b_vals * u.deg, frame='galactic')
            sky_coords_icrs = sky_coords.icrs
            ra_vals = sky_coords_icrs.ra.deg
            dec_vals = sky_coords_icrs.dec.deg
            polygon_coords_icrs = list(zip(ra_vals, dec_vals))
        elif polygon_frame.lower() == 'icrs':
            polygon_coords_icrs = polygon_coords
        else:
            logging.error(f"Unsupported coordinate frame: {polygon_frame}")
            return False

        # Create Shapely Polygon in ICRS frame
        polygon = Polygon(polygon_coords_icrs)

        # Create point in ICRS frame
        point = Point(ra, dec)

        return polygon.contains(point) or polygon.touches(point)

    except Exception as e:
        logging.error(f"Error in is_point_in_convex_polygon: {e}")
        return False


def search_fits_by_point(metadata_df, ra, dec):
    """
    Searches for FITS files in the collection whose coverage includes the given point.
    Assumes the point is given in Equatorial coordinates (ICRS frame).
    """
    matching_indices = []

    for idx, row in metadata_df.iterrows():
        polygon_coords_str = row['Polygon_Coords']
        polygon_frame = row.get('Coordinate_Frame', 'icrs')  # Default to 'icrs' if not specified

        if pd.isnull(polygon_coords_str):
            continue
        try:
            # Deserialize Polygon_Coords
            polygon_coords = json.loads(polygon_coords_str)

            # Check if the point is inside the polygon
            if is_point_in_convex_polygon(ra, dec, polygon_coords, polygon_frame=polygon_frame):
                matching_indices.append(idx)
        except Exception as e:
            logging.error(f"Error processing Polygon_Coords for index {idx}: {e}")

    # Return the subset of metadata_df with matching FITS files
    return metadata_df.loc[matching_indices]


def search_fits_by_region(metadata_df, region):
    """
    Searches for FITS files in the collection whose coverage intersects the given region.
    Assumes the region is given in Equatorial coordinates (ICRS frame).
    """
    matching_indices = []

    # Prepare the region Polygon in ICRS frame
    if region['type'] == 'circle':
        center_ra, center_dec = region['center']
        radius = region['radius']  # In degrees
        # Create a circular region using SkyCoord and astropy's angle units
        center = SkyCoord(ra=center_ra * u.deg, dec=center_dec * u.deg, frame='icrs')
        circle_region = center.directional_offset_by(0 * u.deg, radius * u.deg)
        # Create a polygon approximation of the circle
        num_points = 100
        angles = np.linspace(0, 2 * np.pi, num_points)
        ra_vals = center_ra + radius * np.cos(angles)
        dec_vals = center_dec + radius * np.sin(angles)
        region_polygon = Polygon(zip(ra_vals, dec_vals))
    elif region['type'] == 'polygon':
        coordinates = region['coordinates']
        # Assume coordinates are in ICRS frame
        region_polygon = Polygon(coordinates)
    else:
        raise ValueError("Unsupported region type. Supported types are 'circle' and 'polygon'.")

    for idx, row in metadata_df.iterrows():
        polygon_coords_str = row['Polygon_Coords']
        polygon_frame = row.get('Coordinate_Frame', 'icrs')  # Default to 'icrs' if not specified

        if pd.isnull(polygon_coords_str):
            continue
        try:
            # Deserialize Polygon_Coords
            polygon_coords = json.loads(polygon_coords_str)

            # Convert polygon_coords to ICRS frame if needed
            if polygon_frame.lower() == 'galactic':
                # Convert polygon_coords from Galactic to ICRS
                l_vals, b_vals = zip(*polygon_coords)
                sky_coords = SkyCoord(l=l_vals * u.deg, b=b_vals * u.deg, frame='galactic')
                sky_coords_icrs = sky_coords.icrs
                ra_vals = sky_coords_icrs.ra.deg
                dec_vals = sky_coords_icrs.dec.deg
                polygon_coords_icrs = list(zip(ra_vals, dec_vals))
            elif polygon_frame.lower() == 'icrs':
                polygon_coords_icrs = polygon_coords
            else:
                logging.error(f"Unsupported coordinate frame: {polygon_frame}")
                continue

            fits_polygon = Polygon(polygon_coords_icrs)

            # Check for intersection
            if fits_polygon.intersects(region_polygon):
                matching_indices.append(idx)
        except Exception as e:
            logging.error(f"Error processing Polygon_Coords for index {idx}: {e}")

    # Return the subset of metadata_df with matching FITS files
    return metadata_df.loc[matching_indices]

