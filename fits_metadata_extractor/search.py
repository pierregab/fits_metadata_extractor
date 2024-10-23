# fits_metadata_extractor/search.py

import logging
import pandas as pd
import json
from shapely.geometry import Point, Polygon

def is_point_in_convex_polygon(ra, dec, polygon_coords):
    """
    Determines if a point (ra, dec) is inside a convex polygon defined by polygon_coords.
    """
    try:
        if not polygon_coords or len(polygon_coords) < 3:
            logging.error("Polygon must have at least 3 vertices.")
            return False

        # Create Shapely Polygon
        polygon = Polygon(polygon_coords)
        point = Point(ra, dec)
        return polygon.contains(point) or polygon.touches(point)

    except Exception as e:
        logging.error(f"Error in is_point_in_convex_polygon: {e}")
        return False

def search_fits_by_point(metadata_df, ra, dec):
    """
    Searches for FITS files in the collection whose coverage includes the given point.
    """
    matching_indices = []

    for idx, row in metadata_df.iterrows():
        polygon_coords_str = row['Polygon_Coords']
        if pd.isnull(polygon_coords_str):
            continue
        try:
            # Deserialize Polygon_Coords
            polygon_coords = json.loads(polygon_coords_str)

            # Check if the point is inside the polygon
            if is_point_in_convex_polygon(ra, dec, polygon_coords):
                matching_indices.append(idx)
        except Exception as e:
            logging.error(f"Error processing Polygon_Coords for index {idx}: {e}")

    # Return the subset of metadata_df with matching FITS files
    return metadata_df.loc[matching_indices]

def search_fits_by_region(metadata_df, region):
    """
    Searches for FITS files in the collection whose coverage intersects the given region.
    """
    matching_indices = []

    if region['type'] == 'circle':
        center_ra, center_dec = region['center']
        radius = region['radius']
        # Create a Shapely circle (buffered point)
        center_point = Point(center_ra, center_dec)
        circle = center_point.buffer(radius)
    elif region['type'] == 'polygon':
        coordinates = region['coordinates']
        region_polygon = Polygon(coordinates)
    else:
        raise ValueError("Unsupported region type. Supported types are 'circle' and 'polygon'.")

    for idx, row in metadata_df.iterrows():
        polygon_coords_str = row['Polygon_Coords']
        if pd.isnull(polygon_coords_str):
            continue
        try:
            # Deserialize Polygon_Coords
            polygon_coords = json.loads(polygon_coords_str)
            fits_polygon = Polygon(polygon_coords)

            # Check for intersection
            if region['type'] == 'circle':
                if fits_polygon.intersects(circle):
                    matching_indices.append(idx)
            elif region['type'] == 'polygon':
                if fits_polygon.intersects(region_polygon):
                    matching_indices.append(idx)
        except Exception as e:
            logging.error(f"Error processing Polygon_Coords for index {idx}: {e}")

    # Return the subset of metadata_df with matching FITS files
    return metadata_df.loc[matching_indices]
