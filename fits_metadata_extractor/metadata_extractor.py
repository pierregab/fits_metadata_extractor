# fits_metadata_extractor/metadata_extractor.py

import os
import logging
import re
import json
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, Angle
import astropy.units as u
from astropy.utils.exceptions import AstropyWarning
import warnings
from mocpy import MOC
from regions import PolygonSkyRegion
from shapely.geometry import Polygon
from .resolver import ObjectNameResolver
from .polygon_utils import polygonskyregion_to_shapely
from .logger import log_success

# Remove Astropy warnings
warnings.simplefilter('ignore', AstropyWarning)

# Constants
REQUIRED_KEYWORDS = [
    'NAXIS', 'NAXIS1', 'NAXIS2', 'OBJECT', 'RA', 'DEC',
    'RADESYS', 'DATE-OBS', 'MJD-OBS', 'EXPTIME',
    'INSTRUME', 'TELESCOP', 'JD'
]

class FITSMetadataExtractor:
    """
    Class for extracting and homogenizing metadata from FITS files.
    """

    def __init__(self):
        self.resolver = ObjectNameResolver()

    def infer_object_name(self, fits_file):
        """
        Infers the object name from the FITS file name based on naming conventions.

        Parameters:
            fits_file (str): Path to the FITS file.

        Returns:
            str: Inferred object name or 'Unknown'.
        """
        filename = os.path.basename(fits_file)

        # Enhanced regex to capture more precise object names
        match = re.match(r'^([A-Za-z0-9\.\-]+)', filename)
        if match:
            inferred_name = match.group(1)
            # Further sanitize if needed
            inferred_name = re.sub(r'[^A-Za-z0-9\.\-]', '', inferred_name)
            return inferred_name

        return 'Unknown'

    def extract_fits_metadata(self, fits_file):
        """
        Extracts relevant metadata from a FITS header.

        Parameters:
            fits_file (str): Path to the FITS file.

        Returns:
            dict: Dictionary containing extracted and homogenized metadata.
        """
        # Initialize metadata with all expected keys set to None
        metadata = {
            'FITS_File': os.path.basename(fits_file),
            'Polygon': None,
            'MOC': None,
            'Polygon_Coords': None,
            'RA_WCS': None,      # New column for RA via WCS
            'DEC_WCS': None      # New column for DEC via WCS
        }
        try:
            with fits.open(fits_file) as hdul:
                header = hdul[0].header

                # Extract required keywords
                for key in REQUIRED_KEYWORDS:
                    metadata[key] = header.get(key, None)

                # Handle NAXIS and NAXISi
                naxis = metadata.get('NAXIS', None)
                if naxis is not None:
                    for i in range(1, naxis + 1):
                        key = f'NAXIS{i}'
                        metadata[key] = header.get(key, None)
                else:
                    logging.warning(f"NAXIS keyword missing in {fits_file}.")

                # Handle RADESYS keyword (use RADESYSa if available)
                radesys = header.get('RADESYSa', header.get('RADECSYS', 'ICRS'))
                metadata['RADESYS'] = radesys

                # Check 'OBJECT' keyword
                object_name = metadata.get('OBJECT', 'Unknown')
                if not object_name or object_name.strip().lower() == 'unknown':
                    # Attempt to infer object name from file name
                    inferred_name = self.infer_object_name(fits_file)
                    logging.info(f"Inferred object name '{inferred_name}' for file '{fits_file}'.")
                    object_name = inferred_name

                # Resolve object name using the enhanced resolver with fallbacks
                resolved_name, resolution_method = self.resolver.resolve_object_name(object_name)
                metadata['Resolved_Object'] = resolved_name

                # Compute Polygon and MOC
                try:
                    # Initialize WCS object from header
                    wcs = WCS(header)
                    naxis = header.get('NAXIS', None)
                    naxis1 = header.get('NAXIS1', None)
                    naxis2 = header.get('NAXIS2', None)
                    if naxis is not None and naxis1 is not None and naxis2 is not None:
                        # Determine number of extra axes
                        extra_axes = [1] * (naxis - 2) if naxis > 2 else []

                        # Define the four corners of the image in pixel coordinates
                        corners_pixel = np.array([
                            [0.5, 0.5] + extra_axes,
                            [naxis1 + 0.5, 0.5] + extra_axes,
                            [naxis1 + 0.5, naxis2 + 0.5] + extra_axes,
                            [0.5, naxis2 + 0.5] + extra_axes
                        ])

                        # Extract the coordinate types
                        ctype1 = wcs.wcs.ctype[0]
                        ctype2 = wcs.wcs.ctype[1]

                        # Determine the coordinate frame
                        if 'RA' in ctype1.upper() and 'DEC' in ctype2.upper():
                            coord_frame = 'icrs'
                        elif 'GLON' in ctype1.upper() and 'GLAT' in ctype2.upper():
                            coord_frame = 'galactic'
                        else:
                            logging.error(f"Unsupported coordinate system in {fits_file}: {ctype1}, {ctype2}")
                            coord_frame = 'unknown'

                        # Convert pixel coordinates to world coordinates
                        try:
                            corners_world = wcs.all_pix2world(corners_pixel, 1)  # origin=1 for FITS convention
                            x_coord = corners_world[:, 0]
                            y_coord = corners_world[:, 1]
                        except Exception as e:
                            logging.error(f"WCS transformation failed for {fits_file}: {e}")
                            corners_world = None
                            x_coord = None
                            y_coord = None

                        if x_coord is not None and y_coord is not None:
                            # Create SkyCoord objects with the appropriate frame
                            if coord_frame == 'icrs':
                                sky_coords = SkyCoord(ra=x_coord * u.deg, dec=y_coord * u.deg, frame='icrs')
                            elif coord_frame == 'galactic':
                                sky_coords = SkyCoord(l=x_coord * u.deg, b=y_coord * u.deg, frame='galactic')
                            else:
                                logging.error(f"Cannot create SkyCoord for unknown coordinate frame in {fits_file}.")
                                sky_coords = None

                            if sky_coords is not None:
                                # Create a PolygonSkyRegion
                                polygon_region = PolygonSkyRegion(vertices=sky_coords)

                                # Convert to Shapely polygon manually
                                shapely_polygon = polygonskyregion_to_shapely(polygon_region)

                                if shapely_polygon:
                                    # Store polygon WKT and coordinates
                                    metadata['Polygon'] = shapely_polygon.wkt
                                    metadata['Polygon_Coords'] = list(shapely_polygon.exterior.coords)

                                    # Compute MOC using mocpy
                                    try:
                                        moc = MOC.from_polygon_skycoord(sky_coords, max_depth=10)
                                        moc_str = moc.serialize(format='str')
                                        metadata['MOC'] = moc_str

                                        # Gather MOC details
                                        moc_details = (
                                            f"Sky Fraction: {moc.sky_fraction:.6f}, "
                                            f"Number of Cells: {len(moc.uniq_hpx)}, "
                                            f"Max Order: {moc.max_order}"
                                        )

                                        # Log success message in green with resolution method and MOC details
                                        success_message = (
                                            f"Successfully processed '{fits_file}'. "
                                            f"Resolved using {resolution_method}: '{resolved_name}'. "
                                            f"MOC Details - {moc_details}. "
                                            f"Polygon: {shapely_polygon.wkt}."
                                        )
                                        log_success(success_message)
                                    except Exception as e:
                                        logging.error(f"Failed to compute MOC for {fits_file}: {e}")
                                else:
                                    logging.error(f"Shapely polygon is None for {fits_file}.")
                            else:
                                logging.error(f"Skipping polygon and MOC computation for {fits_file} due to SkyCoord errors.")
                        else:
                            logging.error(f"Skipping polygon and MOC computation for {fits_file} due to WCS errors.")

                        # **New Section: Extract RA and DEC via WCS**
                        try:
                            # Calculate the center pixel
                            center_pixel = [(naxis1 + 1) / 2, (naxis2 + 1) / 2] + extra_axes
                            # Convert to world coordinates
                            center_world = wcs.all_pix2world([center_pixel], 1)[0]
                            center_x = center_world[0] % 360  # Ensure longitude is within [0, 360)
                            center_y = center_world[1]
                            if coord_frame == 'icrs':
                                metadata['RA_WCS'] = center_x
                                metadata['DEC_WCS'] = center_y
                                logging.info(f"Extracted RA_WCS: {center_x}, DEC_WCS: {center_y} for {fits_file}.")
                            elif coord_frame == 'galactic':
                                metadata['GLON_WCS'] = center_x
                                metadata['GLAT_WCS'] = center_y
                                logging.info(f"Extracted GLON_WCS: {center_x}, GLAT_WCS: {center_y} for {fits_file}.")
                        except Exception as e:
                            logging.error(f"Failed to extract center coordinates for {fits_file}: {e}")
                    else:
                        logging.warning(f"NAXIS1 or NAXIS2 is missing in the header of {fits_file}.")
                except Exception as e:
                    logging.error(f"Failed to compute polygon and MOC for {fits_file}: {e}")

        except Exception as e:
            logging.error(f"Failed to extract metadata from {fits_file}: {e}", exc_info=True)

        return metadata
