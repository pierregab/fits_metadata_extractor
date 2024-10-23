#!/usr/bin/env python3
# fits_metadata_extractor.py

"""
FITS Metadata Extractor Module

This module provides functionalities to extract and homogenize metadata from FITS headers.
It processes a collection of FITS files, resolves object names, computes coverage regions,
and compiles the metadata into a structured dataset.

Author: Bibal Sobeaux Pierre Gabriel
Date: October 20, 2024
"""

import os
import glob
import logging
import argparse
import re
from urllib.parse import quote
from functools import lru_cache
import json

import numpy as np
import pandas as pd
import requests
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.simbad import Simbad
from astroquery.ipac.ned import Ned  # Updated import
from astroquery.vizier import Vizier
from astroquery.exceptions import InvalidQueryError, TimeoutError
from filelock import FileLock, Timeout as LockTimeout
import xml.etree.ElementTree as ET
from astropy.wcs import WCS
from mocpy import MOC
from astropy.utils.exceptions import AstropyWarning
import warnings
from astropy.coordinates import Angle
import matplotlib.pyplot as plt
from regions import PolygonSkyRegion
from shapely.geometry import Polygon  # Added import
from reproject import reproject_interp  
from astropy.wcs import FITSFixedWarning

# Remove FITSFixedWarning from warnings
warnings.simplefilter('ignore', FITSFixedWarning)

# Import custom mappings
from fits_metadata_extractor.custom_mapping import CUSTOM_NAME_MAPPING  # Ensure this module exists

# Suppress non-critical warnings
warnings.simplefilter('ignore', category=AstropyWarning)

# Constants
REQUIRED_KEYWORDS = [
    'NAXIS', 'NAXIS1', 'NAXIS2', 'OBJECT', 'RA', 'DEC',
    'RADESYS', 'DATE-OBS', 'MJD-OBS', 'EXPTIME',
    'INSTRUME', 'TELESCOP', 'JD'
]

# Configure custom Simbad and Vizier
custom_simbad = Simbad()
custom_simbad.add_votable_fields('otype', 'ra(d)', 'dec(d)')

custom_vizier = Vizier(columns=['*'], row_limit=1)

def setup_logging():
    """
    Sets up logging to output to both file and console with appropriate formatting.
    """
    # Define ANSI color codes for colored logging
    GREEN = '\033[92m'
    RESET = '\033[0m'

    class ConsoleFormatter(logging.Formatter):
        """Custom logging formatter to add colors to log messages."""
        def format(self, record):
            message = super().format(record)
            if record.levelno == logging.INFO and getattr(record, 'is_success', False):
                message = f"{GREEN}{message}{RESET}"
            return message

    class FileFormatter(logging.Formatter):
        """Standard logging formatter without colors for file logs."""
        def format(self, record):
            return super().format(record)

    # Configure logging to output to both file and console
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # Create handlers
    file_handler = logging.FileHandler("fits_metadata_extractor.log")
    console_handler = logging.StreamHandler()

    # Create formatter instances
    console_formatter = ConsoleFormatter('%(asctime)s - %(levelname)s - %(message)s')
    file_formatter = FileFormatter('%(asctime)s - %(levelname)s - %(message)s')

    # Assign formatters to handlers
    console_handler.setFormatter(console_formatter)
    file_handler.setFormatter(file_formatter)

    # Add handlers to the logger
    logger.addHandler(console_handler)
    logger.addHandler(file_handler)

    return logger

# Initialize logger
logger = setup_logging()

def log_success(message):
    """
    Logs a success message in green.

    Parameters:
        message (str): The message to log.
    """
    extra = {'is_success': True}
    logging.info(message, extra=extra)

def standardize_object_name(object_name):
    """
    Standardizes the astronomical object name to conform to expected formats.

    Parameters:
        object_name (str): The original object name.

    Returns:
        list: A list of potential standardized object names.
    """
    if not object_name or object_name.strip().lower() == 'unknown':
        return ['Unknown']
    
    object_name = object_name.strip()

    # Create a list of potential standardized names
    potential_names = [object_name]

    # Attempt prefix correction: 'SNR' to 'SN' if applicable
    if object_name.upper().startswith('SNR'):
        sn_name = 'SN' + object_name[3:]
        potential_names.append(sn_name)

    # Remove spaces and special characters for alternative attempts
    sanitized_name = re.sub(r'\s+', '', object_name)
    sanitized_name = re.sub(r'[^\w\-]', '', sanitized_name)
    potential_names.append(sanitized_name)

    # Add more standardized forms if needed
    # Example: Append common prefixes/suffixes
    if not object_name.startswith(('SN', 'SNR')):
        potential_names.append('SN' + object_name)
        potential_names.append(object_name + 'SN')

    return potential_names

def resolve_with_sesame(object_name, retries=3):
    """
    Resolves object name using CDS Sesame service with XML output.

    Parameters:
        object_name (str): The object name to resolve.
        retries (int): Number of retry attempts in case of failure.

    Returns:
        str: Resolved object name or 'Unknown'.
    """
    if not object_name or object_name.strip().lower() == 'unknown':
        logging.warning("Object name is unknown or empty.")
        return 'Unknown'

    url = f"http://cdsweb.u-strasbg.fr/cgi-bin/nph-sesame/-oxp/{quote(object_name)}"

    attempt = 0
    while attempt < retries:
        try:
            response = requests.get(url, timeout=10)
            response.raise_for_status()

            # Parse XML response
            root = ET.fromstring(response.content)
            obj = root.find('.//obj')
            if obj is not None and 'name' in obj.attrib:
                resolved_name = obj.attrib['name']
                logging.info(f"Resolved '{object_name}' to '{resolved_name}' using Sesame XML.")
                return resolved_name

            logging.warning(f"No resolved name found in Sesame XML for '{object_name}'.")
            return 'Unknown'

        except ET.ParseError as e:
            logging.error(f"XML parsing failed for '{object_name}': {e}")
            attempt += 1
            logging.info(f"Retrying ({attempt}/{retries})...")
        except requests.RequestException as e:
            logging.error(f"Sesame request failed for '{object_name}': {e}")
            attempt += 1
            logging.info(f"Retrying ({attempt}/{retries})...")
        except Exception as e:
            logging.error(f"Unexpected error during Sesame resolution for '{object_name}': {e}")
            return 'Unknown'

    logging.error(f"Failed to resolve '{object_name}' after {retries} attempts.")
    return 'Unknown'

@lru_cache(maxsize=1000)
def resolve_object_name(object_name):
    """
    Resolves the celestial object name using multiple astronomical databases and custom mappings.

    Parameters:
        object_name (str): The original object name.

    Returns:
        tuple: (resolved_name (str), method_used (str))
    """
    if not object_name or object_name.strip().lower() == 'unknown':
        logging.warning("Object name is unknown or empty.")
        return ('Unknown', 'Unknown')

    # Check if object_name exists in custom mappings
    if object_name in CUSTOM_NAME_MAPPING:
        resolved_name = CUSTOM_NAME_MAPPING[object_name]
        logging.info(f"Resolved '{object_name}' to '{resolved_name}' using custom mapping.")
        return (resolved_name, 'Custom Mapping')

    potential_names = standardize_object_name(object_name)

    for name in potential_names:
        if name == 'Unknown':
            continue

        # Attempt resolution with Simbad
        try:
            result_simbad = custom_simbad.query_object(name)
            if result_simbad:
                resolved_name = result_simbad['MAIN_ID'][0]
                resolved_name = resolved_name.decode('utf-8') if isinstance(resolved_name, bytes) else resolved_name
                logging.info(f"Resolved '{object_name}' to '{resolved_name}' using Simbad with name '{name}'.")
                return (resolved_name, 'Simbad')
        except (InvalidQueryError, TimeoutError) as e:
            logging.error(f"Simbad query failed for '{name}': {e}")
        except Exception as e:
            logging.error(f"Unexpected error during Simbad query for '{name}': {e}")

        # Attempt resolution with NED
        try:
            result_ned = Ned.query_object(name)
            if result_ned:
                resolved_name = result_ned['Object Name'][0]
                logging.info(f"Resolved '{object_name}' to '{resolved_name}' using NED with name '{name}'.")
                return (resolved_name, 'NED')
        except (InvalidQueryError, TimeoutError) as e:
            logging.error(f"NED query failed for '{name}': {e}")
        except Exception as e:
            logging.error(f"Unexpected error during NED query for '{name}': {e}")

        # Attempt resolution with VizieR
        try:
            result_vizier = custom_vizier.query_object(name)
            if result_vizier:
                if 'Identifier' in result_vizier[0].columns:
                    resolved_name = result_vizier[0]['Identifier'][0]
                    resolved_name = resolved_name.decode('utf-8') if isinstance(resolved_name, bytes) else resolved_name
                else:
                    resolved_name = name  # Fallback to standardized name
                logging.info(f"Resolved '{object_name}' to '{resolved_name}' using VizieR with name '{name}'.")
                return (resolved_name, 'VizieR')
        except (InvalidQueryError, TimeoutError) as e:
            logging.error(f"VizieR query failed for '{name}': {e}")
        except Exception as e:
            logging.error(f"Unexpected error during VizieR query for '{name}': {e}")

    # If all resolvers fail, attempt with Sesame as a last resort
    try:
        resolved_name = resolve_with_sesame(object_name)
        if resolved_name != 'Unknown':
            return (resolved_name, 'Sesame')
    except Exception as e:
        logging.error(f"Sesame fallback failed for '{object_name}': {e}")

    # If all resolvers fail
    logging.warning(f"No resolved name found for '{object_name}'.")
    return ('Unknown', 'Unknown')

def infer_object_name(fits_file):
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

@lru_cache(maxsize=1000)
def resolve_object_name_with_fallback(object_name):
    """
    Resolves the celestial object name using multiple astronomical databases, custom mappings, and fallbacks.

    Parameters:
        object_name (str): The original object name.

    Returns:
        tuple: (resolved_name (str), method_used (str))
    """
    return resolve_object_name(object_name)

def polygonskyregion_to_shapely(polygon_region):
    """
    Manually converts a PolygonSkyRegion to a Shapely Polygon.

    Parameters:
        polygon_region (PolygonSkyRegion): The PolygonSkyRegion object.

    Returns:
        shapely.geometry.Polygon: The equivalent Shapely Polygon or None if conversion fails.
    """
    try:
        # Extract longitude and latitude from SkyCoord vertices
        lon = polygon_region.vertices.data.lon.deg
        lat = polygon_region.vertices.data.lat.deg
        coords = list(zip(lon, lat))
        shapely_polygon = Polygon(coords)
        return shapely_polygon
    except Exception as e:
        logging.error(f"Failed to convert PolygonSkyRegion to Shapely Polygon: {e}")
        return None


def extract_fits_metadata(fits_file):
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
                inferred_name = infer_object_name(fits_file)
                logging.info(f"Inferred object name '{inferred_name}' for file '{fits_file}'.")
                object_name = inferred_name

            # Resolve object name using the enhanced resolver with fallbacks
            resolved_name, resolution_method = resolve_object_name_with_fallback(object_name)
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



def serialize_polygon(x):
    """
    Serializes the Polygon_Coords to a JSON string if it's not None.
    Converts NumPy arrays to lists before serialization.
    
    Parameters:
        x (list, tuple, np.ndarray, or None): The polygon coordinates.
    
    Returns:
        str or None: JSON string of coordinates or None.
    """
    if x is not None:
        if isinstance(x, np.ndarray):
            x = x.tolist()  # Convert NumPy array to list
        return json.dumps(x)
    else:
        return None

def homogenize_metadata(df):
    """
    Homogenizes units and cleans the metadata DataFrame.
    
    Parameters:
        df (pandas.DataFrame): DataFrame with extracted metadata.
    
    Returns:
        pandas.DataFrame: Homogenized DataFrame.
    """
    # Convert RA and DEC to float degrees if possible
    for coord in ['RA', 'DEC', 'RA_WCS', 'DEC_WCS']:
        if coord in df.columns:
            if df[coord].dtype == object:
                def convert_coord(coord_str):
                    try:
                        angle = Angle(coord_str)
                        return angle.degree
                    except Exception:
                        return pd.NA
                df[coord] = df[coord].apply(convert_coord)
            else:
                df[coord] = pd.to_numeric(df[coord], errors='coerce')
        else:
            logging.warning(f"Column '{coord}' not found in DataFrame.")

    # Convert DATE-OBS to ISO 8601 format
    if 'DATE-OBS' in df.columns:
        df['DATE-OBS'] = pd.to_datetime(df['DATE-OBS'], errors='coerce').dt.strftime('%Y-%m-%dT%H:%M:%S')
    else:
        logging.warning("Column 'DATE-OBS' not found in DataFrame.")

    # Convert MJD-OBS to float
    if 'MJD-OBS' in df.columns:
        df['MJD-OBS'] = pd.to_numeric(df['MJD-OBS'], errors='coerce')
    else:
        logging.warning("Column 'MJD-OBS' not found in DataFrame.")

    # Convert JD to float
    if 'JD' in df.columns:
        df['JD'] = pd.to_numeric(df['JD'], errors='coerce')
    else:
        logging.warning("Column 'JD' not found in DataFrame.")

    # Convert EXPTIME to float
    if 'EXPTIME' in df.columns:
        df['EXPTIME'] = pd.to_numeric(df['EXPTIME'], errors='coerce')
    else:
        logging.warning("Column 'EXPTIME' not found in DataFrame.")

    # Fill missing RADESYS with 'ICRS'
    if 'RADESYS' in df.columns:
        df['RADESYS'] = df['RADESYS'].fillna('ICRS')
    else:
        logging.warning("Column 'RADESYS' not found in DataFrame.")

    # Fill missing Resolved_Object names with 'Unknown'
    if 'Resolved_Object' in df.columns:
        df['Resolved_Object'] = df['Resolved_Object'].fillna('Unknown')
    else:
        logging.warning("Column 'Resolved_Object' not found in DataFrame.")

    # Ensure Polygon, MOC, and Polygon_Coords are present
    for col in ['Polygon', 'MOC', 'Polygon_Coords']:
        if col not in df.columns:
            df[col] = pd.NA
            logging.warning(f"Column '{col}' was missing and has been added with NA values.")

    # Serialize 'Polygon_Coords' as JSON strings using the helper function
    if 'Polygon_Coords' in df.columns:
        df['Polygon_Coords'] = df['Polygon_Coords'].apply(serialize_polygon)
    else:
        logging.warning("Column 'Polygon_Coords' not found in DataFrame.")

    # Serialize 'Polygon' and 'MOC' if they are not strings
    for col in ['Polygon', 'MOC']:
        if col in df.columns:
            df[col] = df[col].apply(lambda x: x if isinstance(x, str) else (json.dumps(x) if pd.notnull(x) else None))
        else:
            logging.warning(f"Column '{col}' not found in DataFrame.")

    # **New Section: Handle RA_WCS and DEC_WCS Serialization**
    for coord in ['RA_WCS', 'DEC_WCS']:
        if coord in df.columns:
            df[coord] = pd.to_numeric(df[coord], errors='coerce')
        else:
            logging.warning(f"Column '{coord}' not found in DataFrame.")

    # Log the number of missing values
    missing_values = df.isnull().sum()
    logging.info(f"Missing values after homogenization:\n{missing_values}")

    return df


def process_fits_file(fits_file):
    """
    Processes a single FITS file and extracts its metadata.

    Parameters:
        fits_file (str): Path to the FITS file.

    Returns:
        dict: Extracted metadata.
    """
    logging.info(f"Processing file: {fits_file}")
    metadata = extract_fits_metadata(fits_file)
    # 'FITS_File' is already set in extract_fits_metadata
    return metadata

def process_fits_directory_parallel(directory_path, max_workers=10):
    """
    Processes all FITS files in a directory using parallel processing and extracts metadata.

    Parameters:
        directory_path (str): Path to the directory containing FITS files.
        max_workers (int): Maximum number of threads to use.

    Returns:
        pandas.DataFrame: DataFrame containing metadata for all FITS files.
    """
    from concurrent.futures import ThreadPoolExecutor, as_completed

    fits_pattern = os.path.join(directory_path, '*.fits')
    fits_files = glob.glob(fits_pattern)
    metadata_list = []

    if not fits_files:
        logging.warning(f"No FITS files found in directory: {directory_path}")

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_fits = {executor.submit(process_fits_file, fits_file): fits_file for fits_file in fits_files}
        for future in as_completed(future_to_fits):
            fits_file = future_to_fits[future]
            try:
                metadata = future.result()
                metadata_list.append(metadata)
            except Exception as e:
                logging.error(f"Error processing {fits_file}: {e}")

    df = pd.DataFrame(metadata_list)
    df = homogenize_metadata(df)

    # Debugging: Log the columns present in the DataFrame
    logging.debug(f"DataFrame columns: {df.columns.tolist()}")

    return df

def save_metadata_to_csv(df, output_file):
    """
    Saves the metadata DataFrame to a CSV file.

    Parameters:
        df (pandas.DataFrame): DataFrame containing metadata.
        output_file (str): Path to the output CSV file.
    """
    try:
        df.to_csv(output_file, index=False)
        logging.info(f"Metadata successfully saved to {output_file}.")
    except Exception as e:
        logging.error(f"Failed to save metadata to CSV: {e}")

def skycoord_to_unit_vector(ra, dec):
    """
    Converts RA, Dec in degrees to a unit vector in 3D Cartesian coordinates.

    Parameters:
        ra (float): Right Ascension in degrees.
        dec (float): Declination in degrees.

    Returns:
        numpy.ndarray: Unit vector [x, y, z]
    """
    ra_rad = np.deg2rad(ra)
    dec_rad = np.deg2rad(dec)
    x = np.cos(dec_rad) * np.cos(ra_rad)
    y = np.cos(dec_rad) * np.sin(ra_rad)
    z = np.sin(dec_rad)
    return np.array([x, y, z])

def is_point_in_convex_polygon(ra, dec, polygon_coords):
    """
    Determines if a point (ra, dec) is inside a convex polygon defined by polygon_coords.

    Parameters:
        ra (float): Right Ascension of the point in degrees.
        dec (float): Declination of the point in degrees.
        polygon_coords (list of tuples): List of (RA, Dec) tuples representing the polygon vertices in clockwise or counter-clockwise order.

    Returns:
        bool: True if the point is inside the polygon or on its edge, False otherwise.
    """
    try:
        if not polygon_coords or len(polygon_coords) < 3:
            logging.error("Polygon must have at least 3 vertices.")
            return False

        # Validate that each coordinate is a tuple (RA, Dec)
        for idx, coord in enumerate(polygon_coords):
            if not isinstance(coord, (list, tuple)) or len(coord) != 2:
                logging.error(f"Invalid coordinate at index {idx}: {coord}. Each coordinate must be a tuple of (RA, Dec).")
                return False

        # Extract RA and Dec from polygon coordinates
        ras = [coord[0] for coord in polygon_coords]
        decs = [coord[1] for coord in polygon_coords]

        # Validate RA and Dec ranges
        for ra_val, dec_val in zip(ras, decs):
            if not (0.0 <= ra_val < 360.0):
                logging.error(f"Invalid RA value: {ra_val}. Must be in [0, 360).")
                return False
            if not (-90.0 <= dec_val <= 90.0):
                logging.error(f"Invalid Dec value: {dec_val}. Must be in [-90, 90].")
                return False

        # Close the polygon by adding the first point at the end
        ras_closed = ras + [ras[0]]
        decs_closed = decs + [decs[0]]

        # Initialize list to store cross product signs
        cross_signs = []

        for i in range(len(ras)):
            # Current vertex
            x1, y1 = ras_closed[i], decs_closed[i]
            # Next vertex
            x2, y2 = ras_closed[i + 1], decs_closed[i + 1]

            # Compute delta_ra with RA wrapping for the edge
            delta_ra = x2 - x1
            if delta_ra > 180:
                delta_ra -= 360
            elif delta_ra < -180:
                delta_ra += 360

            # Adjusted next RA
            x2_adjusted = x1 + delta_ra

            # Edge vector
            edge = np.array([x2_adjusted - x1, y2 - y1])

            # Compute delta_ra with RA wrapping for the point
            delta_ra_point = ra - x1
            if delta_ra_point > 180:
                delta_ra_point -= 360
            elif delta_ra_point < -180:
                delta_ra_point += 360

            # Adjusted point RA
            adjusted_point_ra = x1 + delta_ra_point

            # Vector from current vertex to point
            vector = np.array([adjusted_point_ra - x1, dec - y1])

            # Compute the 2D cross product (scalar)
            cross = edge[0] * vector[1] - edge[1] * vector[0]

            # Store the sign of the cross product
            cross_signs.append(np.sign(cross))

        cross_signs = np.array(cross_signs)

        # Determine if all signs are non-negative or non-positive
        if np.all(cross_signs >= 0) or np.all(cross_signs <= 0):
            return True
        else:
            return False

    except Exception as e:
        logging.error(f"Error in is_point_in_convex_polygon: {e}")
        return False

def plot_moc_and_polygon(polygon_coords, moc_str, title="MOC_Debug", fits_file=None):
    """
    Plots the MOC, Polygon, and optionally FITS data on the same WCS projection,
    handling both Equatorial (RA/Dec) and Galactic (GLON/GLAT) coordinate systems.

    Parameters:
        polygon_coords (list of lists/tuples or str): 
            - List of [RA, Dec] or [GLON, GLAT] pairs defining the polygon.
            - Can also be a JSON string representing the list of coordinates.
        moc_str (str): 
            - Serialized MOC string.
        title (str): 
            - Title for the plot.
        fits_file (str, optional): 
            - Path to the FITS file to be plotted.
            - If provided, the FITS data will be reprojected and overlaid on the plot.

    Returns:
        None
    """
    # Suppress warnings for FITS header parsing
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', FITSFixedWarning)
        warnings.simplefilter('ignore', AstropyWarning)
        try:
            # Determine the coordinate frame based on the FITS file's header
            if fits_file:
                with fits.open(fits_file) as hdul:
                    fits_header = hdul[0].header
                    fits_wcs = WCS(fits_header)
                    ctype1 = fits_wcs.wcs.ctype[0].upper()
                    ctype2 = fits_wcs.wcs.ctype[1].upper()

                    if 'RA' in ctype1 and 'DEC' in ctype2:
                        coord_frame = 'icrs'  # Equatorial Coordinates
                    elif 'GLON' in ctype1 and 'GLAT' in ctype2:
                        coord_frame = 'galactic'  # Galactic Coordinates
                    else:
                        logging.error(f"Unsupported coordinate system in {fits_file}: {ctype1}, {ctype2}")
                        return
            else:
                # Default to Equatorial if no FITS file is provided
                coord_frame = 'icrs'

            # Process polygon coordinates
            if isinstance(polygon_coords, str):
                try:
                    polygon_coords = json.loads(polygon_coords)
                    logging.debug("Polygon coordinates deserialized from JSON string.")
                except json.JSONDecodeError as e:
                    logging.error(f"Failed to decode Polygon_Coords: {e}")
                    return

            # Handle flat list of numbers by reshaping into list of [RA, Dec] or [GLON, GLAT] pairs
            if isinstance(polygon_coords, list):
                if all(isinstance(item, (int, float)) for item in polygon_coords):
                    # Flat list detected, reshape into list of pairs
                    if len(polygon_coords) % 2 != 0:
                        logging.error(f"Polygon_Coords list length is not even: {polygon_coords}")
                        return
                    else:
                        polygon_coords = [polygon_coords[i:i+2] for i in range(0, len(polygon_coords), 2)]
                        logging.debug("Polygon coordinates reshaped from flat list to [RA/GLON, Dec/GLAT] pairs.")
                elif not all(isinstance(item, (list, tuple)) and len(item) == 2 for item in polygon_coords):
                    logging.error(f"Invalid Polygon_Coords format: {polygon_coords}")
                    return
            else:
                logging.error(f"Polygon_Coords is not a list or JSON string: {polygon_coords}")
                return

            # Extract RA/GLON and Dec/GLAT from polygon coordinates
            try:
                x_coord, y_coord = zip(*polygon_coords)
                x_coord = np.array(x_coord, dtype=np.float64)
                y_coord = np.array(y_coord, dtype=np.float64)
                logging.debug(f"Extracted X and Y coordinates from polygon: {len(x_coord)} points.")
            except ValueError as ve:
                logging.error(f"Invalid polygon coordinates format: {ve}")
                return

            # Check for NaN or infinite values
            if not (np.isfinite(x_coord).all() and np.isfinite(y_coord).all()):
                logging.error("RA/GLON or Dec/GLAT contains non-finite values.")
                return

            # Close the polygon by adding the first point at the end if not already closed
            if (x_coord[0], y_coord[0]) != (x_coord[-1], y_coord[-1]):
                x_coord = np.append(x_coord, x_coord[0])
                y_coord = np.append(y_coord, y_coord[0])
                logging.debug("Closed the polygon by appending the first coordinate to the end.")

            # Create SkyCoord object for the polygon
            if coord_frame == 'icrs':
                sky_coords = SkyCoord(ra=x_coord * u.deg, dec=y_coord * u.deg, frame='icrs')
            elif coord_frame == 'galactic':
                sky_coords = SkyCoord(l=x_coord * u.deg, b=y_coord * u.deg, frame='galactic')
            else:
                logging.error("Unknown coordinate frame. Cannot create SkyCoord object.")
                return

            # Create a PolygonSkyRegion
            polygon_region = PolygonSkyRegion(vertices=sky_coords)

            # Convert to Shapely polygon manually
            shapely_polygon = polygonskyregion_to_shapely(polygon_region)
            if shapely_polygon:
                polygon_wkt = shapely_polygon.wkt
                polygon_coords_list = list(shapely_polygon.exterior.coords)
                logging.debug("Polygon WKT and coordinates extracted successfully.")
            else:
                logging.error("Failed to convert PolygonSkyRegion to Shapely Polygon.")
                return

            # Deserialize the MOC string
            try:
                moc = MOC.from_string(moc_str)
                logging.debug("Successfully deserialized MOC.")
            except Exception as e:
                logging.error(f"Failed to deserialize MOC: {e}")
                return

            # Initialize FITS data variables
            fits_data_reprojected = None
            im = None

            # If FITS file is provided, load the data
            if fits_file:
                try:
                    with fits.open(fits_file) as hdul:
                        fits_data = hdul[0].data
                        fits_header = hdul[0].header
                        fits_wcs = WCS(fits_header)
                    logging.debug("FITS data loaded successfully.")
                except Exception as e:
                    logging.error(f"Failed to load FITS file '{fits_file}': {e}")
                    fits_file = None  # Proceed without FITS data

            # Create a matplotlib figure
            fig = plt.figure(figsize=(10, 5))

            try:
                # Obtain WCS from MOC by passing the figure
                wcs = moc.wcs(fig)

                # Use the WCS for the projection
                ax = fig.add_subplot(111, projection=wcs)
                ax.grid(True)

                # Determine shape_out based on figure size and DPI
                dpi = fig.get_dpi()
                width_px = int(fig.get_figwidth() * dpi)
                height_px = int(fig.get_figheight() * dpi)
                shape_out = (height_px, width_px)
                logging.debug(f"Determined shape_out for reproject: {shape_out}")

                # If FITS data is available, reproject it to match MOC's WCS
                if fits_file:
                    try:
                        fits_data_reprojected, footprint = reproject_interp(
                            (fits_data, fits_wcs), 
                            wcs, 
                            shape_out=shape_out
                        )
                        logging.debug("FITS data reprojected successfully to match MOC's WCS.")

                        # Plot the reprojected FITS data using imshow
                        im = ax.imshow(
                            fits_data_reprojected, 
                            origin='lower', 
                            cmap='gray', 
                            aspect='auto', 
                            interpolation='none', 
                            alpha=0.7, 
                            label='FITS Data'
                        )
                        logging.debug("FITS data plotted successfully.")
                    except Exception as e:
                        logging.error(f"Failed to reproject or plot FITS data: {e}")

                # Plot the MOC
                moc.fill(ax=ax, wcs=wcs, alpha=0.3, color="blue", label="MOC")
                moc.border(ax=ax, wcs=wcs, color="red", linewidth=1.5)
                logging.debug("MOC filled and bordered successfully.")

                # Plot the polygon
                if coord_frame == 'icrs':
                    transform = ax.get_transform('icrs')
                    xlabel = 'Right Ascension (J2000)'
                    ylabel = 'Declination (J2000)'
                elif coord_frame == 'galactic':
                    transform = ax.get_transform('galactic')
                    xlabel = 'Galactic Longitude (GLON)'
                    ylabel = 'Galactic Latitude (GLAT)'
                else:
                    logging.error("Unknown coordinate frame. Cannot plot polygon.")
                    transform = None
                    xlabel = 'X Coordinate'
                    ylabel = 'Y Coordinate'

                if transform:
                    ax.plot(
                        x_coord, 
                        y_coord, 
                        transform=transform, 
                        color='green', 
                        linewidth=2, 
                        label="Polygon"
                    )
                    logging.debug("Polygon plotted successfully.")

                    # Set axis labels based on coordinate frame
                    ax.set_xlabel(xlabel)
                    ax.set_ylabel(ylabel)
                    logging.debug(f"Axis labels set to '{xlabel}' and '{ylabel}'.")

                # Add legend and title
                ax.legend(loc='upper right')
                plt.title(title)

                # Add a colorbar for FITS data if plotted
                if fits_file and im is not None:
                    try:
                        cbar = plt.colorbar(im, ax=ax, orientation='vertical', pad=0.02)
                        cbar.set_label('FITS Data Intensity')
                        logging.debug("Colorbar added successfully.")
                    except Exception as e:
                        logging.error(f"Failed to add colorbar: {e}")

                # Save the plot as a PNG file
                try:
                    filename = f"{title.replace(' ', '_')}.png"
                    plt.savefig(filename, dpi=300, bbox_inches='tight')
                    logging.info(f"Plot saved as {filename}")
                except Exception as e:
                    logging.error(f"Error saving plot: {e}")

                # Optionally display the plot
                # plt.show()

            except Exception as e:
                logging.error(f"Error plotting MOC, polygon, and FITS data: {e}")
                plt.close(fig)
                return

            finally:
                # Close the figure to free memory
                plt.close(fig)

        except Exception as e:
            logging.error(f"Error plotting MOC, polygon, and FITS data: {e}")
            plt.close(fig)
            return
        
def plot_moc_and_polygon_from_dataset_notebook(
    metadata_df, 
    input_dir, 
    output_dir='plots', 
    max_plots=None,
    fits_files=None,
    indices=None,
    filter_func=None
):
    """
    Facilitates the plotting of MOC and Polygon regions for selected FITS files in the metadata dataset.
    
    This function iterates over each selected entry in the provided metadata DataFrame, retrieves the corresponding
    FITS file, and generates a plot of its MOC and Polygon region. The plots are saved as PNG files in
    the specified output directory.
    
    Parameters:
        metadata_df (pandas.DataFrame): DataFrame containing metadata for FITS files.
            Expected columns:
                - 'FITS_File': Basename of the FITS file.
                - 'Polygon_Coords': JSON string of polygon coordinates.
                - 'MOC': Serialized MOC string.
        input_dir (str): Directory path where the FITS files are located.
        output_dir (str, optional): Directory to save the generated plots. Defaults to 'plots'.
        max_plots (int, optional): Maximum number of plots to generate. Useful for testing.
            If None, plots all selected FITS files.
        fits_files (list of str, optional): List of FITS filenames to plot.
            If provided, only these FITS files will be plotted.
        indices (list of int, optional): List of DataFrame indices to plot.
            If provided, only these rows will be plotted.
        filter_func (callable, optional): A function that takes a DataFrame row and returns True
            if the row should be plotted, False otherwise.
            Example: lambda row: row['Resolved_Object'] == 'SN2023abc'

    Returns:
        None
    """
    import os
    import json
    import logging
    from tqdm.notebook import tqdm  # For progress bar in Jupyter Notebook
    
    # Ensure the plotting function is accessible
    if 'plot_moc_and_polygon' not in globals():
        raise ImportError("The function 'plot_moc_and_polygon' must be defined or imported before using this plotting function.")
    
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    logging.info(f"Plots will be saved to: {output_dir}")
    
    # Select rows based on provided parameters
    if fits_files:
        selected_df = metadata_df[metadata_df['FITS_File'].isin(fits_files)]
        if selected_df.empty:
            logging.warning("No FITS files matched the provided filenames.")
    elif indices:
        selected_df = metadata_df.iloc[indices]
        if selected_df.empty:
            logging.warning("No rows matched the provided indices.")
    elif filter_func:
        selected_df = metadata_df[metadata_df.apply(filter_func, axis=1)]
        if selected_df.empty:
            logging.warning("No rows matched the filter criteria.")
    else:
        # If no selection criteria provided, select all
        selected_df = metadata_df
        logging.info("No specific selection criteria provided. All FITS files will be plotted.")
    
    # Determine the number of plots to generate
    total_plots = len(selected_df) if max_plots is None else min(len(selected_df), max_plots)
    logging.info(f"Generating {total_plots} plots.")
    
    # Iterate over the selected DataFrame rows with a progress bar
    for idx, row in tqdm(selected_df.head(total_plots).iterrows(), total=total_plots, desc="Generating Plots"):
        fits_file = row.get('FITS_File', None)
        polygon_coords_str = row.get('Polygon_Coords', None)
        moc_str = row.get('MOC', None)
        
        if pd.isnull(fits_file):
            logging.warning(f"Row {idx}: 'FITS_File' is missing. Skipping.")
            continue
        
        if pd.isnull(polygon_coords_str):
            logging.warning(f"Row {idx}: 'Polygon_Coords' is missing for '{fits_file}'. Skipping.")
            continue
        
        if pd.isnull(moc_str):
            logging.warning(f"Row {idx}: 'MOC' is missing for '{fits_file}'. Skipping.")
            continue
        
        # Construct the full path to the FITS file
        fits_path = os.path.join(input_dir, fits_file)
        if not os.path.isfile(fits_path):
            logging.error(f"FITS file not found: {fits_path}. Skipping.")
            continue
        
        # Define a title for the plot
        title = os.path.splitext(fits_file)[0]
        
        try:
            # Generate the plot using the existing plot_moc_and_polygon function
            plot_moc_and_polygon(
                polygon_coords=polygon_coords_str,
                moc_str=moc_str,
                title=title,
                fits_file=fits_path
            )
            
            # Move the generated plot to the output directory
            plot_filename = f"{title.replace(' ', '_')}.png"
            if os.path.exists(plot_filename):
                destination = os.path.join(output_dir, plot_filename)
                os.rename(plot_filename, destination)
                logging.info(f"Plot saved: {destination}")
            else:
                logging.error(f"Plot file was not created for '{fits_file}'.")
        
        except Exception as e:
            logging.error(f"Failed to plot '{fits_file}': {e}")
            continue

    
def search_fits_by_point(metadata_df, ra, dec):
    """
    Searches for FITS files in the collection whose coverage includes the given point.

    Parameters:
        metadata_df (pandas.DataFrame): DataFrame containing FITS metadata, including 'Polygon_Coords' column.
        ra (float): Right Ascension of the point in degrees.
        dec (float): Declination of the point in degrees.

    Returns:
        pandas.DataFrame: Subset of metadata_df with FITS files that include the point.
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

    Parameters:
        metadata_df (pandas.DataFrame): DataFrame containing FITS metadata, including 'Polygon_Coords' column.
        region (dict): Dictionary defining the region. Expected keys:
            - 'type': 'circle' or 'polygon'
            - For 'circle':
                - 'center': (ra, dec) in degrees
                - 'radius': radius in degrees
            - For 'polygon':
                - 'coordinates': list of (ra, dec) pairs in degrees

    Returns:
        pandas.DataFrame: Subset of metadata_df with FITS files that intersect the region.
    """
    from shapely.geometry import Point, Polygon

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


def load_metadata_from_csv(csv_file):
    """
    Loads the metadata DataFrame from a CSV file.

    Parameters:
        csv_file (str): Path to the CSV file containing the metadata.

    Returns:
        pandas.DataFrame: Loaded metadata DataFrame.
    """
    try:
        df = pd.read_csv(csv_file)
        logging.info(f"Metadata loaded from {csv_file}.")
        return df
    except Exception as e:
        logging.error(f"Failed to load metadata from CSV: {e}")
        return pd.DataFrame()


def main():
    """
    Main function to execute the metadata extraction and homogenization process.
    """
    # Lock file setup
    LOCK_FILE = 'fits_metadata_extractor.lock'
    lock = FileLock(LOCK_FILE, timeout=1)

    try:
        with lock:
            logging.info("Starting FITS metadata extraction process.")

            parser = argparse.ArgumentParser(description="FITS Header Extraction and Metadata Homogenization")
            parser.add_argument(
                '-i', '--input_dir',
                type=str,
                required=True,
                help='Path to the directory containing FITS files.'
            )
            parser.add_argument(
                '-o', '--output_csv',
                type=str,
                default='metadata_dataset.csv',
                help='Path to the output CSV file.'
            )

            args = parser.parse_args()

            input_dir = args.input_dir
            output_csv = args.output_csv

            logging.info(f"Input directory: {input_dir}")
            logging.info(f"Output CSV: {output_csv}")

            # Process FITS files using parallel processing
            metadata_df = process_fits_directory_parallel(input_dir, max_workers=10)

            # Save the metadata to CSV
            save_metadata_to_csv(metadata_df, output_csv)

            logging.info("FITS metadata extraction and homogenization process completed successfully.")

    except LockTimeout:
        logging.error("Another instance of the script is running. Exiting.")
        exit(1)
    except Exception as e:
        logging.error("An unexpected error occurred:", exc_info=True)  # Updated to log traceback
        exit(1)

if __name__ == "__main__":
    main()
