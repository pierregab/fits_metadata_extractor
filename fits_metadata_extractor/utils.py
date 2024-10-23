# fits_metadata_extractor/utils.py

import logging
import pandas as pd
import json
from astropy.coordinates import Angle

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
        if isinstance(x, (list, tuple)):
            return json.dumps(x)
        else:
            return None
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