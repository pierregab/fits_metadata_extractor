#!/usr/bin/env python3
# test_fits_metadata_extractor.py

"""
Test Script for FITS Metadata Extractor Module

This script processes a collection of FITS files, extracts and saves metadata,
plots MOCs and polygons, and tests the search functionalities by point and region.

Author: Your Name
Date: October 23, 2024
"""

import os
import logging
import pandas as pd

# Import the necessary classes and functions from the new structure
from fits_metadata_extractor.processor import FITSProcessor
from fits_metadata_extractor.utils import save_metadata_to_csv, load_metadata_from_csv
from fits_metadata_extractor.search import search_fits_by_point, search_fits_by_region
from fits_metadata_extractor.plotter import plot_moc_and_polygon_from_dataset_notebook

def main():
    """
    Main function for processing FITS files, plotting MOCs and polygons,
    and testing search functionalities.
    """

    # Specify the directory containing FITS files
    fits_directory = "fits_collection"  # <-- CHANGE THIS PATH TO YOUR FITS FILES DIRECTORY
    output_csv = 'metadata.csv'
    plot_output_dir = 'plots'  # Directory where plots will be saved

    # Check if the FITS directory exists
    if not os.path.isdir(fits_directory):
        logging.error(f"The specified FITS directory does not exist: {fits_directory}")
        return

    # Process FITS files using FITSProcessor
    logging.info("Starting processing of FITS files.")
    processor = FITSProcessor(max_workers=5)
    metadata_df = processor.process_fits_directory_parallel(fits_directory)
    logging.info("Finished processing FITS files.")

    # Save the metadata to CSV
    save_metadata_to_csv(metadata_df, output_csv)
    logging.info(f"Metadata saved to {output_csv}.")

    # Load the metadata from CSV for further tests
    metadata_loaded_df = load_metadata_from_csv(output_csv)
    if metadata_loaded_df.empty:
        logging.error("Loaded metadata DataFrame is empty. Exiting tests.")
        return

    # Test plotting using plot_moc_and_polygon_from_dataset_notebook
    logging.info("Testing plotting using plot_moc_and_polygon_from_dataset_notebook.")

    try:
        # Plot all FITS files
        plot_moc_and_polygon_from_dataset_notebook(
            metadata_df=metadata_loaded_df,
            input_dir=fits_directory,
            output_dir=plot_output_dir,
            max_plots=10  # Limit to 10 plots for testing
        )
        logging.info("Finished plotting MOCs and polygons.")
    except Exception as e:
        logging.error(f"Error occurred during plotting: {e}")

    # Test the search functionality
    logging.info("Starting search functionality tests.")

    # Define test inputs for searches

    # 1. Search by Point
    test_point = {
        'ra': 150.0,  # Right Ascension in degrees
        'dec': 2.2    # Declination in degrees
    }
    logging.info(f"Testing search by point: RA={test_point['ra']}, Dec={test_point['dec']}")
    matching_fits_by_point = search_fits_by_point(metadata_loaded_df, test_point['ra'], test_point['dec'])
    if not matching_fits_by_point.empty:
        logging.info(f"FITS files containing the point (RA={test_point['ra']}, Dec={test_point['dec']}):")
        print(matching_fits_by_point[['FITS_File', 'Resolved_Object']])
    else:
        logging.info("No FITS files found containing the specified point.")

    # 2. Search by Circle
    test_circle = {
        'type': 'circle',
        'center': (150.0, 2.2),  # Center RA and Dec in degrees
        'radius': 1.0             # Radius in degrees
    }
    logging.info(f"Testing search by circular region: Center RA={test_circle['center'][0]}, Dec={test_circle['center'][1]}, Radius={test_circle['radius']} degrees")
    matching_fits_by_circle = search_fits_by_region(metadata_loaded_df, test_circle)
    if not matching_fits_by_circle.empty:
        logging.info(f"FITS files intersecting the circular region (RA={test_circle['center'][0]}, Dec={test_circle['center'][1]}, Radius={test_circle['radius']} degrees):")
        print(matching_fits_by_circle[['FITS_File', 'Resolved_Object']])
    else:
        logging.info("No FITS files found intersecting the specified circular region.")

    # 3. Search by Polygon
    test_polygon = {
        'type': 'polygon',
        'coordinates': [
            (149.0, 1.0),
            (151.0, 1.0),
            (151.0, 3.0),
            (149.0, 3.0)
        ]
    }
    logging.info(f"Testing search by polygon region with coordinates: {test_polygon['coordinates']}")
    matching_fits_by_polygon = search_fits_by_region(metadata_loaded_df, test_polygon)
    if not matching_fits_by_polygon.empty:
        logging.info("FITS files intersecting the polygonal region:")
        print(matching_fits_by_polygon[['FITS_File', 'Resolved_Object']])
    else:
        logging.info("No FITS files found intersecting the specified polygonal region.")

    logging.info("Search tests completed.")

if __name__ == "__main__":
    # Configure logging to show messages on the console
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    main()
