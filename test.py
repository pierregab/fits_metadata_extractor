import os
import logging
import pandas as pd
import matplotlib.pyplot as plt

# Import classes and functions from the fits_metadata_extractor package
from fits_metadata_extractor.processor import FITSProcessor
from fits_metadata_extractor.utils import save_metadata_to_csv, load_metadata_from_csv
from fits_metadata_extractor.search import search_fits_by_point, search_fits_by_region
from fits_metadata_extractor.plotter import plot_moc_and_polygon_from_dataset_notebook, test, plot_search_region_and_find_fits
from fits_metadata_extractor.logger import setup_logging

# Define the output CSV file path
output_csv = 'test_metadata.csv'

# Load the metadata from CSV
metadata_df = load_metadata_from_csv(output_csv)

# Define a search region (e.g., a circle)
search_region = {
    'type': 'circle',
    'center': (270, 0),  # RA=150.0°, Dec=2.2°
    'radius': 20           # 5 degrees radius
}

# Call the plotting function
plot_search_region_and_find_fits(
    metadata_df=metadata_df,
    region=search_region,
    input_dir='fits_collection',
    output_dir='search_plots',
    max_plots=10,
    plot_search_region=True
)