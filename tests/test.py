import os
import logging
import pandas as pd
import matplotlib.pyplot as plt

# Import classes and functions from the fits_metadata_extractor package
from fits_metadata_extractor.processor import FITSProcessor
from fits_metadata_extractor.utils import save_metadata_to_csv, load_metadata_from_csv
from fits_metadata_extractor.search import search_fits_by_point, search_fits_by_region
from fits_metadata_extractor.plotter import plot_moc_and_polygon_from_dataset_notebook, plot_search_region_and_find_fits
from fits_metadata_extractor.logger import setup_logging

# Define the output CSV file path
output_csv = 'test_metadata.csv'

# Load the metadata from CSV
metadata_df = load_metadata_from_csv(output_csv)

# Define a search region (e.g., a circle)
search_region = {
    'type': 'circle',
    'center': (220, 0),  
    'radius' : 60         
}


test_polygon = {
    'type': 'polygon',
    'coordinates': [
    (160, -60),  # Bottom-Left Corner
    (160, 60),   # Top-Left Corner
    (290, 60),   # Top-Right Corner
    (290, -60)   # Bottom-Right Corner
    ]
}

# Call the plotting function
plot_search_region_and_find_fits(
    metadata_df=metadata_df,
    region=test_polygon,
    input_dir='fits_collection',
    output_dir='search_plots',
    max_plots=10,
    plot_search_region=True
)