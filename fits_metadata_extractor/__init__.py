# fits_metadata_extractor/__init__.py

"""
FITS Metadata Extractor Package

Provides functionalities to extract and homogenize metadata from FITS headers.
"""

__version__ = "1.0.0"

from .processor import FITSProcessor
from .utils import save_metadata_to_csv, load_metadata_from_csv
from .search import search_fits_by_point, search_fits_by_region
from .plotter import plot_moc_and_polygon_from_dataset_notebook, plot_search_region_and_find_fits
from .logger import setup_logging