# fits_metadata_extractor/processor.py

import os
import glob
import logging
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from .metadata_extractor import FITSMetadataExtractor
from .utils import homogenize_metadata

class FITSProcessor:
    """
    Class for processing directories of FITS files and extracting metadata.
    """

    def __init__(self, max_workers=10, extensions=None):
        self.max_workers = max_workers
        self.extractor = FITSMetadataExtractor()
        # Default extensions if none provided
        self.extensions = extensions if extensions else ['.fits', '.fit']

    def process_fits_file(self, fits_file):
        """
        Processes a single FITS file and extracts its metadata.

        Parameters:
            fits_file (str): Path to the FITS file.

        Returns:
            dict: Extracted metadata.
        """
        logging.info(f"Processing file: {fits_file}")
        metadata = self.extractor.extract_fits_metadata(fits_file)
        return metadata

    def process_fits_directory_parallel(self, directory_path):
        """
        Processes all FITS files in a directory using parallel processing and extracts metadata.

        Parameters:
            directory_path (str): Path to the directory containing FITS files.

        Returns:
            pandas.DataFrame: DataFrame containing metadata for all FITS files.
        """
        fits_files = []
        for ext in self.extensions:
            pattern = os.path.join(directory_path, f'*{ext}')
            fits_files.extend(glob.glob(pattern))
        
        metadata_list = []

        if not fits_files:
            logging.warning(f"No FITS files with extensions {self.extensions} found in directory: {directory_path}")

        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            future_to_fits = {executor.submit(self.process_fits_file, fits_file): fits_file for fits_file in fits_files}
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