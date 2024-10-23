# fits_metadata_extractor/main.py

import logging
import argparse
from filelock import FileLock, Timeout as LockTimeout
from .logger import setup_logging
from .processor import FITSProcessor
from .utils import save_metadata_to_csv

def main():
    """
    Main function to execute the metadata extraction and homogenization process.
    """
    # Initialize logger
    logger = setup_logging()

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
            processor = FITSProcessor(max_workers=10)
            metadata_df = processor.process_fits_directory_parallel(input_dir)

            # Save the metadata to CSV
            save_metadata_to_csv(metadata_df, output_csv)

            logging.info("FITS metadata extraction and homogenization process completed successfully.")

    except LockTimeout:
        logging.error("Another instance of the script is running. Exiting.")
        exit(1)
    except Exception as e:
        logging.error("An unexpected error occurred:", exc_info=True)
        exit(1)

if __name__ == "__main__":
    main()
