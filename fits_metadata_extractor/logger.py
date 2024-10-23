# fits_metadata_extractor/logger.py

import logging

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

def log_success(message):
    """
    Logs a success message in green.

    Parameters:
        message (str): The message to log.
    """
    extra = {'is_success': True}
    logging.info(message, extra=extra)
