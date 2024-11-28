
# FITS Metadata Extractor

*A Python package for processing, visualizing, and searching metadata from FITS files.*

---

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)
    - [Prerequisites](#prerequisites)
    - [Installing from Source](#installing-from-source)
- [Usage](#usage)
    - [Processing FITS Files](#processing-fits-files)
    - [Saving and Loading Metadata](#saving-and-loading-metadata)
    - [Plotting MOCs and Polygons](#plotting-mocs-and-polygons)
    - [Searching FITS Files](#searching-fits-files)
        - [Search by Point](#search-by-point)
        - [Search by Circular Region](#search-by-circular-region)
        - [Search by Polygonal Region](#search-by-polygonal-region)
    - [Plotting Search Results](#plotting-search-results)
- [API Reference](#api-reference)
- [Contributing](#contributing)
- [License](#license)

---

## Introduction

**FITS Metadata Extractor** is a Python package designed to streamline the processing, visualization, and searching of metadata from FITS (Flexible Image Transport System) files. It provides tools to:

- Extract metadata from FITS files in parallel.
- Save and load metadata to and from CSV files.
- Generate and plot Multi-Order Coverage maps (MOCs) and polygons.
- Perform spatial searches (point, circular, and polygonal) on FITS metadata.
- Visualize search results by plotting matching FITS file coverages.

---

## Features

- **Parallel Processing**: Efficiently process large directories of FITS files using multithreading.
- **Metadata Extraction**: Extract key metadata fields, including WCS information, object names, and coordinate frames.
- **Visualization**: Generate MOCs and polygons for FITS files and visualize them using Matplotlib.
- **Spatial Search**: Search for FITS files covering specific points or regions in the sky.
- **Integration with Astronomical Databases**: Resolve object names using Simbad, NED, and VizieR services via `astroquery`.

---

## Installation

### Prerequisites

- **Python 3.7 or higher**: Ensure you have a compatible Python version installed.
- **pip**: Python package installer.

### Installing from Source

1. **Clone the Repository**

   ```bash
   git clone https://github.com/pierregab/fits_metadata_extractor.git
   cd fits_metadata_extractor
   ```

2. **Install Dependencies**

   Install the required packages listed in `requirements.txt`:

   ```bash
   pip install -r requirements.txt
   ```

3. **Install the Package**

   Install the package locally using `setup.py`:

   ```bash
   pip install .
   ```

   *For development purposes, you can install it in editable mode:*

   ```bash
   pip install -e .
   ```

---

## Usage

Below is an example of how to use the `fits_metadata_extractor` package in a Python script or Jupyter Notebook.

### Processing FITS Files

```python
import os
import logging
from fits_metadata_extractor.processor import FITSProcessor
from fits_metadata_extractor.logger import setup_logging

# Initialize and configure logging
logger = setup_logging()

# Define the directory containing FITS files
fits_directory = "fits_collection"  # <-- Replace with your FITS files directory

# Check if the directory exists
if not os.path.isdir(fits_directory):
    logger.error(f"The specified FITS directory does not exist: {fits_directory}")
else:
    logger.info(f"FITS directory found: {fits_directory}")

# Instantiate FITSProcessor with desired number of workers
processor = FITSProcessor(max_workers=5)

# Process the FITS directory
metadata_df = processor.process_fits_directory_parallel(fits_directory)

# Display the first few rows of the metadata DataFrame
print(metadata_df.head())
```

### Saving and Loading Metadata

Save the extracted metadata to a CSV file:

```python
from fits_metadata_extractor.utils import save_metadata_to_csv

# Define the output CSV file path
output_csv = 'metadata.csv'

# Save the metadata DataFrame to CSV
save_metadata_to_csv(metadata_df, output_csv)
```

Load the metadata from the CSV file:

```python
from fits_metadata_extractor.utils import load_metadata_from_csv

# Load the metadata DataFrame from CSV
metadata_loaded_df = load_metadata_from_csv(output_csv)
```

### Plotting MOCs and Polygons

```python
from fits_metadata_extractor.plotter import plot_moc_and_polygon_from_dataset_notebook

# Define the directory to save plots
plot_output_dir = 'plots'

# Plot MOCs and polygons
plot_moc_and_polygon_from_dataset_notebook(
    metadata_df=metadata_loaded_df,
    input_dir=fits_directory,
    output_dir=plot_output_dir,
    max_plots=10,       # Adjust as needed
    plot_directly=True  # Set to True to display plots in a notebook
)
```

### Searching FITS Files

#### Search by Point

```python
from fits_metadata_extractor.search import search_fits_by_point

# Define a test point (RA, Dec in degrees)
test_point = {'ra': 270.0, 'dec': -30.0}

# Perform the search
matching_fits_by_point = search_fits_by_point(metadata_loaded_df, test_point['ra'], test_point['dec'])

# Display results
if not matching_fits_by_point.empty:
    print("FITS files containing the point:")
    print(matching_fits_by_point[['FITS_File', 'Resolved_Object']])
else:
    print("No FITS files found containing the specified point.")
```

#### Search by Circular Region

```python
from fits_metadata_extractor.search import search_fits_by_region

# Define a circular region
test_circle = {
    'type': 'circle',
    'center': (270.0, -29.0),  # RA, Dec in degrees
    'radius': 5.0              # Radius in degrees
}

# Perform the search
matching_fits_by_circle = search_fits_by_region(metadata_loaded_df, test_circle)

# Display results
if not matching_fits_by_circle.empty:
    print("FITS files intersecting the circular region:")
    print(matching_fits_by_circle[['FITS_File', 'Resolved_Object']])
else:
    print("No FITS files found intersecting the specified circular region.")
```

#### Search by Polygonal Region

```python
# Define a polygonal region
test_polygon = {
    'type': 'polygon',
    'coordinates': [
        (180.0, -60.0),  # Bottom-Left Corner
        (180.0, 60.0),   # Top-Left Corner
        (280.0, 60.0),   # Top-Right Corner
        (280.0, -60.0)   # Bottom-Right Corner
    ]
}

# Perform the search
matching_fits_by_polygon = search_fits_by_region(metadata_loaded_df, test_polygon)

# Display results
if not matching_fits_by_polygon.empty:
    print("FITS files intersecting the polygonal region:")
    print(matching_fits_by_polygon[['FITS_File', 'Resolved_Object']])
else:
    print("No FITS files found intersecting the specified polygonal region.")
```

### Plotting Search Results

```python
from fits_metadata_extractor.plotter import plot_search_region_and_find_fits

# Define a search region (e.g., a circle)
search_region = {
    'type': 'circle',
    'center': (220.0, 0.0),
    'radius': 60.0
}

# Plot the search region and matching FITS files
plot_search_region_and_find_fits(
    metadata_df=metadata_loaded_df,
    region=search_region,
    input_dir=fits_directory,
    output_dir='search_plots',
    max_plots=10,
    plot_search_region=True
)
```

---

## API Reference

### `FITSProcessor`

- **Description**: Class for processing FITS files in parallel.

- **Methods**:
    - `process_fits_directory_parallel(fits_directory)`: Processes all FITS files in a directory using multithreading.

### `utils`

- **Functions**:
    - `save_metadata_to_csv(metadata_df, output_csv)`: Saves the metadata DataFrame to a CSV file.
    - `load_metadata_from_csv(input_csv)`: Loads metadata from a CSV file into a DataFrame.

### `search`

- **Functions**:
    - `search_fits_by_point(metadata_df, ra, dec)`: Searches for FITS files covering a specific point.
    - `search_fits_by_region(metadata_df, region)`: Searches for FITS files intersecting a region (circle or polygon).

### `plotter`

- **Functions**:
    - `plot_moc_and_polygon_from_dataset_notebook(...)`: Plots MOCs and polygons from the metadata DataFrame.
    - `plot_search_region_and_find_fits(...)`: Plots the search region and matching FITS file coverages.

---

## Contributing

Contributions are welcome! If you have suggestions, bug reports, or patches, please open an issue or submit a pull request on the [GitHub repository](https://github.com/yourusername/fits_metadata_extractor).

### Steps to Contribute

1. **Fork the Repository**

   Click the "Fork" button on the top right of the repository page.

2. **Clone Your Fork**

   ```bash
   git clone https://github.com/pierregab/fits_metadata_extractor.git
   cd fits_metadata_extractor
   ```

3. **Create a New Branch**

   ```bash
   git checkout -b feature/your-feature-name
   ```

4. **Make Changes**

   Implement your changes or additions.

5. **Commit and Push**

   ```bash
   git add .
   git commit -m "Description of your changes"
   git push origin feature/your-feature-name
   ```

6. **Submit a Pull Request**

   Go to your forked repository on GitHub and click "New pull request".

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Acknowledgments

- **Astropy**: For providing core astronomical utilities.
- **Astroquery**: For facilitating queries to astronomical databases.
- **MOCpy**: For handling Multi-Order Coverage maps.
- **Matplotlib**: For plotting and visualization.
- **Shapely**: For geometric operations.

---

**Disclaimer**: This package is provided as-is without any guarantees. Use it at your own risk.

