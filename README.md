# FITS Metadata Extractor

A Python tool for extracting and homogenizing metadata from FITS files, resolving astronomical object names, computing sky regions, and saving the metadata into a CSV file. It also provides functionalities to plot the MOCs (Multi-Order Coverage maps) and polygons of the FITS files and search FITS files by sky regions.

## Features

- **Metadata Extraction**: Extracts relevant metadata from FITS files, including WCS information, and homogenizes it for easy analysis.
- **Object Name Resolution**: Resolves astronomical object names using multiple services (Simbad, NED, VizieR, and Sesame), with support for custom mappings.
- **Sky Region Computation**: Computes sky regions (Polygon and MOC) from FITS files based on WCS information.
- **Metadata Storage**: Saves extracted metadata into a CSV file for further processing or analysis.
- **Plotting**: Provides functions to plot MOCs and Polygons of FITS files, with optional overlay of FITS data.
- **Search Functionality**: Allows searching FITS files by sky coordinates or regions.

## Installation

### Prerequisites

- **Python**: 3.6 or higher
- **pip**: Package installer for Python

### Install via GitHub

Clone the repository:

```bash
git clone https://github.com/yourusername/fits_metadata_extractor.git
cd fits_metadata_extractor
```

Install the required packages:

```bash
pip install -r requirements.txt
```

Alternatively, you can install the package using pip:

```bash
pip install git+https://github.com/yourusername/fits_metadata_extractor.git
```

### Dependencies

The project depends on the following Python packages:

- astropy
- numpy
- pandas
- mocpy
- shapely
- regions
- astroquery
- matplotlib
- filelock
- reproject
- tqdm

These will be installed automatically when you install the package using the instructions above.

## Usage

### Command Line Interface

Run the main script to extract metadata from FITS files:

```bash
python -m fits_metadata_extractor.main -i /path/to/fits/files -o metadata_dataset.csv
```

#### Options

- `-i`, `--input_dir`: **(Required)** Path to the directory containing FITS files.
- `-o`, `--output_csv`: **(Optional)** Path to the output CSV file. Default is `metadata_dataset.csv`.

### Examples

#### Extract Metadata from FITS Files

Extract metadata from FITS files in a directory and save it to a CSV file:

```bash
python -m fits_metadata_extractor.main -i ./fits_files -o metadata.csv
```

#### Plotting MOC and Polygon

You can use the `plotter.py` module to generate plots of MOC and Polygon regions.

```python
from fits_metadata_extractor.plotter import plot_moc_and_polygon

# Example usage
polygon_coords = '[[10.0, 20.0], [15.0, 25.0], [10.0, 30.0], [5.0, 25.0]]'
moc_str = '1/0-4'

plot_moc_and_polygon(polygon_coords, moc_str, title="Example Plot", fits_file='path/to/fits_file.fits')
```

#### Searching FITS Files by Sky Region

Search for FITS files that cover a given sky coordinate or region.

```python
from fits_metadata_extractor.search import search_fits_by_point, search_fits_by_region
from fits_metadata_extractor.utils import load_metadata_from_csv

# Load metadata from CSV
metadata_df = load_metadata_from_csv('metadata_dataset.csv')

# Search by point
ra = 150.0  # degrees
dec = 2.0   # degrees
matching_fits = search_fits_by_point(metadata_df, ra, dec)
print(matching_fits)

# Search by region (e.g., circle)
region = {
    'type': 'circle',
    'center': (150.0, 2.0),
    'radius': 1.0  # degrees
}
matching_fits = search_fits_by_region(metadata_df, region)
print(matching_fits)
```

## Documentation

Detailed documentation for the FITS Metadata Extractor can be found below.

### Table of Contents

- [Introduction](#introduction)
- [Modules](#modules)
  - [custom_mapping.py](#custom_mappingpy)
  - [logger.py](#loggerpy)
  - [main.py](#mainpy)
  - [metadata_extractor.py](#metadata_extractorpy)
  - [plotter.py](#plotterpy)
  - [polygon_utils.py](#polygon_utilspy)
  - [processor.py](#processorpy)
  - [resolver.py](#resolverpy)
  - [search.py](#searchpy)
  - [utils.py](#utilspy)
- [Examples](#examples)
- [License](#license)
- [Acknowledgments](#acknowledgments)

---

## Introduction

The FITS Metadata Extractor is a Python package designed to extract and homogenize metadata from FITS files. It provides tools to:

- Extract metadata from FITS headers.
- Resolve astronomical object names using various services.
- Compute sky regions (Polygon and MOC) based on WCS information.
- Plot MOCs and Polygons.
- Search FITS files by sky coordinates or regions.

This documentation provides detailed information on the modules, classes, and functions in the package.

## Modules

### custom_mapping.py

#### Overview

Contains custom mappings for astronomical object names. It allows you to define custom name resolutions for objects that may not be resolved by standard services.

#### Variables

- `CUSTOM_NAME_MAPPING`: A dictionary where keys are original object names and values are resolved names.

#### Usage

You can add your custom mappings to the `CUSTOM_NAME_MAPPING` dictionary. For example:

```python
CUSTOM_NAME_MAPPING = {
    'MyCustomObject': 'ResolvedNameForMyObject',
    'OriginalName': 'ResolvedName',
}
```

---

### logger.py

#### Overview

Sets up logging for the package. It configures logging to output to both a file and the console with appropriate formatting, including colored output for success messages.

#### Functions

- `setup_logging()`: Sets up logging handlers and formatters. Returns a logger instance.

  **Usage:**

  ```python
  from fits_metadata_extractor.logger import setup_logging

  logger = setup_logging()
  ```

- `log_success(message)`: Logs a success message in green color.

  **Usage:**

  ```python
  from fits_metadata_extractor.logger import log_success

  log_success("This is a success message.")
  ```

---

### main.py

#### Overview

Contains the entry point for the package when used as a script. It parses command-line arguments, processes FITS files, and saves the extracted metadata to a CSV file.

#### Functions

- `main()`: The main function that orchestrates the metadata extraction process.

#### Usage

Run the main script from the command line:

```bash
python -m fits_metadata_extractor.main -i /path/to/fits/files -o metadata_dataset.csv
```

---

### metadata_extractor.py

#### Overview

Contains the `FITSMetadataExtractor` class, which is responsible for extracting and homogenizing metadata from individual FITS files.

#### Classes

##### `FITSMetadataExtractor`

Extracts and homogenizes metadata from FITS files.

**Methods:**

- `__init__(self)`: Initializes the extractor, including the object name resolver.

- `infer_object_name(self, fits_file)`: Infers the object name from the FITS file name based on naming conventions.

  **Parameters:**

  - `fits_file (str)`: Path to the FITS file.

  **Returns:**

  - `str`: Inferred object name or `'Unknown'`.

- `extract_fits_metadata(self, fits_file)`: Extracts relevant metadata from a FITS header and computes the Polygon and MOC.

  **Parameters:**

  - `fits_file (str)`: Path to the FITS file.

  **Returns:**

  - `dict`: Dictionary containing extracted and homogenized metadata.

---

### plotter.py

#### Overview

Contains functions for plotting MOCs, Polygons, and FITS data.

#### Functions

- `plot_moc_and_polygon(polygon_coords, moc_str, title="MOC_Debug", fits_file=None)`

  **Description:**

  Plots the MOC, Polygon, and optionally FITS data on the same WCS projection, handling both Equatorial (RA/Dec) and Galactic (GLON/GLAT) coordinate systems.

  **Parameters:**

  - `polygon_coords (list or str)`: List of [RA, Dec] or [GLON, GLAT] pairs defining the polygon, or a JSON string representing the list of coordinates.
  - `moc_str (str)`: Serialized MOC string.
  - `title (str)`: Title for the plot.
  - `fits_file (str, optional)`: Path to the FITS file to be plotted.

  **Usage:**

  ```python
  from fits_metadata_extractor.plotter import plot_moc_and_polygon

  polygon_coords = '[[10.0, 20.0], [15.0, 25.0], [10.0, 30.0], [5.0, 25.0]]'
  moc_str = '1/0-4'
  plot_moc_and_polygon(polygon_coords, moc_str, title="Example Plot", fits_file='path/to/fits_file.fits')
  ```

- `plot_moc_and_polygon_from_dataset_notebook(...)`

  **Description:**

  Facilitates plotting MOCs and Polygons for selected FITS files in a metadata dataset.

  **Parameters:**

  - `metadata_df (pandas.DataFrame)`: DataFrame containing metadata for FITS files.
  - `input_dir (str)`: Directory path where the FITS files are located.
  - `output_dir (str, optional)`: Directory to save the generated plots. Defaults to `'plots'`.
  - `max_plots (int, optional)`: Maximum number of plots to generate.
  - `fits_files (list of str, optional)`: List of FITS filenames to plot.
  - `indices (list of int, optional)`: List of DataFrame indices to plot.
  - `filter_func (callable, optional)`: Function that takes a DataFrame row and returns True if the row should be plotted.

  **Usage:**

  ```python
  from fits_metadata_extractor.plotter import plot_moc_and_polygon_from_dataset_notebook
  from fits_metadata_extractor.utils import load_metadata_from_csv

  metadata_df = load_metadata_from_csv('metadata_dataset.csv')
  plot_moc_and_polygon_from_dataset_notebook(metadata_df, input_dir='./fits_files', max_plots=10)
  ```

---

### polygon_utils.py

#### Overview

Provides utility functions for working with polygons.

#### Functions

- `polygonskyregion_to_shapely(polygon_region)`

  **Description:**

  Manually converts a `PolygonSkyRegion` to a Shapely `Polygon`.

  **Parameters:**

  - `polygon_region (PolygonSkyRegion)`: The `PolygonSkyRegion` object.

  **Returns:**

  - `shapely.geometry.Polygon`: The equivalent Shapely `Polygon` or `None` if conversion fails.

  **Usage:**

  Used internally to convert polygon regions for further processing.

---

### processor.py

#### Overview

Contains the `FITSProcessor` class, which processes directories of FITS files and extracts metadata in parallel.

#### Classes

##### `FITSProcessor`

Processes FITS files and extracts metadata.

**Methods:**

- `__init__(self, max_workers=10)`: Initializes the processor with a maximum number of worker threads.

  **Parameters:**

  - `max_workers (int)`: Maximum number of worker threads.

- `process_fits_file(self, fits_file)`: Processes a single FITS file and extracts its metadata.

  **Parameters:**

  - `fits_file (str)`: Path to the FITS file.

  **Returns:**

  - `dict`: Extracted metadata.

- `process_fits_directory_parallel(self, directory_path)`: Processes all FITS files in a directory using parallel processing.

  **Parameters:**

  - `directory_path (str)`: Path to the directory containing FITS files.

  **Returns:**

  - `pandas.DataFrame`: DataFrame containing metadata for all FITS files.

**Usage:**

```python
from fits_metadata_extractor.processor import FITSProcessor

processor = FITSProcessor(max_workers=10)
metadata_df = processor.process_fits_directory_parallel('/path/to/fits/files')
```

---

### resolver.py

#### Overview

Contains the `ObjectNameResolver` class, which resolves astronomical object names using various services.

#### Classes

##### `ObjectNameResolver`

Resolves astronomical object names using Simbad, NED, VizieR, and Sesame.

**Methods:**

- `__init__(self)`: Initializes the resolver with custom configurations for Simbad and Vizier.

- `standardize_object_name(self, object_name)`: Standardizes the object name to conform to expected formats.

  **Parameters:**

  - `object_name (str)`: The original object name.

  **Returns:**

  - `list`: A list of potential standardized object names.

- `resolve_with_sesame(self, object_name, retries=3)`: Resolves object name using CDS Sesame service.

  **Parameters:**

  - `object_name (str)`: The object name to resolve.
  - `retries (int)`: Number of retry attempts in case of failure.

  **Returns:**

  - `str`: Resolved object name or `'Unknown'`.

- `resolve_object_name(self, object_name)`: Resolves the object name using multiple services and custom mappings.

  **Parameters:**

  - `object_name (str)`: The original object name.

  **Returns:**

  - `tuple`: `(resolved_name (str), method_used (str))`

**Usage:**

Used internally by the metadata extractor to resolve object names.

---

### search.py

#### Overview

Provides functions to search FITS files by sky coordinates or regions.

#### Functions

- `is_point_in_convex_polygon(ra, dec, polygon_coords)`

  **Description:**

  Determines if a point (RA, Dec) is inside a convex polygon defined by `polygon_coords`.

  **Parameters:**

  - `ra (float)`: Right Ascension of the point.
  - `dec (float)`: Declination of the point.
  - `polygon_coords (list)`: List of coordinates defining the polygon.

  **Returns:**

  - `bool`: `True` if the point is inside the polygon, `False` otherwise.

- `search_fits_by_point(metadata_df, ra, dec)`

  **Description:**

  Searches for FITS files in the collection whose coverage includes the given point.

  **Parameters:**

  - `metadata_df (pandas.DataFrame)`: DataFrame containing metadata.
  - `ra (float)`: Right Ascension of the point.
  - `dec (float)`: Declination of the point.

  **Returns:**

  - `pandas.DataFrame`: Subset of `metadata_df` with matching FITS files.

- `search_fits_by_region(metadata_df, region)`

  **Description:**

  Searches for FITS files in the collection whose coverage intersects the given region.

  **Parameters:**

  - `metadata_df (pandas.DataFrame)`: DataFrame containing metadata.
  - `region (dict)`: Dictionary defining the region (e.g., circle or polygon).

  **Returns:**

  - `pandas.DataFrame`: Subset of `metadata_df` with matching FITS files.

**Usage:**

```python
from fits_metadata_extractor.search import search_fits_by_point
from fits_metadata_extractor.utils import load_metadata_from_csv

metadata_df = load_metadata_from_csv('metadata_dataset.csv')
ra = 150.0
dec = 2.0
matching_fits = search_fits_by_point(metadata_df, ra, dec)
```

---

### utils.py

#### Overview

Provides utility functions for data homogenization and saving/loading metadata.

#### Functions

- `serialize_polygon(x)`

  **Description:**

  Serializes the `Polygon_Coords` to a JSON string if it's not `None`. Converts NumPy arrays to lists before serialization.

  **Parameters:**

  - `x (list, tuple, np.ndarray, or None)`: The polygon coordinates.

  **Returns:**

  - `str or None`: JSON string of coordinates or `None`.

- `homogenize_metadata(df)`

  **Description:**

  Homogenizes units and cleans the metadata DataFrame.

  **Parameters:**

  - `df (pandas.DataFrame)`: DataFrame with extracted metadata.

  **Returns:**

  - `pandas.DataFrame`: Homogenized DataFrame.

- `save_metadata_to_csv(df, output_file)`

  **Description:**

  Saves the metadata DataFrame to a CSV file.

  **Parameters:**

  - `df (pandas.DataFrame)`: DataFrame containing metadata.
  - `output_file (str)`: Path to the output CSV file.

- `load_metadata_from_csv(csv_file)`

  **Description:**

  Loads the metadata DataFrame from a CSV file.

  **Parameters:**

  - `csv_file (str)`: Path to the CSV file containing the metadata.

  **Returns:**

  - `pandas.DataFrame`: Loaded metadata DataFrame.

**Usage:**

```python
from fits_metadata_extractor.utils import save_metadata_to_csv, load_metadata_from_csv

save_metadata_to_csv(metadata_df, 'metadata_dataset.csv')
metadata_df = load_metadata_from_csv('metadata_dataset.csv')
```

---

## Examples

### Extracting Metadata

```python
from fits_metadata_extractor.processor import FITSProcessor

processor = FITSProcessor(max_workers=5)
metadata_df = processor.process_fits_directory_parallel('/path/to/fits/files')
```

### Plotting MOC and Polygon

```python
from fits_metadata_extractor.plotter import plot_moc_and_polygon

# Assuming you have polygon coordinates and MOC string
polygon_coords = '[[10.0, 20.0], [15.0, 25.0], [10.0, 30.0], [5.0, 25.0]]'
moc_str = '1/0-4'

plot_moc_and_polygon(polygon_coords, moc_str, title="Example Plot", fits_file='path/to/fits_file.fits')
```

### Searching by Sky Coordinate

```python
from fits_metadata_extractor.search import search_fits_by_point
from fits_metadata_extractor.utils import load_metadata_from_csv

metadata_df = load_metadata_from_csv('metadata_dataset.csv')
ra = 150.0
dec = 2.0
matching_fits = search_fits_by_point(metadata_df, ra, dec)
print(matching_fits)
```

---

## Contributing

Contributions are welcome! Please open an issue or submit a pull request for any improvements or bug fixes.

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Acknowledgments

- [Astropy](https://www.astropy.org/): A community-developed core Python package for Astronomy.
- [Astroquery](https://astroquery.readthedocs.io/en/latest/): A package to access online Astronomical data.
- [MOCpy](https://cds-astro.github.io/mocpy/): A Python library for manipulating MOCs.
- [Shapely](https://shapely.readthedocs.io/en/stable/): A Python package for manipulation and analysis of planar geometric objects.

---
