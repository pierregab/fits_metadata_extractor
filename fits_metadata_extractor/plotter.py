# fits_metadata_extractor/plotter.py

import os
import logging
import json
import matplotlib.pyplot as plt
from tqdm import tqdm
from .polygon_utils import polygonskyregion_to_shapely
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import warnings
from astropy.utils.exceptions import AstropyWarning
from mocpy import MOC
from astropy.coordinates import SkyCoord
import astropy.units as u
from reproject import reproject_interp
from astropy.wcs import FITSFixedWarning
from regions import PolygonSkyRegion
import pandas as pd




def plot_moc_and_polygon(polygon_coords, moc_str, title="MOC_Debug", fits_file=None):
    """
    Plots the MOC, Polygon, and optionally FITS data on the same WCS projection,
    handling both Equatorial (RA/Dec) and Galactic (GLON/GLAT) coordinate systems.

    Parameters:
        polygon_coords (list of lists/tuples or str): 
            - List of [RA, Dec] or [GLON, GLAT] pairs defining the polygon.
            - Can also be a JSON string representing the list of coordinates.
        moc_str (str): 
            - Serialized MOC string.
        title (str): 
            - Title for the plot.
        fits_file (str, optional): 
            - Path to the FITS file to be plotted.
            - If provided, the FITS data will be reprojected and overlaid on the plot.

    Returns:
        None
    """
    # Suppress warnings for FITS header parsing
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', FITSFixedWarning)
        warnings.simplefilter('ignore', AstropyWarning)
        try:
            # Determine the coordinate frame based on the FITS file's header
            if fits_file:
                with fits.open(fits_file) as hdul:
                    fits_header = hdul[0].header
                    fits_wcs = WCS(fits_header)
                    ctype1 = fits_wcs.wcs.ctype[0].upper()
                    ctype2 = fits_wcs.wcs.ctype[1].upper()

                    if 'RA' in ctype1 and 'DEC' in ctype2:
                        coord_frame = 'icrs'  # Equatorial Coordinates
                    elif 'GLON' in ctype1 and 'GLAT' in ctype2:
                        coord_frame = 'galactic'  # Galactic Coordinates
                    else:
                        logging.error(f"Unsupported coordinate system in {fits_file}: {ctype1}, {ctype2}")
                        return
            else:
                # Default to Equatorial if no FITS file is provided
                coord_frame = 'icrs'

            # Process polygon coordinates
            if isinstance(polygon_coords, str):
                try:
                    polygon_coords = json.loads(polygon_coords)
                    logging.debug("Polygon coordinates deserialized from JSON string.")
                except json.JSONDecodeError as e:
                    logging.error(f"Failed to decode Polygon_Coords: {e}")
                    return

            # Handle flat list of numbers by reshaping into list of [RA, Dec] or [GLON, GLAT] pairs
            if isinstance(polygon_coords, list):
                if all(isinstance(item, (int, float)) for item in polygon_coords):
                    # Flat list detected, reshape into list of pairs
                    if len(polygon_coords) % 2 != 0:
                        logging.error(f"Polygon_Coords list length is not even: {polygon_coords}")
                        return
                    else:
                        polygon_coords = [polygon_coords[i:i+2] for i in range(0, len(polygon_coords), 2)]
                        logging.debug("Polygon coordinates reshaped from flat list to [RA/GLON, Dec/GLAT] pairs.")
                elif not all(isinstance(item, (list, tuple)) and len(item) == 2 for item in polygon_coords):
                    logging.error(f"Invalid Polygon_Coords format: {polygon_coords}")
                    return
            else:
                logging.error(f"Polygon_Coords is not a list or JSON string: {polygon_coords}")
                return

            # Extract RA/GLON and Dec/GLAT from polygon coordinates
            try:
                x_coord, y_coord = zip(*polygon_coords)
                x_coord = np.array(x_coord, dtype=np.float64)
                y_coord = np.array(y_coord, dtype=np.float64)
                logging.debug(f"Extracted X and Y coordinates from polygon: {len(x_coord)} points.")
            except ValueError as ve:
                logging.error(f"Invalid polygon coordinates format: {ve}")
                return

            # Check for NaN or infinite values
            if not (np.isfinite(x_coord).all() and np.isfinite(y_coord).all()):
                logging.error("RA/GLON or Dec/GLAT contains non-finite values.")
                return

            # Close the polygon by adding the first point at the end if not already closed
            if (x_coord[0], y_coord[0]) != (x_coord[-1], y_coord[-1]):
                x_coord = np.append(x_coord, x_coord[0])
                y_coord = np.append(y_coord, y_coord[0])
                logging.debug("Closed the polygon by appending the first coordinate to the end.")

            # Create SkyCoord object for the polygon
            if coord_frame == 'icrs':
                sky_coords = SkyCoord(ra=x_coord * u.deg, dec=y_coord * u.deg, frame='icrs')
            elif coord_frame == 'galactic':
                sky_coords = SkyCoord(l=x_coord * u.deg, b=y_coord * u.deg, frame='galactic')
            else:
                logging.error("Unknown coordinate frame. Cannot create SkyCoord object.")
                return

            # Create a PolygonSkyRegion
            polygon_region = PolygonSkyRegion(vertices=sky_coords)

            # Convert to Shapely polygon manually
            shapely_polygon = polygonskyregion_to_shapely(polygon_region)
            if shapely_polygon:
                polygon_wkt = shapely_polygon.wkt
                polygon_coords_list = list(shapely_polygon.exterior.coords)
                logging.debug("Polygon WKT and coordinates extracted successfully.")
            else:
                logging.error("Failed to convert PolygonSkyRegion to Shapely Polygon.")
                return

            # Deserialize the MOC string
            try:
                moc = MOC.from_string(moc_str)
                logging.debug("Successfully deserialized MOC.")
            except Exception as e:
                logging.error(f"Failed to deserialize MOC: {e}")
                return

            # Initialize FITS data variables
            fits_data_reprojected = None
            im = None

            # If FITS file is provided, load the data
            if fits_file:
                try:
                    with fits.open(fits_file) as hdul:
                        fits_data = hdul[0].data
                        fits_header = hdul[0].header
                        fits_wcs = WCS(fits_header)
                    logging.debug("FITS data loaded successfully.")
                except Exception as e:
                    logging.error(f"Failed to load FITS file '{fits_file}': {e}")
                    fits_file = None  # Proceed without FITS data

            # Create a matplotlib figure
            fig = plt.figure(figsize=(10, 5))

            try:
                # Obtain WCS from MOC by passing the figure
                wcs = moc.wcs(fig)

                # Use the WCS for the projection
                ax = fig.add_subplot(111, projection=wcs)
                ax.grid(True)

                # Determine shape_out based on figure size and DPI
                dpi = fig.get_dpi()
                width_px = int(fig.get_figwidth() * dpi)
                height_px = int(fig.get_figheight() * dpi)
                shape_out = (height_px, width_px)
                logging.debug(f"Determined shape_out for reproject: {shape_out}")

                # If FITS data is available, reproject it to match MOC's WCS
                if fits_file:
                    try:
                        fits_data_reprojected, footprint = reproject_interp(
                            (fits_data, fits_wcs), 
                            wcs, 
                            shape_out=shape_out
                        )
                        logging.debug("FITS data reprojected successfully to match MOC's WCS.")

                        # Plot the reprojected FITS data using imshow
                        im = ax.imshow(
                            fits_data_reprojected, 
                            origin='lower', 
                            cmap='gray', 
                            aspect='auto', 
                            interpolation='none', 
                            alpha=0.7, 
                            label='FITS Data'
                        )
                        logging.debug("FITS data plotted successfully.")
                    except Exception as e:
                        logging.error(f"Failed to reproject or plot FITS data: {e}")

                # Plot the MOC
                moc.fill(ax=ax, wcs=wcs, alpha=0.3, color="blue", label="MOC")
                moc.border(ax=ax, wcs=wcs, color="red", linewidth=1.5)
                logging.debug("MOC filled and bordered successfully.")

                # Plot the polygon
                if coord_frame == 'icrs':
                    transform = ax.get_transform('icrs')
                    xlabel = 'Right Ascension (J2000)'
                    ylabel = 'Declination (J2000)'
                elif coord_frame == 'galactic':
                    transform = ax.get_transform('galactic')
                    xlabel = 'Galactic Longitude (GLON)'
                    ylabel = 'Galactic Latitude (GLAT)'
                else:
                    logging.error("Unknown coordinate frame. Cannot plot polygon.")
                    transform = None
                    xlabel = 'X Coordinate'
                    ylabel = 'Y Coordinate'

                if transform:
                    ax.plot(
                        x_coord, 
                        y_coord, 
                        transform=transform, 
                        color='green', 
                        linewidth=2, 
                        label="Polygon"
                    )
                    logging.debug("Polygon plotted successfully.")

                    # Set axis labels based on coordinate frame
                    ax.set_xlabel(xlabel)
                    ax.set_ylabel(ylabel)
                    logging.debug(f"Axis labels set to '{xlabel}' and '{ylabel}'.")

                # Add legend and title
                ax.legend(loc='upper right')
                plt.title(title)

                # Add a colorbar for FITS data if plotted
                if fits_file and im is not None:
                    try:
                        cbar = plt.colorbar(im, ax=ax, orientation='vertical', pad=0.02)
                        cbar.set_label('FITS Data Intensity')
                        logging.debug("Colorbar added successfully.")
                    except Exception as e:
                        logging.error(f"Failed to add colorbar: {e}")

                # Save the plot as a PNG file
                try:
                    filename = f"{title.replace(' ', '_')}.png"
                    plt.savefig(filename, dpi=300, bbox_inches='tight')
                    logging.info(f"Plot saved as {filename}")
                except Exception as e:
                    logging.error(f"Error saving plot: {e}")

                # Optionally display the plot
                # plt.show()

            except Exception as e:
                logging.error(f"Error plotting MOC, polygon, and FITS data: {e}")
                plt.close(fig)
                return

            finally:
                # Close the figure to free memory
                plt.close(fig)

        except Exception as e:
            logging.error(f"Error plotting MOC, polygon, and FITS data: {e}")
            plt.close(fig)
            return
        
def plot_moc_and_polygon_from_dataset_notebook(
    metadata_df, 
    input_dir, 
    output_dir='plots', 
    max_plots=None,
    fits_files=None,
    indices=None,
    filter_func=None
):
    """
    Facilitates the plotting of MOC and Polygon regions for selected FITS files in the metadata dataset.
    
    This function iterates over each selected entry in the provided metadata DataFrame, retrieves the corresponding
    FITS file, and generates a plot of its MOC and Polygon region. The plots are saved as PNG files in
    the specified output directory.
    
    Parameters:
        metadata_df (pandas.DataFrame): DataFrame containing metadata for FITS files.
            Expected columns:
                - 'FITS_File': Basename of the FITS file.
                - 'Polygon_Coords': JSON string of polygon coordinates.
                - 'MOC': Serialized MOC string.
        input_dir (str): Directory path where the FITS files are located.
        output_dir (str, optional): Directory to save the generated plots. Defaults to 'plots'.
        max_plots (int, optional): Maximum number of plots to generate. Useful for testing.
            If None, plots all selected FITS files.
        fits_files (list of str, optional): List of FITS filenames to plot.
            If provided, only these FITS files will be plotted.
        indices (list of int, optional): List of DataFrame indices to plot.
            If provided, only these rows will be plotted.
        filter_func (callable, optional): A function that takes a DataFrame row and returns True
            if the row should be plotted, False otherwise.
            Example: lambda row: row['Resolved_Object'] == 'SN2023abc'

    Returns:
        None
    """
    import os
    import json
    import logging
    from tqdm.notebook import tqdm  # For progress bar in Jupyter Notebook
    
    # Ensure the plotting function is accessible
    if 'plot_moc_and_polygon' not in globals():
        raise ImportError("The function 'plot_moc_and_polygon' must be defined or imported before using this plotting function.")
    
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    logging.info(f"Plots will be saved to: {output_dir}")
    
    # Select rows based on provided parameters
    if fits_files:
        selected_df = metadata_df[metadata_df['FITS_File'].isin(fits_files)]
        if selected_df.empty:
            logging.warning("No FITS files matched the provided filenames.")
    elif indices:
        selected_df = metadata_df.iloc[indices]
        if selected_df.empty:
            logging.warning("No rows matched the provided indices.")
    elif filter_func:
        selected_df = metadata_df[metadata_df.apply(filter_func, axis=1)]
        if selected_df.empty:
            logging.warning("No rows matched the filter criteria.")
    else:
        # If no selection criteria provided, select all
        selected_df = metadata_df
        logging.info("No specific selection criteria provided. All FITS files will be plotted.")
    
    # Determine the number of plots to generate
    total_plots = len(selected_df) if max_plots is None else min(len(selected_df), max_plots)
    logging.info(f"Generating {total_plots} plots.")
    
    # Iterate over the selected DataFrame rows with a progress bar
    for idx, row in tqdm(selected_df.head(total_plots).iterrows(), total=total_plots, desc="Generating Plots"):
        fits_file = row.get('FITS_File', None)
        polygon_coords_str = row.get('Polygon_Coords', None)
        moc_str = row.get('MOC', None)
        
        if pd.isnull(fits_file):
            logging.warning(f"Row {idx}: 'FITS_File' is missing. Skipping.")
            continue
        
        if pd.isnull(polygon_coords_str):
            logging.warning(f"Row {idx}: 'Polygon_Coords' is missing for '{fits_file}'. Skipping.")
            continue
        
        if pd.isnull(moc_str):
            logging.warning(f"Row {idx}: 'MOC' is missing for '{fits_file}'. Skipping.")
            continue
        
        # Construct the full path to the FITS file
        fits_path = os.path.join(input_dir, fits_file)
        if not os.path.isfile(fits_path):
            logging.error(f"FITS file not found: {fits_path}. Skipping.")
            continue
        
        # Define a title for the plot
        title = os.path.splitext(fits_file)[0]
        
        try:
            # Generate the plot using the existing plot_moc_and_polygon function
            plot_moc_and_polygon(
                polygon_coords=polygon_coords_str,
                moc_str=moc_str,
                title=title,
                fits_file=fits_path
            )
            
            # Move the generated plot to the output directory
            plot_filename = f"{title.replace(' ', '_')}.png"
            if os.path.exists(plot_filename):
                destination = os.path.join(output_dir, plot_filename)
                os.rename(plot_filename, destination)
                logging.info(f"Plot saved: {destination}")
            else:
                logging.error(f"Plot file was not created for '{fits_file}'.")
        
        except Exception as e:
            logging.error(f"Failed to plot '{fits_file}': {e}")
            continue