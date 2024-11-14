# fits_metadata_extractor/plotter.py

# Import statements at the top of plotter.py
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
from astropy.wcs import FITSFixedWarning
from regions import PolygonSkyRegion
import pandas as pd
from matplotlib.patches import Polygon as MplPolygon
from shapely.geometry import Polygon
from shapely.ops import unary_union
from reproject import reproject_interp
from astropy.coordinates import SkyCoord
import astropy.units as u
from mocpy import MOC
from .search import search_fits_by_point, search_fits_by_region





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

def test():
    return "Test successful."

def plot_search_region_and_find_fits(
    metadata_df, 
    region, 
    input_dir, 
    output_dir='search_plots', 
    max_plots=None, 
    plot_search_region=True
):
    """
    Plots the search region and overlays the coverage of matching FITS files using WCS transformations.

    Parameters:
        metadata_df (pandas.DataFrame):
            DataFrame containing metadata for FITS files.
            Expected columns:
                - 'FITS_File': Basename of the FITS file.
                - 'Polygon_Coords': JSON string of polygon coordinates.
                - 'MOC': Serialized MOC string.
                - 'Coordinate_Frame': Coordinate frame of the polygon ('icrs' or 'galactic').

        region (dict):
            Dictionary defining the search region. Supported types:
                - Point:
                    {
                        'type': 'point',
                        'coordinates': (ra, dec),  # In degrees, ICRS frame
                        'coordinate_frame': 'icrs'  # or 'galactic'
                    }
                - Circle:
                    {
                        'type': 'circle',
                        'center': (ra, dec),  # In degrees, ICRS frame
                        'radius': radius_deg    # Radius in degrees
                    }
                - Polygon:
                    {
                        'type': 'polygon',
                        'coordinates': [(ra1, dec1), (ra2, dec2), ...],  # In degrees, ICRS frame
                    }

        input_dir (str):
            Directory path where the FITS files are located.

        output_dir (str, optional):
            Directory to save the generated plot. Defaults to 'search_plots'.

        max_plots (int, optional):
            Maximum number of FITS file plots to generate. Useful for testing.
            If None, plots all matching FITS files.

        plot_search_region (bool, optional):
            If True, plots the search region on the map. Defaults to True.

    Returns:
        None
    """
    import os
    import logging
    import matplotlib.pyplot as plt
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from tqdm.notebook import tqdm
    from mocpy import MOC
    import numpy as np
    import json
    import pandas as pd
    from astropy.io import fits
    from astropy.wcs import WCS
    from astropy.utils.exceptions import AstropyWarning
    import warnings
    from astropy.visualization.wcsaxes import WCSAxes

    # Ensure necessary functions are imported
    if 'search_fits_by_point' not in globals() and 'search_fits_by_region' not in globals():
        raise ImportError("Search functions from 'search.py' must be imported before using this plotting function.")

    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    logging.info(f"Plots will be saved to: {output_dir}")

    # Determine the type of search and find matching FITS files
    if region['type'] == 'point':
        ra, dec = region['coordinates']
        coordinate_frame = region.get('coordinate_frame', 'icrs').lower()
        if coordinate_frame == 'icrs':
            matching_df = search_fits_by_point(metadata_df, ra, dec)
        elif coordinate_frame == 'galactic':
            # Convert point to ICRS for searching
            sky_coord = SkyCoord(l=ra * u.deg, b=dec * u.deg, frame='galactic').icrs
            matching_df = search_fits_by_point(metadata_df, sky_coord.ra.deg, sky_coord.dec.deg)
        else:
            logging.error(f"Unsupported coordinate frame: {coordinate_frame}")
            return
        plot_title = f"Search Region: Point (RA={ra}, Dec={dec})"
    elif region['type'] == 'circle':
        center_ra, center_dec = region['center']
        radius = region['radius']
        search_region = {
            'type': 'circle',
            'center': (center_ra, center_dec),
            'radius': radius
        }
        matching_df = search_fits_by_region(metadata_df, search_region)
        plot_title = f"Search Region: Circle (RA={center_ra}, Dec={center_dec}, Radius={radius}°)"
    elif region['type'] == 'polygon':
        coordinates = region['coordinates']
        search_region = {
            'type': 'polygon',
            'coordinates': coordinates
        }
        matching_df = search_fits_by_region(metadata_df, search_region)
        plot_title = f"Search Region: Polygon ({len(coordinates)} vertices)"
    else:
        logging.error(f"Unsupported region type: {region['type']}")
        return

    if matching_df.empty:
        logging.warning("No FITS files matched the search criteria.")
        return

    # Limit the number of plots if max_plots is set
    total_plots = len(matching_df) if max_plots is None else min(len(matching_df), max_plots)
    logging.info(f"Generating plots for {total_plots} matching FITS files.")

    # Initialize a matplotlib figure with WCS projection (Mollweide)
    # Create a WCS for the Mollweide projection
    main_wcs = WCS(naxis=2)
    main_wcs.wcs.crval = [0, 0]  # Reference coordinates (RA, Dec)
    main_wcs.wcs.crpix = [180, 90]  # Reference pixel (center of the plot)
    main_wcs.wcs.ctype = ["RA---MOL", "DEC--MOL"]
    main_wcs.wcs.cdelt = [-1, 1]  # Degrees per pixel
    main_wcs.wcs.cunit = ["deg", "deg"]

    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(111, projection=main_wcs)

    ax.grid(True, color='lightgray')
    ax.set_title(plot_title + ' (Mollweide Projection)', pad=20)

    # Set axis labels
    ax.set_xlabel('Right Ascension (degrees)')
    ax.set_ylabel('Declination (degrees)')

    # Set coordinate grid labels and formats
    lon = ax.coords[0]
    lat = ax.coords[1]
    lon.set_major_formatter('hh:mm')
    lat.set_major_formatter('dd')
    lon.set_ticks(spacing=15 * u.deg)
    lat.set_ticks(spacing=15 * u.deg)
    lon.display_minor_ticks(True)
    lat.display_minor_ticks(True)

    # Plot the search region if required
    if plot_search_region:
        if region['type'] == 'point':
            ra, dec = region['coordinates']
            coordinate_frame = region.get('coordinate_frame', 'icrs').lower()
            if coordinate_frame == 'icrs':
                sky_coord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame='icrs')
            elif coordinate_frame == 'galactic':
                sky_coord = SkyCoord(l=ra * u.deg, b=dec * u.deg, frame='galactic').icrs
            ra_deg = sky_coord.ra.wrap_at(180 * u.deg).deg
            dec_deg = sky_coord.dec.deg
            ax.plot(ra_deg, dec_deg, marker='*', color='yellow', markersize=15, label='Search Point', transform=ax.get_transform('world'))
        elif region['type'] == 'circle':
            center_ra, center_dec = region['center']
            radius = region['radius']
            center = SkyCoord(ra=center_ra * u.deg, dec=center_dec * u.deg, frame='icrs')
            # Create circle coordinates
            angles = np.linspace(0, 2 * np.pi, 100)
            ra_circle = center.ra.deg + (radius * np.cos(angles)) / np.cos(center.dec.radian)
            dec_circle = center.dec.deg + radius * np.sin(angles)
            sky_circle = SkyCoord(ra=ra_circle * u.deg, dec=dec_circle * u.deg, frame='icrs')
            ra_deg = sky_circle.ra.wrap_at(180 * u.deg).deg
            dec_deg = sky_circle.dec.deg
            ax.plot(ra_deg, dec_deg, color='yellow', linestyle='--', linewidth=2, label='Search Circle', transform=ax.get_transform('world'))
        elif region['type'] == 'polygon':
            polygon_coords = region['coordinates']
            sky_polygon = SkyCoord(
                ra=[c[0] for c in polygon_coords] * u.deg,
                dec=[c[1] for c in polygon_coords] * u.deg,
                frame='icrs'
            )
            ra_deg = sky_polygon.ra.wrap_at(180 * u.deg).deg
            dec_deg = sky_polygon.dec.deg
            ax.plot(ra_deg, dec_deg, color='yellow', linestyle='-', linewidth=2, label='Search Polygon', transform=ax.get_transform('world'))

    # Iterate over matching FITS files and plot their MOCs and polygons
    for idx, row in tqdm(
        matching_df.head(total_plots).iterrows(), 
        total=total_plots, 
        desc="Plotting FITS Coverages"
    ):
        fits_file = row.get('FITS_File', None)
        polygon_coords_str = row.get('Polygon_Coords', None)
        moc_str = row.get('MOC', None)

        if pd.isnull(fits_file) or pd.isnull(polygon_coords_str) or pd.isnull(moc_str):
            logging.warning(f"Missing data for FITS file in row {idx}. Skipping.")
            continue

        fits_path = os.path.join(input_dir, fits_file)
        if not os.path.isfile(fits_path):
            logging.error(f"FITS file not found: {fits_path}. Skipping.")
            continue

        try:
            # Suppress warnings for FITS header parsing
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', AstropyWarning)
                
                # Determine the coordinate frame based on the FITS file's header
                with fits.open(fits_path) as hdul:
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
                        continue  # Skip this FITS file

                # Deserialize the MOC string
                moc = MOC.from_string(moc_str)

                # If the MOC is in a different coordinate frame, transform it to ICRS
                if coord_frame == 'galactic':
                    moc = moc.to_icrs()

                # Plot the MOC using the main WCS
                try:
                    # Use unique label only once
                    label_moc = 'FITS MOC' if idx == matching_df.index[0] else ""
                    moc.fill(ax=ax, wcs=main_wcs, alpha=0.3, label=label_moc, color='blue')
                    moc.border(ax=ax, wcs=main_wcs, color='blue', linewidth=1.0)
                except Exception as e:
                    logging.error(f"Failed to plot MOC for '{fits_file}': {e}")
                    continue

                # Deserialize and plot the polygon
                try:
                    polygon_coords = json.loads(polygon_coords_str)
                    if isinstance(polygon_coords, list):
                        if all(isinstance(item, (int, float)) for item in polygon_coords):
                            # Flat list detected, reshape into list of pairs
                            if len(polygon_coords) % 2 != 0:
                                logging.error(f"Polygon_Coords list length is not even: {polygon_coords}")
                                continue
                            else:
                                polygon_coords = [polygon_coords[i:i+2] for i in range(0, len(polygon_coords), 2)]
                        elif not all(isinstance(item, (list, tuple)) and len(item) == 2 for item in polygon_coords):
                            logging.error(f"Invalid Polygon_Coords format: {polygon_coords}")
                            continue
                    else:
                        logging.error(f"Polygon_Coords is not a list or JSON string: {polygon_coords}")
                        continue

                    # Extract RA/GLON and Dec/GLAT from polygon coordinates
                    x_coord, y_coord = zip(*polygon_coords)
                    x_coord = np.array(x_coord, dtype=np.float64)
                    y_coord = np.array(y_coord, dtype=np.float64)

                    # Check for NaN or infinite values
                    if not (np.isfinite(x_coord).all() and np.isfinite(y_coord).all()):
                        logging.error("RA/GLON or Dec/GLAT contains non-finite values.")
                        continue

                    # Close the polygon by adding the first point at the end if not already closed
                    if (x_coord[0], y_coord[0]) != (x_coord[-1], y_coord[-1]):
                        x_coord = np.append(x_coord, x_coord[0])
                        y_coord = np.append(y_coord, y_coord[0])

                    # Create SkyCoord object for the polygon
                    if coord_frame == 'icrs':
                        sky_coords = SkyCoord(ra=x_coord * u.deg, dec=y_coord * u.deg, frame='icrs')
                    elif coord_frame == 'galactic':
                        sky_coords = SkyCoord(l=x_coord * u.deg, b=y_coord * u.deg, frame='galactic').icrs
                    else:
                        logging.error("Unknown coordinate frame. Cannot create SkyCoord object.")
                        continue

                    # Convert to degrees for plotting
                    ra_deg = sky_coords.ra.wrap_at(180 * u.deg).deg
                    dec_deg = sky_coords.dec.deg

                    # Use unique label only once
                    label_polygon = 'FITS Polygon' if idx == matching_df.index[0] else ""
                    ax.plot(
                        ra_deg,
                        dec_deg,
                        color='green',
                        linewidth=1,
                        label=label_polygon,
                        transform=ax.get_transform('world')
                    )

                except json.JSONDecodeError as e:
                    logging.error(f"Failed to decode Polygon_Coords for '{fits_file}': {e}")
                    continue
                except Exception as e:
                    logging.error(f"Failed to plot polygon for '{fits_file}': {e}")
                    continue

        except Exception as e:
            logging.error(f"Failed to process FITS file '{fits_file}': {e}")
            continue

    # Finalize the plot
    # To avoid duplicate labels in the legend
    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    ax.legend(unique.values(), unique.keys(), loc='upper right', fontsize='small')

    output_file = os.path.join(output_dir, "search_region_and_fits.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    logging.info(f"Search region and matching FITS coverages plotted and saved to '{output_file}'")
    plt.show()
