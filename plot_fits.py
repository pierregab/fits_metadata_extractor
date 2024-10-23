import argparse
import sys
import warnings

from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt

def get_corner_coordinates(data, wcs):
    """
    Get the pixel and world coordinates of the four corners of the image.

    Parameters:
    - data: 2D numpy array of the image data
    - wcs: WCS object for coordinate transformations

    Returns:
    - corners: List of tuples containing (pixel_x, pixel_y, world_coord)
    """
    ny, nx = data.shape  # Note: FITS images are typically in (y, x) format

    # Define corner pixel coordinates
    pixel_corners = [
        (0, 0),          # Bottom-left
        (nx-1, 0),       # Bottom-right
        (nx-1, ny-1),    # Top-right
        (0, ny-1)        # Top-left
    ]

    # Convert pixel coordinates to world coordinates if WCS is available
    world_corners = []
    for x, y in pixel_corners:
        if wcs is not None:
            try:
                world = wcs.pixel_to_world(x, y)
                world_corners.append(world)
            except Exception as e:
                print(f"Error converting pixel to world coordinates at ({x}, {y}): {e}")
                world_corners.append(None)
        else:
            world_corners.append(None)

    corners = list(zip(pixel_corners, world_corners))
    return corners

def plot_fits_with_corners(fits_path):
    """
    Plot the FITS image and mark its corner coordinates.

    Parameters:
    - fits_path: Path to the FITS file
    """
    # Open the FITS file
    try:
        hdulist = fits.open(fits_path)
    except Exception as e:
        print(f"Error opening FITS file: {e}")
        sys.exit(1)

    # Assume the image data is in the primary HDU
    data = hdulist[0].data
    header = hdulist[0].header
    hdulist.close()

    if data is None:
        print("No image data found in the FITS file.")
        sys.exit(1)

    # If the data has more than 2 dimensions, take the first 2
    if data.ndim > 2:
        data = data[0]

    # Fix deprecated 'RADECSYS' keyword by replacing it with 'RADESYS'
    if 'RADECSYS' in header:
        header['RADESYS'] = header['RADECSYS']
        del header['RADECSYS']
        print("Updated FITS header: Replaced 'RADECSYS' with 'RADESYS'.")

    # Suppress specific WCS warnings related to deprecations
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=fits.verify.VerifyWarning)
        try:
            wcs = WCS(header)
            if not wcs.has_celestial:
                print("WCS does not contain celestial information. Proceeding with pixel coordinates only.")
                wcs = None
            else:
                print("WCS Information:")
                print(wcs)
                print("\nWCS Axes:")
                for i, axis in enumerate(wcs.wcs.ctype):
                    print(f"  Axis {i+1}: {axis}")
        except Exception as e:
            print(f"Error parsing WCS information: {e}")
            wcs = None

    # Get corner coordinates
    corners = get_corner_coordinates(data, wcs)

    # Plot the image
    plt.figure(figsize=(8, 6))
    if wcs is not None:
        ax = plt.subplot(projection=wcs)
        im = ax.imshow(data, origin='lower', cmap='gray')
        ax.set_xlabel('RA')
        ax.set_ylabel('Dec')
    else:
        ax = plt.subplot()
        im = ax.imshow(data, origin='lower', cmap='gray')
        ax.set_xlabel('X Pixel')
        ax.set_ylabel('Y Pixel')


    # Plot and annotate corners
    for idx, ((x, y), world) in enumerate(corners):
        ax.plot(x, y, marker='o', markersize=8, markeredgecolor='red', markerfacecolor='yellow')
        if wcs is not None and world is not None:
            # Check if the world coordinate has RA and Dec
            coord_text = ""
            if hasattr(world, 'ra') and hasattr(world, 'dec'):
                coord_text = (f"Corner {idx+1}:\n"
                              f"RA={world.ra.deg:.4f} deg\n"
                              f"Dec={world.dec.deg:.4f} deg")
            elif hasattr(world, 'spherical') and hasattr(world.spherical, 'lon') and hasattr(world.spherical, 'lat'):
                # Generic spherical coordinates
                coord_text = (f"Corner {idx+1}:\n"
                              f"Lon={world.spherical.lon.deg:.4f} deg\n"
                              f"Lat={world.spherical.lat.deg:.4f} deg")
            else:
                # Attempt to retrieve generic world coordinates
                try:
                    coord = world.to_string('hmsdms', sep=':')
                    coord_text = f"Corner {idx+1}:\nCoord={coord}"
                except Exception:
                    coord_text = f"Corner {idx+1}:\nUnrecognized WCS coordinates"

            if coord_text == "":
                coord_text = f"Corner {idx+1}:\nUnrecognized WCS coordinates"
        else:
            coord_text = f"Corner {idx+1}:\nX={x} px\nY={y} px"

        # Adjust text position slightly to avoid overlapping with the marker
        ax.text(x + 10, y + 10, coord_text, fontsize=9, color='white',
                bbox=dict(facecolor='red', alpha=0.5, boxstyle='round,pad=0.2'))

    plt.title('FITS Image with Corner Coordinates')
    plt.tight_layout()
    plt.show()

def main():
    parser = argparse.ArgumentParser(description='Plot FITS image with corner coordinates.')
    parser.add_argument('fits_image', help='Path to the input FITS image.')
    args = parser.parse_args()

    plot_fits_with_corners(args.fits_image)

if __name__ == '__main__':
    main()
