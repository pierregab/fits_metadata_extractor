# test_fits_metadata_extractor.py

import unittest
from fits_metadata_extractor.extractor import is_point_in_convex_polygon
import matplotlib.pyplot as plt
from mocpy import MOC
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs import WCS
import numpy as np
import logging

import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.units as u
from mocpy import MOC
import logging

def debug_plot_moc_and_polygon(polygon_coords, point_ra, point_dec, title="Polygon and MOC Debug Plot"):
    """
    Plots the polygon, MOC, and the point for debugging purposes.

    Parameters:
        polygon_coords (list of tuples): List of (RA, Dec) tuples representing the polygon vertices.
        point_ra (float): Right Ascension of the point in degrees.
        point_dec (float): Declination of the point in degrees.
        title (str): Title for the plot.
    """
    try:
        # Separate RA and Dec for the polygon
        ra_vals, dec_vals = zip(*polygon_coords)

        # Close the polygon by adding the first point at the end
        ra_vals = list(ra_vals) + [ra_vals[0]]
        dec_vals = list(dec_vals) + [dec_vals[0]]

        # Convert the polygon coordinates to SkyCoord
        sky_coords = SkyCoord(ra=ra_vals, dec=dec_vals, unit='deg', frame='icrs')

        # Create the MOC from the polygon
        moc = MOC.from_polygon_skycoord(sky_coords, max_depth=10)

    except ValueError as e:
        logging.error(f"Invalid polygon coordinates: {e}")
        return
    except Exception as e:
        logging.error(f"Error creating SkyCoord or MOC: {e}")
        return

    # Create a figure and generate the WCS from the MOC
    fig = plt.figure(figsize=(10, 10))
    wcs = moc.wcs(fig)

    # Add a subplot to the figure with the WCS projection
    ax = fig.add_subplot(111, projection=wcs)

    # Plot the MOC using the WCS projection
    try:
        moc.fill(ax=ax, wcs=wcs, alpha=0.5, color="blue", label="MOC")
        moc.border(ax=ax, wcs=wcs, color="red")
    except Exception as e:
        logging.error(f"Error plotting MOC: {e}")
        return

    # Plot the polygon as a green line with markers
    try:
        ax.plot(ra_vals, dec_vals, color="green", marker='o', linestyle='-', label="Polygon", transform=ax.get_transform('world'))
    except Exception as e:
        logging.error(f"Error plotting polygon: {e}")
        return

    # Plot the point (RA, Dec) as a red cross marker
    try:
        ax.plot(point_ra, point_dec, 'rx', markersize=10, label="Point", transform=ax.get_transform('world'))
    except Exception as e:
        logging.error(f"Error plotting point: {e}")
        return

    # Add labels, grid, legend, and title
    ax.set_xlabel('Right Ascension (degrees)')
    ax.set_ylabel('Declination (degrees)')
    ax.set_title(title)
    ax.grid(True)
    ax.legend()

    # Show the plot with title in name
    plt.savefig(f"{title}.png")


class TestIsPointInConvexPolygon(unittest.TestCase):
    def test_point_inside_polygon(self):
        """
        Test a point that is clearly inside the convex polygon.
        """
        polygon = [
            (10.0, 10.0),
            (20.0, 10.0),
            (20.0, 20.0),
            (10.0, 20.0)
        ]
        point_ra, point_dec = 15.0, 15.0
        result = is_point_in_convex_polygon(point_ra, point_dec, polygon)
        debug_plot_moc_and_polygon(polygon, point_ra, point_dec, title="Point Inside Polygon")
        self.assertTrue(result, "Point should be inside the polygon.")

    def test_point_outside_polygon(self):
        """
        Test a point that is clearly outside the convex polygon.
        """
        polygon = [
            (10.0, 10.0),
            (20.0, 10.0),
            (20.0, 20.0),
            (10.0, 20.0)
        ]
        point_ra, point_dec = 25.0, 25.0
        result = is_point_in_convex_polygon(point_ra, point_dec, polygon)
        debug_plot_moc_and_polygon(polygon, point_ra, point_dec, title="Point Outside Polygon")
        self.assertFalse(result, "Point should be outside the polygon.")

    def test_point_on_edge(self):
        """
        Test a point that lies exactly on one edge of the convex polygon.
        """
        polygon = [
            (10.0, 10.0),
            (20.0, 10.0),
            (20.0, 20.0),
            (10.0, 20.0)
        ]
        point_ra, point_dec = 15.0, 10.0
        result = is_point_in_convex_polygon(point_ra, point_dec, polygon)
        debug_plot_moc_and_polygon(polygon, point_ra, point_dec, title="Point On Edge")
        self.assertTrue(result, "Point on the edge should be considered inside the polygon.")

    def test_point_at_vertex(self):
        """
        Test a point that coincides with one of the polygon's vertices.
        """
        polygon = [
            (10.0, 10.0),
            (20.0, 10.0),
            (20.0, 20.0),
            (10.0, 20.0)
        ]
        point_ra, point_dec = 10.0, 10.0
        result = is_point_in_convex_polygon(point_ra, point_dec, polygon)
        debug_plot_moc_and_polygon(polygon, point_ra, point_dec, title="Point At Vertex")
        self.assertTrue(result, "Point at the vertex should be considered inside the polygon.")

    def test_polygon_wrapping_ra(self):
        """
        Test a polygon that wraps around the RA=0/360 boundary.
        """
        polygon = [
            (350.0, 10.0),
            (10.0, 10.0),
            (10.0, 20.0),
            (350.0, 20.0)
        ]
        # Point inside the wrapped polygon
        point_ra_inside, point_dec_inside = 355.0, 15.0
        result_inside = is_point_in_convex_polygon(point_ra_inside, point_dec_inside, polygon)
        debug_plot_moc_and_polygon(polygon, point_ra_inside, point_dec_inside, title="Point Inside Wrapped Polygon")
        self.assertTrue(result_inside, "Point should be inside the RA-wrapped polygon.")

        # Point outside the wrapped polygon
        point_ra_outside, point_dec_outside = 15.0, 15.0
        result_outside = is_point_in_convex_polygon(point_ra_outside, point_dec_outside, polygon)
        debug_plot_moc_and_polygon(polygon, point_ra_outside, point_dec_outside, title="Point Outside Wrapped Polygon")
        self.assertFalse(result_outside, "Point should be outside the RA-wrapped polygon.")

    def test_invalid_polygon(self):
        """
        Test handling of invalid polygon coordinates.
        """
        # Polygon with an invalid coordinate (only RA provided)
        polygon = [
            (10.0, 10.0),
            (20.0, 10.0),
            (20.0,),  # Invalid tuple
            (10.0, 20.0)
        ]
        point_ra, point_dec = 15.0, 15.0
        result = is_point_in_convex_polygon(point_ra, point_dec, polygon)
        self.assertFalse(result, "Invalid polygon coordinates should result in point being considered outside.")

    def test_polygon_with_negative_dec(self):
        """
        Test a polygon that includes negative declinations.
        """
        polygon = [
            (10.0, -10.0),
            (20.0, -10.0),
            (20.0, 0.0),
            (10.0, 0.0)
        ]
        point_ra, point_dec = 15.0, -5.0
        result = is_point_in_convex_polygon(point_ra, point_dec, polygon)
        debug_plot_moc_and_polygon(polygon, point_ra, point_dec, title="Point Inside Polygon with Negative Dec")
        self.assertTrue(result, "Point should be inside the polygon with negative declinations.")

    def test_point_near_boundary(self):
        """
        Test a point that is very close to the polygon boundary.
        """
        polygon = [
            (10.0, 10.0),
            (20.0, 10.0),
            (20.0, 20.0),
            (10.0, 20.0)
        ]
        # Point very close to the edge
        point_ra, point_dec = 19.9999, 15.0
        result = is_point_in_convex_polygon(point_ra, point_dec, polygon)
        debug_plot_moc_and_polygon(polygon, point_ra, point_dec, title="Point Near Boundary")
        self.assertTrue(result, "Point very close to the edge should be considered inside the polygon.")

    def test_large_convex_polygon(self):
        """
        Test a large convex polygon covering a significant area.
        """
        polygon = [
            (0.0, -30.0),
            (0.0, 30.0),
            (180.0, 30.0),
            (180.0, -30.0)
        ]
        # Point inside the large polygon
        point_ra_inside, point_dec_inside = 90.0, 0.0
        result_inside = is_point_in_convex_polygon(point_ra_inside, point_dec_inside, polygon)
        debug_plot_moc_and_polygon(polygon, point_ra_inside, point_dec_inside, title="Point Inside Large Polygon")
        self.assertTrue(result_inside, "Point should be inside the large polygon.")

        # Point outside the large polygon
        point_ra_outside, point_dec_outside = 90.0, 40.0
        result_outside = is_point_in_convex_polygon(point_ra_outside, point_dec_outside, polygon)
        debug_plot_moc_and_polygon(polygon, point_ra_outside, point_dec_outside, title="Point Outside Large Polygon")
        self.assertFalse(result_outside, "Point should be outside the large polygon.")

    def test_small_convex_polygon(self):
        """
        Test a small convex polygon covering a limited area.
        """
        polygon = [
            (0.0, 0.0),
            (0.0, 10.0),
            (10.0, 10.0),
            (10.0, 0.0)
        ]
        # Point inside the small polygon
        point_ra, point_dec = 5.0, 5.0
        result = is_point_in_convex_polygon(point_ra, point_dec, polygon)
        debug_plot_moc_and_polygon(polygon, point_ra, point_dec, title="Point Inside Small Polygon")
        self.assertTrue(result, "Point should be inside the smaller polygon.")

        # Point outside the small polygon
        point_ra, point_dec = 15.0, 5.0
        result = is_point_in_convex_polygon(point_ra, point_dec, polygon)
        debug_plot_moc_and_polygon(polygon, point_ra, point_dec, title="Point Outside Small Polygon")
        self.assertFalse(result, "Point should be outside the smaller polygon.")

if __name__ == '__main__':
    # Configure logging to display debug messages during testing
    logging.basicConfig(level=logging.DEBUG)
    unittest.main()