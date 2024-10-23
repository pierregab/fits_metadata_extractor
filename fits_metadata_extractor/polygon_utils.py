# fits_metadata_extractor/polygon_utils.py

import logging
import numpy as np
from shapely.geometry import Polygon
from regions import PolygonSkyRegion

def polygonskyregion_to_shapely(polygon_region):
    """
    Manually converts a PolygonSkyRegion to a Shapely Polygon.

    Parameters:
        polygon_region (PolygonSkyRegion): The PolygonSkyRegion object.

    Returns:
        shapely.geometry.Polygon: The equivalent Shapely Polygon or None if conversion fails.
    """
    try:
        # Extract longitude and latitude from SkyCoord vertices
        lon = polygon_region.vertices.data.lon.deg
        lat = polygon_region.vertices.data.lat.deg
        coords = list(zip(lon, lat))
        shapely_polygon = Polygon(coords)
        return shapely_polygon
    except Exception as e:
        logging.error(f"Failed to convert PolygonSkyRegion to Shapely Polygon: {e}")
        return None
