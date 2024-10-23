# fits_metadata_extractor/resolver.py

import logging
import re
import requests
import xml.etree.ElementTree as ET
from urllib.parse import quote
from functools import lru_cache
from astroquery.simbad import Simbad
from astroquery.ipac.ned import Ned
from astroquery.vizier import Vizier
from astroquery.exceptions import InvalidQueryError, TimeoutError

from .custom_mapping import CUSTOM_NAME_MAPPING

class ObjectNameResolver:
    """
    Class for resolving astronomical object names using various services.
    """

    def __init__(self):
        # Configure custom Simbad and Vizier
        self.custom_simbad = Simbad()
        self.custom_simbad.add_votable_fields('otype', 'ra(d)', 'dec(d)')

        self.custom_vizier = Vizier(columns=['*'], row_limit=1)

    def standardize_object_name(self, object_name):
        """
        Standardizes the astronomical object name to conform to expected formats.

        Parameters:
            object_name (str): The original object name.

        Returns:
            list: A list of potential standardized object names.
        """
        if not object_name or object_name.strip().lower() == 'unknown':
            return ['Unknown']
        
        object_name = object_name.strip()

        # Create a list of potential standardized names
        potential_names = [object_name]

        # Attempt prefix correction: 'SNR' to 'SN' if applicable
        if object_name.upper().startswith('SNR'):
            sn_name = 'SN' + object_name[3:]
            potential_names.append(sn_name)

        # Remove spaces and special characters for alternative attempts
        sanitized_name = re.sub(r'\s+', '', object_name)
        sanitized_name = re.sub(r'[^\w\-]', '', sanitized_name)
        potential_names.append(sanitized_name)

        # Add more standardized forms if needed
        # Example: Append common prefixes/suffixes
        if not object_name.startswith(('SN', 'SNR')):
            potential_names.append('SN' + object_name)
            potential_names.append(object_name + 'SN')

        return potential_names

    def resolve_with_sesame(self, object_name, retries=3):
        """
        Resolves object name using CDS Sesame service with XML output.

        Parameters:
            object_name (str): The object name to resolve.
            retries (int): Number of retry attempts in case of failure.

        Returns:
            str: Resolved object name or 'Unknown'.
        """
        if not object_name or object_name.strip().lower() == 'unknown':
            logging.warning("Object name is unknown or empty.")
            return 'Unknown'

        url = f"http://cdsweb.u-strasbg.fr/cgi-bin/nph-sesame/-oxp/{quote(object_name)}"

        attempt = 0
        while attempt < retries:
            try:
                response = requests.get(url, timeout=10)
                response.raise_for_status()

                # Parse XML response
                root = ET.fromstring(response.content)
                obj = root.find('.//obj')
                if obj is not None and 'name' in obj.attrib:
                    resolved_name = obj.attrib['name']
                    logging.info(f"Resolved '{object_name}' to '{resolved_name}' using Sesame XML.")
                    return resolved_name

                logging.warning(f"No resolved name found in Sesame XML for '{object_name}'.")
                return 'Unknown'

            except ET.ParseError as e:
                logging.error(f"XML parsing failed for '{object_name}': {e}")
                attempt += 1
                logging.info(f"Retrying ({attempt}/{retries})...")
            except requests.RequestException as e:
                logging.error(f"Sesame request failed for '{object_name}': {e}")
                attempt += 1
                logging.info(f"Retrying ({attempt}/{retries})...")
            except Exception as e:
                logging.error(f"Unexpected error during Sesame resolution for '{object_name}': {e}")
                return 'Unknown'

        logging.error(f"Failed to resolve '{object_name}' after {retries} attempts.")
        return 'Unknown'

    @lru_cache(maxsize=1000)
    def resolve_object_name(self, object_name):
        """
        Resolves the celestial object name using multiple astronomical databases and custom mappings.

        Parameters:
            object_name (str): The original object name.

        Returns:
            tuple: (resolved_name (str), method_used (str))
        """
        if not object_name or object_name.strip().lower() == 'unknown':
            logging.warning("Object name is unknown or empty.")
            return ('Unknown', 'Unknown')

        # Check if object_name exists in custom mappings
        if object_name in CUSTOM_NAME_MAPPING:
            resolved_name = CUSTOM_NAME_MAPPING[object_name]
            logging.info(f"Resolved '{object_name}' to '{resolved_name}' using custom mapping.")
            return (resolved_name, 'Custom Mapping')

        potential_names = self.standardize_object_name(object_name)

        for name in potential_names:
            if name == 'Unknown':
                continue

            # Attempt resolution with Simbad
            try:
                result_simbad = self.custom_simbad.query_object(name)
                if result_simbad:
                    resolved_name = result_simbad['MAIN_ID'][0]
                    resolved_name = resolved_name.decode('utf-8') if isinstance(resolved_name, bytes) else resolved_name
                    logging.info(f"Resolved '{object_name}' to '{resolved_name}' using Simbad with name '{name}'.")
                    return (resolved_name, 'Simbad')
            except (InvalidQueryError, TimeoutError) as e:
                logging.error(f"Simbad query failed for '{name}': {e}")
            except Exception as e:
                logging.error(f"Unexpected error during Simbad query for '{name}': {e}")

            # Attempt resolution with NED
            try:
                result_ned = Ned.query_object(name)
                if result_ned:
                    resolved_name = result_ned['Object Name'][0]
                    logging.info(f"Resolved '{object_name}' to '{resolved_name}' using NED with name '{name}'.")
                    return (resolved_name, 'NED')
            except (InvalidQueryError, TimeoutError) as e:
                logging.error(f"NED query failed for '{name}': {e}")
            except Exception as e:
                logging.error(f"Unexpected error during NED query for '{name}': {e}")

            # Attempt resolution with VizieR
            try:
                result_vizier = self.custom_vizier.query_object(name)
                if result_vizier:
                    if 'Identifier' in result_vizier[0].columns:
                        resolved_name = result_vizier[0]['Identifier'][0]
                        resolved_name = resolved_name.decode('utf-8') if isinstance(resolved_name, bytes) else resolved_name
                    else:
                        resolved_name = name  # Fallback to standardized name
                    logging.info(f"Resolved '{object_name}' to '{resolved_name}' using VizieR with name '{name}'.")
                    return (resolved_name, 'VizieR')
            except (InvalidQueryError, TimeoutError) as e:
                logging.error(f"VizieR query failed for '{name}': {e}")
            except Exception as e:
                logging.error(f"Unexpected error during VizieR query for '{name}': {e}")

        # If all resolvers fail, attempt with Sesame as a last resort
        try:
            resolved_name = self.resolve_with_sesame(object_name)
            if resolved_name != 'Unknown':
                return (resolved_name, 'Sesame')
        except Exception as e:
            logging.error(f"Sesame fallback failed for '{object_name}': {e}")

        # If all resolvers fail
        logging.warning(f"No resolved name found for '{object_name}'.")
        return ('Unknown', 'Unknown')
