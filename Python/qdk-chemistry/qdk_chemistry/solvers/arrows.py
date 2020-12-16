import requests
import re
import yaml
import logging

from yaml.scanner import ScannerError
from urllib.parse import urljoin

_log = logging.getLogger(__name__)

# <molecule> theory{qsharp_chem} qsharp_chem_filling{<num_orbitals> <num_electrons_alpha> <num_electrons_beta>} qsharp_chem_nroots{<nroots>}
def save_broombridge(molecule: str, file_path: str, url: str = "https://arrows.emsl.pnnl.gov/api/broombridge/") -> dict:
    _log.info(f"Submitting molecule '{molecule}' to EMSL Arrows at {url}...")
    target = urljoin(url, molecule)
    response = requests.get(target)
    assert response.status_code == requests.status_codes.codes.ALL_GOOD, f"Received invalid response at {target}: {response}, {response.text}"

    if response.text.startswith("<!DOCTYPE html>"):
        _log.error(f"Could not get Broombridge file: molecule {molecule} was submitted to EMSL Arrows for simulation. Please try again later.")
    else:
        try:
            yaml.full_load(response.text)
        except ScannerError:
            _log.error(f"Could not parse Broombridge file. Please check if molecule name or SMILES string '{molecule}' is valid.")
        else:
            _log.debug(f"Saving to file {file_path}")
            with open(file_path, "w") as fp:
                fp.write(response.text)
