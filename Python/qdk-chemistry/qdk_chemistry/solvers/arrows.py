import requests
import re
import yaml
import logging

_log = logging.getLogger(__name__)

HTTPS_ = "https?:\/\/"
URL_DOMAIN = "[-a-zA-Z0-9@:%._\+~#=]{1,256}\.[a-zA-Z0-9()]{1,6}"
URL_DIR = "[-a-zA-Z0-9()@:%_\+.~#?&//=]"
DOWNLOAD_DIR = "nwdatafile_download"
FILE_PATTERN = "microsoft_qsharp_chem.yaml-[0-9,-]+[0-9,:,%]+"
BROOMBRIDGE_DL_PATTERN = f"{HTTPS_}{URL_DOMAIN}{URL_DIR}*{DOWNLOAD_DIR}{URL_DIR}* {FILE_PATTERN}"


def get_broombridge_data(text) -> str:
    download_url = re.findall(BROOMBRIDGE_DL_PATTERN, text)[0]
    if download_url:
        response = requests.get(download_url)
        assert response.status_code == requests.status_codes.codes.ALL_GOOD, f"Received invalid response: {response}, {response.text}"
        return response.text


def save_broombridge(molecule: str, file_path: str, url: str) -> dict:
    myobj = {"smi": molecule + " theory{qsharp_chem}"}

    _log.info(f"Submitting molecule '{molecule}' to EMSL Arrows at {url}...")
    response = requests.post(url, data = myobj)
    assert response.status_code == requests.status_codes.codes.ALL_GOOD, f"Received invalid response: {response}, {response.text}"
    data = get_broombridge_data(response.text)

    _log.debug(f"Saving to file {file_path}")
    with open(file_path, "w") as fp:
        fp.write(data)
