from .tiling_functions import read_limits
import requests 
import time
from dnora import msg
import zipfile
from pathlib import Path
import os
def find_tiles(lon, lat):
    tiles, lons, lats = read_limits("kartverket50m_tile_coords.txt")
    needed_tiles = []

    for t, lo, la in zip(tiles, lons, lats):
        lon_match = lo[0]<lon[1] and lo[1]>lon[0]
        lat_match = la[0]<lat[1] and la[1]>lat[0]
        if lon_match and lat_match:
            needed_tiles.append(t)

    return needed_tiles

def get_geonorge_email():
    # Define the path to the configuration file in the user's home directory
    rcfile = "~/.geonorgeapirc"
    config_file = os.path.expanduser(rcfile)

    # Check if the configuration file exists
    if os.path.exists(config_file):
        # Read the email from the file
        with open(config_file, "r") as file:
            email = file.read().strip()
            if email:
                msg.get_value(email, rcfile, 'Reading email')
                return email
    
    # If the file doesn't exist or is empty, prompt the user for their email
    email = input("Give email for Geonorge API: ").strip()

    # Save the email to the configuration file
    while True:
        save_in_rc = input(f"Save email in {rcfile} for future reference? [Y/n]: ")
        if not save_in_rc or save_in_rc.lower() == 'y':
            with open(config_file, "w") as file:
                file.write(email)
            msg.to_file(rcfile)
            return email
        elif save_in_rc.lower() == 'n':
            return email
    

    


def download_tile(tile:str, folder:str, email: str) -> None:
    """Downloads a KArtverket50m tile using the GeoNorge-API"""
    order_url = "https://nedlasting.geonorge.no/api/order"
    metadata_uuid = "bbd687d0-d34f-4d95-9e60-27e330e0f76e"
    
    
    # Define the payload for the order
    order_payload = {
        "downloadAsBundle": False,   # Set to True if you want all data in a single bundle
        "email": email,             # Add your email address
        "orderLines": [
            {
                "areas": [          # Specify the tile code
                    {
                        "code": tile,  # Use the tile code
                        "name": tile,  # Use the tile name
                        "type": "celle"           # Type of area (e.g., 'celle' for tiles)
                    }
                ],
                "formats": [        # Specify the format
                    {
                        "name": "XYZ"
                    }
                ],
                "metadataUuid": metadata_uuid,  # Specify the metadata UUID for the dataset
                "projections": [    # Specify the projection
                    {
                        "code": "25833",
                        "name": "EUREF89 UTM sone 33, 2d",
                        "codespace": "http://www.opengis.net/def/crs/EPSG/0/25833"
                    }
                ]
            }
        ]
    }

    # Submit the order
    response = requests.post(order_url, json=order_payload)
    if response.status_code == 200:
        order_details = response.json()  # Order confirmation with download link
        msg.plain(f"Order submitted: {order_details}")
    else:
        msg.plain(f"Order failed with status code {response.status_code}")
        msg.plain("Response content:", response.text)


    for n in range(10):
        files = order_details.get("files", [])
        if files:
            break
        if n == 0:
            msg.plain('Waiting for order to be ready')
        else:
            msg.plain('.')
        time.sleep(1)
        

    if not files:
        msg.error('No files available in the API order!')
        return

    # Download the file
    file_info = files[0]  
    download_url = file_info.get("downloadUrl")
    file_name = file_info.get("name")
    response = requests.get(download_url)
    filepath = Path(folder)/Path(file_name)
    if response.status_code == 200:
        # Save the file locally
        with open(filepath, "wb") as file:
            file.write(response.content)
    msg.blank()
    with zipfile.ZipFile(filepath, 'r') as zip_ref:
        file_list = zip_ref.namelist()
        msg.to_multifile([Path(folder)/Path(fn) for fn in file_list])
        zip_ref.extractall(folder)
    
    os.remove(filepath)