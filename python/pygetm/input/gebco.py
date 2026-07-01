import shutil
import tempfile
import os
import zipfile
import logging
import requests
import tqdm
from typing import Optional, Union
from pathlib import Path

URL = "https://dap.ceda.ac.uk/thredds/dodsC/bodc/gebco/global/gebco_2025/ice_surface_elevation/netcdf/GEBCO_2025.nc"
URL = "https://dap.ceda.ac.uk/bodc/gebco/global/gebco_2025/ice_surface_elevation/netcdf/gebco_2025.zip?download=1"
URL = "https://dap.ceda.ac.uk/bodc/gebco/global/gebco_2026/sub_ice_topography_bathymetry/netcdf/GEBCO_2026_sub_ice.zip?download=1"


def download(
    outfile: Union[str, os.PathLike, None] = None,
    logger: Optional[logging.Logger] = None,
    source: str = URL,
) -> Path:
    logging.basicConfig(level=logging.INFO)
    if logger is None:
        logger = logging.getLogger()

    with tempfile.TemporaryFile(prefix="gebco_", suffix=".zip") as fout:
        logger.info(f"Downloading {source}...")
        with requests.get(source, stream=True) as response:
            response.raise_for_status()
            total = int(response.headers.get("content-length", 0)) or None
            with tqdm.tqdm(
                total=total,
                unit="B",
                unit_scale=True,
                unit_divisor=1024,
                desc="Downloading",
            ) as bar:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    fout.write(chunk)
                    bar.update(len(chunk))
        fout.flush()
        fout.seek(0)

        with zipfile.ZipFile(fout) as zipf:
            for name in zipf.namelist():
                if name.endswith(".nc"):
                    break
            else:
                raise ValueError("No .nc file found in the zip archive")
            if outfile is None:
                outfile = name
            outfile = Path(outfile)
            info = zipf.getinfo(name)
            logger.info(f"Extracting {name} to {outfile}...")
            with zipf.open(name) as ncfile, open(outfile, "wb") as outnc:
                with tqdm.tqdm(
                    total=info.file_size,
                    unit="B",
                    unit_scale=True,
                    unit_divisor=1024,
                    desc="Extracting",
                ) as bar:
                    while True:
                        chunk = ncfile.read(1024 * 1024)
                        if not chunk:
                            break
                        outnc.write(chunk)
                        bar.update(len(chunk))
    return outfile


if __name__ == "__main__":
    download()
