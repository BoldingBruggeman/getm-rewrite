import shutil
import tempfile
import os
import zipfile
import logging
import urllib.request
from typing import Optional, Union
from pathlib import Path


URL = "https://dap.ceda.ac.uk/thredds/dodsC/bodc/gebco/global/gebco_2025/ice_surface_elevation/netcdf/GEBCO_2025.nc"
URL = "https://dap.ceda.ac.uk/bodc/gebco/global/gebco_2025/ice_surface_elevation/netcdf/gebco_2025.zip?download=1"


def download(
    outfile: Union[str, os.PathLike, None] = None,
    logger: Optional[logging.Logger] = None,
    source: str = URL,
) -> Path:
    logging.basicConfig(level=logging.INFO)
    if logger is None:
        logger = logging.getLogger()

    with tempfile.TemporaryFile(prefix="gebco_", suffix=".zip") as fout:
        logger.info(f"Downloading {source} to {fout.name}...")
        with urllib.request.urlopen(source) as response:
            shutil.copyfileobj(response, fout)
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
            logger.info(f"Extracting {name} from {fout.name} to {outfile}...")
            with zipf.open(name) as ncfile, open(outfile, "wb") as outnc:
                shutil.copyfileobj(ncfile, outnc)
    return outfile


if __name__ == "__main__":
    download()
