from ast import literal_eval
import contextlib
from cyvcf2 import VCF, Writer
from datetime import datetime, timezone
from hashlib import sha256
import json
from loguru import logger
import numpy as np
import os
import subprocess
from pathlib import Path
from pandas import DataFrame, notna
import warnings
import _pickle as cPickle

# ===========================================================================================
# I/O
# ===========================================================================================

def runcmd(cmd: list, stdout: str = None) -> subprocess.CompletedProcess:

    logger.debug(f"Running command: {cmd}")

    if stdout:

        with open(stdout, mode = "wb") as out:

            return subprocess.run(cmd, check=True, stdout=out)
    
    else:

        return subprocess.run(cmd, check=True, capture_output=True)

def save(obj: DataFrame, path: Path, format: str = "pickle", lookup: dict = None, out: str = os.getcwd()) -> int:
    """ Serialize DataFrame to file of specified format """

    assert (isinstance(obj,DataFrame)), "Input provided is not a Dataframe instance"

    # Map file formats to respective file extensions, write modes, and functions
    FILES = {
            "json": {"ext": "json", "mode": "w", "func": json.dump},
            "pickle": {"ext": "pkl", "mode": "wb", "func": cPickle.dump},
            "vcf": {"ext": "vcf", "mode": "wz", "func": None},
            "vcf.gz": {"ext": "vcf.gz", "mode":"wz", "func": None}
        }
    
    assert (
        format in FILES
    ), "File format is not supported"

    if format in FILES:

        outf = str(Path(out, path.stem))

        if format in ["json", "pickle"]:
            # Open data stream
            with open(
                file=f"{outf}_delta.{FILES[format]['ext']}",
                mode=FILES[format]["mode"],
            ) as f:
                obj = obj.to_dict(orient='list')
                FILES[format]["func"](obj, f)
        # Write a VCF
        else:
            # Open the initial VCF file to get the template
            vcf = VCF(f"{str(path)}.gz")
            # Add match data in INFO column
            vcf.add_info_to_header({'ID': 'match', 'Description': 'overlapping variant', 'Type': 'String', 'Number': '1'})
            # Open the output VCF file for writing, open the data stream
            w = Writer(f"{outf}_delta.{FILES[format]['ext']}",vcf)
            # Iterate over the variants in the VCF file
            for i, v in enumerate(vcf):
                # Hash the variant to allow fast lookup
                hash = sha256(
                    string=f"{(v.CHROM).removeprefix('chr')}:{v.POS}:{v.REF}:{'|'.join(v.ALT)}".encode()
                ).hexdigest()
                # O(1) lookup
                if hash in lookup:
                    # Get the corresponding row from the DataFrame, O(1) lookup by index
                    series = obj.iloc[lookup[hash]]
                    # Add the match information to the INFO dictionary
                    v.INFO["match"] = int((notna(series["Variant.L"])) & (notna(series["Variant.R"])))
                # Write the variant to the output VCF file
                w.write_record(v)
            # Close data streams
            w.close(); vcf.close()

        assert os.path.isfile(
            f"{outf}_delta.{FILES[format]['ext']}"
        ), "File was not created"

        logger.success(f"Results are seralized to {out}")

        return 1
    else:
        raise ValueError(f"Error: The file format {format} is not supported.")

def file_infos(path: str) -> dict:
    """ Get file stats """
    statinfo = os.stat(path, follow_symlinks=True)

    return {"basename": os.path.basename(path),
            "path": os.path.dirname(path),
            "size": round(statinfo.st_size / pow(1024,2),2),
            "mtime": datetime.fromtimestamp(statinfo.st_mtime, tz=timezone.utc)}

@contextlib.contextmanager
def suppress_warnings():
    """Suppress warnings from being printed to the console"""
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="no intervals found for*")
        warnings.filterwarnings("ignore", message="Mean of empty slice.")
        yield

# ===========================================================================================
# SETS
# ===========================================================================================


def intersect(a: set[str], b: set[str]) -> set:
    """ Return the intersection of two sets """
    return a & b


def difference(a: set[str], b: set[str]) -> set:
    """ Return the difference of two sets """
    return a - b

def jaccard_index(shared: int, total: dict) -> float|None:
    """ Calculate the Jaccard index """
    try:
        return shared / total
    except ZeroDivisionError:
        return None

# ===========================================================================================
# Variables
# ===========================================================================================

def convert(a: object) -> object:
    """ Convert variable to appropriate type """
    try:
        # If the variable contains a '/' or '|' character, it is a genotype information, return the variable as is
        # Else return the variable as an evaluated expression
        return a if sum(list(map(lambda x: x in a, ('/','|')))) else literal_eval(a)
    except Exception:
        # If the variable cannot be evaluated, return the variable as is
        return a


# ===========================================================================================
# Variants
# ===========================================================================================

def hamming_distance(a: np.ndarray, b: np.ndarray) -> float:
    """ Calculate the Hamming distance """
    try:
        return 1 - ((np.sum(np.not_equal(a, b))) / (a.size + b.size))
    except ZeroDivisionError:
        return None