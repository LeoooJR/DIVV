import contextlib
from cyvcf2 import VCF, Writer
from datetime import datetime, timezone
import errors
import filetype
import gzip
from hashlib import sha256
import json
from loguru import logger
import numpy as np
import os
import subprocess
from pathlib import Path
from pandas import Series, DataFrame, notna, isna
import warnings
import _pickle as cPickle

# I/O

def runcmd(cmd: list, stdout: str = None) -> subprocess.CompletedProcess:

    logger.debug(f"Running command: {cmd}")

    if stdout:

        with open(stdout, mode = "wb") as out:

            return subprocess.run(cmd, check=True, stdout=out)
    
    else:

        return subprocess.run(cmd, check=True, capture_output=True)

def compressing(file: str) -> str|None:

    logger.debug(f"Compressing file {file}.")

    EXT: str = "gz"

    archive: str = f"{file}.{EXT}"

    CMDS: dict = {"bcftools": ["bcftools", "view", "-O", "z", "-o", archive, file],
                  "bgzip": ["bgzip", "-c", file]}

    outcode: int = 1
    retry: int = 0

    # Call a process to compress the file
    while retry < len(CMDS) and outcode != 0:
        bin: str = list(CMDS.keys())[retry]
        logger.debug(f"Compressing with {bin} binary.")
        try:
            process: subprocess.CompletedProcess = runcmd(CMDS[bin], stdout=archive) if bin == "bgzip" else runcmd(CMDS[bin])
            outcode: int = process.returncode
            logger.success(f"{file} has been successfully compressed.")
        except subprocess.CalledProcessError as e:
            logger.warning(f"Compressing {file} with {bin} binary did not succeed with exit code {outcode}.")
            logger.warning(f"Standard output: {e.stdout}")
            logger.warning(f"Standard error: {e.stderr}")
            retry += 1

    # Return archive path ONLY if exit code is 0
    return None if outcode else archive

def indexing(file: str) -> str|None:

    logger.debug(f"Indexing file {file}.")

    EXT: str = "tbi"

    index: str = f"{file}.{EXT}"

    CMDS: dict = {"bcftools": ["bcftools", "index", "-t", file],
                  "tabix": ["tabix", "-p", "vcf", file]}

    # cmd = ["./src/indexing.sh", file]

    outcode: int = 1
    retry: int = 0

    # Call a process to index the file
    while retry < len(CMDS) and outcode != 0:
        bin: str = list(CMDS.keys())[retry]
        logger.debug(f"Compressing {file} with {bin} binary.")
        try:
            process: subprocess.CompletedProcess = runcmd(CMDS[bin])
            outcode: int = process.returncode
            logger.success(f"{file} has been successfully indexed.")
        except subprocess.CalledProcessError as e:
            logger.warning(f"Compressing {file} with {bin} binary did not succeed with exit code {outcode}.")
            logger.warning(f"Standard output: {e.stdout}")
            logger.warning(f"Standard error: {e.stderr}")
            retry += 1

    # Verify that the index file has been created
    # verify_files(file=file, index=index)

    # Return index path ONLY if exit code is 0
    return None if outcode else index

def save(obj: DataFrame, path: Path, format: str = "pickle", target: str = "L", lookup: dict = None, out: str = os.getcwd()) -> int:
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

        outf = str(Path(out,path.stem))

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


def is_file(path: str) -> bool:
    """ Check if path is a file """
    return os.path.isfile(path)


def is_empty(path: str) -> bool:
    """ Check if file is empty """
    return os.path.getsize(path) == 0


def is_indexed(path: str) -> bool:
    """ Check if file is indexed """
    pass

def is_compressed(path: str, type: str) -> bool:
    """ Checks if a file is a compressed archive type """
    return filetype.archive_match(path) == type

def verify_files(file: str, index: str = None) -> dict|None:
    """ Verify files exists and are not empty """

    files = {}

    def verify_VCF(file: str) -> str:

        assert isinstance(file, str), "VCF must be of string instance"

        TYPE = {"compression": "Gz"}

        if not is_file(file):

            raise errors.VCFError(f"Error: The file {file} does not exist.")

        if is_empty(file):

            raise errors.VCFError(f"Error: The file {file} is empty.")
        
        type = filetype.archive_match(file)

        # Is a archive
        if type:

            if type != TYPE["compression"]:

                raise errors.VCFError(f"Error: The compression of {file} is not supported.")
            
            # The file is a Gz archive
            else:

                try:

                    with open(file, mode='rb') as f:

                        header = f.read(3)

                        # Gzip magic number and flag byte (3rd byte)
                        # If 3rd bit (0x04) is set, header has extra field.
                        if not header[0:2] != b'\x1f\x8b' and header[3] & 0x04:

                            raise errors.VCFError(f"Error: {file} is not a BGZF archive")
                        
                        else:

                            # Check for BC extra sub-field

                            pass

                    with gzip.open(file,mode='rt') as f:
                                
                        line = f.readline()

                        # Check if first line is empty
                        if not line:

                            raise errors.VCFError(f"Error: First line of {file} is empty.")
                                
                        else:
                            # Check if first line start with "#"
                            if line[0] != '#':

                                raise errors.VCFError(f"Error: First line inconsistent with VCF header format")
                            
                except FileNotFoundError:

                    raise errors.VCFError(f"{file} is not a valid path")
                
                except IOError:

                    raise errors.VCFError(f"An error occurred while reading {file}")
                
            return "archive"
        
        # Is not a archive
        else:

            try:

                with open(file, mode='r') as f:

                    line = f.readline()

                    if not line:

                        raise errors.VCFError(f"Error: First line of {file} is empty.")
                        
                    else:
                        # Check if first line start with "#"
                        if line[0] != '#':

                            raise errors.VCFError(f"Error: First line inconsistent with VCF header format")
                        
            except FileNotFoundError:

                raise errors.VCFError(f"{file} is not a valid path")
                
            except IOError:

                raise errors.VCFError(f"An error occurred while reading {file}")
                    
            return "file"

    def verify_index(vcf: str, index: str) -> str:

        assert isinstance(index, str), "VCF must be of string instance"

        if not is_file(index):

            raise errors.IndexError(f"Error: The file {index} does not exist.")

        if is_empty(index):

            raise errors.IndexError(f"Error: The file {index} is empty.")
        
        # Index should be newer than VCF
        if os.path.getmtime(index) < os.path.getmtime(vcf):

            raise errors.IndexError(f"Error: {index} is older than {vcf}.")

        else:
            # Try using the index file
            cmd = ["tabix", "-l", vcf]

            try:

                runcmd(cmd)

            except subprocess.CalledProcessError:

                raise errors.IndexError(f"Error: {index} does not allow fast lookup for {vcf}")
            
        return "index"
    
    files[verify_VCF(file)] = file

    # Index is considered ONLY if file is an archive
    if "archive" in files and index:

        files[verify_index(index)] = index

    return files
    
def file_stats(path: str) -> dict:
    """ Get file stats """
    statinfo = os.stat(path)

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

# SETS


def intersect(a: set[str], b: set[str]) -> set:
    """ Return the intersection of two sets """
    return a & b


def difference(a: set[str], b: set[str]) -> set:
    """ Return the difference of two sets """
    return a - b

def jaccard_index(shared: int, total: dict) -> float:
    """ Calculate the Jaccard index """
    try:
        return shared / total
    except ZeroDivisionError:
        return None

# Variables

def convert(a: object) -> object:
    """ Convert variable to appropriate type """
    try:
        # If the variable contains a '/' or '|' character, it is a genotype information, return the variable as is
        # Else return the variable as an evaluated expression
        return a if any(list(map(lambda x: x in a,('/','|')))) else eval(a)
    except Exception:
        # If the variable cannot be evaluated, return the variable as is
        return a


# Variants

def evaluate(df: DataFrame) -> DataFrame:
    """ Evaluate the performance of the variant caller based on the truth set """

    assert "Filter.L" in df.columns, "Missing truth filter column"
    assert "Filter.R" in df.columns, "Missing query filter column"
    assert "Type" in df.columns, "Missing variant type column"

    # Create a mask of variants that passed the filter in one or both of the truth and query sets
    pass_mask = (df["Filter.L"] == "PASS") | (df["Filter.R"] == "PASS")

    series = []
    
    # Compute the performance metrics for SNPs and INDELs
    for v in ['snp','indel']:

        # Create a mask of variants that are of the specified type
        type_mask = df["Type"] == v

        # Filter the DataFrame based on the variant type and the filter mask
        filtered_df = df[type_mask & pass_mask]

        # If FILTER.L and FILTER.R are equal, the variant is a true positive
        is_match = filtered_df["Filter.L"] == filtered_df["Filter.R"]
        # If FILTER.L is PASS and FILTER.R is either missing or FAIL, the variant is a false negative
        is_truth_only = (filtered_df["Filter.L"] == "PASS") & ((isna(filtered_df["Filter.R"])) | (filtered_df["Filter.R"] == "FAIL"))
        # If FILTER.R is PASS and FILTER.L is either missing or FAIL, the variant is a false positive
        is_query_only = ((isna(filtered_df["Filter.L"])) | (filtered_df["Filter.R"] == "FAIL")) & (filtered_df["Filter.R"] == "PASS")

        # Create a Series of the performance metrics
        summary = Series([v, 
                          'PASS', 
                          (notna(filtered_df["Filter.L"])).sum(), 
                          is_match.sum(), 
                          is_truth_only.sum(), 
                          (notna(filtered_df["Filter.R"])).sum(), 
                          is_query_only.sum()], 
                          index=["TYPE","FILTER","TRUTH.TOTAL","TRUTH.TP","TRUTH.FN","QUERY.TOTAL","QUERY.FP"])
        
        # Compute the recall, precision, and F1 score
        try:
            summary["RECALL"] = summary["TRUTH.TP"] / (summary["TRUTH.TP"] + summary["TRUTH.FN"])
        except ZeroDivisionError:
            summary["RECALL"] = 0.00

        try:
            summary["PRECISION"] = summary["TRUTH.TP"] / (summary["TRUTH.TP"] + summary["QUERY.FP"])
        except ZeroDivisionError:
            summary["PRECISION"] = 0.00
            
        try:
            summary["F1"] = 2 * (summary["PRECISION"] * summary["RECALL"]) / (summary["PRECISION"] + summary["RECALL"])
        except ZeroDivisionError:
            summary["F1"] = 0.00

        series.append(summary)

    # Return the performance metrics as a DataFrame, column types are specified for better memory management
    return DataFrame(series).astype({"TYPE": "category", 
                                    "FILTER": "category", 
                                    "TRUTH.TOTAL": "int64", 
                                    "TRUTH.TP": "int64", 
                                    "TRUTH.FN": "int64", 
                                    "QUERY.TOTAL": "int64", 
                                    "QUERY.FP": "int64",
                                    "RECALL": "float64",
                                    "PRECISION": "float64",
                                    "F1": "float64"})

def hamming_distance(a: np.ndarray, b: np.ndarray) -> float:
    """ Calculate the Hamming distance """
    try:
        return 1 - ((np.sum(np.not_equal(a, b))) / (a.size + b.size))
    except ZeroDivisionError:
        return None

def is_homozygous(GT: str):
    """ Check if variant is homozygous """
    if '/' in GT:
        alleles = GT.split('/')
    elif '|' in GT:
        alleles = GT.split('|')
    
    return alleles[0] == alleles[1]

def is_heterozygous(GT: str):
    """ Check if variant is heterozygous """
    if '/' in GT:
        alleles = GT.split('/')
    elif '|' in GT:
        alleles = GT.split('|')
    
    return alleles[0] != alleles[1]

def exclude(v: object, filters: dict = None) -> bool:
    """ Check if variant should be excluded from analysis """
    return (
            v.is_indel and filters["exclude"]["exclude_indels"]
        ) or (
            v.is_snp and filters["exclude"]["exclude_snps"]
        ) or (
            v.is_mnp and filters["exclude"]["exclude_mnps"]
        ) or (
            v.is_sv and filters["exclude"]["exclude_svs"]
        ) or (
            v.is_transition and filters["exclude"]["exclude_transitions"]
        ) or (
            (v.FILTER != None) and filters["exclude"]["pass_only"]
        ) if filters else False

def format_to_values(format: str, values: str|list[str]) -> dict:
    """ map FORMAT string to respective SAMPLE values """

    # Split the format string into a list of fields
    format = format.split(":")

    # VCF is composed of multiples samples
    if isinstance(values,list):
        # Split the values string into a lists of list of values,
        # Each list of values represents a sample
        values = list(map(lambda x: x.split(":"), values))
        # Return a dictionary of the format and values mapped to each sample
        return {f"sample{s}": {f: convert(v)} for s in range(len(values)) for f, v in zip(format,values[s])}

    # VCF is composed of a unique sample
    else:
        values = values.split(":")

        return {f: convert(v) for f, v in zip(format,values)}