import contextlib
from cyvcf2 import VCF, Writer
from datetime import datetime, timezone
from hashlib import sha256
import json
from loguru import logger
import numpy as np
import os
from pathlib import Path
from pandas import Series, DataFrame, notna, isna
import warnings
import _pickle as cPickle

# I/O

def save(obj: DataFrame, path: Path, format: str = "pickle", target: str = "L", lookup: dict = None, out: str = os.getcwd()) -> int:

    assert (
        format == "json" or format == "pickle" or format == "vcf"
    ), "File format is not supported"

    assert (isinstance(obj,DataFrame)), "Input provided is not a Dataframe instance"

    FILES = {
            "json": {"ext": "json", "mode": "w", "func": json.dump},
            "pickle": {"ext": "pkl", "mode": "wb", "func": cPickle.dump},
            "vcf": {"ext": "vcf", "mode": "wz", "func": None}
        }

    if format in FILES:

        outf = str(Path(out,path.stem))

        if format in ["json", "pickle"]:
            with open(
                file=f"{outf}_delta.{FILES[format]['ext']}",
                mode=FILES[format]["mode"],
            ) as f:
                obj = obj.to_dict(orient='list')
                FILES[format]["func"](obj, f)
        else:
            vcf = VCF(f"{str(path)}.gz")
            vcf.add_info_to_header({'ID': 'match', 'Description': 'overlapping variant', 'Type': 'String', 'Number': '1'})
            w = Writer(f"{outf}_delta.{FILES[format]['ext']}",vcf)
            for i, v in enumerate(vcf):
                hash = sha256(
                    string=f"{v.CHROM}:{v.POS}:{v.REF}:{'|'.join(v.ALT)}".encode()
                ).hexdigest()
                series = obj.iloc[lookup[hash]]
                v.INFO["match"] = int((notna(series["Variant.L"])) & (notna(series["Variant.R"])))
                w.write_record(v)
            w.close(); vcf.close()

        assert os.path.isfile(
            f"{outf}_delta.{FILES[format]['ext']}"
        ), "File was not created"

        logger.debug(f"Results are seralized to {out}")

        return 1
    else:
        raise ValueError(f"Error: The file format {format} is not supported.")


def is_file(path: str) -> bool:
    return os.path.isfile(path)


def is_empty(path: str) -> bool:
    return os.path.getsize(path) == 0


def is_indexed(path: str) -> bool:
    return verify_file(path)


def verify_file(file: str):

    assert isinstance(file, str), "Values must be of string instance"

    if not is_file(file):

        raise FileNotFoundError(f"Error: The file '{file}' does not exist.")

    if is_empty(file):

        raise ValueError(f"Error: The file '{file}' is empty.")
    
def file_stats(path: str) -> dict:
    
    statinfo = os.stat(path)

    return {"basename": os.path.basename(path),
            "path": os.path.dirname(path),
            "size": round(statinfo.st_size / pow(1024,2),2),
            "mtime": datetime.fromtimestamp(statinfo.st_mtime, tz=timezone.utc)}

@contextlib.contextmanager
def suppress_warnings():
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="no intervals found for*")
        warnings.filterwarnings("ignore", message="Mean of empty slice.")
        yield

# SETS


def intersect(a: set[str], b: set[str]) -> set:
    return a & b


def difference(a: set[str], b: set[str]) -> set:
    return a - b

def jaccard_index(shared: int, total: dict) -> float:
    try:
        return shared / total
    except ZeroDivisionError:
        return None

# Variables

def convert(a: object) -> object:

    try:
        return a if any(list(map(lambda x: x in a,('/','|')))) else eval(a)
    except Exception:
        return a


# Variants

def evaluate(df: DataFrame) -> DataFrame:

    assert "Filter.L" in df.columns, "Missing truth filter column"
    assert "Filter.R" in df.columns, "Missing query filter column"
    assert "Type" in df.columns, "Missing variant type column"

    pass_mask = (df["Filter.L"] == "PASS") | (df["Filter.R"] == "PASS")

    series = []
        
    for v in ['snp','indel']:

        type_mask = df["Type"] == v

        filtered_df = df[type_mask & pass_mask]

        is_match = filtered_df["Filter.L"] == filtered_df["Filter.R"]
        is_truth_only = (filtered_df["Filter.L"] == "PASS") & ((isna(filtered_df["Filter.R"])) | (filtered_df["Filter.R"] == "FAIL"))
        is_query_only = ((isna(filtered_df["Filter.L"])) | (filtered_df["Filter.R"] == "FAIL")) & (filtered_df["Filter.R"] == "PASS")

        summary = Series([v, 
                          'PASS', 
                          (notna(filtered_df["Filter.L"])).sum(), 
                          is_match.sum(), 
                          is_truth_only.sum(), 
                          (notna(filtered_df["Filter.R"])).sum(), 
                          is_query_only.sum()], 
                          index=["TYPE","FILTER","TRUTH.TOTAL","TRUTH.TP","TRUTH.FN","QUERY.TOTAL","QUERY.FP"])
        
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
    try:
        return 1 - ((np.sum(np.not_equal(a, b))) / (a.size + b.size))
    except ZeroDivisionError:
        return None

def is_homozygous(GT: str):
    if '/' in GT:
        alleles = GT.split('/')
    elif '|' in GT:
        alleles = GT.split('|')
    
    return alleles[0] == alleles[1]

def is_heterozygous(GT: str):
    if '/' in GT:
        alleles = GT.split('/')
    elif '|' in GT:
        alleles = GT.split('|')
    
    return alleles[0] != alleles[1]

def exclude(v: object, filters: dict = None) -> bool:

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

    format = format.split(":")

    if isinstance(values,list):

        values = list(map(lambda x: x.split(":"), values))

        return {f"sample{s}": {f: convert(v)} for s in range(len(values)) for f, v in zip(format,values[s])}

    else:

        values = values.split(":")

        return {f: convert(v) for f, v in zip(format,values)}