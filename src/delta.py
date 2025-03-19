import concurrent.futures
from cyvcf2 import VCF
from hashlib import sha256
from itertools import chain, repeat
from loguru import logger
from memory_profiler import profile
import numpy as np
from operator import itemgetter
from os import getcwd
from os.path import basename
from pandas import Series, DataFrame, Index, concat
from pathlib import Path
from plot import visualization
import subprocess
from sys import argv
from tabulate import tabulate
from template import Report
import pprint
import utils


@logger.catch
def process_chromosome(
    chrom: str,
    samples: list,
    header: dict,
    FILES: dict[str:str],
    filters: dict = None,
    compute: bool = False,
) -> dict:

    logger.debug(
        f"Processing chromosome {chrom} for file {FILES['compression']}"
    )

    try:
        vcf = VCF(FILES["compression"], lazy=True)
        vcf.set_index(index_path=FILES["index"])
    except FileNotFoundError as e:
        logger.error(e)

    FORMAT = {
        "genotype": ["GT"],
        "genotype_quality": ["GQ"],
        "depth": ["DP", "TRC"],
    }

    try:

        exclude: bool = False

        variants, filtered = {}, {
            "snp": 0,
            "mnp": 0,
            "indel": 0,
            "sv": 0,
            "transition": 0,
        }

        if compute:
            stats = {
                "variant": {
                    "snp": {
                        "transition": 0,
                        "transversion": 0,
                        "A": {"A": 0, "T": 0, "C": 0, "G": 0},
                        "T": {"A": 0, "T": 0, "C": 0, "G": 0},
                        "C": {"A": 0, "T": 0, "C": 0, "G": 0},
                        "G": {"A": 0, "T": 0, "C": 0, "G": 0},
                    },
                    "mnp": 0,
                    "indel": {"insertion": 0, "deletion": 0},
                    "sv": 0,
                },
                "depth": [],
                "quality": [],
                "GQ": [],
                "ref": 0,
                "het": 0,
                "hom": 0,
            }

        with utils.suppress_warnings():

            for i, v in enumerate(vcf(f"{chrom}")):

                parts = str(v).split("\t")

                parts[header["INFO"]] = '.'

                if not i:

                    format = parts[header["FORMAT"]]

                    logger.debug(f"FORMAT for chromosome {chrom}: {format}")

                samples_data = {
                    s: utils.format_to_values(
                        format=format, values=parts[header[s]]
                    )
                    for s in samples
                }

                if filters:

                    exclude: bool = utils.exclude(v, filters)

                if exclude:

                    filtered[v.var_type] += 1

                else:
                    
                    hash = sha256(
                        string=f"{v.CHROM}:{v.POS}:{v.REF}:{'|'.join(v.ALT)}".encode()
                    ).hexdigest()

                    if not hash in variants:

                        variants[hash] = [
                            v.CHROM,
                            v.POS,
                            v.var_type,
                            "FAIL" if v.FILTER else "PASS",
                            '\t'.join(parts),
                        ]
                        if compute:

                            if v.var_type == "snp":
                                if v.is_transition:
                                    mutation = "transition"
                                else:
                                    mutation = "transversion"
                                stats["variant"][v.var_type][mutation] += 1
                                stats["variant"][v.var_type][v.REF][v.ALT[0]] += 1
                            elif v.var_type == "indel":
                                if v.is_deletion:
                                    mutation = "deletion"
                                else:
                                    mutation = "insertion"
                                stats["variant"][v.var_type][mutation] += 1
                            else:
                                stats["variant"][v.var_type] += 1

                            if (
                                FORMAT["genotype_quality"][0]
                                in samples_data[samples[0]]
                            ):
                                stats[FORMAT["genotype_quality"][0]].append(
                                    [
                                        samples_data[s][
                                            FORMAT["genotype_quality"][0]
                                        ]
                                        for s in samples
                                    ]
                                )

                            stats["hom"] += sum(
                                list(
                                    map(
                                        lambda x: utils.is_homozygous(
                                            GT=samples_data[x][
                                                FORMAT["genotype"][0]
                                            ]
                                        ),
                                        samples,
                                    )
                                )
                            )
                            stats["het"] += sum(
                                list(
                                    map(
                                        lambda x: utils.is_heterozygous(
                                            GT=samples_data[x][
                                                FORMAT["genotype"][0]
                                            ]
                                        ),
                                        samples,
                                    )
                                )
                            )

                            if v.QUAL:
                                stats["quality"].append(v.QUAL)

                            stats["depth"].append(
                                [
                                    (
                                        samples_data[s][FORMAT["depth"][0]]
                                        if FORMAT["depth"][0] in samples_data[s]
                                        else [samples_data[s][FORMAT["depth"][1]]]
                                    )
                                    for s in samples
                                ]
                            )

    except UserWarning as e:
        logger.warning(e)

    variants = DataFrame.from_dict(
        variants,
        orient="index",
        columns=["Chromosome", "Position", "Type", "Filter", "Variant"],
    ).astype(
        {
            "Chromosome": "category",
            "Position": "int",
            "Type": "category",
            "Filter": "category",
            "Variant": "string[pyarrow]",
        }
    )

    variants.index = Index(variants.index.values, dtype="string[pyarrow]")

    logger.debug(
        f"Filtered: {filtered['snp']} SNP(s), {filtered['indel']} INDEL(s), {filtered['sv']} structural variant(s) variant(s) for chromosome {chrom} in file {FILES['compression']}"
    )

    if compute:
        stats["depth"], stats["quality"], stats["GQ"] = (
            np.array(stats["depth"], dtype=np.uint16),
            np.array(stats["quality"], dtype=np.float16),
            np.array(stats["GQ"], dtype=np.uint8),
        )

    return (variants, filtered, stats) if compute else (variants, filtered, {})


@logger.catch
def process_files(
    file: str, index: str = None, filters: dict = None, compute: bool = False
) -> dict:

    logger.debug(f"Processing file: {file}")

    try:
        utils.verify_file(file=file)
    except (FileNotFoundError, ValueError) as e:
        logger.error(e)

    FILES = {"compression": f"{file}.gz", "index": f"{file}.gz.tbi"}

    if index:

        try:
            utils.verify_file(file=index)
        except (FileNotFoundError, ValueError) as e:
            logger.error(e)

    else:
        # Try to look for indexing files
        try:
            utils.verify_file(FILES["compression"])
            utils.is_indexed(FILES["index"])
        except (FileNotFoundError, ValueError) as e:
            logger.warning(e)

            # Indexing
            logger.debug(f"Indexing file: {file}")

            try:
                code = subprocess.run(
                    ["./src/indexing.sh", file],
                    capture_output=True,
                    text=True,
                    check=True,
                )
            except (subprocess.CalledProcessError, FileNotFoundError) as e:
                logger.error(e.stderr)

            try:
                utils.verify_file(file=FILES["compression"])
            except (FileNotFoundError, ValueError) as e:
                logger.error(e)

    try:
        vcf = VCF(FILES["compression"])
        vcf.set_index(index_path=FILES["index"])
    except FileNotFoundError as e:
        logger.error(e)

    chromosomes: list = vcf.seqnames

    samples: list = vcf.samples

    HEADER = {
        "CHROM": 0,
        "POS": 1,
        "ID": 2,
        "REF": 3,
        "ALT": 4,
        "QUAL": 5,
        "FILTER": 6,
        "INFO": 7,
        "FORMAT": 8,
    }

    HEADER.update(
        {
            s: i
            for i, s in zip(
                range(
                    (HEADER["FORMAT"] + 1),
                    (HEADER["FORMAT"] + 1) + len(samples),
                ),
                samples,
            )
        }
    )

    logger.debug(
        f"{len(samples)} samples have been found in {file}: {samples}"
    )

    logger.debug(f"File {file} is composed of {chromosomes} chromosomes")

    logger.debug(
        f"Header for {file} has such format: {' '.join(HEADER.keys())}"
    )

    variants, filtered, stats = {}, {}, {}

    with concurrent.futures.ProcessPoolExecutor(
        max_workers=2
    ) as chrom_executor:

        futures_to_chrom = {
            chrom_executor.submit(
                process_chromosome,
                chrom,
                samples,
                HEADER,
                FILES,
                filters,
                compute,
            ): chrom
            for chrom in chromosomes
        }

        for future in concurrent.futures.as_completed(futures_to_chrom):

            logger.success(
                f"Process {future} for chromosome {futures_to_chrom[future]} in file {FILES['compression']} has completed."
            )

            try:
                (
                    variants[futures_to_chrom[future]],
                    filtered[futures_to_chrom[future]],
                    stats[futures_to_chrom[future]],
                ) = future.result()

                if compute:

                    stats[futures_to_chrom[future]]["length"] = vcf.seqlens[
                        chromosomes.index(futures_to_chrom[future])
                    ]

            except Exception as e:
                logger.warning(
                    f"Chromosome {futures_to_chrom[future]} generated an exception: {e}"
                )

    if compute:
        library = visualization(file=basename(file), stats=stats)

    return {
        "info": utils.file_stats(file),
        "header": "\t".join(list(HEADER.keys())),
        "variants": concat(
            list(itemgetter(*list(sorted(variants.keys())))(variants))
        ),
        "filter": filtered,
        "plots": library if compute else None,
    }


def delta(params: object) -> int:

    if not params.verbosity:

        logger.remove(0)

        logger.add("VCFDelta.log")

    logger.debug(f"VCFS: {params.vcfs}")

    logger.debug(
        f"{params.vcfs[0]} is set as truth"
        if params.truth
        else "No VCF has been set as truth"
    )

    logger.debug(f"Indexes: {params.indexes}")

    logger.debug(f"Serialize output: {params.serialize}")

    logger.debug(f"Output a report: {params.report}")

    assert len(params.vcfs) == 2, "Two VCF files are required."

    assert isinstance(params.vcfs[0], str) and isinstance(
        params.vcfs[1], str
    ), "Input vcf should be string instance"

    assert isinstance(params.process,int) and params.process >= 0, "Number of processes available must be an positive unsigned integer."

    if params.process <= 0:
        raise ValueError("Number of processes available must be an postive unsigned integer.")

    result = {}

    FILTERS = (
        {
            "threshold": params.threshold,
            "exclude": {
                "exclude_snps": params.exclude_snps,
                "exclude_indels": params.exclude_indels,
                "exclude_vars": params.exclude_vars,
                "exclude_mnps": params.exclude_mnps,
                "exclude_transitions": params.exclude_trans,
                "exclude_svs": params.exclude_svs,
                "pass_only": params.pass_only,
            },
        }
        if any(
            [
                params.threshold,
                params.exclude_snps,
                params.exclude_indels,
                params.exclude_vars,
                params.exclude_mnps,
                params.exclude_trans,
                params.exclude_svs,
                params.pass_only,
            ]
        )
        else None
    )

    logger.debug(f"Filters used: {FILTERS}")

    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as files_pool:

        futures_to_vcf = {
            files_pool.submit(
                process_files, vcf, index, FILTERS, params.report
            ): vcf
            for vcf, index in zip(params.vcfs, params.indexes)
        }

        for future in concurrent.futures.as_completed(futures_to_vcf):

            logger.success(
                f"Process {future} for file {futures_to_vcf[future]} has completed."
            )

            try:
                (result[futures_to_vcf[future]]) = future.result()
            except Exception as e:
                logger.error(
                    f"File {futures_to_vcf[future]} generated an exception: {e}"
                )

    result["delta"] = {
        "common": {},
        "unique": {params.vcfs[0]: {}, params.vcfs[1]: {}},
        "jaccard": {},
    }

    uniqueL = utils.difference(
        a=frozenset(result[params.vcfs[0]]["variants"].index),
        b=frozenset(result[params.vcfs[1]]["variants"].index),
    )

    uniqueR = utils.difference(
        a=frozenset(result[params.vcfs[1]]["variants"].index),
        b=frozenset(result[params.vcfs[0]]["variants"].index),
    )

    common = utils.intersect(
        a=frozenset(result[params.vcfs[0]]["variants"].index),
        b=frozenset(result[params.vcfs[1]]["variants"].index),
    )

    (
        result["delta"]["common"],
        result["delta"]["unique"][params.vcfs[0]],
        result["delta"]["unique"][params.vcfs[1]],
    ) = (
        len(common),
        len(uniqueL),
        len(uniqueR),
    )

    logger.debug(
        f"{result['delta']['common']} variant(s) is/are commom in both files"
    )

    logger.debug(
        f"{result['delta']['unique'][params.vcfs[0]]} variant(s) is/are unique in files {params.vcfs[0]}"
    )

    logger.debug(
        f"{result['delta']['unique'][params.vcfs[1]]} variant(s) is/are unique in files {params.vcfs[1]}"
    )

    result["delta"]["jaccard"] = utils.jaccard_index(
        shared=result["delta"]["common"],
        total=(
            result["delta"]["common"]
            + result["delta"]["unique"][params.vcfs[0]]
            + result["delta"]["unique"][params.vcfs[1]]
        ),
    )

    logger.debug(f"Jaccard index: {result['delta']['jaccard']}")

    list(
        map(
            (
                lambda x, n: x.rename(
                    columns={
                        c: f"{c}.{n}"
                        for c in x.columns
                        if not (
                            c in ["Chromosome", "Position", "Type"]
                            and n == "L"
                        )
                    },
                    inplace=True,
                )
            ),
            [
                result[params.vcfs[0]]["variants"],
                result[params.vcfs[1]]["variants"],
            ],
            ["L", "R"],
        )
    )

    df: DataFrame = concat(
        [
            result[params.vcfs[0]]["variants"],
            result[params.vcfs[1]]["variants"],
        ],
        axis=1,
        join="outer",
        sort=False,
    )

    df["Chromosome"] = df["Chromosome"].fillna(df["Chromosome.R"])
    df["Position"] = df["Position"].fillna(df["Position.R"])
    df["Type"] = df["Type"].fillna(df["Type.R"])

    df.drop(columns=["Chromosome.R", "Position.R", "Type.R"], inplace=True)

    df.reset_index(drop=False, names="Hash", inplace=True)

    df = df.astype(
        {
            "Hash": "string[pyarrow]",
            "Chromosome": "category",
            "Position": "int64",
            "Type": "category",
            "Filter.L": "category",
            "Filter.R": "category",
        }
    )

    df.sort_values(
        by=["Chromosome", "Position"],
        axis=0,
        ascending=True,
        inplace=True,
        kind="mergesort",
    )

    del result[params.vcfs[0]]["variants"]
    del result[params.vcfs[1]]["variants"]

    # Key is a references to the Hash string object in the Dataframe
    # Allow a O(1) lookup while keeping memory footprint low
    lookup = {hash: row for row, hash in enumerate(df["Hash"])}

    if params.truth:

        summary = utils.evaluate(df)

    if params.serialize:

        with concurrent.futures.ProcessPoolExecutor(max_workers=2) as files_pool:

            futures_to_vcf = {
                files_pool.submit(
                    utils.save, df, vcf, params.serialize, t, lookup
                ): vcf for vcf, t in zip([Path(params.vcfs[0]),Path(params.vcfs[1])],['L','R'])
            }

            for future in concurrent.futures.as_completed(futures_to_vcf):

                try:
                   out: int = future.result()
                   logger.success(f"Process {future} has successfully serialized {futures_to_vcf[future]}") if out else logger.error(f"Process {future} did not successfully serialized with a {out} output code.")
                except ValueError as e:
                    logger.error(e)

    if params.report:

        Report(
            vcfs=params.vcfs,
            prefix=params.out,
            cmd=" ".join(argv),
            infos={
                params.vcfs[0]: result[params.vcfs[0]]["info"],
                params.vcfs[1]: result[params.vcfs[1]]["info"],
            },
            view={
                "headers": {
                    params.vcfs[0]: result[params.vcfs[0]]["header"],
                    params.vcfs[1]: result[params.vcfs[1]]["header"],
                },
                "variants": df,
                "stats": result["delta"]
            },
            plots={
                params.vcfs[0]: result[params.vcfs[0]]["plots"],
                params.vcfs[1]: result[params.vcfs[1]]["plots"],
            },
            summary=summary if params.truth else None,
        ).create()

    else:
        print(f"{params.vcfs[0]}: [{result['delta']['unique'][params.vcfs[0]]} unique]────[{result['delta']['common']} common]────[{result['delta']['unique'][params.vcfs[1]]} unique] :{params.vcfs[1]}")
        print(f"Jaccard index: {result['delta']['jaccard']}")
        if params.truth:
            print(tabulate(summary,headers='keys',tablefmt='grid',numalign='center', stralign='center'))

    return 1
