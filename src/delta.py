import concurrent.futures
from cyvcf2 import VCF
from hashlib import sha256
from itertools import chain, repeat
from loguru import logger
from memory_profiler import profile
import numpy as np
from operator import itemgetter
from os.path import basename, getsize
from pandas import DataFrame, Index, concat
from pathlib import Path
from plot import visualization, PlotLibrary
import subprocess
from sys import argv
from tabulate import tabulate
from template import Report
import pprint
from psutil import cpu_count
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
    """
    Process a chromosome from a VCF file
        chrom: containing the chromosome to process
        samples: containing the samples to process
        header: containing the header of the VCF file
        FILES: containing the path to the VCF file and the index file
        filters: containing the filters to apply
        compute: to compute the statistics if a report is wanted
    """

    logger.debug(
        f"Processing chromosome {chrom} for file {FILES['compression']}"
    )

    # Open the VCF file and set the index for faster lookup
    try:
        vcf = VCF(FILES["compression"], lazy=True)
        vcf.set_index(index_path=FILES["index"])
    except FileNotFoundError as e:
        logger.error(e)

    # Map the FORMAT values to there respectives informations
    FORMAT = {
        "genotype": ["GT"],
        "genotype_quality": ["GQ"],
        "depth": ["DP", "TRC"],
    }

    try:
        # At first, set the filters to False
        exclude: bool = False
        # Save the filtered variants number of operations
        variants, filtered = {}, {
            "snp": 0,
            "mnp": 0,
            "indel": 0,
            "sv": 0,
            "transition": 0,
        }
        # If a report is wanted, set the statistics to 0
        # If no report is wanted, the dictionary is not created reducing memory footprint
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

        # Suppress warnings from cyvcf2 in case of missing values
        with utils.suppress_warnings():

            for i, v in enumerate(vcf(f"{chrom}")):
                # Get the values from the VCF file
                parts = str(v).split("\t")
                # Set the INFO values as a single character to reduce memory footprint
                parts[header["INFO"]] = '.'
                # First iteration, get the FORMAT values
                if not i:

                    format = parts[header["FORMAT"]]

                    logger.debug(f"FORMAT for chromosome {chrom}: {format}")
                # From FORMAT get the values for each sample
                samples_values = {
                    s: utils.format_to_values(
                        format=format, values=parts[header[s]]
                    )
                    for s in samples
                }
                # Should variant be filtered ?
                if filters:

                    exclude: bool = utils.exclude(v, filters)
                # Should the variant be excluded ?
                if exclude:

                    filtered[v.var_type] += 1
                # Variants pass the filters
                else:
                    # Hash the variant to avoid duplicates and enhance lookup
                    hash = sha256(
                        string=f"{(v.CHROM).removeprefix('chr')}:{v.POS}:{v.REF}:{'|'.join(v.ALT)}".encode()
                    ).hexdigest()

                    # Is it a new variant ?
                    if not hash in variants:
                        # Save the variant
                        variants[hash] = [
                            (v.CHROM).removeprefix('chr'),
                            v.POS,
                            v.var_type,
                            "FAIL" if v.FILTER else "PASS",
                            '\t'.join(parts),
                        ]
                        # Should the statistics be computed ?
                        if compute:
                            # The variant type is common to all samples
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

                            # Statistics unique to each samples
                            if (
                                FORMAT["genotype_quality"][0]
                                in samples_values[samples[0]]
                            ):
                                stats[FORMAT["genotype_quality"][0]].append(
                                    [
                                        samples_values[s][
                                            FORMAT["genotype_quality"][0]
                                        ]
                                        for s in samples
                                    ]
                                )

                            stats["hom"] += sum(
                                list(
                                    map(
                                        lambda x: utils.is_homozygous(
                                            GT=samples_values[x][
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
                                            GT=samples_values[x][
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
                                        samples_values[s][FORMAT["depth"][0]]
                                        if FORMAT["depth"][0] in samples_values[s]
                                        else [samples_values[s][FORMAT["depth"][1]]]
                                    ) # Try to get the depth value from the first FORMAT value, if not found, get it from the second FORMAT value
                                    for s in samples
                                ]
                            )

    except UserWarning as e:
        logger.warning(e)

    # Close the data stream, avoid memory leaks
    vcf.close()

    # Create a DataFrame from the variants,
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
    ) # Set the columns to the right type for better memory management

    # Set the hash values as the index
    variants.index = Index(variants.index.values, dtype="string[pyarrow]")

    logger.debug(
        f"Filtered: {filtered['snp']} SNP(s), {filtered['indel']} INDEL(s), {filtered['sv']} structural variant(s) variant(s) for chromosome {chrom} in file {FILES['compression']}"
    )

    # Set the statistics computed as the right type for better memory management
    if compute:
        stats["depth"], stats["quality"], stats["GQ"] = (
            np.array(stats["depth"], dtype=np.uint16),
            np.array(stats["quality"], dtype=np.float16),
            np.array(stats["GQ"], dtype=np.uint8),
        )

    return (variants, filtered, stats) if compute else (variants, filtered, {})


@logger.catch
def process_files(
    file: str, index: str = None, filters: dict = None, compute: bool = False, pool: int = 1
) -> dict:
    """ 
    Main function to process a VCF file
        file: containing the path to the VCF file
        index: containing the path to the index file
        filters: containing the filters to apply
        compute: to compute the statistics if a report is wanted
        pool: to set the number of process available
    """

    logger.debug(f"Processing file: {file}")

    # Check if the file exists and is not empty
    try:
        utils.verify_file(file=file)
    except (FileNotFoundError, ValueError) as e:
        logger.error(e)

    # Map files to there respectives extensions
    FILES = {"compression": f"{file}.gz", "index": f"{file}.gz.tbi"}

    # Is a index file provided ?
    if index:
        # Check if the file exists and is not empty
        try:
            utils.verify_file(file=index)
        except (FileNotFoundError, ValueError) as e:
            logger.error(e)
    # No index file provided
    else:
        # Try to look for indexing files
        try:
            utils.verify_file(FILES["compression"])
            utils.is_indexed(FILES["index"])
        except (FileNotFoundError, ValueError) as e:
            logger.warning(e)

            # Indexing
            logger.debug(f"Indexing file: {file}")

            # Call a process to index the file
            try:
                code = subprocess.run(
                    ["./src/indexing.sh", file],
                    capture_output=True,
                    text=True,
                    check=True,
                )
            except (subprocess.CalledProcessError, FileNotFoundError) as e:
                logger.error(e.stderr)

            # Verify that the index file has been created
            try:
                utils.verify_file(file=FILES["compression"])
            except (FileNotFoundError, ValueError) as e:
                logger.error(e)
    # From this point, the VCF file is supposed indexed
    try:
        # Open the VCF file and set the index for faster lookup
        vcf = VCF(FILES["compression"])
        vcf.set_index(index_path=FILES["index"])
    except FileNotFoundError as e:
        logger.error(e)

    # Get the chromosomes from the VCF file
    chromosomes: list = vcf.seqnames

    logger.debug(f"File {file} is composed of {chromosomes} chromosomes")

    # Get the samples from the VCF file
    samples: list = vcf.samples

    logger.debug(
        f"{len(samples)} samples have been found in {file}: {samples}"
    )

    # Map the header values to there respectives indexes
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

    # Add the samples to the mapper
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
        f"Header for {file} has such format: {' '.join(HEADER.keys())}"
    )

    variants, filtered, stats = {}, {}, {}

    # Process the chromosomes concurrently, using as much process as available
    # use a context manager to avoid memory leaks
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=pool
    ) as chrom_executor:
        # Submit the processes
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
            for chrom in chromosomes # VCF object cannot be pickled thus cannot be passed to a process
        }
        # Check if a process is completed
        for future in concurrent.futures.as_completed(futures_to_chrom):

            logger.success(
                f"Process {future} for chromosome {futures_to_chrom[future]} in file {FILES['compression']} has completed."
            )
            # Get the returned result
            try:
                (
                    variants[futures_to_chrom[future]],
                    filtered[futures_to_chrom[future]],
                    stats[futures_to_chrom[future]],
                ) = future.result()
                # If a report is wanted, add the length of the chromosome
                if compute:

                    stats[futures_to_chrom[future]]["length"] = vcf.seqlens[
                        chromosomes.index(futures_to_chrom[future])
                    ]

            except Exception as e:
                logger.warning(
                    f"Chromosome {futures_to_chrom[future]} generated an exception: {e}"
                )
    # Close the data stream
    vcf.close()

    # If a report is wanted, create plots
    if compute:
        library = visualization(file=basename(file), stats=stats)

    return {
        "info": utils.file_stats(file), # Get the information about the VCF file
        "header": "\t".join(list(HEADER.keys())), # Get the header of the VCF file
        "variants": concat(
            list(itemgetter(*list(sorted(variants.keys())))(variants))
        ), # Concatenate the chromosomes DataFrames
        "filter": filtered, # Get the stats about filtered variants
        "plots": library if compute else None, # Get the plots
    }


def delta(params: object) -> int:
    """ 
    Main function to compute the delta between two VCF files 
            params: object containing the parameters parsed by the CLI
    """

    # Should the log be printed to CLI or saved in a file ?
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
    else:
        if params.process > cpu_count(logical=True):
            params.process = cpu_count(logical=True)

    # Number of files to process
    PROCESS_FILE = 2

    result = {}

    # Filters to apply during the parsing
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
        # If no filters are set, used None for better memory management
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

    # Process the files concurrently, using as much process as VCF files inputed
    # use a context manager to avoid memory leaks
    with concurrent.futures.ProcessPoolExecutor(max_workers=PROCESS_FILE) as files_pool:
        # Check if the number of process available is odd or even
        pavailable = (params.process - PROCESS_FILE) % 2

        # Allocate to each file executor the number of process available for further parallelization by chromosomes
        palloc = [(params.process - PROCESS_FILE)//2]*2

        # If the number of process available is odd, allocate the remaining process to the file with the biggest size
        if pavailable != 0:
            
            fsizes = list(map(getsize,params.vcfs))

            maxid = 0 if fsizes[0] > fsizes[1] else 1

            palloc[maxid] += 1
        # Submit the process for each file mapped to the vcf path
        futures_to_vcf = {
            files_pool.submit(
                process_files, vcf, index, FILTERS, params.report, proc
            ): vcf
            for vcf, index, proc in zip(params.vcfs, params.indexes, palloc)
        }
        # Check if a process is completed
        for future in concurrent.futures.as_completed(futures_to_vcf):

            logger.success(
                f"Process {future} for file {futures_to_vcf[future]} has completed."
            )
            # Get the returned result
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
    # Compute unique variants in first VCF file
    uniqueL: set = utils.difference(
        a=frozenset(result[params.vcfs[0]]["variants"].index),
        b=frozenset(result[params.vcfs[1]]["variants"].index),
    )
    # Compute unique variants in second VCF file
    uniqueR: set = utils.difference(
        a=frozenset(result[params.vcfs[1]]["variants"].index),
        b=frozenset(result[params.vcfs[0]]["variants"].index),
    )
    # Compute common variants between both VCF files
    common: set = utils.intersect(
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

    # Compute the Jaccard Index
    result["delta"]["jaccard"] = utils.jaccard_index(
        shared=result["delta"]["common"],
        total=(
            result["delta"]["common"]
            + result["delta"]["unique"][params.vcfs[0]]
            + result["delta"]["unique"][params.vcfs[1]]
        ),
    )

    logger.debug(f"Jaccard index: {result['delta']['jaccard']}")

    # Rename columns to avoid conflicts, inplace for better memory management
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
    # Merge the two DataFrames with a outer join algorithm,
    # Merge is made on the hash index
    df: DataFrame = concat(
        [
            result[params.vcfs[0]]["variants"],
            result[params.vcfs[1]]["variants"],
        ],
        axis=1,
        join="outer",
        sort=False,
    )
    # Fill missing values with the values from the other VCF file
    # This is done to avoid NaN values in the DataFrame
    # These specific columns are used to identify the variants and are common to both VCF files
    df["Chromosome"] = df["Chromosome"].fillna(df["Chromosome.R"])
    df["Position"] = df["Position"].fillna(df["Position.R"])
    df["Type"] = df["Type"].fillna(df["Type.R"])

    # Drop redondant columns to reduce memory footprint
    df.drop(columns=["Chromosome.R", "Position.R", "Type.R"], inplace=True)

    # Reset the index for later use
    df.reset_index(drop=False, names="Hash", inplace=True)

    # Convert the DataFrame columns for better memory management,
    # Make use of PyArrow for better performance to store string values
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
    # Sort values with a stable algorithm for better view in the report
    df.sort_values(
        by=["Chromosome", "Position"],
        axis=0,
        ascending=True,
        inplace=True,
        kind="mergesort",
    )

    # Delete DataFrames not of use anymore, to reduce memory footprint
    del result[params.vcfs[0]]["variants"]
    del result[params.vcfs[1]]["variants"]

    # Key is a references to the Hash string object in the Dataframe
    # Allow a O(1) lookup while keeping memory footprint low
    lookup = {hash: row for row, hash in enumerate(df["Hash"])}

    # Should the truth be computed ?
    if params.truth:

        summary = utils.evaluate(df)

    # Should the output be serialized ?
    if params.serialize:

        # Parallelize the serialization process, with as much process as VCF files inputed
        # use a context manager to avoid memory leaks
        with concurrent.futures.ProcessPoolExecutor(max_workers=PROCESS_FILE) as files_pool:
            # Submit the process for each file mapped to the vcf path
            futures_to_vcf = {
                files_pool.submit(
                    utils.save, df, vcf, params.serialize, t, lookup
                ): vcf for vcf, t in zip([Path(params.vcfs[0]),Path(params.vcfs[1])],['L','R'])
            }
            # Check if a process is completed
            for future in concurrent.futures.as_completed(futures_to_vcf):
                # Get the returned result
                try:
                   out: int = future.result()
                   logger.success(f"Process {future} has successfully serialized {futures_to_vcf[future]}") if out else logger.error(f"Process {future} did not successfully serialized with a {out} output code.")
                except ValueError as e:
                    logger.error(e)

    # Should the output be reported ?
    if params.report:

        # Library of common plots between the two VCF files
        pcommon: PlotLibrary = PlotLibrary()

        # Create a Venn diagram to display the common, unique variants between the two VCF files
        pcommon.venn((result["delta"]["unique"][params.vcfs[0]], result["delta"]["unique"][params.vcfs[1]], result["delta"]["common"]), ['L','R'])

        # Create a report with the results
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
                "common": pcommon
            },
            summary=summary if params.truth else None,
        ).create()
    # Print the results to the CLI
    else:
        print(f"{params.vcfs[0]}: [{result['delta']['unique'][params.vcfs[0]]} unique]────[{result['delta']['common']} common]────[{result['delta']['unique'][params.vcfs[1]]} unique] :{params.vcfs[1]}")
        print(f"Jaccard index: {result['delta']['jaccard']}")
        if params.truth:
            print(tabulate(summary,headers='keys',tablefmt='grid',numalign='center', stralign='center'))

    return 1
