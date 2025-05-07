import concurrent.futures
import errors
import files
from hashlib import sha256
from itertools import chain, repeat
from loguru import logger
from memory_profiler import profile
import numpy as np
from operator import itemgetter
from os.path import basename, getsize
from pandas import DataFrame, Index, concat
from plots import visualization, PlotLibrary
from sys import argv
from tabulate import tabulate
from psutil import cpu_count
import utils

def process_chromosome(
    chrom: str,
    file: files.VCF,
    compute: bool = False,
) -> dict:
    """
    Process a chromosome from a VCF file
        chrom: containing the chromosome to process
        file: containing the path to the VCF file and the index file
        compute: to compute the statistics if a report is wanted
    """

    logger.debug(
        f"Processing chromosome {chrom} for file {file}"
    )

    # Open the VCF file and set the index for faster lookup
    try:
        file.open(context=False)
    except (FileNotFoundError, errors.VCFError) as e:
        logger.error(e)

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
            stats: dict = {
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

            # Record warnings for caller dependent FORMAT field
            warnings: dict = {}

        # Suppress warnings from cyvcf2 in case of missing values
        with utils.suppress_warnings():
            # Iterate over the VCF file
            for i, v in enumerate(file.stdin(f"{chrom}")):
                # Get the values from the VCF file
                parts: list[str] = str(v).split("\t")
                # Set the INFO values as a single character to reduce memory footprint
                parts[file.header["INFO"]] = '.'
                # First iteration, get the FORMAT values
                if not i:

                    format: str = parts[file.header["FORMAT"]]

                    logger.debug(f"FORMAT for chromosome {chrom}: {format}")
                # From FORMAT get the values for each sample
                samples_values: dict[str:dict] = {
                    s: utils.format_to_values(
                        format=format, values=parts[file.header[s]]
                    )
                    for s in file.samples
                }
                # Should variant be filtered ?
                if file.variants.filters:

                    exclude: bool = utils.exclude(v, file.variants.filters)
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
                                    mutation: str = "transition"
                                else:
                                    mutation: str = "transversion"
                                stats["variant"][v.var_type][mutation] += 1
                                stats["variant"][v.var_type][v.REF][v.ALT[0]] += 1
                            elif v.var_type == "indel":
                                if v.is_deletion:
                                    mutation: str = "deletion"
                                else:
                                    mutation: str = "insertion"
                                stats["variant"][v.var_type][mutation] += 1
                            else:
                                stats["variant"][v.var_type] += 1

                            # Statistics unique to each samples
                            # If previous passes have raised a warning, do not search for genotype quality score.
                            if not "genotype_quality" in warnings:
                                try:
                                    stats[file.format["genotype_quality"][0]].append(
                                        [
                                            samples_values[s][
                                                file.format["genotype_quality"][0]
                                            ]
                                            for s in file.samples
                                        ]
                                    )
                                except KeyError:
                                    logger.warning(
                                        f"Genotype quality value cannot be retrieved with key(s): {file.format['genotype_quality']}"
                                    )
                                    # Keep record of exception
                                    warnings["genotype_quality"] = True

                            # If previous passes have raised a warning, do not search for genotype.
                            if not "genotype" in warnings:
                                try:
                                    # Sum homzygous genotypes for each sample
                                    stats["hom"] += sum(
                                        list(
                                            map(
                                                lambda sample: utils.is_homozygous(
                                                    GT=samples_values[sample][
                                                        file.format["genotype"][0]
                                                    ]
                                                ),
                                                file.samples,
                                            )
                                        )
                                    )
                                    # Sum heterozygous genotypes for each sample
                                    stats["het"] += sum(
                                        list(
                                            map(
                                                lambda sample: utils.is_heterozygous(
                                                    GT=samples_values[sample][
                                                        file.format["genotype"][0]
                                                    ]
                                                ),
                                                file.samples,
                                            )
                                        )
                                    )
                                except KeyError:
                                    logger.warning(f"Genotype type cannot be retrieved with key(s): {file.format['genotype']}")
                                    # Keep record of exception
                                    warnings["genotype"] = True

                            # Save the call quality
                            if v.QUAL:
                                stats["quality"].append(v.QUAL)

                            # If previous passes have raised a warning, do not search for depth metric.
                            if not "depth" in warnings:
                                try:
                                    stats["depth"].append(
                                        [
                                            (
                                                samples_values[s][file.format["depth"][0]]
                                                if file.format["depth"][0] in samples_values[s]
                                                else [samples_values[s][file.format["depth"][1]]]
                                            ) # Try to get the depth value from the first FORMAT value, if not found, get it from the second FORMAT value
                                            for s in file.samples
                                        ]
                                    )
                                except KeyError:
                                    logger.warning(f"Sequencing depth value cannot be retrieved with key(s): {file.format['depth']}")
                                    # Keep record of exception
                                    warnings["depth"] = True
    except UserWarning as e:
        logger.warning(e)

    # Close the data stream, avoid memory leaks
    file.close()

    # Create a DataFrame from the variants,
    variants: DataFrame = DataFrame.from_dict(
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
        f"Filtered: {filtered['snp']} SNP(s), {filtered['indel']} INDEL(s), {filtered['sv']} structural variant(s) variant(s) for chromosome {chrom} in file {file}"
    )

    # Set the statistics computed as the right type for better memory management
    if compute:
        stats["depth"], stats["quality"], stats["GQ"] = (
            np.array(stats["depth"], dtype=np.uint16),
            np.array(stats["quality"], dtype=np.float16),
            np.array(stats["GQ"], dtype=np.uint16),
        )

    return (variants, filtered, stats) if compute else (variants, filtered, {})


def process_vcf(
    file: files.VCF, compute: bool = False, pool: int = 1
) -> dict:
    """ 
    Main function to process a VCF file
        file: containing the path to the VCF file
        index: containing the path to the index file
        compute: to compute the statistics if a report is wanted
        pool: to set the number of process available
    """

    logger.debug(f"Processing file: {file}")

    assert pool >= 0, "Number of processes alocated must be an positive unsigned integer."

    if pool < 0:
        raise ValueError("Number of processes alocated must be an positive unsigned integer.")

    logger.debug(f"Computation will be {'parallelized' if pool else 'made sequentially'} by chromosome(s).")

    # Check if the file exists and is not empty
    try:

        file.verify()

        if file.is_indexed():

            file.index.verify()

    except (errors.VCFError, errors.IndexError) as e:

        logger.error(e)

        # Error raised by VCF file are critical
        if isinstance(e,errors.VCFError):
            logger.error(f"{file} is not a valid VCF file.")
            raise errors.ProcessError(f"{file} is not a valid VCF file.")
        
        # Errors caused by the provided index are treated as warnings.
        elif isinstance(e,errors.IndexError):
            logger.warning(f"{file.index} is not a valid VCF index file.")
            file.index = None
            
    # File must be indexed
    if not file.is_indexed():
        
        # Compression
        if not file.archive:

            file.compressing()

        # Call a process to index the file
        file.indexing()
    
    if file.archive and file.index:
        # From this point, the VCF file is supposed indexed
        try:
            # Open the VCF with CyVCF2 file and set the index for faster lookup
            file.open(context=True)
        except (FileNotFoundError, errors.VCFError) as e:
            logger.error(e)

        variants, filtered, stats = {}, {}, {}

        if pool:
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
                        file,
                        compute,
                    ): chrom
                    for chrom in file.chromosomes # VCF object cannot be pickled thus cannot be passed to a process
                }
                # Check if a process is completed
                for future in concurrent.futures.as_completed(futures_to_chrom):

                    logger.success(
                        f"Process {future} for chromosome {futures_to_chrom[future]} in file {file} has completed."
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

                            stats[futures_to_chrom[future]]["length"] = file.seqlens[
                                file.chromosomes.index(futures_to_chrom[future])
                            ]

                    except Exception as e:
                        logger.warning(
                            f"Chromosome {futures_to_chrom[future]} generated an exception: {e}"
                        )
        # Computation is carried out sequentially.                
        else:
            for chrom in file.chromosomes:
                variants[chrom], filtered[chrom], stats[chrom] = process_chromosome(chrom=chrom, 
                                                                                    file=file, 
                                                                                    compute=compute)

        # If a report is wanted, create plots
        if compute:
            library = visualization(file=basename(file.path), stats=stats)

        return {
            "info": utils.file_stats(file.path), # Get the information about the VCF file
            "header": "\t".join(list(file.header.keys())), # Get the header of the VCF file
            "variants": concat(
                list(itemgetter(*list(sorted(variants.keys())))(variants))
            ), # Concatenate the chromosomes DataFrames
            "filter": filtered, # Get the stats about filtered variants
            "plots": library if compute else None, # Get the plots
        }
    
    else:

        raise errors.CompressionIndexError(f"Failed to compress and index {file}")


def delta(params: object) -> int:
    """ 
    Main function to compute the delta between two VCF files 
            params: object containing the command line parameters parsed.
    """

    logger.debug(f"VCFS: {params.vcfs}")

    logger.debug(
        f"{params.vcfs[0]} is set as truth."
        if params.benchmark
        else "No VCF has been set as truth."
    )

    logger.debug(f"Indexes: {params.indexes}")

    logger.debug(f"Serialize output: {params.serialize}")

    logger.debug(f"Output a report: {params.report}")

    assert len(params.vcfs) == 2, "Two VCF files are required."

    assert isinstance(params.vcfs[0], str) and isinstance(
        params.vcfs[1], str
    ), "Input VCF should be string instance"

    # Number of files to process
    PROCESS_FILE: int = 2

    # Number of CPUs available
    CPUS = cpu_count(logical=True)

    logger.debug(f"Number of logical(s) CPU(s) detected: {CPUS}.")

    assert isinstance(params.process,int) and params.process >= 0, "Number of processes alocated must be an positive unsigned integer."

    if params.process < 0:
        raise ValueError("Number of processes alocated must be an postive unsigned integer.")
    else:
        if CPUS > 1:
            logger.debug("Computation will be parallelized by the number of VCF files.")
            params.process = min(params.process, cpu_count(logical=True))
        else:
            logger.debug("Computation will be carried out sequentially.")

    filters: dict[str:bool] = {"SNP": params.exclude_snps,
                               "TRANSITION": params.exclude_trans,
                               "MNP": params.exclude_mnps,
                               "INDELS": params.exclude_indels,
                               "SV": params.exclude_svs,
                               "VARS": params.exclude_vars,
                               "PASS_ONLY": params.pass_only}
    
    logger.debug(f"Filters used: {filters}")

    vcfs: list[files.VCF] = [files.VCF(path=params.vcfs[0], reference=True if params.benchmark else False, index=params.indexes[0], filters=filters), 
                             files.VCF(path=params.vcfs[1], index=params.indexes[1], filters=filters)]

    results: dict = {}

    if params.process:

        # Process the files concurrently, using as much process as VCF files inputed
        # use a context manager to avoid memory leaks
        with concurrent.futures.ProcessPoolExecutor(max_workers=min(params.process,PROCESS_FILE)) as files_pool:

            if params.process >= PROCESS_FILE:

                # Check if the number of process available is odd or even
                pavailable = (params.process - PROCESS_FILE) % 2

                # Allocate to each file executor the number of process available for further parallelization by chromosomes
                palloc = [(params.process - PROCESS_FILE)//2]*2

                # If the number of process available is odd, allocate the remaining process to the file with the biggest size
                if pavailable != 0:
                    
                    fsizes = list(map(getsize, params.vcfs))

                    maxid = 0 if fsizes[0] > fsizes[1] else 1

                    palloc[maxid] += 1

                iterable: zip = zip(vcfs, palloc)
            
            else:

                iterable: tuple = [(vcfs[1], 0)]

            # Submit the process for each file mapped to the vcf path
            futures_to_vcf = {
                files_pool.submit(
                    process_vcf, vcf, params.report, proc
                ): vcf
                for vcf, proc in iterable
            }

            # Should the main process compute a VCF, contributing to the workload ?
            if PROCESS_FILE > params.process:   
                results[vcfs] = process_vcf(file=vcfs[0], compute=params.report, pool=0)

            # Check if a process is completed
            for future in concurrent.futures.as_completed(futures_to_vcf):

                # Get the returned results
                try:
                    (results[futures_to_vcf[future]]) = future.result()
                    logger.success(
                    f"Process {future} for file {futures_to_vcf[future]} has completed."
                    )
                except errors.ProcessError as e:
                    logger.error(
                        f"File {futures_to_vcf[future]} generated an exception: {e}"
                    )
                    if len(futures_to_vcf) == 2:
                        # Retrieve all futures
                        futures: list[concurrent.futures.Future] = list(futures_to_vcf.keys())
                        if future is futures[0]:
                            # Attempt to cancel the other call.
                            # If the call is currently being executed or finished running and cannot be cancelled then this attempt will fail.
                                futures[1].cancel()
                        else:
                            # Attempt to cancel the other call.
                            # If the call is currently being executed or finished running and cannot be cancelled then this attempt will fail.
                            futures[0].cancel()
                    
    # Computation is carried out sequentially.
    else:
        for vcf in vcfs:
            results[vcf] = process_vcf(file=vcf, compute=params.report, pool=0)

    results["delta"] = {
        "common": {},
        "unique": {vcfs[0]: {}, vcfs[1]: {}},
        "jaccard": {},
    }
    # Compute unique variants in first VCF file
    uniqueL: set = utils.difference(
        a=frozenset(results[vcfs[0]]["variants"].index),
        b=frozenset(results[vcfs[1]]["variants"].index),
    )
    # Compute unique variants in second VCF file
    uniqueR: set = utils.difference(
        a=frozenset(results[vcfs[1]]["variants"].index),
        b=frozenset(results[vcfs[0]]["variants"].index),
    )
    # Compute common variants between both VCF files
    common: set = utils.intersect(
        a=frozenset(results[vcfs[0]]["variants"].index),
        b=frozenset(results[vcfs[1]]["variants"].index),
    )

    (
        results["delta"]["common"],
        results["delta"]["unique"][vcfs[0]],
        results["delta"]["unique"][vcfs[1]],
    ) = (
        len(common),
        len(uniqueL),
        len(uniqueR),
    )

    logger.debug(
        f"{results['delta']['common']} variant(s) is/are commom in both files"
    )

    logger.debug(
        f"{results['delta']['unique'][vcfs[0]]} variant(s) is/are unique in files {vcfs[0]}"
    )

    logger.debug(
        f"{results['delta']['unique'][vcfs[1]]} variant(s) is/are unique in files {vcfs[1]}"
    )

    # Compute the Jaccard Index
    results["delta"]["jaccard"] = utils.jaccard_index(
        shared=results["delta"]["common"],
        total=(
            results["delta"]["common"]
            + results["delta"]["unique"][vcfs[0]]
            + results["delta"]["unique"][vcfs[1]]
        ),
    )

    logger.debug(f"Jaccard index: {results['delta']['jaccard']}")

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
                results[vcfs[0]]["variants"],
                results[vcfs[1]]["variants"],
            ],
            ["L", "R"],
        )
    )
    # Merge the two DataFrames with a outer join algorithm,
    # Merge is made on the hash index
    df: DataFrame = concat(
        [
            results[vcfs[0]]["variants"],
            results[vcfs[1]]["variants"],
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
    del results[vcfs[0]]["variants"]
    del results[vcfs[1]]["variants"]

    # Key is a references to the Hash string object in the Dataframe
    # Allow a O(1) lookup while keeping memory footprint low
    lookup = {hash: row for row, hash in enumerate(df["Hash"])}

    # Should a benchmark be computed ?
    if params.benchmark:

        logger.debug(f"Computing benchmark metrics from {vcfs[0]}.")

        table: DataFrame = utils.evaluate(df)

    # Should the output be serialized ?
    if params.serialize:

        logger.debug(f"Serializing output as {params.serialize}.")

        if params.process:
            # Parallelize the serialization process, with as much process as VCF files inputed
            # use a context manager to avoid memory leaks
            with concurrent.futures.ProcessPoolExecutor(max_workers=min(params.process, PROCESS_FILE)) as files_pool:
                # Submit the process for each file mapped to the vcf path
                futures_to_vcf = {
                    files_pool.submit(
                        utils.save, df, vcf, params.serialize, t, lookup
                    ): vcf for vcf, t in zip([vcfs[0].path, vcfs[1].path], ['L','R'])
                }
                # Check if a process is completed
                for future in concurrent.futures.as_completed(futures_to_vcf):
                    # Get the returned result
                    try:
                        out: int = future.result()
                        logger.success(f"Process {future} has successfully serialized {futures_to_vcf[future]}") if out else logger.error(f"Process {future} did not successfully serialized with a {out} output code.")
                    except ValueError as e:
                       logger.error(e)
        # Computation is carried out sequentially.
        else:
            for vcf, t in zip(vcfs, ['L','R']):
                utils.save(obj=df, path=vcf.path, format=params.serialize, target=t, lookup=lookup)

    # Should the output be reported ?
    if params.report:

        logger.debug("Generating a HTML report.")

        # Library of common plots between the two VCF files
        pcommon: PlotLibrary = PlotLibrary()

        # Create a Venn diagram to display the common, unique variants between the two VCF files
        pcommon.venn((results["delta"]["unique"][vcfs[0]], results["delta"]["unique"][vcfs[1]], results["delta"]["common"]), ['L','R'])

        # Create a report with the results
        files.Report(
            vcfs=vcfs,
            prefix=params.out,
            cmd=" ".join(argv),
            infos={
                vcfs[0]: results[vcfs[0]]["info"],
                vcfs[1]: results[vcfs[1]]["info"],
            },
            view={
                "headers": {
                    vcfs[0]: results[vcfs[0]]["header"],
                    vcfs[1]: results[vcfs[1]]["header"],
                },
                "variants": df,
                "stats": results["delta"]
            },
            plots={
                vcfs[0]: results[vcfs[0]]["plots"],
                vcfs[1]: results[vcfs[1]]["plots"],
                "common": pcommon
            },
            table=table if params.benchmark else None,
        ).create()
    # Print the results to the CLI
    else:
        print(f"{vcfs[0]}: [{results['delta']['unique'][vcfs[0]]} unique]────[{results['delta']['common']} common]────[{results['delta']['unique'][vcfs[1]]} unique] :{vcfs[1]}")
        print(f"Jaccard index: {results['delta']['jaccard']}")
        if params.benchmark:
            print(tabulate(table,headers='keys',tablefmt='grid',numalign='center', stralign='center'))

    return 1
