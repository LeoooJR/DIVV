import concurrent.futures
import errors
import files
from loguru import logger
from memory_profiler import profile
from operator import itemgetter
from pandas import DataFrame, concat
from plots import PlotLibrary
import processes
from sys import argv
from tabulate import tabulate
import utils


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

    filters: dict[str:bool] = {"SNP": params.exclude_snps,
                               "TRANSITION": params.exclude_trans,
                               "MNP": params.exclude_mnps,
                               "INDELS": params.exclude_indels,
                               "SV": params.exclude_svs,
                               "VARS": params.exclude_vars,
                               "PASS_ONLY": params.pass_only}
    
    logger.debug(f"Filters used: {filters}")

    try:

        vcfs: list[files.VCF] = [files.VCF(path=params.vcfs[0], reference=True if params.benchmark else False, index=params.indexes[0], filters=filters, lazy=False), 
                                files.VCF(path=params.vcfs[1], index=params.indexes[1], filters=filters, lazy=False)]
        
    except errors.VCFError as e: # Error raised by VCF file are critical

        logger.error(f"Error: {e}")

        raise SystemExit()
    
    processor: processes.VCFProcessor = processes.VCFProcessor()

    for vcf in vcfs:
        try:
            processor.preprocessing(vcf)
        except (errors.CompressionIndexError, errors.VCFError) as e:
            logger.error(e)
            raise SystemExit()

    manager: processes.TasksManager = processes.TasksManager(vcfs, params.process)

    manager.scheduling(tasks=[vcfs[0].chromosomes, vcfs[1].chromosomes])

    try:

        manager.commit(job=processor.process_chromosome, jobargs=[True])

    except errors.ProcessError as e:

        logger.error(e)

        raise SystemExit()
    
    for vcf in vcfs:

        for task in manager.tasks[vcf]:
    
            vcf.variants.update_repository(task[1], *manager.results[task])

    results: dict = {vcfs[0]: {"variants": vcfs[0].variants.collapse()},
                     vcfs[1]: {"variants": vcfs[1].variants.collapse()}}

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

        vcfs[0].variants.visualization()

        vcfs[1].variants.visualization()

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
                vcfs[0]: vcfs[0].informations(),
                vcfs[1]: vcfs[1].informations(),
            },
            view={
                "headers": {
                    vcfs[0]: "\t".join(list(vcfs[0].header.keys())),
                    vcfs[1]: "\t".join(list(vcfs[1].header.keys())),
                },
                "variants": df,
                "stats": results["delta"]
            },
            plots={
                vcfs[0]: vcfs[0].variants.plots,
                vcfs[1]: vcfs[1].variants.plots,
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
