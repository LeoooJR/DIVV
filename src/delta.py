from console import stdout_console
import contextlib
import errors
import files
from itertools import repeat
from loguru import logger
from memory_profiler import profile
import os
from pathlib import Path
import processes
from rich.table import Table
from rich.errors import NotRenderableError
from sys import argv


def delta(params: object) -> int:
    """ 
    Main function to compute the delta between two VCF files 
            params: Namespace containing the command line parameters parsed.
    """

    assert len(params.vcfs) == 2, "Two VCF files are required."

    assert isinstance(params.vcfs[0], str) and isinstance(
        params.vcfs[1], str
    ), "Input VCF should be string instance"

    if not (os.path.isdir(params.output)):

        logger.error(f"No such directory: '{params.output}'")

        raise SystemExit(f"No such directory: '{params.output}'")
        
    else:

        if not os.access(params.output, os.W_OK):

            logger.error(f"Write permissions are not granted for the directory: {params.output}")

            raise SystemExit(f"Write permissions are not granted for the directory: {params.output}")
            
    params.output = Path(params.output)

    logger.debug(f"VCFS: {params.vcfs}")

    logger.debug(
        f"{params.vcfs[0]} is set as truth."
        if params.benchmark
        else "No VCF has been set as truth."
    )

    logger.debug(f"Indexes: {params.indexes}")

    logger.debug(f"Serialize output: {params.serialize}")

    logger.debug(f"Output a report: {params.report}")

    filters: dict[str:bool] = {"SNP": params.exclude_snps,
                               "TRANSITION": params.exclude_trans,
                               "MNP": params.exclude_mnps,
                               "INDELS": params.exclude_indels,
                               "SV": params.exclude_svs,
                               "VARS": params.exclude_vars,
                               "PASS_ONLY": params.pass_only}
    
    logger.debug(f"Filters used: {filters}")

    try:

        vcfs: files.VCFRepository = files.VCFRepository(vcfs=params.vcfs, index=params.indexes, reference=params.benchmark, filters=filters)
        
    except errors.VCFError as e: # Error raised by VCF file are critical

        logger.error(f"Error: {e}")

        raise SystemExit(e)

    for vcf in vcfs.repository:

        try:

            vcfs.processor.preprocessing(vcf, bins="env" if params.env_binaries else "project")

        except (errors.CompressionIndexError, errors.VCFError) as e:

            logger.error(e)

            raise SystemExit(e)

    manager: processes.TasksManager = processes.TasksManager(vcfs, params.process)

    manager.scheduling(tasks=[vcfs.repository[0].variants.chromosomes, vcfs.repository[1].variants.chromosomes])

    try:
        
        manager.commit(job=vcfs.processor.process_chromosome, jobargs=[True])

    except errors.ProcessError as e:

        logger.error(e)

        raise SystemExit(e)
    
    for vcf in vcfs.repository:

        for task in manager.tasks[vcf]:
    
            vcf.variants.update_repository(task[1], *manager.results[task])

    comparaisons: dict = vcfs.compare()

    # Should the output be serialized ?
    if params.serialize:

        logger.debug(f"Serializing output as {params.serialize} to {params.output}.")

        if params.serialize in ["vcf", "vcf.gz"]:

            manager.scheduling(tasks=[[params.output], [params.output]])

            try:

                manager.commit(job=vcfs.processor.serialize, jobargs=[comparaisons[(vcfs.repository[0],vcfs.repository[1])], params.serialize])

            except errors.ProcessError as e:

                logger.error(e)

                raise SystemExit(e)
            
        else:

            pass
            # files.VCFRepository.processor.serialize([None, os.getcwd()], comparaisons[(vcfs.repository[0],vcfs.repository[1])]["variants"], params.serialize)

    # Should the output be reported ?
    if params.report:

        logger.debug(f"Generating a HTML report to {params.output}.")

        vcfs.repository[0].variants.visualization()

        vcfs.repository[1].variants.visualization()

        tags = list(repeat(None, 2))

        if isinstance(params.tags, list):

            for i in range(len(params.tags[:2])):

                if isinstance(params.tags[i], str):

                    tags[i] = params.tags[i].split(',')

        # Create a report with the results
        files.Report(
            vcfs=vcfs,
            tags=tags,
            cmd=" ".join(argv),
            view=comparaisons[(vcfs.repository[0],vcfs.repository[1])],
            table=comparaisons[(vcfs.repository[0],vcfs.repository[1])]["benchmark"] if params.benchmark else None,
        ).create(output=params.output)

    # Print the results to the CLI
    else:
        stdout_console.rule("[bold sky_blue3] Results")
        stdout_console.print(f":vs: Comparaison: {vcfs.repository[0]} [{comparaisons[(vcfs.repository[0],vcfs.repository[1])]['unique'][vcfs.repository[0]]} unique]────[{comparaisons[(vcfs.repository[0],vcfs.repository[1])]['common']} common]────[{comparaisons[(vcfs.repository[0],vcfs.repository[1])]['unique'][vcfs.repository[1]]} unique] {vcfs.repository[1]}", style="result")
        stdout_console.print(f":heavy_large_circle: Jaccard index: {comparaisons[(vcfs.repository[0],vcfs.repository[1])]['jaccard']}", style="result")
        if params.benchmark:
            stdout_console.rule("[bold sky_blue3] Benchmark")
            stdout_console.print(f":bookmark: Reference: {vcfs.repository[0]}", style="info")
            stdout_console.print(f":dart: Query: {vcfs.repository[1]}", style="info")            
            table = Table(title="Benchmark table", style="bold")
            df = comparaisons[(vcfs.repository[0],vcfs.repository[1])]["benchmark"].astype(str)
            for col in df.columns:
                table.add_column(col, justify="center", vertical="middle")
            for row in df.values:
                with contextlib.suppress(NotRenderableError):
                    table.add_row(*row)
            stdout_console.print(table)