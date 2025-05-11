import concurrent.futures
import errors
import files
from loguru import logger
from memory_profiler import profile
import processes
from sys import argv
from tabulate import tabulate
import utils


def delta(params: object) -> int:
    """ 
    Main function to compute the delta between two VCF files 
            params: object containing the command line parameters parsed.
    """

    assert len(params.vcfs) == 2, "Two VCF files are required."

    assert isinstance(params.vcfs[0], str) and isinstance(
        params.vcfs[1], str
    ), "Input VCF should be string instance"

    logger.debug(f"VCFS: {params.vcfs}")

    logger.debug(
        f"{params.vcfs[0]} is set as truth."
        if params.benchmark
        else "No VCF has been set as truth."
    )

    logger.debug(f"Indexes: {params.indexes}")

    logger.debug(f"Serialize output: {params.serialize}")

    logger.debug(f"Output a report: {params.report}")

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

        vcfs: files.VCFRepository = files.VCFRepository(vcfs=params.vcfs, index=params.indexes, reference=params.benchmark, filters=filters)
        
    except errors.VCFError as e: # Error raised by VCF file are critical

        logger.error(f"Error: {e}")

        raise SystemExit(e)

    for vcf in vcfs.repository:

        try:

            vcfs.processor.preprocessing(vcf)

        except (errors.CompressionIndexError, errors.VCFError) as e:

            logger.error(e)

            raise SystemExit(e)

    manager: processes.TasksManager = processes.TasksManager(vcfs, params.process)

    manager.scheduling(tasks=[vcfs.repository[0].chromosomes, vcfs.repository[1].chromosomes])

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

        logger.debug(f"Serializing output as {params.serialize}.")

        if params.process:
            # Parallelize the serialization process, with as much process as VCF files inputed
            # use a context manager to avoid memory leaks
            with concurrent.futures.ProcessPoolExecutor(max_workers=min(params.process, PROCESS_FILE)) as files_pool:
                # Submit the process for each file mapped to the vcf path
                futures_to_vcf = {
                    files_pool.submit(
                        utils.save, comparaisons[(vcfs.repository[0],vcfs.repository[1])]["variants"], vcf, params.serialize, comparaisons[(vcfs.repository[0],vcfs.repository[1])]["index"]
                    ): vcf for vcf in [vcfs.repository[0].path, vcfs.repository[1].path]
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
            for vcf, t in zip(vcfs.repository, ['L','R']):
                utils.save(obj=comparaisons[(vcfs.repository[0],vcfs.repository[1])]["variants"], path=vcf.path, format=params.serialize, target=t, lookup=comparaisons[(vcfs.repository[0],vcfs.repository[1])]["index"])

    # Should the output be reported ?
    if params.report:

        logger.debug("Generating a HTML report.")

        vcfs.repository[0].variants.visualization()

        vcfs.repository[1].variants.visualization()

        # Create a report with the results
        files.Report(
            vcfs=vcfs,
            prefix=params.out,
            cmd=" ".join(argv),
            view=comparaisons[(vcfs.repository[0],vcfs.repository[1])],
            table=comparaisons[(vcfs.repository[0],vcfs.repository[1])]["benchmark"] if params.benchmark else None,
        ).create()

    # Print the results to the CLI
    else:
        print(f"{vcfs.repository[0]}: [{comparaisons[(vcfs.repository[0],vcfs.repository[1])]['unique'][vcfs.repository[0]]} unique]────[{comparaisons[(vcfs.repository[0],vcfs.repository[1])]['common']} common]────[{comparaisons[(vcfs.repository[0],vcfs.repository[1])]['unique'][vcfs.repository[1]]} unique] :{vcfs.repository[1]}")
        print(f"Jaccard index: {comparaisons[(vcfs.repository[0],vcfs.repository[1])]['jaccard']}")
        if params.benchmark:
            print(tabulate(comparaisons[(vcfs.repository[0],vcfs.repository[1])]["benchmark"],headers='keys',tablefmt='grid',numalign='center', stralign='center'))