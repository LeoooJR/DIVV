from console import stdout_console
import contextlib
import enum
import errors
import files
from itertools import repeat
from loguru import logger
import os
from pathlib import Path
import processes
from rich.table import Table
from rich.errors import NotRenderableError
from rich.panel import Panel
from sys import argv

class RunningMode(enum.Enum):
    """Running mode enum."""
    # Production mode
    production = (1, ('prod',), ['false', '0', 'no', ''])
    # Development mode
    development = (0, ('dev',), ['true', '1', 'yes'])

    def __init__(self, num: int, aliases: tuple[str], env_values: list[str]):
        """Initialize the RunningMode object."""
        # Id of the running mode
        self.num: int = num
        # Aliases of the running mode
        self.aliases: tuple[str] = aliases
        # Environment values of the running mode
        self.env_values: list[str] = env_values

    @classmethod
    def from_env(cls) -> "RunningMode":
        """Determine the running mode from the DIVV_DEV environment variable."""
        # Get the value of the DIVV_DEV environment variable
        value = os.getenv("DIVV_DEV", '').lower()
        # If the env value matches any of the development values, return development
        if value in cls.development.env_values:
            return cls.development
        # Otherwise, treat as production (default)
        return cls.production

    def is_production(self) -> bool:
        """Check if the running mode is production."""
        return self is RunningMode.production


def supervisor(params: object) -> int:
    """ 
    Main function to compute the delta between two VCF files 
            params: Namespace containing the command line parameters parsed.
    """

    running_mode = RunningMode.from_env()
    in_production: bool = running_mode.is_production()

    if in_production:
        logger.debug("Running in production mode.")
    else:
        logger.debug("Running in development mode.")

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

    # Filters used to filter the VCF files
    filters: dict[str:bool] = {"SNP": params.exclude_snps,
                               "TRANSITION": params.exclude_trans,
                               "MNP": params.exclude_mnps,
                               "INDELS": params.exclude_indels,
                               "SV": params.exclude_svs,
                               "VARS": params.exclude_vars,
                               "PASS_ONLY": params.pass_only}
    
    logger.debug(f"Filters used: {filters}")

    # If a report is requested, check output path
    if params.report:

        logger.debug(f"Output a report: {params.report}")

        # If a archive | file (prod) is requested, the output path must link to a single file
        if params.archive or in_production:

            if params.archive:
                logger.debug("Generated report will be archived as ZIP file.")
            else:
                logger.debug("Generated report is a single HTML file.")

            # If the output path is a directory, raise an error
            if os.path.isdir(params.output):

                logger.error(f"Output path '{params.output}' is a directory, expected a file {'with --archive option' if params.archive else ''}.")

                raise SystemExit(f"Output path '{params.output}' is a directory, expected a file {'with --archive option' if params.archive else ''}.")

            # If the output path is a file, check the parent directory
            else:

                parent_dir = os.path.dirname(params.output)

                # If there is no parent directory, it is a relative path in the current directory
                if not parent_dir:

                    parent_dir = os.getcwd()

                # If the parent directory does not exist, raise an error
                if not os.path.isdir(parent_dir):

                    logger.error(f"No such parent directory: '{parent_dir}'")

                    raise SystemExit(f"No such parent directory: '{parent_dir}'")
                
                else:

                    # If the parent directory is not writable, raise an error
                    if not os.access(parent_dir, os.W_OK):

                        logger.error(f"Write permissions are not granted for the parent directory: {parent_dir}")

                        raise SystemExit(f"Write permissions are not granted for the parent directory: {parent_dir}")
        
        # A plain directory (dev) is expected 
        else:

            # If the output path is a directory, check if it exists
            if not (os.path.isdir(params.output)):
                # Log the error
                logger.error(f"No such output directory: '{params.output}'")
                # Raise an error
                raise SystemExit(f"No such output directory: '{params.output}'")
                    
            else:

                # If the output directory is not writable, raise an error
                if not os.access(params.output, os.W_OK):
                    # Log the error
                    logger.error(f"Write permissions are not granted for the output directory: {params.output}")
                    # Raise an error
                    raise SystemExit(f"Write permissions are not granted for the output directory: {params.output}")
            
    # Convert the output path to a Path object allowing easier manipulation
    params.output = Path(params.output)

    # Trace
    logger.debug(f"Serialize output: {params.serialize}")

    # Try to create a VCFRepository object from the input VCF files and indexes
    try:

        # Create a VCFRepository object from the input VCF files and indexes
        vcfs: files.VCFRepository = files.VCFRepository(vcfs=params.vcfs, index=params.indexes, reference=params.benchmark, filters=filters)
        
    except errors.VCFError as e: # Error raised by VCF file are critical

        logger.error(e)

        raise SystemExit(e)

    # Preprocess the VCF files
    for vcf in vcfs.repository:

        # Try to preprocess the VCF file
        try:

            vcfs.processor.preprocessing(vcf, bins="env" if params.env_binaries else "project")

        # If an error occurs while preprocessing the VCF file, raise an error
        except (errors.CompressionIndexError, errors.VCFError) as e:

            logger.error(e)

            raise SystemExit(e)

    # Create a TasksManager object to manage the tasks
    manager: processes.TasksManager = processes.TasksManager(vcfs, params.process)

    # Schedule the chromosomes processing tasks
    manager.scheduling(tasks=[vcf.variants.chromosomes for vcf in vcfs.repository])

    # Commit the tasks to the TasksManager
    try:
        
        manager.commit(job=vcfs.processor.process_chromosome, jobargs=[True])

    except errors.ProcessError as e:

        logger.error(e)

        raise SystemExit(e)
    
    # Update the VCFRepository with the results
    for vcf in vcfs.repository:

        for task in manager.tasks[vcf]:
    
            vcf.variants.update_repository(task[1], *manager.results[task])

    # Compare the VCF files
    comparaisons: dict = vcfs.compare()

    # If a serialization is requested, serialize the output
    if params.serialize:
        # Trace
        logger.debug(f"Serializing output as {params.serialize} to {params.output}.")
        # If the serialization format is VCF or VCF.GZ
        if params.serialize in ["vcf", "vcf.gz"]:
            # Schedule the tasks
            manager.scheduling(tasks=[[params.output], [params.output]])
            # Try to serialize the output
            try:
                # Commit the tasks
                manager.commit(job=vcfs.processor.serialize, jobargs=[comparaisons[(vcfs.repository[0],vcfs.repository[1])], params.serialize])
            # If an error occurs while serializing the output, raise an error
            except errors.ProcessError as e:
                # Log the error
                logger.error(e)
                # Raise an error
                raise SystemExit(e)
            
        else:

            pass
            # files.VCFRepository.processor.serialize([None, os.getcwd()], comparaisons[(vcfs.repository[0],vcfs.repository[1])]["variants"], params.serialize)

    # If a report is requested, generate the report
    if params.report:
        # Trace
        logger.debug(f"Generating a HTML report to {params.output}.")
        # Create figures for the first VCF
        vcfs.repository[0].variants.visualization()
        # Create figures for the second VCF
        vcfs.repository[1].variants.visualization()

        # Initialize a list of tags
        tags = list(repeat(None, 2))

        # If tags are provided, parse them
        if isinstance(params.tags, list):

            for i in range(len(params.tags[:2])):

                if isinstance(params.tags[i], str):

                    tags[i] = params.tags[i].split(',')

        try:
            # Create a report with the results
            dest: Path = files.Report(
                vcfs=vcfs,
                tags=tags,
                cmd=" ".join(argv),
                view=comparaisons[(vcfs.repository[0],vcfs.repository[1])],
                table=comparaisons[(vcfs.repository[0],vcfs.repository[1])]["benchmark"] if params.benchmark else None,
                archive=params.archive
            ).create(output=params.output, bundle=in_production)
            # Print a success message
            stdout_console.print(Panel.fit(f"Report successfully generated at '{dest}'.", title="Success", highlight=True), style="result")
        # If an error occurs while generating the report, raise an error
        except errors.ReportError as e:
            raise SystemExit(e)

    # Print the results to the CLI
    else:

        # If filters have been applied, print the summary of the filtering steps
        if vcfs.repository[0].variants.is_filtered() or vcfs.repository[1].variants.is_filtered():
            # Divider
            stdout_console.rule("[bold sky_blue3] Filters")
            # Print the summary of the filtering steps for the first VCF
            stdout_console.print(f":magnifying_glass_tilted_right: {vcfs.repository[0]}({vcfs.repository[0].variants.filtered.total()}) {dict(vcfs.repository[0].variants.filtered)}", style="info")
            # Print the summary of the filtering steps for the second VCF
            stdout_console.print(f":magnifying_glass_tilted_right: {vcfs.repository[1]}({vcfs.repository[1].variants.filtered.total()}) {dict(vcfs.repository[1].variants.filtered)}", style="info")
        
        # Print the results of the comparison
        # Divider
        stdout_console.rule("[bold sky_blue3] Results")
        # Print the comparison results
        stdout_console.print(f":vs: Comparaison: {vcfs.repository[0]} [{comparaisons[(vcfs.repository[0],vcfs.repository[1])]['unique'][vcfs.repository[0]]} unique]────[{comparaisons[(vcfs.repository[0],vcfs.repository[1])]['common']} common]────[{comparaisons[(vcfs.repository[0],vcfs.repository[1])]['unique'][vcfs.repository[1]]} unique] {vcfs.repository[1]}", style="result")
        # Print the Jaccard index
        stdout_console.print(f":heavy_large_circle: Jaccard index: {comparaisons[(vcfs.repository[0],vcfs.repository[1])]['jaccard']}", style="result")
        
        # If a benchmark is requested, print the benchmark results
        if params.benchmark:
            # Divider
            stdout_console.rule("[bold sky_blue3] Benchmark")
            # Print the reference VCF
            stdout_console.print(f":bookmark: Reference: {vcfs.repository[0]}", style="info")
            # Print the query VCF
            stdout_console.print(f":dart: Query: {vcfs.repository[1]}", style="info")            
            # Create a table to print the benchmark results
            table = Table(title="Benchmark table", style="bold")
            df = comparaisons[(vcfs.repository[0],vcfs.repository[1])]["benchmark"].astype(str)
            # Add the columns to the table
            for col in df.columns:
                table.add_column(col, justify="center", vertical="middle")
            # Add the rows to the table
            for row in df.values:
                with contextlib.suppress(NotRenderableError):
                    table.add_row(*row)
            stdout_console.print(table)