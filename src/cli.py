from __init__ import __version__
import argparse
from supervisor import supervisor
from loguru import logger
import os
from sys import argv
import validation

class EntryPoint:
    """
    A class to parse the command line arguments and launch the program.
    """
    # Dictionary of functions to call
    FUNC = {"call": supervisor}

    def __init__(self):

        # Create the parser
        self.parser = argparse.ArgumentParser(
            prog="DIVV",
            description="Compare VCF files.",
        )
        # Add the version argument
        self.parser.add_argument(
            "-v",
            "--version",
            action="version",
            version=f"DIVV v{__version__}",
        )
        # Add the vcfs argument
        self.parser.add_argument(
            "--vcfs",
            dest="vcfs",
            metavar="VCF",
            nargs=2,
            type=str,
            required=True,
            help="Paths to the inputs VCF files [.vcf[.gz]]",
        )
        # Add the indexes argument
        self.parser.add_argument(
            "-i",
            "--indexes",
            dest="indexes",
            metavar="INDEX",
            nargs=2,
            type=str,
            required=False,
            default=[None, None],
            help="Paths to the index of the VCF files",
        )
        # Add the output argument
        self.parser.add_argument(
            "-o",
            "--output",
            dest="output",
            metavar="PATH",
            type=str,
            required=False,
            default=os.getcwd(),
            help="Path to the output",
        )
        # Add the serialize argument
        self.parser.add_argument(
            "-s",
            "--serialize",
            dest="serialize",
            type=str,
            choices=["json", "pickle", "vcf", "vcf.gz"],
            required=False,
            help="Should the result be seralized [.json, .pickle, .vcf[.gz]]?",
        )
        # Add the process argument
        self.parser.add_argument(
            "-p",
            "--process",
            dest="process",
            metavar="INT",
            type=int,
            required=False,
            default=0,
            action=validation.ValidateProcessAction,
            help="Number of processes to be used in addition to the main process.",
        )
        # Add the env_binaries argument
        self.parser.add_argument(
            "-e",
            "--env-binaries",
            dest="env_binaries",
            required=False,
            default=False,
            action="store_true",
            help="Binaries to be used for processes such as compressing and indexing VCF. If this option flag is specified, the program will attempt to call binaries such as bgzip and tabix from the user's local environment."
        )
        # Add the exclude_snps argument
        self.parser.add_argument(
            "--exclude-snps",
            dest="exclude_snps",
            action="store_true",
            default=False,
            help="Exclude single nucleotide polymorphisme calls. A heterozygous call with both snp and indel is not excluded unless both snps and indels are excluded.",
        )
        # Add the exclude_mnps argument
        self.parser.add_argument(
            "--exclude-mnps",
            dest="exclude_mnps",
            action="store_true",
            default=False,
            help="Exclude mutliple nucleotide polymorphisme calls.",
        )
        # Add the exclude_indels argument
        self.parser.add_argument(
            "--exclude-indels",
            dest="exclude_indels",
            action="store_true",
            default=False,
            help="Exclude insertions and deletions. A heterozygous call with both snp and indel is not excluded unless both snps and indels are excluded.",
        )
        # Add the exclude_vars argument
        self.parser.add_argument(
            "--exclude-vars",
            dest="exclude_vars",
            action="store_true",
            default=False,
            help="Exclude variants other than snps and indels.",
        )
        # Add the exclude_svs argument
        self.parser.add_argument(
            "--exclude-svs",
            dest="exclude_svs",
            action="store_true",
            default=False,
            help="Exclude structural variant calls.",
        )
        # Add the exclude_transitions argument
        self.parser.add_argument(
            "--exclude-transitions",
            dest="exclude_trans",
            action="store_true",
            default=False,
            help="Exclude transition calls.",
        )
        # Add the exclude_refs argument
        self.parser.add_argument(
            "--exclude-refs",
            dest="exclude_refs",
            action="store_true",
            default=False,
            help="Exclude reference calls.",
        )
        # Add the exclude_hetero argument
        self.parser.add_argument(
            "--exclude-hetero",
            dest="exclude_hetero",
            action="store_true",
            default=False,
            help="Exclude heterozygous calls.",
        )
        # Add the exclude_filtered argument
        self.parser.add_argument(
            "--exclude-filtered",
            dest="exclude_filtered",
            action="store_true",
            default=False,
            help="Exclude filtered calls (FILTER value is not PASS).",
        )
        # Add the exclude_missing argument
        self.parser.add_argument(
            "--exclude-missing",
            dest="exclude_missing",
            action="store_true",
            default=False,
            help="Exclude calls with all data elements missing.",
        )
        # Add the pass_only argument
        self.parser.add_argument(
            "--pass-only",
            dest="pass_only",
            action="store_true",
            default=False,
            help="Keep only PASS calls.",
        )
        # Add the benchmark argument
        self.parser.add_argument(
            "--benchmark",
            dest="benchmark",
            action="store_true",
            default=False,
            help="Additional metrics are generated assuming the first VCF file is the truth.",
        )
        # Add the report argument
        self.parser.add_argument(
            "-r",
            "--report",
            dest="report",
            action="store_true",
            help="Should a HTML report be generated?",
            default=False,
        )
        # Add the archive argument
        self.parser.add_argument(
            "-a",
            "--archive",
            dest="archive",
            action="store_true",
            help="Should the report be archived?",
            default=False
        )
        # Add the tags argument
        self.parser.add_argument(
            "-t",
            "--tags",
            dest="tags",
            type=str,
            nargs='+',
            required=False,
            help="For the final report, add more information about VCFs with tags. Values must be separated by commas, e.g. [tag,tag,...]."
        )
        # Add the debug argument
        self.parser.add_argument(
            "-d",
            "--debug",
            dest="debug",
            action="store_true",
            help="Should logs be saved?",
            default=False,
        )

        # Set the default function to call
        self.parser.set_defaults(func=self.FUNC["call"])

    # Launch the program with command line arguments
    def launch(self) -> int:
        """Launch the program with command line arguments.
        
        Returns:
            int: Exit code of the program.
        """        
        
        # Parse the command line arguments
        cmd = self.parser.parse_args(args=argv[1:])

        # Remove the default handler
        logger.remove(0)

        # Should the log be saved in a file ?
        if cmd.debug:

            logger.add("DIVV.log")

        # Call the function with the command line arguments
        return cmd.func(params=cmd)

    def __str__(self):
        """String representation of the program.
        
        Returns:
            str: The name of the program.
        """

        return "DIVV"

    def __repr__(self):
        """String representation of the program for debugging.
        
        Returns:
            str: The name of the program.
        """        
        
        return "DIVV"
