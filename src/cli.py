from __init__ import __version__
import argparse
from delta import delta
from loguru import logger
from sys import argv

class Program:

    FUNC = {"delta": delta}

    def __init__(self):

        self.parser = argparse.ArgumentParser(
            prog="VCFDelta",
            description="Compare VCF files.",
        )
        self.parser.add_argument(
            "-v",
            "--version",
            action="version",
            version=f"VCFDelta v{__version__}",
        )
        self.parser.add_argument(
            "--vcfs",
            dest="vcfs",
            nargs=2,
            type=str,
            metavar="",
            required=True,
            help="Paths to the inputs VCF files [.vcf[.gz]]",
        )
        self.parser.add_argument(
            "-i",
            "--indexes",
            dest="indexes",
            nargs=2,
            type=str,
            required=False,
            default=[None, None],
            help="Paths to the index of the VCF files",
        )
        self.parser.add_argument(
            "-o",
            "--output",
            dest="out",
            type=str,
            metavar="",
            required=False,
            default="delta",
            help="Prefix of the outputs files",
        )
        self.parser.add_argument(
            "-s",
            "--serialize",
            dest="serialize",
            type=str,
            choices=["json", "pickle", "vcf", "vcf.gz"],
            required=False,
            help="Should the result be seralized [.json, .pickle, .vcf[.gz]]",
        )
        self.parser.add_argument(
            "-p",
            "--process",
            dest="process",
            type=int,
            required=False,
            default=0,
            help="Number of processes to be used in addition to the main process.",
        )
        self.parser.add_argument(
            "-e",
            "--env-binaries",
            dest="env_binaries",
            required=False,
            default=False,
            action="store_true",
            help="Binaries to be used for processes such as compressing and indexing VCF. If this option flag is specified, the program will attempt to call binaries such as bgzip and tabix from the user's local environment."
        )
        self.parser.add_argument(
            "--exclude-snps",
            dest="exclude_snps",
            action="store_true",
            default=False,
            help="Exclude single nucleotide polymorphisme calls. A heterozygous call with both snp and indel is not excluded unless both snps and indels are excluded.",
        )
        self.parser.add_argument(
            "--exclude-mnps",
            dest="exclude_mnps",
            action="store_true",
            default=False,
            help="Exclude mutliple nucleotide polymorphisme calls.",
        )
        self.parser.add_argument(
            "--exclude-indels",
            dest="exclude_indels",
            action="store_true",
            default=False,
            help="Exclude insertions and deletions. A heterozygous call with both snp and indel is not excluded unless both snps and indels are excluded.",
        )
        self.parser.add_argument(
            "--exclude-vars",
            dest="exclude_vars",
            action="store_true",
            default=False,
            help="Exclude variants other than snps and indels.",
        )
        self.parser.add_argument(
            "--exclude-svs",
            dest="exclude_svs",
            action="store_true",
            default=False,
            help="Exclude structural variant calls.",
        )
        self.parser.add_argument(
            "--exclude-transitions",
            dest="exclude_trans",
            action="store_true",
            default=False,
            help="Exclude transition calls.",
        )
        self.parser.add_argument(
            "--exclude-refs",
            dest="exclude_refs",
            action="store_true",
            default=False,
            help="Exclude reference calls.",
        )
        self.parser.add_argument(
            "--exclude-hetero",
            dest="exclude_hetero",
            action="store_true",
            default=False,
            help="Exclude heterozygous calls.",
        )
        self.parser.add_argument(
            "--exclude-filtered",
            dest="exclude_filtered",
            action="store_true",
            default=False,
            help="Exclude filtered calls (FILTER value is not PASS).",
        )
        self.parser.add_argument(
            "--exclude-missing",
            dest="exclude_missing",
            action="store_true",
            default=False,
            help="Exclude calls with all data elements missing.",
        )
        self.parser.add_argument(
            "--pass-only",
            dest="pass_only",
            action="store_true",
            default=False,
            help="Keep only PASS calls.",
        )
        self.parser.add_argument(
            "--benchmark",
            dest="benchmark",
            action="store_true",
            default=False,
            help="Additional metrics are generated assuming the first VCF file is the truth.",
        )
        self.parser.add_argument(
            "-r",
            "--report",
            dest="report",
            action="store_true",
            help="Should a HTML report be generated.",
            default=False,
        )
        self.parser.add_argument(
            "-t",
            "--tags",
            dest="tags",
            type=str,
            nargs=2,
            required=False,
            help="For the final report, add more information about VCFs with tags. Values must be separated by commas, e.g. [tag,tag,...]."
        )
        self.parser.add_argument(
            "--verbosity",
            dest="verbosity",
            action="store_true",
            help="Should logs be printed to the shell.",
            default=False,
        )

        self.parser.set_defaults(func=self.FUNC["delta"])

    def launch(self) -> int:

        cmd = self.parser.parse_args(args=argv[1:])

        # Should the log be printed to CLI or saved in a file ?
        if not cmd.verbosity:

            logger.remove(0)

            logger.add("VCFDelta.log")

        return cmd.func(params=cmd)

    def __str__(self):

        return "VCFDelta"

    def __repr__(self):
        
        return "VCFDelta"
