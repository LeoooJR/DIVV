from __init__ import __version__
import argparse
from delta import delta
from psutil import cpu_count
import sys
import utils


class Programm:

    FUNC = {"delta": delta}

    def __init__(self):

        self.parser = argparse.ArgumentParser(
            prog="VCFDelta",
            description="Program for comparing Variant Call Format files.",
        )
        self.parser.add_argument(
            "-v",
            "--version",
            action="version",
            version="VCFDelta v" + __version__,
        )
        self.parser.add_argument(
            "--vcfs",
            dest="vcfs",
            nargs=2,
            type=str,
            metavar="",
            required=True,
            help="Paths to the inputs VCF files",
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
            help="prefix of the outputs files",
        )
        self.parser.add_argument(
            "-s",
            "--serialize",
            dest="serialize",
            type=str,
            choices=["json", "pickle"],
            required=False,
            help="Should the result be seralized",
        )
        self.parser.add_argument(
            "-p",
            "--process",
            dest="process",
            type=int,
            required=False,
            default=0,
            help="Number of process to use in addition to the current one",
        )
        self.parser.add_argument(
            "--threshold",
            dest="threshold",
            type=float,
            default=0.01,
            metavar="",
            required=False,
            help="Threshold at which differences in ARR and other raise error",
        )
        self.parser.add_argument(
            "--skip_perl_error_pvalue",
            action="store_true",
            required=False,
            help="Skip pvalue comparison if perl script return an error score (-1). Set to 1 to activate.",
        )
        self.parser.add_argument(
            "--specific_threshold_pvalue",
            type=float,
            default=0.0,
            metavar="",
            required=False,
            help="Take a float. Instead of checking if difference of pvalue is lower than the value given"
            " by option --threshold. It check if both pvalue are lower or higher"
            " than the value given here.",
        )
        self.parser.add_argument(
            "--soft_check_bkg",
            action="store_true",
            required=False,
            help="Instead of comparing exact BKG value, check only if it is noisy or clean and skip the"
            "the likely or probably status.",
        )
        self.parser.add_argument(
            "--exclude_snps",
            dest="exclude_snps",
            action="store_true",
            default=False,
            help="Exclude snp calls. A heterozygous call with both snp and indel is not excluded unless both snps and indels are excluded.",
        )
        self.parser.add_argument(
            "--exclude_indels",
            dest="exclude_indels",
            action="store_true",
            default=False,
            help="Exclude insertions and deletions. A heterozygous call with both snp and indel is not excluded unless both snps and indels are excluded.",
        )
        self.parser.add_argument(
            "--exclude_vars",
            dest="exclude_vars",
            action="store_true",
            default=False,
            help="Exclude variants other than snps and indels.",
        )
        self.parser.add_argument(
            "--exclude_refs",
            dest="exclude_refs",
            action="store_true",
            default=False,
            help="Exclude reference calls.",
        )
        self.parser.add_argument(
            "--exclude_hetero",
            dest="exclude_hetero",
            action="store_true",
            default=False,
            help="Exclude heterozygous calls.",
        )
        self.parser.add_argument(
            "--exclude_filtered",
            dest="exclude_filtered",
            action="store_true",
            default=False,
            help="Exclude filtered calls (FT or FILTER is not PASS).",
        )
        self.parser.add_argument(
            "--exclude_missing",
            dest="exclude_missing",
            action="store_true",
            default=False,
            help="Exclude calls with all data elements missing.",
        )
        self.parser.add_argument(
            "--truth",
            dest="truth_flag",
            action="store_true",
            default=False,
            help="Additional metrics are generated assuming the first VCF file is the truth. This also triggers extra analysis of filtered calls.",
        )
        self.parser.add_argument(
            "--stats",
            dest="stats",
            action="store_true",
            help="Should statistiques be computed.",
            default=False,
        )
        self.parser.add_argument(
            "-r",
            "--report",
            dest="report",
            action="store_true",
            help="Should a report be generated.",
            default=False,
        )

        self.parser.set_defaults(func=self.FUNC["delta"])

    def launch(self) -> int:
        cmd = self.parser.parse_args(args=sys.argv[1:])
        return cmd.func(params=cmd)

    def __str__(self):
        pass
