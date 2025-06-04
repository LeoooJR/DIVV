import copy as cp
import datetime
import json
import cyvcf2
import errors
import filetype
from glob import glob
import gzip
from hashlib import sha256, md5
from itertools import pairwise
from jinja2 import Environment, FileSystemLoader
from loguru import logger
import numpy as np
import os
from pandas import isna, concat, Index, DataFrame, Series, notna
import pathlib
from plots import PlotLibrary
from shutil import copy, copytree
import stat
import subprocess
import utils
from variants import VariantRepository
import webbrowser
import zipfile
import _pickle as cPickle
from __init__ import __version__

class GenomicFile():

    """Base class for genomic files"""

    def __init__(self, path: str):

        # Path to the file
        self.path: pathlib.Path = pathlib.Path(path)

        self._hash = None

    def is_empty(self) -> bool:
        """Check if file is empty"""

        return os.path.getsize(self.path) == 0

    def is_file(self) -> bool:
        """Check if path is a file"""

        return os.path.isfile(self.path)

    def get_path(self) -> pathlib.Path:

        return self.path
    
    def basename(self) -> str:

        return os.path.basename(self.path)
    
    def informations(self) -> dict:

        return utils.file_infos(self.path)
    
    def __str__(self):
        
        return str(self.path)
    
    def __repr__(self):
        
        return repr(self.path)
    
    def __hash__(self):

        if self._hash is None:

            megabyte: int = 1_048_576

            buffer: int = 10 * megabyte

            if self.path.exists():

                hash_md5: object  = md5()

                with open(self.path, "rb") as f:

                    for chunk in iter(lambda: f.read(buffer), b""):

                        hash_md5.update(chunk)

                self._hash = int(hash_md5.hexdigest(), 16)

            else:

                raise errors.FileError(f"File {self.path} cannot be hashed.")

        return self._hash
        
    def __eq__(self, value):

        if isinstance(value, GenomicFile):

            if self.__hash__():
            
                return self.__hash__() == hash(value)

        return False

class VCF(GenomicFile):

    """Class for VCF files"""

    # Map the FORMAT values to there respectives informations
    FORMAT: dict = {
        "genotype": ["GT"],
        "genotype_quality": ["GQ"],
        "depth": ["DP", "TRC"],
    }

    def __init__(self, path: str, reference: bool = False, index: str = None, filters: dict[str:bool] = None,  lazy: bool = True):
        
        super().__init__(path)

        if not lazy:

            self.verify()

        self.reference: bool = reference

        self.archive: str = None

        self.stdin = None

        # Map the header values to there respectives indexes
        self.HEADER = {
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

        self.SAMPLES = None

        self.variants: VariantRepository = VariantRepository(filters)

        if self.archive and index:

            try:

                self._index: VCFIndex = VCFIndex(path=index, vcf=self, lazy=False)

            # Errors caused by the provided index are treated as warnings.
            except errors.IndexError:

                logger.warning(f"{index} is not a valid VCF index file.")

                self._index: VCFIndex = None

        else:

            self._index: VCFIndex = None

    @property
    def index(self):

        return getattr(self, "_index", None)
    
    @index.setter
    def index(self, value: str):

        if isinstance(value, str):

            try:

                self._index: VCFIndex = VCFIndex(path=value, vcf=self, lazy=False)

            except errors.IndexError:

                self._index = None

        elif isinstance(value, VCFIndex):

            self._index: VCFIndex = value

    def is_indexed(self) -> bool:

        return isinstance(self._index, VCFIndex)
    
    def is_compressed(path: str, type: str) -> bool:
        """ Checks if a file is a compressed archive type """
        return filetype.archive_match(path) == type
    
    @property
    def header(self):

        return getattr(self, "HEADER", None)
    
    @property
    def format(self):

        return getattr(self, "FORMAT", None)
    
    @property
    def samples(self):

        return getattr(self, "SAMPLES", None)
    
    @samples.setter
    def samples(self, values: list[str]):

        self.SAMPLES = values

        self.update_header(samples=values)

    def is_reference(self):

        return self.reference

    def verify(self):

        TYPE = {"compression": "Gz"}

        if not self.is_file():

            raise errors.VCFError(f"The file {self} does not exist.")

        if self.is_empty():

            raise errors.VCFError(f"The file {self} is empty.")
        
        archive = filetype.archive_match(self.path)
        
        # Is a archive
        if archive:

            if not (archive.__class__.__name__ == TYPE["compression"] and archive.EXTENSION == TYPE["compression"].lower()):

                raise errors.VCFError(f"The compression of {self} is not supported.")
            
            # The file is a Gz archive
            else:

                try:

                    with open(self.path, mode='rb') as vcf:

                        header = vcf.read(4)

                        # Gzip magic number and flag byte (3rd byte)
                        # If 3rd bit (0x04) is set, header has extra field.
                        if not (header[0:2] == b'\x1f\x8b' and header[3] & 0x04):

                            raise errors.VCFError(f"{self} is not a BGZF archive")
                        
                        else:

                            # Check for BC extra sub-field

                            pass

                    with gzip.open(self.path,mode='rt') as vcf:
                                
                        line = vcf.readline()

                        # Check if first line is empty
                        if not line:

                            raise errors.VCFError(f"First line of {self} is empty.")
                                
                        else:
                            # Check if first line start with "#"
                            if line[0] != '#':

                                raise errors.VCFError(f"First line inconsistent with VCF header format")
                            
                    self.archive: str = self.path
                            
                except FileNotFoundError:

                    raise errors.VCFError(f"{self} is not a valid path")
                
                except IOError:

                    raise errors.VCFError(f"An error occurred while reading {self}")
        
        # Is not a archive
        else:

            try:

                with open(self.path, mode='r') as vcf:

                    line = vcf.readline()

                    if not line:

                        raise errors.VCFError(f"First line of {self} is empty.")
                        
                    else:
                        # Check if first line start with "#"
                        if line[0] != '#':

                            raise errors.VCFError(f"First line inconsistent with VCF header format")
                        
            except FileNotFoundError:

                raise errors.VCFError(f"{self} is not a valid path")
                
            except IOError:

                raise errors.VCFError(f"An error occurred while reading {self}")
            
    def compressing(self, bin: str = "project"):

        logger.debug(f"Compressing file {self}.")

        EXT: str = "gz"

        archive: str = f"{self.path}.{EXT}"

        if bin == "env":
            CMDS: dict = {"bcftools": ["bcftools", "view", "-O", "z", "-o", archive, str(self.path)],
                        "bgzip": ["bgzip", "-c", "-f", str(self.path)]}
        else:
            binary = f"{os.path.dirname(os.path.abspath(__file__))}/htslib/bin/bgzip"
            os.chmod(binary, os.stat(binary).st_mode | stat.S_IEXEC)
            CMDS: dict = {"bgzip": [binary, "-c", "-f", str(self.path)]}

        outcode: int = 1
        retry: int = 0

        # Call a process to compress the file
        while retry < len(CMDS) and outcode != 0:
            bin: str = list(CMDS.keys())[retry]
            logger.debug(f"Compressing with {bin} binary.")
            try:
                process: subprocess.CompletedProcess = utils.runcmd(CMDS[bin], stdout=archive) if bin == "bgzip" else utils.runcmd(CMDS[bin])
                outcode: int = process.returncode
                logger.success(f"{self} has been successfully compressed.")
            except subprocess.CalledProcessError as e:
                logger.warning(f"Compressing {self} with {bin} binary did not succeed with exit code {outcode}.")
                logger.warning(f"Standard output: {e.stdout}")
                logger.warning(f"Standard error: {e.stderr}")
                retry += 1

        # Set archive path ONLY if exit code is 0
        self.archive = None if outcode else archive

    def indexing(self, bin: str = "project"):

        logger.debug(f"Indexing file {self}.")

        EXT: str = "tbi"

        index: str = f"{self.archive}.{EXT}"

        if bin == "env":
            CMDS: dict = {"bcftools": ["bcftools", "index", "-t", self.archive],
                        "tabix": ["tabix", "-f", "-p", "vcf", self.archive]}
        else:
            binary = f"{os.path.dirname(os.path.abspath(__file__))}/htslib/bin/tabix"
            os.chmod(binary, os.stat(binary).st_mode | stat.S_IEXEC)
            CMDS: dict = {"tabix": [binary, "-f", "-p", "vcf", self.archive]}

        outcode: int = 1
        retry: int = 0

        # Call a process to index the file
        while retry < len(CMDS) and outcode != 0:
            bin: str = list(CMDS.keys())[retry]
            logger.debug(f"Indexing {self} with {bin} binary.")
            try:
                process: subprocess.CompletedProcess = utils.runcmd(CMDS[bin])
                outcode: int = process.returncode
                logger.success(f"{self} has been successfully indexed.")
            except subprocess.CalledProcessError as e:
                logger.warning(f"Indexing {self} with {bin} binary did not succeed with exit code {outcode}.")
                logger.warning(f"Standard output: {e.stdout}")
                logger.warning(f"Standard error: {e.stderr}")
                retry += 1

        # Verify that the index file has been created
        # verify_files(file=file, index=index)

        # Return index path ONLY if exit code is 0
        self._index = None if outcode else VCFIndex(path=index, vcf=self)

    def open(self, lookup: str = None, context: bool = True):

        # Open the VCF with CyVCF2 and set the index for faster lookup
        self.stdin: cyvcf2.VCF = cyvcf2.VCF(self.archive, lazy=True)

        self.stdin.set_index(index_path=str(self._index.path))

        if not self.SAMPLES:

            samples = self.stdin.samples

            if len(samples) == 0:

                raise errors.VCFError(f"No sample found in {self}")
            
            else :

                logger.debug(
                f"{len(samples)} samples have been found in {self}: {samples}"
                )

                # Check if the samples are unique
                if len(samples) != len(set(samples)):

                    raise errors.VCFError(f"Duplicated sample names found in {self}")

            self.samples = samples

        if not self.variants.chromosomes:

            try:

                chromosomes = set(self.stdin.seqnames)

            except AttributeError as e:

                raise errors.VCFError(e)

            self.variants.chromosomes = chromosomes

            logger.debug(f"File {self} is composed of {self.variants.chromosomes} chromosomes")

        if not self.variants.seqlens:

            try:

                self.variants.seqlens = self.stdin.seqlens

            except AttributeError as e:

                logger.warning(e)

        if not self.variants.format:

            try:

                self.variants.format = str(next(self.stdin)).split()[self.HEADER["FORMAT"]]

            except Exception as e:

                logger.warning(e)

        if context:

            # Close the data stream
            self.close()

    def close(self):

        if isinstance(self.stdin, cyvcf2.VCF):

            self.stdin.close()

            self.stdin = None

    def update_header(self, samples: list):

        # Add the samples to the header
        self.HEADER.update(
            {
                s: i
                for i, s in zip(
                    range(
                        (self.HEADER["FORMAT"] + 1),
                        (self.HEADER["FORMAT"] + 1) + len(samples),
                    ),
                    samples,
                )
            }
        )

        logger.debug(
            f"Header for {self} has such format: {' '.join(self.HEADER.keys())}"
        )

    def package(self):

        return cp.deepcopy(self)

class VCFIndex(GenomicFile):

    """Class for VCF index files"""

    def __init__(self, path: str, vcf: VCF, lazy: bool = True):
        
        super().__init__(path)

        if not lazy:

            self.verify()

        self.vcf: VCF = vcf

    def verify(self):

        if not self.is_file():

            raise errors.IndexError(f"The file {self.path} does not exist.")

        if self.is_empty():

            raise errors.IndexError(f"The file {self.path} is empty.")
        
        # Index should be newer than VCF
        if os.path.getmtime(self.path) < os.path.getmtime(self.vcf.path):

            raise errors.IndexError(f"{self.path} is older than {self.vcf.path}.")

        else:
            # Try using the index file
            cmd = ["tabix", "-l", self.vcf.path]

            try:

                utils.runcmd(cmd)

            except subprocess.CalledProcessError:

                raise errors.IndexError(f"{self.path} does not allow fast lookup for {self.vcf.path}")
            
class VCFProcessor:

    def __init__(self):

        pass

    def is_supported():

        pass
    
    @staticmethod
    def preprocessing(vcf: VCF, bins: str = "project"):

        # File must be indexed
        if not vcf.is_indexed():
            
            # Compression
            if not vcf.archive:

                try:

                    vcf.compressing(bin=bins)

                except Exception as e:

                    raise errors.CompressionIndexError(f"Failed to compress {vcf}: {e}")

            try:

                # Call a process to index the file
                vcf.indexing(bin=bins)

            except Exception as e:

                raise errors.CompressionIndexError(f"Failed to index {vcf}: {e}")
            
        vcf.open(context=True)

    @staticmethod
    def process_chromosome(
        task,
        profile: bool = False,
    ) -> dict:
        """
        Process a chromosome from a VCF file
            chrom: containing the chromosome to process
            file: containing the path to the VCF file and the index file
            profile: to compute the statistics if a report is wanted
        """

        logger.debug(
            f"Processing chromosome {task[1]} for file {task[0]}"
        )

        # Open the VCF file and set the index for faster lookup
        try:
            task[0].open(context=False)
        except (FileNotFoundError, errors.VCFError) as e:
            logger.error(e)
            raise errors.ProcessError()

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
            if profile:
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
                        "inv": 0,
                        "csv": 0,
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
                for i, v in enumerate(task[0].stdin(f"{task[1]}")):
                    # Get the values from the VCF file
                    parts: list[str] = str(v).split("\t")
                    # Set the INFO values as a single character to reduce memory footprint
                    parts[task[0].header["INFO"]] = '.'
                    # First iteration, get the FORMAT values
                    if not i:

                        format: str = parts[task[0].header["FORMAT"]]

                        logger.debug(f"FORMAT for chromosome {task[1]}: {format}")

                    # From FORMAT get the values for each sample
                    samples_values: dict[str:dict] = {
                        s: VariantRepository.format_to_values(
                            format=format, values=parts[task[0].header[s]]
                        )
                        for s in task[0].samples
                    }
                    vts: list = VariantRepository.get_variant_type(v.REF, tuple(v.ALT))
                    # Should variant be filtered ?
                    if task[0].variants.filters:

                        exclude: bool = VariantRepository.exclude(v, task[0].variants.filters)
                    # Should the variant be excluded ?
                    if exclude:

                        filtered[vts[0]] += 1
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
                                vts[0],
                                "FAIL" if v.FILTER else "PASS",
                                '\t'.join(parts),
                            ]
                            # Should the statistics be computed ?
                            if profile:
                                # The variant type is common to all samples
                                # Iterare over uniq variant type
                                for vt in set(vts):
                                    if vt == "snp":
                                        if v.is_transition:
                                            mutation: str = "transition"
                                        else:
                                            mutation: str = "transversion"
                                        stats["variant"][vt][mutation] += 1
                                        stats["variant"][vt][v.REF][v.ALT[0]] += 1
                                    elif vt in ["del", "ins"]:
                                        mutation: str = "deletion" if vt == "del" else "insertion"
                                        stats["variant"]["indel"][mutation] += 1
                                    else:
                                        stats["variant"][vt] += 1

                                # Statistics unique to each samples
                                # If previous passes have raised a warning, do not search for genotype quality score.
                                if not "genotype_quality" in warnings:
                                    try:
                                        stats[task[0].format["genotype_quality"][0]].append(
                                            [
                                                samples_values[s][
                                                    task[0].format["genotype_quality"][0]
                                                ]
                                                for s in task[0].samples
                                            ]
                                        )
                                    except KeyError:
                                        logger.warning(
                                            f"Genotype quality value cannot be retrieved with key(s): {task[0].format['genotype_quality']}"
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
                                                    lambda sample: VariantRepository.is_homozygous(
                                                        GT=samples_values[sample][
                                                            task[0].format["genotype"][0]
                                                        ]
                                                    ),
                                                    task[0].samples,
                                                )
                                            )
                                        )
                                        # Sum heterozygous genotypes for each sample
                                        stats["het"] += sum(
                                            list(
                                                map(
                                                    lambda sample: VariantRepository.is_heterozygous(
                                                        GT=samples_values[sample][
                                                            task[0].format["genotype"][0]
                                                        ]
                                                    ),
                                                    task[0].samples,
                                                )
                                            )
                                        )
                                    except KeyError:
                                        logger.warning(f"Genotype type cannot be retrieved with key(s): {task[0].format['genotype']}")
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
                                                    samples_values[s][task[0].format["depth"][0]]
                                                    if task[0].format["depth"][0] in samples_values[s]
                                                    else [samples_values[s][task[0].format["depth"][1]]]
                                                ) # Try to get the depth value from the first FORMAT value, if not found, get it from the second FORMAT value
                                                for s in task[0].samples
                                            ]
                                        )
                                    except KeyError:
                                        logger.warning(f"Sequencing depth value cannot be retrieved with key(s): {task[0].format['depth']}")
                                        # Keep record of exception
                                        warnings["depth"] = True
        except Exception as e:

            logger.error(e)

            raise errors.ProcessError(e)

        # Close the data stream, avoid memory leaks
        task[0].close()

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

        variants["Type"].cat.set_categories(VariantRepository.VT)

        # Set the hash values as the index
        variants.index = Index(variants.index.values, dtype="string[pyarrow]")

        logger.debug(
            f"Filtered: {filtered['snp']} SNP(s), {filtered['indel']} INDEL(s), {filtered['sv']} structural variant(s) variant(s) for chromosome {task[1]} in file {task[0]}"
        )

        # Set the statistics computed as the right type for better memory management
        if profile:
            stats["depth"], stats["quality"], stats["GQ"] = (
                np.array(stats["depth"], dtype=np.uint16),
                np.array(stats["quality"], dtype=np.float16),
                np.array(stats["GQ"], dtype=np.uint16),
            )

        return (variants, filtered, stats) if profile else (variants, filtered, {})
    
    @staticmethod
    def serialize(task, obj, format: str = "pickle") -> int:
        """ Serialize DataFrame to file of specified format """

        def write(task, obj, prefix, ext):

            # Open the initial VCF file to get the template
            vcf = cyvcf2.VCF(f"{task[0].archive}")
            # Add match data in INFO column
            vcf.add_info_to_header({'ID': 'match', 'Description': 'overlapping variant', 'Type': 'String', 'Number': '1'})
            # Open the output VCF file for writing, open the data stream
            w = cyvcf2.Writer(f"{prefix}_delta.{ext}", vcf)
            # Iterate over the variants in the VCF file
            for i, v in enumerate(vcf):
                # Hash the variant to allow fast lookup
                hash = sha256(
                    string=f"{(v.CHROM).removeprefix('chr')}:{v.POS}:{v.REF}:{'|'.join(v.ALT)}".encode()
                ).hexdigest()
                # O(1) lookup
                if hash in obj["index"]:
                    # Get the corresponding row from the DataFrame, O(1) lookup by index
                    series = obj["variants"].iloc[obj["index"][hash]]
                    # Add the match information to the INFO dictionary
                    v.INFO["match"] = int((notna(series["Variant.L"])) & (notna(series["Variant.R"])))
                # Write the variant to the output VCF file
                w.write_record(v)
            # Close data streams
            w.close(); vcf.close()

        # Map file formats to respective file extensions, write modes, and functions
        FILES = {
                "json": {"ext": "json", "mode": "w", "func": json.dump},
                "pickle": {"ext": "pkl", "mode": "wb", "func": cPickle.dump},
                "vcf": {"ext": "vcf", "mode": "wz", "func": None},
                "vcf.gz": {"ext": "vcf.gz", "mode":"wz", "func": None}
            }
        
        assert (
            format in FILES
        ), "File format is not supported"

        if format in FILES:

            prefix = pathlib.Path(task[1], task[0].path.stem)

            if format in ["json", "pickle"]:
                # Open data stream
                with open(
                    file=f"{prefix}_delta.{FILES[format]['ext']}",
                    mode=FILES[format]["mode"],
                ) as f:
                    obj = obj.to_dict(orient='list')
                    FILES[format]["func"](obj, f)
            # Write a VCF
            else:
                write(task, obj, prefix, FILES[format]['ext'])

            assert os.path.isfile(
                f"{prefix}_delta.{FILES[format]['ext']}"
            ), "File was not created"

            logger.success(f"Results are seralized to {task[1]}")

            return 1
        else:
            raise ValueError(f"Error: The file format {format} is not supported.")
    
class VCFRepository():

    processor: VCFProcessor = VCFProcessor()

    def __init__(self, vcfs: list[str], index: list[str], reference: bool = False, filters: dict = None):

        assert len(vcfs) == len(index), "VFC(s) and Index(s) collections must be the same size" 

        self.populate(vcfs, index, reference, filters)

    def populate(self, vcfs: list[str], index: list[str], reference: bool = False, filters: dict = None):

        assert len(vcfs) == len(index), "VFC(s) and Index(s) collections must be the same size"

        self.REF = 0 if reference else None

        self.QUERY = 1 if reference else None

        self.repository: list[VCF] = [VCF(path=vcf, reference=(reference and not i), index=index[i], filters=filters, lazy=False) for i, vcf in enumerate(vcfs)]

    def add(self, vcfs: list[str], index: list[str]):

        assert len(vcfs) == len(index), "VFC(s) and Index(s) collections must be the same size"

        self.repository.extend([VCF(path=vcf, index=index[i], lazy=False) for i, vcf in enumerate(vcfs)])

    def get_reference(self) -> VCF | None:

        try:
            vcf = self.repository[self.REF]
        except TypeError:
            vcf = None

        return vcf

    def get_query(self) -> VCF | None:

        try:
            vcf = self.repository[self.QUERY]
        except TypeError:
            vcf = None

        return vcf

    def compare(self):

        results: dict = {}

        for pair in pairwise(self.repository):

            results.setdefault("pairs", []).append(pair)

            results[pair] = {
                "common": 0,
                "unique": {pair[0]: 0, pair[1]: 0},
                "jaccard": 0,
            }

            variantsL: DataFrame = pair[0].variants.collapse()
            
            variantsR: DataFrame = pair[1].variants.collapse()
            
            chromosomes = sorted(list(set(variantsL['Chromosome'].dropna().unique()) | set(variantsR['Chromosome'].dropna().unique())))

            variantsL["Chromosome"] = variantsL["Chromosome"].cat.set_categories(chromosomes)
            
            variantsR["Chromosome"] = variantsR["Chromosome"].cat.set_categories(chromosomes)

            filters = list(set(variantsL["Filter"].dropna()) | set(variantsR["Filter"].dropna()))

            variantsL["Filter"] = variantsL["Filter"].cat.set_categories(filters)
            
            variantsR["Filter"] = variantsR["Filter"].cat.set_categories(filters)

            # Compute unique variants in first VCF file
            results[pair]["unique"][pair[0]] = len(utils.difference(
                a=frozenset(variantsL.index),
                b=frozenset(variantsR.index),
            ))
            # Compute unique variants in second VCF file
            results[pair]["unique"][pair[1]] = len(utils.difference(
                a=frozenset(variantsR.index),
                b=frozenset(variantsL.index),
            ))
            # Compute common variants between both VCF files
            results[pair]["common"] = len(utils.intersect(
                a=frozenset(variantsL.index),
                b=frozenset(variantsR.index),
            ))

            logger.debug(
                f"{results[pair]['common']} variant(s) is/are commom in both files"
            )

            logger.debug(
                f"{results[pair]['unique'][pair[0]]} variant(s) is/are unique in files {pair[0]}"
            )

            logger.debug(
                f"{results[pair]['unique'][pair[1]]} variant(s) is/are unique in files {pair[1]}"
            )

            # Compute the Jaccard Index
            results[pair]["jaccard"] = utils.jaccard_index(
                shared=results[pair]["common"],
                total=(
                    results[pair]["common"]
                    + results[pair]["unique"][pair[0]]
                    + results[pair]["unique"][pair[1]]
                ),
            )

            logger.debug(f"Jaccard index: {results[pair]['jaccard']}")

            results[pair]["headers"] = {pair[0]: "\t".join(list(pair[0].header.keys())),
                                        pair[1]: "\t".join(list(pair[1].header.keys()))}

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
                        variantsL,
                        variantsR,
                    ],
                    ["L", "R"],
                )
            )
            # Merge the two DataFrames with a outer join algorithm,
            # Merge is made on the hash index
            df: DataFrame = concat(
                [
                    variantsL,
                    variantsR,
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
                    "Position": "uint64",
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
                key=lambda serie: serie.astype(str).map(VariantRepository.chromosome_sort_key) if serie.name == "Chromosome" else serie
            )

            results[pair]["variants"] = df

            # Key is a references to the Hash string object in the Dataframe
            # Allow a O(1) lookup while keeping memory footprint low
            lookup = {hash: row for row, hash in enumerate(df["Hash"])}

            results[pair]["index"] = lookup

            results[pair]["plots"] = PlotLibrary()

            # Create a Venn diagram to display the common, unique variants between the two VCF files
            results[pair]["plots"].venn((results[pair]["unique"][pair[0]], results[pair]["unique"][pair[1]], results[pair]["common"]), [pair[0].basename(), pair[1].basename()])

            # Should a benchmark be computed ?
            if pair[0].is_reference():

                logger.debug(f"Computing benchmark metrics from {pair[0]}.")

                # Create a mask of variants that passed the filter in one or both of the truth and query sets
                pass_mask = (df["Filter.L"] == "PASS") | (df["Filter.R"] == "PASS")

                series = []
                
                # Compute the performance metrics for SNPs and INDELs
                for v in ['snp','indel']:

                    # Create a mask of variants that are of the specified type
                    type_mask = df["Type"] == v

                    # Filter the DataFrame based on the variant type and the filter mask
                    filtered_df = df[type_mask & pass_mask]

                    # If FILTER.L and FILTER.R are equal, the variant is a true positive
                    is_match = filtered_df["Filter.L"] == filtered_df["Filter.R"]
                    # If FILTER.L is PASS and FILTER.R is either missing or FAIL, the variant is a false negative
                    is_truth_only = (filtered_df["Filter.L"] == "PASS") & ((isna(filtered_df["Filter.R"])) | (filtered_df["Filter.R"] == "FAIL"))
                    # If FILTER.R is PASS and FILTER.L is either missing or FAIL, the variant is a false positive
                    is_query_only = ((isna(filtered_df["Filter.L"])) | (filtered_df["Filter.R"] == "FAIL")) & (filtered_df["Filter.R"] == "PASS")

                    # Create a Series of the performance metrics
                    summary = Series([v, 
                                    'PASS', 
                                    (notna(filtered_df["Filter.L"])).sum(), 
                                    is_match.sum(), 
                                    is_truth_only.sum(), 
                                    (notna(filtered_df["Filter.R"])).sum(), 
                                    is_query_only.sum()], 
                                    index=["TYPE","FILTER","TRUTH.TOTAL","TRUTH.TP","TRUTH.FN","QUERY.TOTAL","QUERY.FP"])
                    
                    # Compute the recall, precision, and F1 score
                    try:
                        summary["RECALL"] = summary["TRUTH.TP"] / (summary["TRUTH.TP"] + summary["TRUTH.FN"])
                    except ZeroDivisionError:
                        summary["RECALL"] = 0.00

                    try:
                        summary["PRECISION"] = summary["TRUTH.TP"] / (summary["TRUTH.TP"] + summary["QUERY.FP"])
                    except ZeroDivisionError:
                        summary["PRECISION"] = 0.00
                        
                    try:
                        summary["F1"] = 2 * (summary["PRECISION"] * summary["RECALL"]) / (summary["PRECISION"] + summary["RECALL"])
                    except ZeroDivisionError:
                        summary["F1"] = 0.00

                    series.append(summary)

                # Return the performance metrics as a DataFrame, column types are specified for better memory management
                results[pair]["benchmark"] = DataFrame(series).astype({"TYPE": "category", 
                                                                        "FILTER": "category", 
                                                                        "TRUTH.TOTAL": "uint32", 
                                                                        "TRUTH.TP": "uint32", 
                                                                        "TRUTH.FN": "uint32", 
                                                                        "QUERY.TOTAL": "uint32", 
                                                                        "QUERY.FP": "uint32",
                                                                        "RECALL": "float64",
                                                                        "PRECISION": "float64",
                                                                        "F1": "float64"})

        return results

    def __len__(self):

        return len(self.repository)
            
class Report:
    """
    A class to represent a report
        vcfs: The list of VCF files
        prefix: The prefix of the report
        cmd: The command used to generate the report
        view: The view object as a DataFrame
        table: The benchmark table (truth metrics) object as a DataFrame
    """
    # Make use of __slots__ to avoid the creation of __dict__ and __weakref__ for each instance, reducing the memory footprint
    __slots__ = ('vcfs', 'tags', 'cmd', 'view', 'table', 'archive')

    def __init__(self, vcfs: VCFRepository, tags: list[str], cmd: str, view: dict, table: DataFrame = None, archive: bool = False):

        # The list of VCF files
        self.vcfs = vcfs

        # Tags provided by user about VCFs
        self.tags = tags

        # The command used to generate the report
        self.cmd = cmd
        
        # The view object as a DataFrame
        self.view = view

        # Set one library of plots to dark mode
        vcfs.repository[1].variants.plots.dark()

        # The benchmark table (truth metrics) object as a DataFrame
        self.table = table

        # How to save the outputs
        self.archive = archive

    def create(self, output: pathlib.Path = pathlib.Path.cwd()):

        # Custom filter to check if a value is NaN with Jinja2
        def is_nan(value):
            return isna(value)
        
        def format_variant_type(value:str):

            if value in ["ins", "del"]:

                vt = "indel"
            
            else:

                vt = value

            return vt
        
        def format_to_html(value:str):

            html_tags = {}

            if value != "<NA>":

                variant: list[str] = value.split()

                format: str = variant[VariantRepository.INDEX["FORMAT"]]

                values: str = variant[VariantRepository.INDEX["FORMAT"] + 1]

                html_tags: dict = VariantRepository.format_to_values(format, values)

            return html_tags

        def info_to_html(value:str):

            pass
        
        # Path to look for the template
        ressources = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'templates')

        # Path to look for the assets
        assets = os.path.join(ressources, "assets")

        # Path to look for the CSS stylesheets
        stylesheets = os.path.join(ressources, "*.css")

        # Path to look for the statics
        statics = os.path.join(ressources, 'statics', "*.png")
        
        # Create the environment
        env = Environment(loader=FileSystemLoader(ressources))

        # Add the custom filter to the environment
        env.filters['is_nan'] = is_nan

        env.filters['format_to_html'] = format_to_html

        env.filters['format_variant_type'] = format_variant_type

        # Load the template
        template = env.get_template("template.html")

        # Render the template
        html = template.render(vcfs=self.vcfs, tags=self.tags, cmd=self.cmd, view=self.view, table=self.table, date=datetime.datetime.now().strftime("%Y-%m-%d, %X"), version=__version__)

        if self.archive:

            path = output.joinpath("delta.zip")

            # Open data stream with context manager
            with zipfile.ZipFile(path, mode='w') as zip:
                zip.writestr("delta.html", html)
                for stylesheet in glob(stylesheets):
                    zip.write(stylesheet, arcname=os.path.relpath(path=stylesheet, start=ressources))
                for static in glob(statics):
                    zip.write(static, arcname=os.path.relpath(path=static, start=ressources))
                for root, dirs, files in os.walk(assets):
                    for file in files:
                        file_path = os.path.join(root, file)
                        zip.write(file_path, arcname=os.path.relpath(path=file_path, start=ressources))

        else:

            path = output.joinpath("delta")

            # Open data stream
            try:
                path.mkdir(exist_ok=True)
            except (PermissionError, FileNotFoundError, FileExistsError) as e:
                raise errors.ReportError(f"Report cannot be created: {e}")

            with open(path.joinpath("delta.html"),'w') as f:
                f.writelines(html)

            # Copy the assets and statics files next to the report
            for f in glob(stylesheets):
                copy(f, path)

            copytree(assets, os.path.join(path, "assets"), dirs_exist_ok=True)

            copytree(os.path.dirname(statics), os.path.join(path,'statics'), dirs_exist_ok=True)

            if os.path.exists(path.joinpath("delta.html")):
                webbrowser.open(str(path.joinpath("delta.html")))
            else:
                raise errors.ReportError("Report not found on filesystem.")

    def __str__(self):
        
        pass

    def __repr__(self):
        
        pass