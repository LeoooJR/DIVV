import cyvcf2
import errors
import filetype
import gzip
from loguru import logger
import os
import subprocess
from utils import runcmd

class GenomicFile():

    """Base class for genomic files"""

    def __init__(self, path: str):

        # Path to the file
        self.path = path

    def is_empty(self) -> bool:
        """Check if file is empty"""

        return os.path.getsize(self.path) == 0

    def is_file(self) -> bool:
        """Check if path is a file"""

        return os.path.isfile(self.path)

    def get_path(self) -> str:

        return self.path
    
    def __str__(self):
        
        return self.path
    
    def __repr__(self):
        
        return self.path
    
    def __hash__(self):
        
        return hash(self.path)
    
class VCF(GenomicFile):

    """Class for VCF files"""

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

    def __init__(self, path: str, reference: bool = False, index: str = None,  lazy: bool = True):
        
        super().__init__(path)

        self.reference: bool = reference

        self.archive: str = None

        self.SAMPLES = None

        self.chromosomes = None

        if index:

            self._index: VCFIndex = VCFIndex(path=index, vcf=self)

        if not lazy:

            self.verify()

    @property
    def index(self):

        return self._index
    
    @index.setter
    def index(self, value: str):

        self._index: VCFIndex = VCFIndex(path=value)

    def is_indexed(self) -> bool:

        return isinstance(self._index, VCFIndex)
    
    @property
    def header(self):

        return self.HEADER
    
    @property
    def samples(self):

        return self.SAMPLES
    
    @samples.setter
    def samples(self, values: list[str]):

        self.SAMPLES = values

        self.update_header(samples=values)

    def verify(self):

        TYPE = {"compression": "Gz"}

        if not self.is_file():

            raise errors.VCFError(f"The file {self} does not exist.")

        if not self.is_empty():

            raise errors.VCFError(f"The file {self} is empty.")
        
        archive = filetype.archive_match(self.path)
        
        # Is a archive
        if archive:

            if not (archive.__class__.__name__ == TYPE["compression"] and archive.EXTENSION == TYPE["compression"].lower()):

                raise errors.VCFError(f"The compression of {self} is not supported.")
            
            # The file is a Gz archive
            else:

                self.archive: str = self.path

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
                            
                except FileNotFoundError:

                    raise errors.VCFError(f"{self} is not a valid path")
                
                except IOError:

                    raise errors.VCFError(f"An error occurred while reading {self}")
                
                if self._index:

                    self._index.verify()
        
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
            
    def compressing(self):

        logger.debug(f"Compressing file {self}.")

        EXT: str = "gz"

        archive: str = f"{self.path}.{EXT}"

        CMDS: dict = {"bcftools": ["bcftools", "view", "-O", "z", "-o", archive, self.path],
                    "bgzip": ["bgzip", "-c", self.path]}

        outcode: int = 1
        retry: int = 0

        # Call a process to compress the file
        while retry < len(CMDS) and outcode != 0:
            bin: str = list(CMDS.keys())[retry]
            logger.debug(f"Compressing with {bin} binary.")
            try:
                process: subprocess.CompletedProcess = runcmd(CMDS[bin], stdout=archive) if bin == "bgzip" else runcmd(CMDS[bin])
                outcode: int = process.returncode
                logger.success(f"{self} has been successfully compressed.")
            except subprocess.CalledProcessError as e:
                logger.warning(f"Compressing {self} with {bin} binary did not succeed with exit code {outcode}.")
                logger.warning(f"Standard output: {e.stdout}")
                logger.warning(f"Standard error: {e.stderr}")
                retry += 1

        # Set archive path ONLY if exit code is 0
        self.archive = None if outcode else archive

    def indexing(self):

        logger.debug(f"Indexing file {self}.")

        EXT: str = "tbi"

        index: str = f"{self.path}.{EXT}"

        CMDS: dict = {"bcftools": ["bcftools", "index", "-t", self.path],
                    "tabix": ["tabix", "-f", "-p", "vcf", self.path]}

        outcode: int = 1
        retry: int = 0

        # Call a process to index the file
        while retry < len(CMDS) and outcode != 0:
            bin: str = list(CMDS.keys())[retry]
            logger.debug(f"Indexing {self} with {bin} binary.")
            try:
                process: subprocess.CompletedProcess = runcmd(CMDS[bin])
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

    def open(self, lookup: str, context: bool = True):

        self.stdin: cyvcf2.VCF = cyvcf2.VCF(self.archive)

        self.stdin.set_index(index_path=self._index.path)

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

        if not self.chromosomes:

            chromosomes = self.stdin.seqnames

            if len(chromosomes) == 0:

                raise errors.VCFError(f"No chromosome found in {self}")
            
            else:

                logger.debug(f"File {self} is composed of {chromosomes} chromosomes")

            self.chromosomes = chromosomes

        if context:

            # Close the data stream
            self.close()

    def close(self):

        if isinstance(self.stdin, cyvcf2.VCF):

            self.stdin.close()

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

class VCFIndex(GenomicFile):

    """Class for VCF index files"""

    def __init__(self, path: str, vcf: VCF, lazy: bool = True):
        
        super().__init__(path)

        self.vcf: VCF = vcf

        if not lazy:

            self.verify()

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

                runcmd(cmd)

            except subprocess.CalledProcessError:

                raise errors.IndexError(f"{self.path} does not allow fast lookup for {self.vcf.path}")