import concurrent.futures
from cyvcf2 import VCF
from hashlib import sha256
from loguru import logger
from multiprocessing import Manager
import utils
import subprocess

class VCFLibrary:

    def __init__(self, **kwargs):

        self.library = []
        
        for key, value in kwargs.items():
            setattr(self, key, value)

@logger.catch
def process_chromosome(chrom: str, FILES: dict[str:str]) -> dict:
    
    logger.debug(f"Processing chromosome {chrom} for file {FILES["compression"]}")

    try:
        vcf = VCF(FILES["compression"], lazy=True)
        vcf.set_index(index_path=FILES["index"])
    except FileNotFoundError as e:
        logger.error(e)

    result = {}

    try:

        for v in vcf(f'{chrom}'):

            hash = sha256(string=f"{v.CHROM}:{v.POS}:{v.REF}:{"|".join(v.ALT)}".encode()).hexdigest()
            result[hash] = str(v)
    
    except UserWarning as e:
        logger.warning(e)

    return result

@logger.catch
def process_files(file: str, index: str = None) -> dict:

    logger.debug(f"Processing file: {file}")

    try:
        utils.verify_file(file=file)
    except (FileNotFoundError, ValueError) as e:
        logger.error(e)

    if(index):

        try:
            utils.verify_file(file=index)
        except (FileNotFoundError, ValueError) as e:
            logger.error(e)
            
    else:
        # Indexing
        logger.debug(f"Indexing file: {file}")

        try:
            code = subprocess.run(["./src/indexing.sh", file], capture_output=True, text=True, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            logger.error(e.stderr)

    FILES = {"compression": f"{file}.gz", 
            "index": f"{file}.gz.tbi"}
    
    try:
        utils.verify_file(file=FILES["compression"])
    except (FileNotFoundError, ValueError) as e:
        logger.error(e)

    try:
        vcf = VCF(FILES["compression"])
        vcf.set_index(index_path=FILES["index"])
    except FileNotFoundError as e:
        logger.error(e)

    chromosomes = vcf.seqnames

    logger.debug(f"File {file} is composed of {chromosomes} chromosomes")

    result = {}

    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as chrom_executor:

        futures_to_chrom = {chrom_executor.submit(process_chromosome, chrom, FILES): chrom for chrom in chromosomes}

        for future in concurrent.futures.as_completed(futures_to_chrom):

            try:
                result[futures_to_chrom[future]] = future.result()
            except Exception as e:
                logger.warning(f"Chromosome {futures_to_chrom[future]} generated an exception: {e}")

    return result

def intersect(a: set[str], b: set[str]) -> set:
    return a & b

def difference(a: set[str], b: set[str]) -> set:
    return a - b

def delta(params: object) -> int:

    logger.debug(f"VCFS: {params.vcfs}")

    logger.debug(f"Indexes: {params.indexes}")

    assert len(params.vcfs) == 2, "Two VCF files are required"

    assert isinstance(params.vcfs[0],str) and isinstance(params.vcfs[1],str), "Input vcf should be string instance"

    result = {}

    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as files_pool:

        futures_to_vcf = {files_pool.submit(process_files, vcf, index): vcf for vcf, index in zip(params.vcfs,params.indexes)}

        for future in concurrent.futures.as_completed(futures_to_vcf):

            try:
                result[futures_to_vcf[future]] = future.result()
            except Exception as e:
                logger.error(f"File {futures_to_vcf[future]} generated an exception: {e}")

    for chrom in result[params.vcfs[0]]:
                
        print(difference(a=set(result[params.vcfs[0]][chrom].keys()), b=set(result[params.vcfs[1]][chrom].keys())))

        print(intersect(a=set(result[params.vcfs[0]][chrom].keys()), b=set(result[params.vcfs[1]][chrom].keys())))
    
    return 1