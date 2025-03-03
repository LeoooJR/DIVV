import concurrent.futures
from cyvcf2 import VCF
from loguru import logger
import utils
import subprocess

@logger.catch
def process_chromosome():
    pass

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
    
    try:
        utils.verify_file(file=file+'.gz')
    except (FileNotFoundError, ValueError) as e:
        logger.error(e)

    try:
        vcf = VCF(file+'.gz')
    except FileNotFoundError as e:
        logger.error(e)

    chromosomes = vcf.seqnames

    logger.debug(f"File {file} is composed of {chromosomes} chromosomes")

    return {}

def delta(params: object) -> int:

    logger.debug(f"VCFS: {params.vcfs}")

    logger.debug(f"Indexes: {params.indexes}")

    assert len(params.vcfs) == 2, "Two VCF files are required"

    result = []

    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as pool:

        result = [pool.submit(process_files, vcf, index) for vcf, index in zip(params.vcfs,params.indexes)]

    print(result[0])
    
    return 1