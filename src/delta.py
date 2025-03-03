import concurrent.futures
from cyvcf2 import VCF
from loguru import logger
import utils
import subprocess

@logger.catch
def process_chromosome(chrom: str, file: VCF) -> dict:
    
    logger.debug(f"Processing chromosome {chrom} for file {file}")

    for v in file(f'{chrom}'):
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

    result = {}

    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as pool:

        futures = {pool.submit(process_chromosome, chrom, vcf):chrom for chrom in chromosomes}

        for future in concurrent.futures.as_completed(futures):

            try:
                result[futures[future]] = future.result()
            except Exception as e:
                logger.warning(f"Chromosome {futures[future]} generated an exception: {e}")

    return result

def delta(params: object) -> int:

    logger.debug(f"VCFS: {params.vcfs}")

    logger.debug(f"Indexes: {params.indexes}")

    assert len(params.vcfs) == 2, "Two VCF files are required"

    assert isinstance(params.vcfs[0],str) and isinstance(params.vcfs[1],str), "Input vcf should be string instance"

    result = {}

    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as pool:

        futures = {pool.submit(process_files, vcf, index):vcf for vcf, index in zip(params.vcfs,params.indexes)}

        for future in concurrent.futures.as_completed(futures):

            try:
                result[futures[future]] = future.result
            except Exception as e:
                logger.error(f"File {futures[future]} generated an exception: {e}")
                
    print(result[0])
    
    return 1