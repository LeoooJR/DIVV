import concurrent.futures
from cyvcf2 import VCF
from hashlib import sha256
from loguru import logger
import numpy as np
import utils
import subprocess


class VCFLibrary:

    def __init__(self, **kwargs):

        self.library = []

        for key, value in kwargs.items():
            setattr(self, key, value)


@logger.catch
def process_chromosome(
    chrom: str, FILES: dict[str:str], filters: dict = None, stats: bool = False
) -> dict:

    logger.debug(
        f"Processing chromosome {chrom} for file {FILES['compression']}"
    )

    try:
        vcf = VCF(FILES["compression"], lazy=True)
        vcf.set_index(index_path=FILES["index"])
    except FileNotFoundError as e:
        logger.error(e)

    try:

        exclude: bool = False

        result, filtered = {}, {"snp": 0,
                                "mnp": 0, 
                                "indel": 0, 
                                "sv": 0,
                                "transition": 0}

        if stats:
            values = {
                "variant": 0,
                "depth": [],
                "quality": [],
                "GQ": [],
                "ref": 0,
                "het": 0,
                "hom": 0,
            }

        for v in vcf(f"{chrom}"):

            if filters:

                exclude: bool = utils.exclude(v,filters)

            if exclude:

                filtered[
                    v.var_type
                ] += 1

            else:

                hash = sha256(
                    string=f"{v.CHROM}:{v.POS}:{v.REF}:{'|'.join(v.ALT)}".encode()
                ).hexdigest()
                
                result[hash] = str(v)

                if stats:
                    values["quality"].append(v.QUAL)

    except UserWarning as e:
        logger.warning(e)

    logger.debug(
        f"Filtered: {filtered['snp']} SNP(s), {filtered['indel']} INDEL(s), {filtered['sv']} structural variant(s) variant(s) for chromosome {chrom} in file {FILES['compression']}"
    )

    if stats:
        values["depth"], values["quality"], values["GQ"] = (
            np.array(values["depth"], dtype=np.uint8),
            np.array(values["quality"], dtype=np.float32),
            np.array(values["GQ"]),
        )

    return result, filtered


@logger.catch
def process_files(
    file: str, index: str = None, filters: dict = None, stats: bool = False
) -> dict:

    logger.debug(f"Processing file: {file}")

    try:
        utils.verify_file(file=file)
    except (FileNotFoundError, ValueError) as e:
        logger.error(e)

    FILES = {"compression": f"{file}.gz", "index": f"{file}.gz.tbi"}

    if index:

        try:
            utils.verify_file(file=index)
        except (FileNotFoundError, ValueError) as e:
            logger.error(e)

    else:
        # Try to look for indexing files
        try:
            utils.verify_file(FILES["compression"])
            utils.is_indexed(FILES["index"])
        except (FileNotFoundError, ValueError) as e:
            logger.warning(e)

            # Indexing
            logger.debug(f"Indexing file: {file}")

            try:
                code = subprocess.run(
                    ["./src/indexing.sh", file],
                    capture_output=True,
                    text=True,
                    check=True,
                )
            except (subprocess.CalledProcessError, FileNotFoundError) as e:
                logger.error(e.stderr)

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

    result, filtered = {}, {}

    with concurrent.futures.ProcessPoolExecutor(
        max_workers=2
    ) as chrom_executor:

        futures_to_chrom = {
            chrom_executor.submit(
                process_chromosome, chrom, FILES, filters, stats
            ): chrom
            for chrom in chromosomes
        }

        for future in concurrent.futures.as_completed(futures_to_chrom):

            logger.success(
                f"Process {future} for chromosome {futures_to_chrom[future]} in file {FILES['compression']} has completed."
            )

            try:
                (
                    result[futures_to_chrom[future]],
                    filtered[futures_to_chrom[future]],
                ) = future.result()
            except Exception as e:
                logger.warning(
                    f"Chromosome {futures_to_chrom[future]} generated an exception: {e}"
                )

    return {"info": utils.file_stats(file),
            "data": result, 
            "filter": filtered}


def delta(params: object) -> int:

    logger.debug(f"VCFS: {params.vcfs}")

    logger.debug(f"Indexes: {params.indexes}")

    logger.debug(f"Serialize output: {params.serialize}")

    logger.debug(f"Compute statistics: {params.stats}")

    assert len(params.vcfs) == 2, "Two VCF files are required"

    assert isinstance(params.vcfs[0], str) and isinstance(
        params.vcfs[1], str
    ), "Input vcf should be string instance"

    result = {}

    FILTERS = (
        {
            "threshold": params.threshold,
            "exclude": {
                "exclude_snps": params.exclude_snps,
                "exclude_indels": params.exclude_indels,
                "exclude_vars": params.exclude_vars,
                "exclude_mnps": params.exclude_mnps,
                "exclude_transitions": params.exclude_trans,
                "exclude_svs": params.exclude_svs
            }
        }
        if any(
            [
                params.threshold,
                params.exclude_snps,
                params.exclude_indels,
                params.exclude_vars,
                params.exclude_mnps,
                params.exclude_trans,
                params.exclude_svs
            ]
        )
        else None
    )

    logger.debug(f"Filters used: {FILTERS}")

    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as files_pool:

        futures_to_vcf = {
            files_pool.submit(
                process_files, vcf, index, FILTERS, params.stats
            ): vcf
            for vcf, index in zip(params.vcfs, params.indexes)
        }

        for future in concurrent.futures.as_completed(futures_to_vcf):

            logger.success(
                f"Process {future} for file {futures_to_vcf[future]} has completed."
            )

            try:
                (
                    result[futures_to_vcf[future]]
                ) = future.result()
            except Exception as e:
                logger.error(
                    f"File {futures_to_vcf[future]} generated an exception: {e}"
                )

    common_variants, unique_variants_to_left, unique_variants_to_right = (
        {},
        {},
        {},
    )

    print(result[params.vcfs[0]])

    for chrom in result[params.vcfs[0]]["data"]:

        unique_variants_to_left[chrom] = {
            k: result[params.vcfs[0]]["data"][chrom][k]
            for k in utils.difference(
                a=set(result[params.vcfs[0]]["data"][chrom].keys()),
                b=set(result[params.vcfs[1]]["data"][chrom].keys()),
            )
        }

        unique_variants_to_right[chrom] = {
            k: result[params.vcfs[1]]["data"][chrom][k]
            for k in utils.difference(
                a=set(result[params.vcfs[1]]["data"][chrom].keys()),
                b=set(result[params.vcfs[0]]["data"][chrom].keys()),
            )
        }

        common_variants[chrom] = {
            k: result[params.vcfs[0]]["data"][chrom][k]
            for k in utils.intersect(
                a=set(result[params.vcfs[0]]["data"][chrom].keys()),
                b=set(result[params.vcfs[1]]["data"][chrom].keys()),
            )
        }

    if params.serialize:
        path: str = "/".join(params.vcfs[0].split("/")[:-1])
        logger.debug(f"Results are seralized to {path}")
        try:
            utils.save(
                obj=common_variants,
                prefixe=f"{path}/common",
                format=params.serialize,
            )
        except ValueError as e:
            logger.error(e)

    return 1
