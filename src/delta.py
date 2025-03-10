import concurrent.futures
from cyvcf2 import VCF
from hashlib import sha256
from itertools import chain, repeat
from loguru import logger
import numpy as np
from operator import itemgetter
from os.path import basename
from pandas import DataFrame, concat
import plot
import subprocess
from template import Report
import pprint
import utils


class VCFLibrary:

    def __init__(self, **kwargs):

        self.library = []

        for key, value in kwargs.items():
            setattr(self, key, value)


@logger.catch
def process_chromosome(
    chrom: str,
    samples: list,
    header: dict,
    FILES: dict[str:str],
    filters: dict = None,
    compute: bool = False,
) -> dict:

    logger.debug(
        f"Processing chromosome {chrom} for file {FILES['compression']}"
    )

    try:
        vcf = VCF(FILES["compression"], lazy=True)
        vcf.set_index(index_path=FILES["index"])
    except FileNotFoundError as e:
        logger.error(e)

    FORMAT = {
        "genotype": ["GT"],
        "genotype_quality": ["GQ"],
        "depth": ["DP", "TRC"],
    }

    try:

        exclude: bool = False

        variants, filtered = {}, {
            "snp": 0,
            "mnp": 0,
            "indel": 0,
            "sv": 0,
            "transition": 0,
        }

        if compute:
            stats = {
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
                },
                "depth": [],
                "quality": [],
                "GQ": [],
                "ref": 0,
                "het": 0,
                "hom": 0,
            }

        for i, v in enumerate(vcf(f"{chrom}")):

            v_list = str(v).split("\t")

            if not i:

                format = v_list[header["FORMAT"]]

                logger.debug(f"FORMAT for chromosome {chrom}: {format}")

            samples_data = {
                s: utils.format_to_values(
                    format=format, values=v_list[header[s]]
                )
                for s in samples
            }

            if filters:

                exclude: bool = utils.exclude(v, filters)

            if exclude:

                filtered[v.var_type] += 1

            else:

                hash = sha256(
                    string=f"{v.CHROM}:{v.POS}:{v.REF}:{'|'.join(v.ALT)}".encode()
                ).hexdigest()

                variants[hash] = [v.CHROM,
                                  v.POS,
                                  str(v)
]
                if compute:

                    if v.var_type == "snp":
                        if v.is_transition:
                            mutation = "transition"
                        else:
                            mutation = "transversion"
                        stats["variant"][v.var_type][mutation] += 1
                        stats["variant"][v.var_type][v.REF][v.ALT[0]] += 1
                    elif v.var_type == "indel":
                        if v.is_deletion:
                            mutation = "deletion"
                        else:
                            mutation = "insertion"
                        stats["variant"][v.var_type][mutation] += 1
                    else:
                        stats["variant"][v.var_type] += 1

                    if (
                        FORMAT["genotype_quality"][0]
                        in samples_data[samples[0]].keys()
                    ):
                        stats[FORMAT["genotype_quality"][0]].append(
                            [
                                samples_data[s][FORMAT["genotype_quality"][0]]
                                for s in samples
                            ]
                        )

                    stats["hom"] += sum(
                        list(
                            map(
                                lambda x: utils.is_homozygous(
                                    GT=samples_data[x][FORMAT["genotype"][0]]
                                ),
                                samples,
                            )
                        )
                    )
                    stats["het"] += sum(
                        list(
                            map(
                                lambda x: utils.is_heterozygous(
                                    GT=samples_data[x][FORMAT["genotype"][0]]
                                ),
                                samples,
                            )
                        )
                    )

                    if v.QUAL:
                        stats["quality"].append(v.QUAL)

                    stats["depth"].append(
                        [
                            (
                                samples_data[s][FORMAT["depth"][0]]
                                if FORMAT["depth"][0] in samples_data[s]
                                else [samples_data[s][FORMAT["depth"][1]]]
                            )
                            for s in samples
                        ]
                    )

    except UserWarning as e:
        logger.warning(e)

    variants = DataFrame.from_dict(variants, orient='index', columns=["Chromosome","Position","Variant"])

    logger.debug(
        f"Filtered: {filtered['snp']} SNP(s), {filtered['indel']} INDEL(s), {filtered['sv']} structural variant(s) variant(s) for chromosome {chrom} in file {FILES['compression']}"
    )

    if compute:
        stats["depth"], stats["quality"], stats["GQ"] = (
            np.array(stats["depth"], dtype=np.uint16),
            np.array(stats["quality"], dtype=np.float32),
            np.array(stats["GQ"]),
        )

    return (variants, filtered, stats) if compute else (variants, filtered, {})


@logger.catch
def process_files(
    file: str, index: str = None, filters: dict = None, compute: bool = False
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

    samples = vcf.samples

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

    HEADER.update(
        {
            s: i
            for i, s in zip(
                range(
                    (HEADER["FORMAT"] + 1),
                    (HEADER["FORMAT"] + 1) + len(samples),
                ),
                samples,
            )
        }
    )

    logger.debug(
        f"{len(samples)} samples have been found in {file}: {samples}"
    )

    logger.debug(f"File {file} is composed of {chromosomes} chromosomes")

    logger.debug(
        f"Header for {file} has such format: {' '.join(HEADER.keys())}"
    )

    variants, filtered, stats = {}, {}, {}

    with concurrent.futures.ProcessPoolExecutor(
        max_workers=2
    ) as chrom_executor:

        futures_to_chrom = {
            chrom_executor.submit(
                process_chromosome,
                chrom,
                samples,
                HEADER,
                FILES,
                filters,
                compute,
            ): chrom
            for chrom in chromosomes
        }

        for future in concurrent.futures.as_completed(futures_to_chrom):

            logger.success(
                f"Process {future} for chromosome {futures_to_chrom[future]} in file {FILES['compression']} has completed."
            )

            try:
                (
                    variants[futures_to_chrom[future]],
                    filtered[futures_to_chrom[future]],
                    stats[futures_to_chrom[future]],
                ) = future.result()

                if compute:

                    stats[futures_to_chrom[future]]["length"] = vcf.seqlens[chromosomes.index(futures_to_chrom[future])]
                
            except Exception as e:
                logger.warning(
                    f"Chromosome {futures_to_chrom[future]} generated an exception: {e}"
                )

    if compute:
        library = plot.PlotLibrary(file=file)
        plot.visualization(file=basename(file), stats=stats, library=library)

    return {
        "info": utils.file_stats(file),
        "variants": variants,
        "filter": filtered,
        "plots": library,
    }


def delta(params: object) -> int:

    if not params.verbosity:

        logger.remove(0)

        logger.add("VCFDelta.log")

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
                "exclude_svs": params.exclude_svs,
                "pass_only": params.pass_only
            },
        }
        if any(
            [
                params.threshold,
                params.exclude_snps,
                params.exclude_indels,
                params.exclude_vars,
                params.exclude_mnps,
                params.exclude_trans,
                params.exclude_svs,
                params.pass_only
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
                (result[futures_to_vcf[future]]) = future.result()
            except Exception as e:
                logger.error(
                    f"File {futures_to_vcf[future]} generated an exception: {e}"
                )

    common, unique_vcf0, unique_vcf1 = (
        {},
        {},
        {},
    )

    result["delta"] = {
        "common": {},
        "unique": {params.vcfs[0]: {}, 
                   params.vcfs[1]: {}},
        "jaccard": {}
    }

    for chrom in result[params.vcfs[0]]["variants"]:

        unique_vcf0[chrom] = utils.difference(a=set(result[params.vcfs[0]]["variants"][chrom].index),
                                              b=set(result[params.vcfs[1]]["variants"][chrom].index))

        unique_vcf1[chrom] = utils.difference(a=set(result[params.vcfs[1]]["variants"][chrom].index),
                                              b=set(result[params.vcfs[0]]["variants"][chrom].index))

        common[chrom] = utils.intersect(a=set(result[params.vcfs[0]]["variants"][chrom].index),
                                        b=set(result[params.vcfs[1]]["variants"][chrom].index))

        (
            result["delta"]["common"][chrom],
            result["delta"]["unique"][params.vcfs[0]][chrom],
            result["delta"]["unique"][params.vcfs[1]][chrom],
        ) = (
            len(common[chrom]),
            len(unique_vcf0[chrom]),
            len(unique_vcf1[chrom]),
        )

        logger.debug(
            f"{result['delta']['common'][chrom]} variant(s) is/are commom for chromosome {chrom} in both files"
        )

        result["delta"]["jaccard"][chrom] = utils.jaccard_index(shared=result["delta"]["common"][chrom], 
                                                                total=(result["delta"]["common"][chrom]+
                                                                        result["delta"]["unique"][params.vcfs[0]][chrom]+
                                                                        result["delta"]["unique"][params.vcfs[1]][chrom]))

        logger.debug(f"Jaccard index for chromosome {chrom}: {result['delta']['jaccard'][chrom]}")

        logger.debug(
            f"{result['delta']['unique'][params.vcfs[0]][chrom]} variant(s) is/are unique for chromosme {chrom} in files {params.vcfs[0]}"
        )

        logger.debug(
            f"{result['delta']['unique'][params.vcfs[1]][chrom]} variant(s) is/are unique for chromosome {chrom} in files {params.vcfs[1]}"
        )

    dfs_chroms: list[DataFrame] = list(map(lambda vcf: list(itemgetter(*list(sorted(result[vcf]["variants"].keys())))(result[vcf]["variants"])), params.vcfs))

    dfs_files: list[DataFrame] = list(map(lambda x: concat(x), dfs_chroms))

    list(map((lambda x, n: x.rename(columns={c: f'{c}_vcf{n}' for c in x.columns}, inplace=True)), dfs_files, [1,2]))

    df: DataFrame = concat(dfs_files, axis=1, join='outer', sort=False)

    if params.serialize:
        path: str = "/".join(params.vcfs[0].split("/")[:-1])
        logger.debug(f"Results are seralized to {path}")
        try:
            utils.save(
                obj=df,
                prefixe=f"{path}/common",
                format=params.serialize,
            )
        except ValueError as e:
            logger.error(e)

    if params.report:

        Report(vcfs=params.vcfs, df=df, plots={params.vcfs[0]:result[params.vcfs[0]]["plots"].as_html(),
                                              params.vcfs[1]:result[params.vcfs[1]]["plots"].as_html()}).create()

    return 1
