from collections.abc import Callable
import concurrent.futures
import errors
import files
from hashlib import sha256
from itertools import chain
from loguru import logger
import numpy as np
from pandas import Index, DataFrame
from psutil import cpu_count
import utils

class TasksManager():

    def __init__(self, vcfs: list[files.VCF], processes: int = 0):

        assert processes >= 0, "Number of processes alocated must be an positive unsigned integer."
        
        if processes:

            if processes < 0:

                raise ValueError("Number of processes alocated must be an positive unsigned integer.")

            # Number of CPUs available
            CPUS = cpu_count(logical=True)

            logger.debug(f"Number of logical(s) CPU(s) detected: {CPUS}.")

            if CPUS > 1:
                self.processes = min(processes, CPUS)
                logger.debug(f"Tasks will be parallelized on {self.processes} CPU(s).")
            else:
                self.processes = 0
                logger.debug("Computation will be carried out sequentially.")

        else:
        
            self.processes: int = processes

        self.vcfs: list = vcfs

        self.tasks = None

        self.results: dict = {}

    def scheduling(self, tasks: list[list[object]]):

        if len(tasks) != len(self.vcfs):

            raise ValueError(f"{len(self.vcfs)} collection(s) of tasks is/are required.")

        self.tasks: dict[files.VCF:list] = {vcf: [(vcf, task) for task in tasks[i]] for i, vcf in enumerate(self.vcfs)}

        logger.debug(f"Task(s) scheduled {self.flatten()}")

    def flatten(self):

        return list(chain.from_iterable(self.tasks.values()))

    def commit(self, job: Callable, jobargs: dict[str:object], initializer: Callable = None, initargs: dict[str:object] = None):

        if self.results:

            self.results.clear()

        if initializer:

            for vcf in self.vcfs:

                try:

                    initializer(vcf, *initargs)

                except Exception as e:

                    raise errors.ProcessError(e)

        if self.processes:
            # Process the chromosomes concurrently, using as much process as available
            # use a context manager to avoid memory leaks
            with concurrent.futures.ProcessPoolExecutor(
                max_workers=self.processes
            ) as chrom_executor:
                # Submit the processes
                futures = {
                    chrom_executor.submit(
                        job,
                        task,
                        *jobargs
                    ): task
                    for task in self.flatten()
                }
                # Check if a process is completed
                for future in concurrent.futures.as_completed(futures):

                    # Get the returned result
                    try:

                        self.results[futures[future]] = future.result()
                    
                        logger.success(
                            f"Process {future} for task {futures[future]} has completed."
                        )

                    except errors.ProcessError as e:

                        logger.warning(
                            f"Process {future} for task {futures[future]} generated an exception: {e}"
                        )
                        
        # Computation is carried out sequentially.                
        else:

            for task in self.flatten():

                self.results[futures[future]] = job(task, *jobargs)

class VCFProcessor:

    def __init__(self):

        pass

    def is_supported():

        pass
    
    @staticmethod
    def preprocessing(vcf: files.VCF):

        # File must be indexed
        if not vcf.is_indexed():

            try:
            
                # Compression
                if not vcf.archive:

                    vcf.compressing()

                # Call a process to index the file
                vcf.indexing()

            except Exception as e:

                raise errors.CompressionIndexError(f"Failed to compress or/and index {vcf}")
            
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
                        s: utils.format_to_values(
                            format=format, values=parts[task[0].header[s]]
                        )
                        for s in task[0].samples
                    }
                    # Should variant be filtered ?
                    if task[0].variants.filters:

                        exclude: bool = utils.exclude(v, task[0].variants.filters)
                    # Should the variant be excluded ?
                    if exclude:

                        filtered[v.var_type] += 1
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
                                v.var_type,
                                "FAIL" if v.FILTER else "PASS",
                                '\t'.join(parts),
                            ]
                            # Should the statistics be computed ?
                            if profile:
                                # The variant type is common to all samples
                                if v.var_type == "snp":
                                    if v.is_transition:
                                        mutation: str = "transition"
                                    else:
                                        mutation: str = "transversion"
                                    stats["variant"][v.var_type][mutation] += 1
                                    stats["variant"][v.var_type][v.REF][v.ALT[0]] += 1
                                elif v.var_type == "indel":
                                    if v.is_deletion:
                                        mutation: str = "deletion"
                                    else:
                                        mutation: str = "insertion"
                                    stats["variant"][v.var_type][mutation] += 1
                                else:
                                    stats["variant"][v.var_type] += 1

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
                                                    lambda sample: utils.is_homozygous(
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
                                                    lambda sample: utils.is_heterozygous(
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