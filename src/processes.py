from collections.abc import Callable
import concurrent.futures
import errors
import files
from itertools import chain
from loguru import logger
from psutil import cpu_count
import re

class TasksManager():
    """ 
    A class to manage tasks

    Args:
        vcfs: The VCF repository
        processes: The number of processes to use
    """

    def __init__(self, vcfs: files.VCFRepository, processes: int = 0):

        assert processes >= 0, "Number of processes alocated must be an positive unsigned integer."
        
        if processes:

            if processes < 0:

                raise ValueError("Number of processes alocated must be an positive unsigned integer.")

            try:

                with open("/proc/self/status", mode="r") as f:

                    m = re.search(r"(?m)^Cpus_allowed:\s*(.*)$", f.read())

                    if m:

                        cpuset = bin(int(m.group(1).replace(",", ""), base=16)).count("1")

                        if cpuset > 0:

                            CPUS = min(cpuset, cpu_count(logical=True))
            
            except OSError:

                CPUS = cpu_count(logical=True)

            logger.debug(f"Number of logical(s) CPU(s) available: {CPUS}.")

            if CPUS > 1:
                self.processes = min(processes, CPUS)
                logger.debug(f"Tasks will be parallelized on {self.processes} CPU(s).")
            else:
                self.processes = 0
                logger.debug("Computation will be carried out sequentially.")

        else:
        
            self.processes: int = processes

        self.vcfs: files.VCFRepository = vcfs

        self.tasks = None

        self.results: dict = {}

    def scheduling(self, tasks: list[list[object]]):

        if len(tasks) != len(self.vcfs):

            raise ValueError(f"{len(self.vcfs)} collection(s) of tasks is/are required.")

        self.tasks: dict[files.VCF:list] = {vcf: [(vcf, task) for task in tasks[i]] for i, vcf in enumerate(self.vcfs.repository)}

        logger.debug(f"Task(s) scheduled {self.flatten()}")

    def flatten(self):

        return list(chain.from_iterable(self.tasks.values()))

    def commit(self, job: Callable, jobargs: dict[str:object], initializer: Callable = None, initargs: dict[str:object] = None):

        if self.results:

            self.results.clear()

        if initializer:

            for vcf in self.vcfs.repository:

                try:

                    initializer(vcf, *initargs)

                except Exception as e:

                    raise errors.ProcessError(e)

        if self.processes:
            # Process the chromosomes concurrently, using as much process as available
            # use a context manager to avoid memory leaks
            with concurrent.futures.ProcessPoolExecutor(
                max_workers=self.processes
            ) as executor:
                # Submit the processes
                futures = {
                    executor.submit(
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

                self.results[task] = job(task, *jobargs)