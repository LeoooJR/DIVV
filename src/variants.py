from collections import deque, Counter
import errors
from functools import lru_cache
from loguru import logger
import numpy
from pandas import Series, DataFrame, concat
from plots import PlotLibrary
from sortedcontainers import SortedSet
from utils import suppress_warnings, convert

class VariantRepository():
    """ A class to store variants and perform operations on them """

    INDEX: dict[str:int] = {
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
    
    VT: list[str] = ["snp", "ins", "del", "mnp", "inv", "csv"]

    def __init__(self, filters: dict = None):
        
        if filters and any(filters):

            # Filters to apply during the parsing
            self.FILTERS: dict[str:bool] = {
                                    "exclude_snps": filters.get("SNP", False),
                                    "exclude_indels": filters.get("INDELS", False),
                                    "exclude_vars": filters.get("VARS", False),
                                    "exclude_mnps": filters.get("MNP", False),
                                    "exclude_transitions": filters.get("TRANSITION", False),
                                    "exclude_svs": filters.get("SV", False),
                                    "pass_only": filters.get("PASS_ONLY", False),
                                }
        else: 

            self.FILTERS: dict = None

        self.repository: dict = {}

        self.profile: dict = {}

        self._filtered: dict[str:Counter] = {}

        self._chromosomes: SortedSet = None

        self.seqlens: list = None

        self._format = None

    @staticmethod
    def chromosome_sort_key(item: str):
        """Key to be used for sorting chromosomes."""

        # If the item is a digit, return a tuple with 0 and the integer value
        # If the item is not a digit, return a tuple with 1 and the item itself
        # This will ensure that digits are sorted before non-digits

        item: str = item.removeprefix('chr')

        if item.isdigit():

            return (0, int(item))
        
        else:

            return (1, item)
        
    @property
    def chromosomes(self):

        return self._chromosomes
    
    @chromosomes.setter
    def chromosomes(self, values: set):

        self._chromosomes: SortedSet = SortedSet(iterable=values, key=VariantRepository.chromosome_sort_key)

    @property
    def filters(self):

        return getattr(self, "FILTERS", None)
    
    @filters.setter
    def filters(self, value: dict):

        if any(value):

            self.FILTERS: dict = {
                                    "exclude_snps": value.get("SNP", False),
                                    "exclude_indels": value.get("INDELS", False),
                                    "exclude_vars": value.get("VARS", False),
                                    "exclude_mnps": value.get("MNP", False),
                                    "exclude_transitions": value.get("TRANSITION", False),
                                    "exclude_svs": value.get("SV", False),
                                    "pass_only": value.get("PASS_ONLY", False),
                                }
            
    def is_filtered(self) -> bool:
        """ Check if any filter is applied """

        return any([self.FILTERS[filter] for filter in self.FILTERS])

    @property    
    def filtered(self):

        return sum((self._filtered[chrom] for chrom in self.chromosomes if chrom in self._filtered), Counter())
    
    @property
    def filtered_by_chromosome(self):

        return getattr(self, "_filtered", None)
            
    @property
    def format(self):

        return getattr(self, "_format", None)
    
    @format.setter
    def format(self, value: str|list[str]):

        if isinstance(value, str):
            self._format = value.split(':')
        elif isinstance(value, list):
            self._format = value
    
    @staticmethod
    def exclude(v: object, filters: dict = None) -> bool:
        """ Check if variant should be excluded from analysis """
        return (
                v.is_indel and filters["exclude_indels"]
            ) or (
                v.is_snp and filters["exclude_snps"]
            ) or (
                v.is_mnp and filters["exclude_mnps"]
            ) or (
                v.is_sv and filters["exclude_svs"]
            ) or (
                v.is_transition and filters["exclude_transitions"]
            ) or (
                (v.FILTER != None) and filters["pass_only"]
            ) if filters else False
    
    @staticmethod
    def is_homozygous(GT: str):
        """ Check if variant is homozygous """
        if isinstance(GT, str):
            if GT:
                if '/' in GT:
                    alleles = GT.split('/')
                elif '|' in GT:
                    alleles = GT.split('|')
                else:
                    return False
            
                return alleles[0] == alleles[1]
            
            return False
        else:
            raise errors.VariantError(f"The genotype must be a string instance. The provided genotype is an {type(GT)}")

    @staticmethod
    def is_heterozygous(GT: str):
        """ Check if variant is heterozygous """
        if isinstance(GT, str):
            if GT:
                if '/' in GT:
                    alleles = GT.split('/')
                elif '|' in GT:
                    alleles = GT.split('|')
                else:
                    return False
                
                return alleles[0] != alleles[1]
            
            return False
        else:
            raise errors.VariantError(f"The genotype must be a string instance. The provided genotype is an {type(GT)}")
        
    @staticmethod
    def is_composed_variant(alleles: list[str]) -> bool:
        """ Check if variant is composed """
        return len(alleles) > 1
        
    @staticmethod
    @lru_cache(maxsize=1000) # Use Least Recently Used (LRU) cache to store results, SNP are often repeated
    def get_variant_type(ref: str, alts: tuple[str]) -> list[str]:

        vts = []

        def is_snp(ref: str, alt: str) -> str:
            """Check if the variant is a SNP (Single Nucleotide Polymorphism)."""

            return "snp" if len(ref) == 1 and len(alt) == 1 else ''

        def is_ins(ref: str, alt: str) -> str:
            """Check if the variant is an INS (Insertion)."""
            
            return "ins" if len(ref) == 1 and len(alt) > 1 else ''

        def is_del(ref: str, alt: str) -> str:
            """Check if the variant is a DEL (Deletion)."""

            return "del" if len(ref) > 1 and len(alt) == 1 else ''

        def is_inv(ref: str, alt: str) -> str:
            """Check if the variant is an INV (Inversion)."""
            OLD_CHARS: str = "ACGTacgt"
            REPLACE_CHARS: str = "TGCAtgca"
            rev: str = alt.translate(str.maketrans(OLD_CHARS,REPLACE_CHARS))[::-1]
            return "inv" if len(ref) == len(alt) and ref == rev else ''

        def is_mnv(ref: str, alt: str) -> str:
            """Check if the variant is a MNV (Multi Nucleotide Variant)."""

            OLD_CHARS: str = "ACGTacgt"
            REPLACE_CHARS: str = "TGCAtgca"

            rev: str = alt.translate(str.maketrans(OLD_CHARS,REPLACE_CHARS))[::-1]

            return "mnp" if len(ref) == len(alt) and ref != rev else ''
        
        # List of functions to check variant type
        # Deque is used to pop functions from the left at each iteration
        # until one of them returns a value
        # If no function returns a value, the variant type is set to CSV (complex structural variant)
        funcs: deque = deque([is_snp, is_ins, is_del, is_inv, is_mnv])
        
        for alt in alts:

            # Empty string to store variant type
            # Empty string return False
            variant_type: str = ''

            try:
                
                # Loop until one of the functions returns a value
                # Empty string equivalent to False
                while not variant_type:

                    variant_type: str = funcs.popleft()(ref, alt)

            # IndexError is raised when all functions have been popped from the deque
            # and none of them returned a value
            # In this case, the variant type is set to CSV (complex structural variant)
            # This is a fallback to avoid infinite loop
            except IndexError:

                logger.warning(f"Could not determine variant type from {ref}:{alt}")

                variant_type: str = "csv"
            
            vts.append(variant_type)

        return vts
        
    @staticmethod
    def format_to_values(format: str, values: str|list[str]) -> dict:
        """ map FORMAT string to respective SAMPLE values """

        # Split the format string into a list of fields
        format: list[str] = format.split(":")

        # VCF is composed of multiples samples
        if isinstance(values,list):
            # Split the values string into a lists of list of values,
            # Each list of values represents a sample
            values: list[list[str]] = list(map(lambda x: x.split(":"), values))
            # Return a dictionary of the format and values mapped to each sample
            return {f"sample{s}": {f: convert(v)} for s in range(len(values)) for f, v in zip(format,values[s])}

        # VCF is composed of a unique sample
        else:
            values: list[str] = values.split(":")

            return {f: convert(v) for f, v in zip(format,values)}
            
    def update_repository(self, chromosome: str, variants: DataFrame = None, filtered: dict = None, profile: dict = None):
        """ Update the repository with the new variants """
        if not (variants.empty):

            self.repository[chromosome] = variants

        if filtered:

            self._filtered[chromosome] = filtered

        if profile:

            self.profile[chromosome] = profile

    def collapse(self) -> DataFrame:
        """ Collapse variants into a single DataFrame """

        df: DataFrame = concat([self.repository[chrom] for chrom in self.chromosomes if chrom in self.repository]).astype(
                                                                                                                            {
                                                                                                                                "Chromosome": "category",
                                                                                                                                "Position": "int",
                                                                                                                                "Type": "category",
                                                                                                                                "Filter": "category",
                                                                                                                                "Variant": "string[pyarrow]",
                                                                                                                            }
                                                                                                                        )
        df["Type"] = df["Type"].cat.set_categories(VariantRepository.VT)

        return df

    def visualization(self):
        """ Create plots from metrics """

        self.plots: PlotLibrary = PlotLibrary()

        data: list[dict] = []

        chromosomes: list[str] = list(self.profile.keys())

        for k in chromosomes:

            data.append({"Chromosome": k, "Type": "Indel", "Count": self.profile[k]["variant"]["indel"]["deletion"] + self.profile[k]["variant"]["indel"]["insertion"]})
            data.append({"Chromosome": k, "Type": "SNP", "Count": self.profile[k]["variant"]["snp"]["transition"] + self.profile[k]["variant"]["snp"]["transversion"]})
            data.append({"Chromosome": k, "Type": "CSV", "Count": self.profile[k]["variant"]["csv"]})

        self.plots.barplot(data, "Chromosome","Count", "Type", "Variant by Chromosome", "VariantByChromosome")

        data: DataFrame = DataFrame([{"Chromosome": k, "Genotype": genotype, "Count": self.profile[k][code]} 
                        for k in chromosomes 
                        for genotype, code in [("Homozygous", "hom"), ("Heterozygous", "het")]])

        self.plots.barplot(data, "Chromosome", "Count", "Genotype", "Genotype by Chromosome","GenotypeByChromosome")

        data: DataFrame = DataFrame([{"Chromosome": k, "SNP": snp, "Count": self.profile[k]["variant"]["snp"][snp]} for k in chromosomes for snp in ["transition", "transversion"]])

        # Are there any SNPs in the VCF?
        if data["Count"].values.any():

            self.plots.barplot(data, "Chromosome", "Count", "SNP", "SNP Type by Chromosome", "SNPTypeByChrom")

        chromosome: str = list(chromosomes)[0]
        
        # Are there any GQ collected values?
        if self.profile[chromosome]["GQ"].size:

            pass

        # Are there any depth collected values?
        if self.profile[chromosome]["depth"].size:

            with suppress_warnings():       
                data: DataFrame = DataFrame(list(map(lambda k: Series([k, numpy.mean(self.profile[k]["depth"])], index=["Chromosome", "Depth"]), chromosomes)))
            
            self.plots.barplot(data, "Chromosome", "Depth", color="Chromosome", title="Mean Depth by Chromosome", prefix="DepthByChromosomeBarPlot")

            data: list[DataFrame] = list(map(lambda k: DataFrame({"Chromosome": [k] * self.profile[k]["depth"].size,
                                    "Depth": self.profile[k]["depth"].flatten()}), chromosomes))

            df: DataFrame = concat(data, ignore_index=True).astype({"Chromosome": "category", "Depth": "int"})

            self.plots.boxplot(df, "Chromosome", "Depth", "Chromosome", "Depth by Chromosome", "DepthByChromosomeBoxPlot")

            self.plots.histogram(df, "Depth", None, "Depth distribution", "DepthHist")

            data.clear()

        if self.profile[chromosome]["quality"].size:

            pass

    def __len__(self):

        return len(self.repository)
    
class Chromosome():

    def __init__(self, name: str, length: int):

        self.name: str = name
        
        self._length: int = length

    @property
    def length(self):

        return self._length
    
    @length.setter
    def length(self, value: int):

        self._length = value