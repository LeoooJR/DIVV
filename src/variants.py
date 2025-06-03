import errors
import numpy
from pandas import Series, DataFrame, concat
from plots import PlotLibrary
from sortedcontainers import SortedSet
from utils import suppress_warnings, convert

class VariantRepository():

    INDEX = {
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

    def __init__(self, filters: dict = None):
        
        if filters and any(filters):

            # Filters to apply during the parsing
            self.FILTERS: dict = {
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

        self.filtered: dict = {}

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
    def format_to_values(format: str, values: str|list[str]) -> dict:
        """ map FORMAT string to respective SAMPLE values """

        # Split the format string into a list of fields
        format = format.split(":")

        # VCF is composed of multiples samples
        if isinstance(values,list):
            # Split the values string into a lists of list of values,
            # Each list of values represents a sample
            values = list(map(lambda x: x.split(":"), values))
            # Return a dictionary of the format and values mapped to each sample
            return {f"sample{s}": {f: convert(v)} for s in range(len(values)) for f, v in zip(format,values[s])}

        # VCF is composed of a unique sample
        else:
            values = values.split(":")

            return {f: convert(v) for f, v in zip(format,values)}
            
    def update_repository(self, chromosome: str, variants: DataFrame = None, filtered: dict = None, profile: dict = None):

        if not (variants.empty):

            self.repository[chromosome] = variants

        if filtered:

            self.filtered[chromosome] =  filtered

        if profile:

            self.profile[chromosome] = profile

    def collapse(self) -> DataFrame:

        df: DataFrame = concat([self.repository[chrom] for chrom in self.chromosomes if chrom in self.repository]).astype(
                                                                                                                            {
                                                                                                                                "Chromosome": "category",
                                                                                                                                "Position": "int",
                                                                                                                                "Type": "category",
                                                                                                                                "Filter": "category",
                                                                                                                                "Variant": "string[pyarrow]",
                                                                                                                            }
                                                                                                                        )
        df["Type"] = df["Type"].cat.set_categories(["snp", "indel", "sv"])

        return df

    def visualization(self):

        self.plots = PlotLibrary()

        data = []

        chromosomes = self.profile.keys()

        for k in chromosomes:

            data.append({"Chromosome": k, "Type": "Indel", "Count": self.profile[k]["variant"]["indel"]["deletion"] + self.profile[k]["variant"]["indel"]["insertion"]})
            data.append({"Chromosome": k, "Type": "SNP", "Count": self.profile[k]["variant"]["snp"]["transition"] + self.profile[k]["variant"]["snp"]["transversion"]})
            data.append({"Chromosome": k, "Type": "SV", "Count": self.profile[k]["variant"]["sv"]})

        self.plots.barplot(data, "Chromosome","Count", "Type", "Variant by Chromosome", "VariantByChromosome")

        data = DataFrame([{"Chromosome": k, "Genotype": genotype, "Count": self.profile[k][code]} 
                        for k in chromosomes 
                        for genotype, code in [("Homozygous", "hom"), ("Heterozygous", "het")]])

        self.plots.barplot(data, "Chromosome", "Count", "Genotype", "Genotype by Chromosome","GenotypeByChromosome")

        data = DataFrame([{"Chromosome": k, "SNP": snp, "Count": self.profile[k]["variant"]["snp"][snp]} for k in chromosomes for snp in ["transition","transversion"]])

        self.plots.barplot(data, "Chromosome", "Count", "SNP", "SNP Type by Chromosome", "SNPTypeByChrom")

        chromosome = list(chromosomes)[0]

        if self.profile[chromosome]["GQ"].size:

            pass

        if self.profile[chromosome]["depth"].size:

            with suppress_warnings():       
                data = DataFrame(list(map(lambda k: Series([k, numpy.mean(self.profile[k]["depth"])], index=["Chromosome", "Depth"]), chromosomes)))
            
            self.plots.barplot(data, "Chromosome", "Depth", color="Chromosome", title="Mean Depth by Chromosome", prefix="DepthByChromosomeBarPlot")

            data = list(map(lambda k: DataFrame({"Chromosome": [k] * self.profile[k]["depth"].size,
                                    "Depth": self.profile[k]["depth"].flatten()}), chromosomes))

            df = concat(data, ignore_index=True).astype({"Chromosome": "category", "Depth": "int"})

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