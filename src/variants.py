from operator import itemgetter
import numpy
from pandas import Series, DataFrame, concat
from plots import PlotLibrary
import utils

class VariantRepository():

    def __init__(self, filters: dict = None):
        
        if filters and any(filters):

            # Filters to apply during the parsing
            self.FILTERS: dict = (
                                    {
                                        "exclude": {
                                            "exclude_snps": filters.get("SNP", False),
                                            "exclude_indels": filters.get("INDELS", False),
                                            "exclude_vars": filters.get("VARS", False),
                                            "exclude_mnps": filters.get("MNP", False),
                                            "exclude_transitions": filters.get("TRANSITION", False),
                                            "exclude_svs": filters.get("SV", False),
                                            "pass_only": filters.get("PASS_ONLY", False),
                                        },
                                    }
                                )
        else: 

            self.FILTERS: dict = None

        self.repository: dict = {}

        self.profile: dict = {}

        self.filtered: dict = {}

    @property
    def filters(self):

        return getattr(self, "FILTERS", None)
    
    @filters.setter
    def filters(self, value: dict):

        if any(value):

            self.FILTERS: dict = (
                                    {
                                        "exclude": {
                                            "exclude_snps": value.get("SNP", False),
                                            "exclude_indels": value.get("INDELS", False),
                                            "exclude_vars": value.get("VARS", False),
                                            "exclude_mnps": value.get("MNP", False),
                                            "exclude_transitions": value.get("TRANSITION", False),
                                            "exclude_svs": value.get("SV", False),
                                            "pass_only": value.get("PASS_ONLY", False),
                                        },
                                    }
                                )
            
    def update_repository(self, chromosome: str, variants: DataFrame = None, filtered: dict = None, profile: dict = None):

        if not (variants.empty):

            self.repository[chromosome] = variants

        if filtered:

            self.filtered[chromosome] =  filtered

        if profile:

            self.profile[chromosome] = profile

    def collapse(self) -> DataFrame:

        df: DataFrame = concat(list(itemgetter(*list(sorted(self.repository.keys())))(self.repository))).astype(
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

        self.plots.barplot(data, "Chromosome","Count", "Type", "Variant by chromosome", "VariantByChromosome")

        data = DataFrame([{"Chromosome": k, "Genotype": genotype, "Count": self.profile[k][code]} 
                        for k in chromosomes 
                        for genotype, code in [("Homozygous", "hom"), ("Heterozygous", "het")]])

        self.plots.barplot(data, "Chromosome", "Count", "Genotype", "Genotype by chromosome","GenotypeByChromosome")

        data = DataFrame([{"Chromosome": k, "SNP": snp, "Count": self.profile[k]["variant"]["snp"][snp]} for k in chromosomes for snp in ["transition","transversion"]])

        self.plots.barplot(data, "Chromosome", "Count", "SNP", "SNP Type by Chromosome", "SNPTypeByChrom")

        chromosome = list(chromosomes)[0]

        if self.profile[chromosome]["GQ"].size:

            pass

        if self.profile[chromosome]["depth"].size:

            with utils.suppress_warnings():       
                data = DataFrame(list(map(lambda k: Series([k, numpy.mean(self.profile[k]["depth"])], index=["Chromosome", "Depth"]), chromosomes)))
            
            self.plots.barplot(data, "Chromosome", "Depth", color="Chromosome", title="Mean depth by chromosome", prefix="DepthByChromosomeBarPlot")

            data = list(map(lambda k: DataFrame({"Chromosome": [k] * self.profile[k]["depth"].size,
                                    "Depth": self.profile[k]["depth"].flatten()}), chromosomes))

            df = concat(data, ignore_index=True).astype({"Chromosome": "category", "Depth": "int"})

            self.plots.boxplot(df, "Chromosome", "Depth", "Chromosome", "Depth by chromosome", "DepthByChromosomeBoxPlot")

            self.plots.histogram(df, "Depth", None, "Depth", "DepthHist")

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