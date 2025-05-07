class VariantRepository():

    def __init__(self, filters: dict = None):
        
        if any(filters):

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
class Chromosome():

    def __init__(self, length):
        
        self._length: int = length

    @property
    def length(self):

        return self._length
    
    @length.setter
    def length(self, value: int):

        self._length = value