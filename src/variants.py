class VariantRepository():

    pass


class Chromosome():

    def __init__(self, length):
        
        self._length: int = length

    @property
    def length(self):

        return self._length
    
    @length.setter
    def length(self, value: int):

        self._length = value