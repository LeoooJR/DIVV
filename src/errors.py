class CompressionIndexError(OSError):
    
    pass

class FileError(ValueError):

    pass

class VCFError(FileError):

    pass

class IndexError(FileError):

    pass

class ProcessError(RuntimeError):

    pass