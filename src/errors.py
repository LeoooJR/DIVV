class CompressionIndexError(OSError):
    
    pass

class FileError(Exception):

    pass

class VCFError(FileError):

    pass

class IndexError(FileError):

    pass

class ProcessError(RuntimeError):

    pass

class ReportError(FileError):

    pass

class VariantError(Exception):

    pass