class CompressionIndexError(OSError):
    """
    Exception raised when compressing or indexing a VCF file.
    This error is raised when the VCF file is not successfully compressed or indexed.
    """
    pass

class FileError(Exception):
    """
    Exception raised when an error about a file occurs.
    This error is raised when a file is not valid.
    """
    pass

class VCFError(FileError):
    """
    Exception raised when an error about a VCF file occurs.
    This error is raised when the VCF file is not valid.
    """
    pass

class IndexError(FileError):
    """
    Exception raised when an error about a VCF index occurs.
    This error is raised when the VCF index is not valid.
    """
    pass

class ProcessError(RuntimeError):
    """
    Exception raised when a error about a process occurs.
    This error is raised when a process fails.
    """
    pass

class ReportError(FileError):
    """
    Exception raised when a error about the HTML report occurs.
    This error is raised when the HTML report is not successfully created.
    """
    pass

class VariantError(Exception):
    """
    Exception raised when a error about a variant occurs.
    This error is raised when a variant is not valid.
    """
    pass