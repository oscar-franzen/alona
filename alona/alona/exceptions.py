class AlonaError(Exception):
    """Base exception for alona errors"""

class FileEmptyError(AlonaError):
    pass

class InvalidFileFormatError(AlonaError):
    pass

class FileCorruptError(AlonaError):
    pass

class InputNotPlainTextError(AlonaError):
    pass

class IrregularColumnCountError(AlonaError):
    pass

class IrregularGeneCountError(AlonaError):
    pass

class GeneDuplicatesError(AlonaError):
    pass
    
class TooFewGenesMappableError(AlonaError):
    pass
