class alona_error(Exception):
    """Base exception for alona errors"""
    
    pass

class file_empty_error(alona_error):
    pass
    
class invalid_file_format(alona_error):
    pass

class file_corrupt(alona_error):
    pass

class input_not_plain_text(alona_error):
    pass

class irregular_column_count(alona_error):
    pass
