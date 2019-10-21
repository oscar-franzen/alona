import sys

name = 'alona'

__version__ = '0.3.4'
__license__ = 'GPLv3'
__author__ = 'Oscar Franzén'
__email__ = 'p.oscar.franzen@gmail.com'

# Check Python version
if sys.version_info < (3, 6):
    print('alona requires at least Python 3.6.')
    sys.exit(1)
