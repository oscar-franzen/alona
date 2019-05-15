import sys
from .utils import is_platform_linux

NAME = 'alona'

# Check Python version
if sys.version_info < (3, 0):
    print('alona requires at least Python 3.0.')
    sys.exit(1)

is_platform_linux()
