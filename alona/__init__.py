import sys
from .utils import is_platform_linux

name = 'alona'

# Check Python version
if sys.version_info < (3, 6):
    print('alona requires at least Python 3.6.')
    sys.exit(1)

is_platform_linux()
