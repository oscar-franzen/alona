import sys

name = 'alona'

# Check Python version
if sys.version_info < (3, 6):
    print('alona requires at least Python 3.6.')
    sys.exit(1)
