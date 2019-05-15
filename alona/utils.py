import os
import sys

from ._version import __version__
from .log import log_error

def is_inside_container():
    """ Checks that alona.py is running inside the Docker container. """
    if not os.environ.get('ALONA_INSIDE_DOCKER', False): return False

    return True

def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    print('alona version %s ' % __version__)
    ctx.exit()

def is_platform_linux():
    if sys.platform != 'linux':
        print('alona is developed on Linux and all other systems are untested.')
        sys.exit(1)
