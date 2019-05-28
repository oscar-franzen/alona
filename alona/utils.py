import os
import sys
import random
import inspect

from ._version import __version__
from .log import log_error

def is_inside_container():
    """ Checks that alona.py is running inside the Docker container. """
    if os.environ.get('ALONA_INSIDE_DOCKER', False):
        return True
    else:
        return False

def print_version(ctx, param, value):
    """ Prints alona version number. """
    if not value or ctx.resilient_parsing:
        return
    print('alona version %s\n\nhttps://alona.panglaodb.se\
    \nhttps://github.com/oscar-franzen/alona' % __version__)
    ctx.exit()

def is_platform_linux():
    """ Checks if we run Linux """
    if sys.platform != 'linux':
        print('alona is developed on Linux and all other systems are untested.')
        sys.exit(1)

def get_alona_dir():
    """ Returns the alona base directory. """
    return os.path.dirname(inspect.getfile(log_error)) + '/'

def get_random_color(pastel_fac=0.5):
    r = random.uniform
    return [(x+pastel_fac)/(1.0+pastel_fac) for x in [r(0, 1.0) for i in [1, 2, 3]]]

def color_distance(c1, c2):
    return sum([abs(x[0]-x[1]) for x in zip(c1, c2)])

def generate_new_color(existing_colors,pastel_fac = 0.5):
    max_distance = None
    best_color = None
    for i in range(0, 100):
        color = get_random_color(pastel_fac = pastel_fac)
        if not existing_colors:
            return color
        best_distance = min([color_distance(color, c) for c in existing_colors])
        if not max_distance or best_distance > max_distance:
            max_distance = best_distance
            best_color = color
    return best_color
