"""
 alona
 
 Description:
 An analysis pipeline for scRNA-seq data.
 
 How to use:
 https://github.com/oscar-franzen/alona/
 
 Details:
 https://alona.panglaodb.se/
 
 Contact:
 Oscar Franzen, <p.oscar.franzen@gmail.com>
"""

import os
import sys
import glob
import click
import logging

from .logo import *
from .alona_base import alona_base
from .exceptions import *

from ._version import __version__

def init_logging(loglevel):
    if loglevel == 'regular':
        _ll = logging.INFO
    elif loglevel == 'debug':
        _ll = logging.DEBUG
    
    logging.basicConfig(filename='alona.log',
                        level=_ll,
                        format='%(asctime)s %(message)s')

    return

def is_inside_container():
    """ Checks that alona.py is running inside the Docker container. """
    if not os.environ.get('ALONA_INSIDE_DOCKER',False): return False
    
    return

def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
        
    print('alona version %s ' % __version__)
    ctx.exit()
    
def is_platform_ok():
    assert(not 'linux' in sys.platform), 'alona is developed on Linux and everything else \
is untested.'
        

@click.command()
@click.argument('filename', type=click.Path(exists=True))
@click.option('--output', help='Specify name of output directory')
@click.option('--loglevel', help='Set how much runtime information is written to \
               the log file.', type=click.Choice(['regular', 'debug']), default='regular')
@click.option('--nologo', help='Hide the logo.', is_flag=True)
@click.option('--version', help='Display version number.', is_flag=True, callback=print_version)

def run(filename, output, loglevel, nologo, version):
    init_logging(loglevel)
    
    logging.info('starting alona with %s' % filename)
    show_logo(nologo)
    
    #if not is_inside_container():
        #sys.exit("alona requires several dependencies and should be run through a \
#Docker container. See https://raw.githubusercontent.com/oscar-franzen/alona/master/\
#README.md")
    
    alona_opts = {
        'input_filename' : filename,
        'output_directory' : output
    }
    
    with alona_base(alona_opts) as data:
        try:
            data.validate()
        except file_empty_error:
            print('ERROR: Input file is empty')
            raise
        
        data.create_work_dir()
    
        try:
            data.unpack_data()
        except (invalid_file_format,file_corrupt) as e:
            print(e)
            raise
        except Exception as e:
            print(e)
            raise
    
    logging.info('finishing alona')
