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

from .log import *
from .logo import *
from .exceptions import *
from .alona_base import alona_base
from ._version import __version__

def is_inside_container():
    """ Checks that alona.py is running inside the Docker container. """
    if not os.environ.get('ALONA_INSIDE_DOCKER',False): return False

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
@click.option('--delimiter', help='Data delimiter.',
              type=click.Choice(['auto', 'tab','space']), default='auto')
@click.option('--header', help='Data has a header line.',
              type=click.Choice(['auto', 'yes','no']), default='auto')
@click.option('--nomito', help='Exclude mitochondrial genes from analysis.',is_flag=True)
@click.option('--species', help='Species your data comes from.',
              type=click.Choice(['human', 'mouse']),default='mouse')
@click.option('--loglevel', help='Set how much runtime information is written to \
               the log file.', type=click.Choice(['regular', 'debug']), default='regular')
@click.option('--nologo', help='Hide the logo.', is_flag=True)
@click.option('--version', help='Display version number.', is_flag=True,
              callback=print_version)

def run(filename, output, delimiter, header, nomito, species, loglevel, nologo, version):
    init_logging(loglevel)
    
    logging.debug('starting alona with %s' % filename)
    show_logo(nologo)
    
    #if not is_inside_container():
        #sys.exit("alona requires several dependencies and should be run through a \
#Docker container. See https://raw.githubusercontent.com/oscar-franzen/alona/master/\
#README.md")
    
    alona_opts = {
        'input_filename' : filename,
        'output_directory' : output,
        'species' : species,
        'delimiter' : delimiter,
        'loglevel' : loglevel,
        'header' : header,
        'nomito' : nomito
    }
    
    ab = alona_base(alona_opts)
    
    try:
        ab.is_file_empty()
    except file_empty_error:
        log_error(None,'Input file is empty.')
    
    ab.create_work_dir()

    try:
        ab.unpack_data()
    except (invalid_file_format,
            file_corrupt,
            input_not_plain_text) as e:
        log_error(None,e)
    except Exception as e:
        logging.error(e)
        raise
    
    ab.get_delimiter()
    ab.has_header()
    
    try:
        ab.sanity_check_columns()
    except (irregular_column_count) as e:
        log_error(None,e)
        
    try:
        ab.sanity_check_genes()
    except (irregular_gene_count) as e:
        log_error(None,e)
        
    if species == 'human':
        ab.ortholog_mapper()
        
    try:
        ab.sanity_check_gene_dups()
    except (gene_duplicates) as e:
        log_error(None,e)
        
    ab.cleanup()
    
    logging.debug('finishing alona')
