"""
 alona

 Description:
 An analysis pipeline for scRNA-seq data.

 How to use:
 https://github.com/oscar-franzen/alona/

 Details:
 https://alona.panglaodb.se/

 Contact:
 Oscar Franzen <p.oscar.franzen@gmail.com>
"""

import logging
import click
import time

from .utils import *
from .log import (init_logging, log_error)
from .logo import show_logo
from .exceptions import *
from .alonabase import AlonaBase
from .cell import Cell

@click.command()
@click.argument('filename', type=click.Path(exists=True))
@click.option('-o','--output', help='Specify name of output directory')
@click.option('-df','--dataformat', help='Data format.\n(raw read counts, rpkm, log2 \
normalized data). Default: raw',
              type=click.Choice(['raw', 'rpkm', 'log2']), default='raw')
              
@click.option('-mr','--minreads', help='Minimum number of reads per cell to keep the cell. Default: 1000',
              default=1000)

@click.option('-d','--delimiter', help='Data delimiter. The character used to separate data\
values. Cannot be a mix. Default: auto',
              type=click.Choice(['auto', 'tab', 'space']), default='auto')
@click.option('-h','--header', help='Data has a header line. Default: auto',
              type=click.Choice(['auto', 'yes', 'no']), default='auto')
@click.option('-m','--nomito', help='Exclude mitochondrial genes from analysis.', is_flag=True)
@click.option('-s','--species', help='Species your data comes from. Default: mouse',
              type=click.Choice(['human', 'mouse']), default='mouse')
@click.option('-lf','--logfile', help='Name of log file. Set to /dev/null if you want to \
disable logging to a file. Default: alona.log', default='alona.log')
@click.option('-ll','--loglevel', help='Set how much runtime information is written to \
the log file. Default: regular', type=click.Choice(['regular', 'debug']), default='regular')
@click.option('-n','--nologo', help='Hide the logo.', is_flag=True)
@click.option('--version', help='Display version number.', is_flag=True,
              callback=print_version)
def run(filename, output, dataformat, minreads, delimiter, header, nomito, species,
        logfile, loglevel, nologo, version):
    
    time_start = time.time()
    init_logging(loglevel, logfile)

    logging.debug('starting alona with %s', filename)
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
        'nomito' : nomito,
        'dataformat' : dataformat,
        'minreads' : minreads,
    }

    alonabase = AlonaBase(alona_opts)

    try:
        alonabase.is_file_empty()
    except FileEmptyError:
        log_error(None, 'Input file is empty.')

    alonabase.create_work_dir()

    try:
        alonabase.unpack_data()
    except (InvalidFileFormatError,
            FileCorruptError,
            InputNotPlainTextError) as err:
        log_error(None, err)
    except Exception as err:
        logging.error(err)
        raise

    alonabase.get_delimiter()
    alonabase.has_header()

    try:
        alonabase.sanity_check_columns()
    except (IrregularColumnCountError) as err:
        log_error(None, err)

    try:
        alonabase.sanity_check_genes()
    except (IrregularGeneCountError) as err:
        log_error(None, err)

    if species == 'human':
        alonabase.ortholog_mapper()

    try:
        alonabase.sanity_check_gene_dups()
    except (GeneDuplicatesError) as err:
        log_error(None, err)

    alonabase.load_mouse_gene_symbols()

    try:
        alonabase.map_input_genes()
    except (TooFewGenesMappableError) as err:
        log_error(None, err)

    alonabase.check_gene_name_column_id_present()

    cell = Cell(alonabase)
    cell.load_data()

    #alonabase.cleanup()
    time_end = time.time()

    logging.debug('alona finished in %.2f minutes' % ((time_end - time_start)/60))
