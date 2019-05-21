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
import time
import inspect
import click

from .utils import *
from .log import (init_logging, log_debug, log_error)
from .logo import show_logo
from .exceptions import *
from .alonabase import AlonaBase
from .cell import AlonaCell
from .constants import GENOME

@click.command()
@click.argument('filename', type=click.Path(exists=True))
@click.option('-o', '--output', help='Specify name of output directory')
@click.option('-df', '--dataformat', help='Data format.\n(raw read counts, rpkm, log2 \
normalized data). Default: raw', type=click.Choice(['raw', 'rpkm', 'log2']), default='raw')

@click.option('-mr', '--minreads', help='Minimum number of reads per cell to keep the \
cell. Default: 1000', default=1000)
@click.option('-mg', '--minexpgenes', help='Minimum number of expressed genes as \
percent of all cells, i.e. genes expressed in fewer cells than this are removed. \
Default: 0.01', default=0.01)

@click.option('--mrnafull', help='Data come from a full-length protocol, such as \
SMART-seq2.', is_flag=True)

@click.option('-d', '--delimiter', help='Data delimiter. The character used to separate data\
values. Cannot be a mix. Default: auto',
              type=click.Choice(['auto', 'tab', 'space']), default='auto')
@click.option('-h', '--header', help='Data has a header line. Default: auto',
              type=click.Choice(['auto', 'yes', 'no']), default='auto')
@click.option('-m', '--nomito', help='Exclude mitochondrial genes from analysis.',
              is_flag=True)
              
@click.option('--nn_k', help='Minimum number of reads per cell to keep the \
cell. Default: 1000', default=10)

@click.option('-s', '--species', help='Species your data comes from. Default: mouse',
              type=click.Choice(['human', 'mouse']), default='mouse')
@click.option('--cleanup', help='Perform cleanup of temporary files.',
              is_flag=True)
@click.option('-lf', '--logfile', help='Name of log file. Set to /dev/null if you want to \
disable logging to a file. Default: alona.log', default='alona.log')
@click.option('-ll', '--loglevel', help='Set how much runtime information is written to \
the log file. Default: regular', type=click.Choice(['regular', 'debug']), default='regular')
@click.option('-n', '--nologo', help='Hide the logo.', is_flag=True)
@click.option('--version', help='Display version number.', is_flag=True,
              callback=print_version)
def run(filename, output, dataformat, minreads, minexpgenes, mrnafull, delimiter, header,
        nomito, species, cleanup, logfile, loglevel, nologo, version):

    # confirm the genome reference files can be found
    for item in GENOME:
        path = os.path.dirname(inspect.getfile(AlonaBase)) + '/' + GENOME[item]

        if not os.path.exists(path):
            print('genome directory is not complete. Cannot find file: "%s". I \
tried this path: %s' % (GENOME[item], path))
            sys.exit(1)

    time_start = time.time()
    init_logging(loglevel, logfile)

    log_debug('starting alona with %s' % filename)
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
        'minexpgenes' : minexpgenes,
        'mrnafull' : mrnafull,
        'cleanup' : cleanup
    }

    alonabase = AlonaBase(alona_opts)
    alonabase.prepare()

    alonacell = AlonaCell(alonabase)
    alonacell.load_data()
    alonacell.analysis()

    alonabase.cleanup()
    time_end = time.time()

    log_debug('alona finished in %.2f minutes' % ((time_end - time_start)/60))
