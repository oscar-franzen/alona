"""
 alona

 Description:
 An analysis pipeline for scRNA-seq data.

 How to use:
 https://github.com/oscar-franzen/alona/

 Contact:
 Oscar Franzen <p.oscar.franzen@gmail.com>
"""

import logging
import time
import click

from .utils import *
from .log import (init_logging, log_debug, log_error, log_info)
from .logo import show_logo
from .exceptions import *
from .alonabase import AlonaBase
from .cell import AlonaCell
from .constants import GENOME

@click.command()
@click.argument('filename', type=click.Path(exists=True))
@click.option('-o', '--output', help='Specify name of output directory')
@click.option('-df', '--dataformat', help='How the input data have been processed. For \
example, if data are raw read counts, then select "raw".',
              type=click.Choice(['raw', 'rpkm', 'log2']), default='raw', show_default=True)

@click.option('-mr', '--minreads', help='Minimum number of reads per cell to keep the \
cell.', default=1000, show_default=True)
@click.option('-mg', '--minexpgenes', help='Minimum number of expressed genes as \
percent of all cells, i.e. genes expressed in fewer cells than this are removed.',
              default=0.01, show_default=True)

@click.option('--qc_auto', help='Automatic filtering of low quality cells.',
              type=click.Choice(['yes', 'no']), default='yes', show_default=True)

@click.option('--mrnafull', help='Data come from a full-length protocol, such as \
SMART-seq2.', is_flag=True, show_default=True)

@click.option('-d', '--delimiter', help='Data delimiter. The character used to separate \
data values. The default setting is to autodetect this character.',
              type=click.Choice(['auto', 'tab', 'space']), default='auto',
              show_default=True)

@click.option('-h', '--header', help='Data has a header line. The default setting is to \
autodetect if a header is present or not.', type=click.Choice(['auto', 'yes', 'no']),
              default='auto', show_default=True)

@click.option('-m', '--nomito', help='Exclude mitochondrial genes from analysis.',
              is_flag=True, show_default=True)

@click.option('--hvg', help='Method to use for identifying highly variable genes.',
              type=click.Choice(['seurat', 'Brennecke2013', 'scran', 'Chen2016',
              'M3Drop_smartseq2', 'M3Drop_UMI']), default='seurat', show_default=True)

@click.option('--hvg_n', help='Number of top highly variable genes to use.',
              default=1000, show_default=True)

@click.option('--nn_k', help='k in the nearest neighbour search.',
              default=10, show_default=True)

@click.option('--prune_snn', help='Threshold for pruning the SNN graph, i.e. the edges \
with lower value (Jaccard index) than this will be removed. Set to 0 to disable \
pruning. Increasing this value will result in fewer edges in the graph.',
              default=0.067, show_default=True)

@click.option('--leiden_partition', help='Partitioning algorithm to use. Can be \
RBERVertexPartition or ModularityVertexPartition.',
              default='RBERVertexPartition',
              type=click.Choice(['RBERVertexPartition', 'ModularityVertexPartition']),
              show_default=True)

@click.option('--leiden_res', help='Resolution parameter for the Leiden algorithm\
 (0-1).', default=0.8, show_default=True)

@click.option('--ignore_small_clusters', help='Ignore clusters with fewer or equal to N \
cells.', default=10, show_default=True)

@click.option('--embedding', help='Method used for data projection. Can be either tSNE or \
UMAP.', default='tSNE', type=click.Choice(['tSNE', 'UMAP']), show_default=True)

@click.option('--perplexity', help='The perplexity parameter in the t-SNE algorithm.',
              default=30, show_default=True)

@click.option('-s', '--species', help='Species your data comes from.',
              type=click.Choice(['human', 'mouse']), default='mouse', show_default=True)

@click.option('--dark_bg', help='Use dark background in scatter plots.', is_flag=True,
              default=False, show_default=True)

@click.option('-lf', '--logfile', help='Name of log file. Set to /dev/null if you want to \
disable logging to a file.', default='alona.log', show_default=True)

@click.option('-ll', '--loglevel', help='Set how much runtime information is written to \
the log file.', type=click.Choice(['regular', 'debug']), default='regular',
show_default=True)

@click.option('-n', '--nologo', help='Hide the logo.', is_flag=True,
              default=False)

@click.option('--seed', help='Set seed to get reproducible results.', type=int,
show_default=True)

@click.option('--version', help='Display version number.', is_flag=True,
              callback=print_version)

def run(filename, output, dataformat, minreads, minexpgenes, qc_auto, mrnafull, delimiter,
        header, nomito, hvg, hvg_n, nn_k, prune_snn, leiden_partition, leiden_res,
        ignore_small_clusters, embedding, perplexity, species, dark_bg,
        logfile, loglevel, nologo, seed, version):

    # confirm the genome reference files can be found
    for item in GENOME:
        path = get_alona_dir() + GENOME[item]

        if not os.path.exists(path):
            print('genome directory is not complete. Cannot find file: "%s". Path \
tried: %s' % (GENOME[item], path))
            sys.exit(1)

    time_start = time.time()
    init_logging(loglevel, logfile)

    log_debug('starting alona with %s' % filename)
    show_logo(nologo)

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
        'nn_k' : nn_k,
        'prune_snn' : prune_snn,
        'dark_bg' : dark_bg,
        'perplexity' : perplexity,
        'leiden_partition' : leiden_partition,
        'leiden_res' : leiden_res,
        'ignore_small_clusters' : ignore_small_clusters,
        'hvg_method' : hvg,
        'hvg_cutoff' : hvg_n,
        'qc_auto' : qc_auto,
        'embedding' : embedding,
        'seed' : seed
    }

    alonabase = AlonaBase(alona_opts)
    alonabase.prepare()

    alonacell = AlonaCell(alonabase)
    alonacell.load_data()
    alonacell.analysis()

    time_end = time.time()

    log_info('alona finished in %.2f minutes' % ((time_end - time_start)/60))
