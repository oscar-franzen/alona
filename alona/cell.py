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

import sys
import os
import re
import logging
import pandas as pd
import numpy as np
from joblib import dump, load

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.pyplot import figure

from .log import (log_info, log_debug, log_error)
from .constants import (OUTPUT, ORANGE, GENOME)
from .clustering import AlonaClustering
from .celltypes import AlonaCellTypePred
from .utils import get_alona_dir

class AlonaCell():
    """
    Cell data class.
    """

    def __init__(self, alonabase):
        self.alonabase = alonabase
        self.data = None
        self.data_norm = None
        self.data_ERCC = None

        self._alona_clustering = AlonaClustering(self, alonabase.params)

        # make matplotlib more quiet
        logging.getLogger('matplotlib').setLevel(logging.WARNING)

    def _validate_counts(self):
        """ Basic data validation of gene expression values. """
        log_debug('Running _validate_counts()')
        data_format = self.alonabase.params['dataformat']
        
        if (np.any(self.data.dtypes != 'int64') and data_format == 'raw'):
            log_error(msg='Non-count values detected in data matrix while data format is \
set to raw read counts.')
        elif (np.any(self.data.dtypes == 'int64') and data_format == 'log2'):
            log_error(msg='Count values detected in data matrix while data format is \
set to log2.')
        elif self.alonabase.params['dataformat'] == 'log2':
            if np.any(self.data > 1000):
                log_error(msg='Data do not appear to be log2 transformed.')
        else:
            log_debug('_validate_counts() finished without a problem')
            
        log_debug('Finished _validate_counts()')

    def _remove_empty(self):
        """ Removes empty cells and genes """
        data_zero = self.data == 0
        
        # Axis 0 will act on all the ROWS in each COLUMN
        # Axis 1 will act on all the COLUMNS in each ROW
        cells = data_zero.sum(axis=0)
        genes = data_zero.sum(axis=1)

        total_genes = self.data.shape[0]
        total_cells = self.data.shape[1]

        if np.sum(cells == total_cells) > 0:
            log_info('%s empty cells will be removed' % (np.sum(cells == total_cells)))
            self.data = self.data[self.data.columns[np.logical_not(cells == total_cells)]]
        if np.sum(genes == total_genes) > 0:
            log_info('%s empty genes will be removed' % (np.sum(genes == total_genes)))
            self.data = self.data.loc[self.data.index[np.logical_not(genes == total_genes)]]
            
        log_info('Current dimensions: %s' % str(self.data.shape))

    def _remove_mito(self):
        """ Remove mitochondrial genes. """
        if self.alonabase.params['nomito']:
            log_debug('removing mitochondrial genes')

            mt_count = self.data.index.str.contains('^mt-', regex=True, case=False)

            if np.sum(mt_count) > 0:
                log_info('detected and removed %s mitochondrial gene(s)' % np.sum(mt_count))
                self.data = self.data.loc[self.data.index[np.logical_not(mt_count)]]

    def _read_counts_per_cell_filter(self):
        """ Generates a bar plot of read counts per cell. """
        if self.alonabase.params['dataformat'] == 'raw':
            min_reads = self.alonabase.params['minreads']
            cell_counts = self.data.sum(axis=0)

            plt.clf()

            figure(num=None, figsize=(5, 5))

            no_cells_remove = np.sum(np.array(sorted(cell_counts, reverse=True)) < \
                 min_reads)
            colors = [ORANGE]*(len(cell_counts)-no_cells_remove) + \
                     ['#999999']*no_cells_remove

            plt.bar(np.arange(len(cell_counts)), sorted(cell_counts, reverse=True),
                    color=colors)
            plt.ylabel('raw read counts')
            plt.xlabel('cells (sorted on highest to lowest)')
            plt.axhline(linewidth=1, color='r', y=min_reads)
            plt.legend(handles=[
                mpatches.Patch(color='#E69F00', label='Passed'),
                mpatches.Patch(color='#999999', label='Removed')
            ])
            plt.title('sequencing reads')
            plt.savefig(self.alonabase.get_wd() + \
                OUTPUT['FILENAME_BARPLOT_RAW_READ_COUNTS'], bbox_inches='tight')
            plt.close()

    def _genes_expressed_per_cell_barplot(self):
        """ Generates a bar plot of number of expressed genes per cell. """
        genes_expressed = self.data.apply(lambda x: sum(x > 0), axis=0)

        plt.clf()
        figure(num=None, figsize=(5, 5))
        plt.bar(np.arange(len(genes_expressed)), sorted(genes_expressed, reverse=True),
                color=[ORANGE]*len(genes_expressed))
        plt.ylabel('number of genes')
        plt.xlabel('cells (sorted on highest to lowest)')
        plt.title('expressed genes')

        plt.savefig(self.alonabase.get_wd() + \
            OUTPUT['FILENAME_BARPLOT_GENES_EXPRESSED'], bbox_inches='tight')
        plt.close()

    def _quality_filter(self):
        """
        Removes cells with too few reads (set by `--minreads`).
        Removes "underexpressed" genes (set by `--minexpgenes`).
        """
        if self.alonabase.params['dataformat'] == 'raw':
            min_reads = self.alonabase.params['minreads']
            cell_counts = self.data.sum(axis=0)
            log_info('Keeping %s out of %s cells' % (
                np.sum(cell_counts > min_reads), len(cell_counts)))

            self.data = self.data[self.data.columns[cell_counts > min_reads]]

            if self.data.shape[1] < 100:
                log_error(self.alonabase, msg='After removing cells with < %s reads, \
                   less than 100 reads remain. Please adjust --minreads' % min_reads)
        if self.alonabase.params['minexpgenes'] > 0:
            log_debug('Filtering genes based on --minexpgenes')

            genes_expressed = self.data.apply(lambda x: sum(x > 0)/len(x), axis=1)
            log_info('Removing %s genes' % (self.data.shape[0] - np.sum(genes_expressed > \
                self.alonabase.params['minexpgenes'])))
            self.data = self.data[genes_expressed > self.alonabase.params['minexpgenes']]

    def _print_dimensions(self):
        """ Prints the new dimensions after quality filtering. """
        log_info('Data dimensions after filtering: %s genes and %s cells' % \
        (self.data.shape[0], self.data.shape[1]))

    def _rpkm(self, data):
        """ Normalize expression values as RPKM. """
        """ doi:10.1186/s13059-016-0881-8 """
        log_debug('Normalizing data to RPKM')
        if self.alonabase.params['species'] == 'human':
            # TODO: Implement RPKM for human. This should be executed _before_ mapping
            # to mouse orthologs (16-May-2019).
            log_info('RPKM for "--species human" is not implemented at the moment.')
            raise NotImplementedError('RPKM for human is not implemented at the moment.')

        # the file contains _combined_ lengths of all exons of the particular gene
        exon_lengths = pd.read_csv(get_alona_dir()+GENOME['MOUSE_EXON_LENGTHS'], delimiter=' ',
                                   header=None)
        exon_lengths.columns = ['gene', 'length']

        # intersects and sorts
        temp = pd.merge(data, exon_lengths, how='inner', left_on=data.index,
                        right_on=exon_lengths['gene'].str.extract('^(.+)\.[0-9]+')[0])

        temp.index = temp['gene']
        exon_lengths = temp['length']
        temp = temp.drop(['gene', 'length'], axis=1)
        temp = temp[temp.columns[1:]]

        # gene length in kilobases
        kb = exon_lengths/1000

        def _foo(x):
            s = sum(x)/1000000
            rpm = x/s
            rpkm = rpm/kb

            return rpkm

        _q = temp.apply(_foo, axis=0) # 0 applying to each column

        data_norm = np.log2(_q+1)
        return data_norm

    @staticmethod
    def _dump(d, fn='foo.joblib'):
        """ Serialize data. """
        log_debug('Writing dump file %s' % fn)
        dump(d, fn)

    def _normalization(self, data, fn_out):
        """ Performs normalization of the gene expression values. """
        log_debug('Inside _normalization()')
        norm_mat_path = self.alonabase.get_wd() + fn_out

        if os.path.exists(norm_mat_path):
            log_debug('Loading data matrix from file')
            return load(norm_mat_path)

        if not self.alonabase.params['mrnafull'] and \
           self.alonabase.params['dataformat'] == 'raw':
            # Axis 0 will act on all the ROWS in each COLUMN
            # Axis 1 will act on all the COLUMNS in each ROW
            col_sums = data.apply(lambda x: sum(x), axis=0)
            data_norm = (data / col_sums) * 10000
            data_norm = np.log2(data_norm+1)
        elif self.alonabase.params['mrnafull'] and \
             self.alonabase.params['dataformat'] == 'raw':
            data_norm = self._rpkm(data)
        elif self.alonabase.params['mrnafull'] and \
             self.alonabase.params['dataformat'] == 'rpkm':
            log_debug('_normalization() Running log2')
            data_norm = np.log2(data+1)
        else:
            data_norm = data
            log_debug('Normalization is not needed.')

        self._dump(data_norm, norm_mat_path)
        
        log_debug('Finished _normalization()')
        
        return data_norm
        
    def _lift_ERCC(self):
        """ Moves ERCC (if present) to a separate matrix. """
        log_debug('Inside _lift_ERCC()')
        self.data_ERCC = self.data[self.data.index.str.contains('^ERCC-[0-9]+$')]
        self.data = self.data[np.logical_not(self.data.index.str.contains('^ERCC-[0-9]+$'))]
        log_debug('Finishing _lift_ERCC()')
        
    def _load_rRNA_genes(self):
        """ Detect rRNA genes. """
        log_debug('Inside _load_rRNA_genes()')
        # TODO: Add support for human rRNA genes (12 Jun 2019).
        rRNA_genes = {}
        with open(get_alona_dir() + GENOME['MOUSE_GENOME_ANNOTATIONS'], 'r') as fh:
            for line in fh:
                if re.search('"rRNA"',line):
                    gene = re.search('^gene_id "(.*?)";',line.split('\t')[8]).group(1)
                    rRNA_genes[gene] = 1
        self.rRNA_genes = rRNA_genes
        log_debug('Finished _load_rRNA_genes()')

    def load_data(self):
        """ Load expression matrix. """
        norm_mat_path = self.alonabase.get_wd() + '/normdata.joblib'
        if os.path.exists(norm_mat_path):
            self.data_norm = load(norm_mat_path)
            return

        log_debug('loading expression matrix')
        self.data = pd.read_csv(self.alonabase.get_matrix_file(),
                                delimiter=self.alonabase._delimiter,
                                header=0 if self.alonabase._has_header else None)

        if self.alonabase._has_gene_id_column_id or not self.alonabase._has_header:
            self.data.index = self.data[self.data.columns[0]]
            self.data = self.data.drop(self.data.columns[0], axis=1)

        # remove duplicate genes (another check)
        if np.any(self.data.index.duplicated(False)):
            log_info('%s duplicated genes still detected, removing them.' % np.sum(
                self.data.index.duplicated(False)))
            self.data = self.data.iloc[np.logical_not(self.data.index.duplicated(False))]

        # initialize data
        self._validate_counts()
        self._remove_empty()
        self._remove_mito()
        self._read_counts_per_cell_filter()
        self._genes_expressed_per_cell_barplot()
        self._quality_filter()
        self._lift_ERCC()
        self._load_rRNA_genes()
        self._print_dimensions()

        # normalize gene expression values
        #self._normalization()
        
        self.data_norm = self._normalization(self.data, '/normdata.joblib')
        self.data_ERCC = self._normalization(self.data_ERCC, '/normdata_ERCC.joblib')

        log_debug('done loading expression matrix')

    def analysis(self):
        """ Runs the analysis pipeline. """
        log_debug('Running analysis...')
        tsne_path = self.alonabase.get_wd() + OUTPUT['FILENAME_EMBEDDINGS']
        pca_path = self.alonabase.get_wd() + OUTPUT['FILENAME_PCA']

        if os.path.exists(tsne_path) and os.path.exists(pca_path):
            log_debug('Loading embeddings from file')
            self._alona_clustering.embeddings = pd.read_csv(tsne_path, header=None)
            self._alona_clustering.pca_components = pd.read_csv(pca_path, header=None)
        else:
            self._alona_clustering.find_variable_genes()
            self._alona_clustering.PCA(pca_path)
            self._alona_clustering.tSNE(tsne_path)

        fn = self.alonabase.get_wd() + OUTPUT['FILENAME_CELL_SCATTER_PLOT']

        self._alona_clustering.cluster()
        self._alona_clustering.cell_scatter_plot(fn, cell_type_obj=None)

        self.pred = AlonaCellTypePred(self.alonabase.get_wd(),
                                      self._alona_clustering,
                                      self)
        self.pred.median_exp()
        self.pred.load_markers()
        self.pred.CTA_RANK_F()

        fn = self.alonabase.get_wd() + OUTPUT['FILENAME_CELL_SCATTER_PLOT_W_CT_LABELS']
        self._alona_clustering.cell_scatter_plot(fn, cell_type_obj=self.pred)
