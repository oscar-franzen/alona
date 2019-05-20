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

import os
import logging
import pandas as pd
import numpy as np
from joblib import dump, load

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.pyplot import figure

from .log import (log_info, log_debug, log_error)
from .constants import *
from .analysis import AlonaAnalysis

class AlonaCell():
    """
    Cell data class.
    """

    def __init__(self, alonabase):
        self.alonabase = alonabase
        self.data = None
        self.data_norm = None

        self._alona_analysis = AlonaAnalysis(self)

        # make matplotlib more quiet
        logging.getLogger('matplotlib').setLevel(logging.WARNING)

    def _validate_counts(self):
        """ Basic data validation of gene expression values. """
        if (np.any(self.data.dtypes != 'int64') and
                self.alonabase.params['dataformat'] == 'raw'):
            log_error('Non-count values detected in data matrix while data format is \
set to raw read counts.')
        elif self.alonabase.params['dataformat'] == 'log2':
            if np.any(self.data > 1000):
                log_error('Data do not appear to be log2 transformed.')
        else:
            log_debug('_validate_counts() finished without a problem')

    def _remove_empty(self):
        """ Removes empty cells and genes """
        # TODO: fix variable names for cells and genes (16-May-2019)
        data_zero = self.data == 0
        cells = data_zero.sum(axis=1) # 1 = rows
        genes = data_zero.sum(axis=0) # 0 = columns

        total_genes = self.data.shape[0]
        total_cells = self.data.shape[1]

        if np.sum(cells == total_cells) > 0:
            log_info('%s empty cells will be removed' % (np.sum(cells == total_cells)))
            self.data = self.data.loc[self.data.index[np.logical_not(cells == total_cells)]]
        if np.sum(genes == total_genes) > 0:
            log_info('%s empty genes will be removed' % (np.sum(genes == total_genes)))
            self.data = self.data[self.data.columns[np.logical_not(genes == total_genes)]]

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
            plt.savefig(self.alonabase.get_working_dir() + '/plots/' + \
                 FILENAME_BARPLOT_RAW_READ_COUNTS, bbox_inches='tight')
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

        plt.savefig(self.alonabase.get_working_dir() + '/plots/' + \
                 FILENAME_BARPLOT_GENES_EXPRESSED, bbox_inches='tight')
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
                log_error(self.alonabase, 'After removing cells with < %s reads, \
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

    def _rpkm(self):
        if self.alonabase.params['species'] == 'human':
            # TODO: Implement RPKM for human. This should be executed _before_ mapping
            # to mouse orthologs (16-May-2019).
            log_info('RPKM for "--species human" is not implemented at the moment.')
            raise NotImplementedError('RPKM for human is not implemented at the moment.')
            
        log_debug('Normalizing data to RPKM')

        # the file contains _combined_ lengths of all exons of the particular gene
        exon_lengths = pd.read_csv(GENOME['MOUSE_EXON_LENGTHS'], delimiter=' ',
                                   header=None)
        exon_lengths.columns = ['gene', 'length']

        # intersects and sorts
        temp = pd.merge(self.data, exon_lengths, how='inner', left_on=self.data.index,
                        right_on=exon_lengths['gene'])
        temp.index = temp['gene']
        exon_lengths = temp['length']
        temp = temp.drop(['gene', 'length'], axis=1)
        temp = temp[temp.columns[1:]]

        # gene length in kilobases
        kb = exon_lengths/1000

        index = 0

        def _foo(x):
            global index
            index = index + 1

            s = sum(x)/1000000
            rpm = x/s
            rpkm = rpm/kb

            return rpkm

        _q = temp.apply(_foo, axis=0) # 0 applying to each column

        self.data_norm = np.log2(_q+1)
        
    @staticmethod
    def _dump(d, fn='foo.joblib'):
        log_debug('Writing dump file %s' % fn)
        dump(d, fn)

    def _normalization(self):
        """ Performs normalization of the gene expression values. """
        norm_mat_path = self.alonabase.get_working_dir() + '/normdata.joblib'
        
        if os.path.exists(norm_mat_path):
            log_debug('Loading data matrix from file')
            self.data_norm = load(norm_mat_path)
            return
        
        if not self.alonabase.params['mrnafull'] and \
           self.alonabase.params['dataformat'] == 'raw':
            log_debug('Starting normalization...')
            # Axis 0 will act on all the ROWS in each COLUMN
            # Axis 1 will act on all the COLUMNS in each ROW
            
            col_sums = self.data.apply(lambda x: sum(x), axis=0)
            self.data_norm = (self.data / col_sums) * 10000
            self.data_norm = np.log2(self.data_norm+1)
        elif self.alonabase.params['mrnafull'] and \
             self.alonabase.params['dataformat'] == 'raw':
            self._rpkm()
        elif self.alonabase.params['mrnafull'] and \
             self.alonabase.params['dataformat'] == 'rpkm':
            log_debug('_normalization() Running log2')
            self.data_norm = np.log2(self.data+1)
        else:
            self.data_norm = self.data
            log_debug('Normalization is not needed.')
        
        self._dump(self.data_norm, norm_mat_path)

    def load_data(self):
        """ Load expression matrix. """
        norm_mat_path = self.alonabase.get_working_dir() + '/normdata.joblib'
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
        self._print_dimensions()

        # normalize gene expression values
        self._normalization()

        log_debug('done loading expression matrix')

    def analysis(self):
        log_debug('Running analysis...')
        self._alona_analysis.find_variable_genes()
        self._alona_analysis.PCA()
        self._alona_analysis.tSNE()
