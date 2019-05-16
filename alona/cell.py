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
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.pyplot import figure

from .log import (log_info, log_error)
from .constants import *

class Cell():
    """
    Cell data class.
    """

    def __init__(self, alonabase):
        self.alonabase = alonabase
        self.data = None

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
            logging.debug('_validate_counts() finished without a problem')

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
            logging.debug('removing mitochondrial genes')

            mt_count = self.data.index.str.contains('^mt-', regex=True, case=False)

            if np.sum(mt_count) > 0:
                log_info('detected and removed %s mitochondrial gene(s)' % np.sum(mt_count))
                self.data = self.data.loc[self.data.index[np.logical_not(mt_count)]]

    def _read_counts_per_cell_filter(self):
        """ Generates a bar plot of read counts per cell. """
        if self.alonabase.params['dataformat'] == 'raw':
            cell_counts = self.data.sum(axis=0)

            plt.clf()

            figure(num=None, figsize=(5, 5))

            no_cells_remove = np.sum(np.array(sorted(cell_counts, reverse=True)) < 1000)
            colors = ['#E69F00']*(len(cell_counts)-no_cells_remove) + \
                     ['#999999']*no_cells_remove

            plt.bar(np.arange(len(cell_counts)), sorted(cell_counts, reverse=True),
                    color=colors)
            plt.ylabel('raw read counts')
            plt.xlabel('cells (sorted on highest to lowest)')
            plt.axhline(linewidth=1, color='r', y=1000)
            plt.legend(handles=[
                mpatches.Patch(color='#E69F00', label='Passed'),
                mpatches.Patch(color='#999999', label='Removed')
            ])

            plt.savefig(self.alonabase.get_working_dir() + '/plots/' + \
                 FILENAME_BARPLOT_RAW_READ_COUNTS, bbox_inches='tight')
            plt.close()

            #log_info('Keeping %s out of %s cells' % (
            #     np.sum(cell_counts>1000), len(cell_counts)))

            #self.data = self.data[self.data.columns[cell_counts>1000]]

    def load_data(self):
        """ Load expression matrix. """
        logging.debug('loading expression matrix')
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

        self._validate_counts()
        self._remove_empty()
        self._remove_mito()
        self._read_counts_per_cell_filter()

        logging.debug('done loading expression matrix')
