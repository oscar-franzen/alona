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
import numpy  as np

from .log import (log_info, log_error)

class Cell():
    """
    Cell data class.
    """

    def __init__(self, alonabase):
        self.alonabase = alonabase
        self.data = None

    def _validate_counts(self):
        """ Basic data validation of gene expression values. """
        if (np.any(self.data.dtypes != 'int64') and
           self.alonabase.params['dataformat'] == 'raw'):
            log_error('Non-count values detected in data matrix while data format is \
set to raw read counts.')
        elif self.alonabase.params['dataformat'] == 'log2':
            if np.any(data>1000):
                log_error('Data do not appear to be log2 transformed.')
        else:
            logging.debug('_validate_counts() finished without a problem')
            
    def _remove_empty(self):
        """ Removes empty cells and genes """
        data_zero = self.data == 0
        cells = data_zero.sum(axis=1) # 1 = columns
        genes = data_zero.sum(axis=0) # 0 = rows

        total_genes = self.data.shape[0]
        total_cells = self.data.shape[1]

        if np.sum(cells == total_cells) > 0:
            log_info('%s empty cells will be removed' % (np.sum(cells==total_cells)))
            self.data = self.data.loc[self.data.index[np.logical_not(cells==total_cells)]]
        if np.sum(genes == total_genes) > 0:
            log_info('%s empty genes will be removed' % (np.sum(genes==total_genes)))
            self.data = self.data[self.data.columns[np.logical_not(genes==total_genes)]]
            
    def _remove_mito(self):
        """ Remove mitochondrial genes. """
        if self.alonabase.params['nomito']:
            logging.debug('removing mitochondrial genes')

            mt_count = self.data.index.str.contains('^mt-',regex=True, case=False)
            
            if np.sum(mt_count) > 0:
                log_info('detected and removed %s mitochondrial genes' % np.sum(mt_count))
                self.data = self.data.loc[self.data.index[np.logical_not(mt_count)]]

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

        logging.debug('done loading expression matrix')
