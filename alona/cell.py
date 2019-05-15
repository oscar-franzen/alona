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
