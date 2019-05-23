"""
 This file contains cell type prediction methods used by alona.

 How to use alona:
 https://github.com/oscar-franzen/alona/
 https://alona.panglaodb.se/

 Contact: Oscar Franzen <p.oscar.franzen@gmail.com>
"""

import pandas as pd
import numpy as np

from .constants import OUTPUT
from .log import (log_info, log_debug, log_error)

class AlonaCellTypePred():
    """
    Cell type prediction.
    """

    def __init__(self, wd, alonaclustering, alonacell):
        self.alonaclustering = alonaclustering
        self.alonacell = alonacell
        self.wd = wd

    def median_exp(self):
        """ Represent each cluster with median gene expression. """
        
        log_debug('median_exp() Computing median expression per cluster...')
        
        clust = self.alonaclustering.leiden_cl
        data = self.alonacell.data_norm
        
        # Axis 0 will act on all the ROWS in each COLUMN
        # Axis 1 will act on all the COLUMNS in each ROW
        
        fn = self.wd + OUTPUT['FILENAME_MEDIAN_EXP']
        
        print(fn)
        
        ret = data.groupby(clust, axis=1).aggregate(np.median)
        ret.to_csv(fn,header=True)

        log_debug('median_exp() finished')
