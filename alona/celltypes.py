"""
 This file contains cell type prediction methods used by alona.

 How to use alona:
 https://github.com/oscar-franzen/alona/
 https://alona.panglaodb.se/

 Contact: Oscar Franzen <p.oscar.franzen@gmail.com>
"""

import pandas as pd
import numpy as np

from .constants import (OUTPUT, GENOME, MARKERS)
from .log import (log_info, log_debug, log_error)
from .utils import get_alona_dir

class AlonaCellTypePred():
    """
    Cell type prediction.
    """

    def __init__(self, wd, alonaclustering, alonacell):
        self.alonaclustering = alonaclustering
        self.alonacell = alonacell
        self.wd = wd
        self.median_expr = None

    def median_exp(self):
        """ Represent each cluster with median gene expression. """
        log_debug('median_exp() Computing median expression per cluster...')

        clust = self.alonaclustering.leiden_cl
        data = self.alonacell.data_norm

        fn = self.wd + OUTPUT['FILENAME_MEDIAN_EXP']        

        # Axis 0 will act on all the ROWS in each COLUMN
        # Axis 1 will act on all the COLUMNS in each ROW
        ret = data.groupby(clust, axis=1).aggregate(np.median)
        ret.to_csv(fn, header=True)
        
        self.median_expr = ret

        log_debug('median_exp() finished')

    def CTA_RANK_F(self):
        """
        Cell Type Activity and Rank-based annotation with a one-sided Fisher's Exact test
        """
        # reference symbols
        fn = get_alona_dir() + GENOME['MOUSE_GENE_SYMBOLS']
        mgs = pd.read_csv(fn, header=None)
        mgs = mgs[0].str.upper()
        
        input_symbols = self.median_expr.index.str.extract('^(.+)_.+')
        input_symbols = input_symbols[0].str.upper()

    def load_markers(self):
        """ Load gene to cell type markers. """
        log_debug('Loading markers...')
        ma = pd.read_csv(get_alona_dir() + MARKERS['PANGLAODB'], sep='\t')
        
        # restrict to mouse
        ma = ma[ma.species.str.match('Mm')]
        self.markers = ma
        
        ui = ma.iloc[:,ma.columns=='ubiquitousness index']
        ma = ma[np.array(ui).flatten()<0.05]
        log_debug('Markers loaded')
