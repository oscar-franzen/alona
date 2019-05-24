"""
 This file contains cell type prediction methods used by alona.

 How to use alona:
 https://github.com/oscar-franzen/alona/
 https://alona.panglaodb.se/

 Contact: Oscar Franzen <p.oscar.franzen@gmail.com>
"""

from collections import defaultdict

import pandas as pd
import numpy as np
from sklearn.preprocessing import scale
from scipy.stats import fisher_exact

from .constants import (OUTPUT, GENOME, MARKERS)
from .log import (log_info, log_debug, log_error)
from .utils import get_alona_dir
from .stats import p_adjust_bh

class AlonaCellTypePred():
    """
    Cell type prediction.
    """

    def __init__(self, wd, alonaclustering, alonacell):
        self.alonaclustering = alonaclustering
        self.alonacell = alonacell
        self.wd = wd
        self.median_expr = None
        self.markers = None
        self.marker_freq = None

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

        log_debug('CTA_RANK_F() starting')

        median_expr = self.median_expr
        markers = self.markers
        marker_freq = self.marker_freq

        input_symbols = median_expr.index.str.extract('^(.+)_.+')[0].str.upper()
        median_expr.index = input_symbols

        # (1) centering is done by subtracting the column means
        # (2) scaling is done by dividing the (centered) by their standard deviations
        median_expr_Z = pd.DataFrame(scale(median_expr, with_mean=True, axis=0))
        median_expr_Z.index = median_expr.index
        median_expr_Z.columns = median_expr.columns

        # reference symbols
        fn = get_alona_dir() + GENOME['MOUSE_GENE_SYMBOLS']
        mgs = pd.read_csv(fn, header=None)
        mgs = mgs[0].str.upper()

        markers = markers[markers[markers.columns[0]].isin(mgs)]

        dd = defaultdict(list)
        for item in markers.groupby('cell type'):
            dd[item[0]] = set(item[1][item[1].columns[0]])

        # Following this reasoning:
        # Down-weighting overlapping genes improves gene set analysis
        # Tarca AL, Draghici S, Bhatti G, Romero R
        # BMC Bioinformatics 2012 13:136
        weights = dict()
        for ct in dd:
            s = dd[ct]
            s_freqs = marker_freq[marker_freq.index.isin(s)]
            weight = 1+np.sqrt(((max(marker_freq)-s_freqs)/(max(marker_freq)-min(marker_freq))))
            weights[ct] = weight

        # Axis 0 will act on all the ROWS in each COLUMN
        # Axis 1 will act on all the COLUMNS in each ROW
        def _guess_cell_type(x):
            rr = median_expr.loc[:, median_expr.columns == x.name].values.flatten()
            # genes expressed in this cell cluster
            genes_exp = set(x.index[rr > 0])
            # genes _not_ expressed in this cell cluster
            genes_not_exp = set(x.index[rr == 0])
            res = list()

            for ct in dd:
                s = dd[ct]
                x_ss = x[x.index.isin(s)]
                if len(x_ss) == 0: continue

                gene_weights = weights[ct]
                gene_weights = gene_weights[gene_weights.index.isin(x_ss.index)]
                gene_weights = pd.Series(gene_weights, x_ss.index)

                activity_score = sum(x_ss * gene_weights)/len(x_ss)**0.3

                # how many expressed genesets are found in the geneset?
                ct_exp = len(genes_exp&s)
                # how many _non_ expressed genes are found in the geneset?
                ct_non_exp = len(genes_not_exp&s)
                # how many expressed genes are NOT found in the geneset?
                ct_exp_not_found = len(genes_exp-s)
                # how many _non_ expressed genes are NOT found in the geneset?
                not_exp_not_found_in_geneset = len(genes_not_exp-s)
                # one sided fisher
                contigency_tbl = [[ct_exp, ct_non_exp],
                                  [ct_exp_not_found, not_exp_not_found_in_geneset]]

                odds_ratio, pval = fisher_exact(contigency_tbl, alternative='greater')
                markers_found = ','.join(list(genes_exp&s))
                if markers_found == '': markers_found = 'NA'
                res.append({'activity_score' : activity_score,
                            'ct' : ct,
                            'pvalue' : pval,
                            'markers' : markers_found})

            res = sorted(res, key=lambda k: k['activity_score'], reverse=True)
            return res

        ret = median_expr_Z.apply(func=_guess_cell_type, axis=0)

        # restructure
        lin = []

        bucket = [pd.DataFrame(k) for i, k in enumerate(ret)]
        final_tbl = pd.concat(bucket)

        padj = p_adjust_bh(final_tbl['pvalue'])
        final_tbl['padj_BH'] = padj
        final_tbl.columns = ['activity score',
                             'cell type',
                             'markers',
                             'p-value',
                             'adjusted p-value BH']

        fn = self.wd + OUTPUT['FILENAME_CTA_RANK_F']
        final_tbl.to_csv(fn, sep='\t')

        log_debug('CTA_RANK_F() finished')

    def load_markers(self):
        """ Load gene to cell type markers. """
        log_debug('Loading markers...')
        ma = pd.read_csv(get_alona_dir() + MARKERS['PANGLAODB'], sep='\t')

        # restrict to mouse
        ma = ma[ma.species.str.match('Mm')]
        self.markers = ma

        ui = ma.iloc[:, ma.columns == 'ubiquitousness index']
        ma = ma[np.array(ui).flatten() < 0.05]
        log_debug('Markers loaded')

        # marker frequency across the cell types
        ma_ss = ma.iloc[:, ma.columns.isin(['official gene symbol', 'cell type'])]
        self.marker_freq = ma_ss[ma_ss.columns[0]].value_counts()
        self.markers = ma_ss
