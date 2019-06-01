"""
 alona

 Description:
 Methods for identifying highly variable genes.

 How to use:
 https://github.com/oscar-franzen/alona/

 Details:
 https://alona.panglaodb.se/

 Contact:
 Oscar Franzen <p.oscar.franzen@gmail.com>
"""

import pandas as pd
import numpy as np
import scipy.stats

from .glm.glm import GLM
from .glm.families import Gamma

from .log import (log_info, log_debug, log_error, log_warning)
from .stats import p_adjust_bh

class AlonaHighlyVariableGenes():
    """
    HVG class.
    """

    def __init__(self, hvg_method, hvg_n, data_norm):
        self.hvg_method = hvg_method
        self.hvg_n = hvg_n
        self.data_norm = data_norm
        
    def find(self):
        """ Finds HVG. Returns an array of HVG. """
        if self.hvg_method == 'seurat':
            hvg = self.hvg_seurat()
        elif self.hvg_method == 'brennecke':
            hvg = self.hvg_brennecke()
        return hvg
        
    @staticmethod
    def _exp_mean(mat):
        # Axis 0 will act on all the ROWS in each COLUMN
        # Axis 1 will act on all the COLUMNS in each ROW
        return mat.mean(axis=1)

    def hvg_brennecke(self, norm_ERCC=None, fdr=0.1, minBiolDisp=0.5):
        """
        Implements the method of Brennecke et al. (2013) to identify highly variable genes.
        Largely follows the function Brennecke_getVariableGenes from the R package M3Drop.

        The below code fits data using GLM with Fisher Scoring. GLM code copied from
        (credits to @madrury for this code): https://github.com/madrury/py-glm

        Brennecke et al. (2013) Accounting for technical noise in single-cell RNA-seq
        experiments. Nature Methods 10.1038/nmeth.2645
        """
        
        data_norm = self.data_norm

        if norm_ERCC == None:
            norm_ERCC = data_norm

        data_norm = 2**data_norm-1
        norm_ERCC = 2**norm_ERCC-1

        norm_ERCC = norm_ERCC.dropna(axis=1, how='all')

        # technical gene (spikes)
        meansSp = norm_ERCC.mean(axis=1)
        varsSp = norm_ERCC.var(axis=1)
        cv2Sp = varsSp/meansSp**2

        # biological genes
        meansGenes = data_norm.mean(axis=1)
        varsGenes = data_norm.var(axis=1)
        cv2Genes = varsGenes/meansGenes**2

        minMeanForFit = np.quantile(meansSp[cv2Sp > 0.3], 0.8)
        useForFit = meansSp >= minMeanForFit

        if np.sum(useForFit) < 20:
            meansAll = data_norm.mean(axis=1)
            cv2All = data_norm.var(axis=1)
            minMeanForFit = np.quantile(meansAll[cv2All > 0.3], 0.8)
            useForFit = meansSp >= minMeanForFit

        gamma_model = GLM(family=Gamma())

        x = pd.DataFrame({'a0' : [1]*len(meansSp[useForFit]), 'a1tilde' : 1/meansSp[useForFit]})

        # modified to use the identity link function
        gamma_model.fit(np.array(x), y=np.array(cv2Sp[useForFit]))
        a0 = gamma_model.coef_[0]
        a1 = gamma_model.coef_[1]

        psia1theta = a1
        minBiolDisp = minBiolDisp**2

        m = norm_ERCC.shape[1]
        cv2th = a0+minBiolDisp+a0*minBiolDisp

        testDenom = (meansGenes*psia1theta+(meansGenes**2)*cv2th)/(1+cv2th/m)
        q = varsGenes * (m - 1)/testDenom

        p = 1-scipy.stats.chi2.cdf(q, m-1)
        padj = p_adjust_bh(p)
        res = pd.DataFrame({'gene': meansGenes.index, 'pvalue' : p, 'padj' : padj})
        filt = res[res['padj'] < fdr]['gene']

        return np.array(filt.head(self.hvg_n))

    def hvg_seurat(self):
        """
        Retrieves a list of highly variable genes. Mimics Seurat's `FindVariableGenes`.
        The function bins the genes according to average expression, then calculates
        dispersion for each bin as variance to mean ratio. Within each bin, Z-scores are
        calculated and returned. Z-scores are ranked and the top 1000 are selected.
        """
        
        data_norm = self.data_norm
        
        log_debug('Entering hvg_seurat()')

        # number of bins
        num_bin = 20

        gene_mean = self._exp_mean(data_norm)
        # equal width (not size) of bins
        bins = pd.cut(gene_mean, num_bin)

        ret = []

        for _, sliced in data_norm.groupby(bins):
            # Axis 0 will act on all the ROWS in each COLUMN
            # Axis 1 will act on all the COLUMNS in each ROW
            dispersion = sliced.var(axis=1)/sliced.mean(axis=1)
            zscores = (dispersion-dispersion.mean())/dispersion.std()
            ret.append(zscores)

        ret = pd.concat(ret)
        ret = ret.sort_values(ascending=False)
        self.top_hvg = ret.head(self.hvg_n)
        
        ret = np.array(self.top_hvg.index)
        log_debug('Finishing hvg_seurat()')
        return ret
