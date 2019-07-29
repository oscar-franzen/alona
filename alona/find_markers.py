"""
 alona

 Description:
 An analysis pipeline for scRNA-seq data.

 How to use:
 https://github.com/oscar-franzen/alona/

 Contact:
 Oscar Franzen <p.oscar.franzen@gmail.com>
"""

import sys

import numpy as np
import statsmodels.api as sm
import scipy.linalg
import pandas as pd
import patsy

from .log import (log_info, log_debug, log_error, log_warning)
from .celltypes import AlonaCellTypePred

class AlonaFindmarkers(AlonaCellTypePred):
    """
    Find markers class.
    """

    def __init__(self):
        super().__init__()

    def findMarkers(self):
        """
        Finds differentially expressed genes between clusters by fitting a linear model
        to gene expression (response variables) and clusters (explanatory variables).
        Model coefficients are estimated via the ordinary least squares method and
        p-values are calculated using t-statistics. One of the benefits of using LM for
        DE is that computations are vectorized and therefore very fast.
        
        The ideas behind using linear models to explore differential expression have been
        extensively covered by the limma R package.
        
        Reference:
        https://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf
        """

        log_debug('Entering findMarkers()')
        data_norm = self.data_norm.transpose()
        leiden_cl = self.leiden_cl

        # full design matrix
        dm_full = patsy.dmatrix('~ 0 + C(cl)', pd.DataFrame({'cl' : leiden_cl}))
        resid_df = dm_full.shape[0] - dm_full.shape[1]

        # gene expression should be the response
        lm = sm.regression.linear_model.OLS(endog=data_norm, # response
                                            exog=dm_full, # design matrix of clusters
                                           )
        res = lm.fit()
        coef = res.params # coefficients

        pred = lm.predict(res.params)

        # computing standard errors
        # https://stats.stackexchange.com/questions/44838/how-are-the-standard-errors-of-coefficients-calculated-in-a-regression
        # http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/lm.series.html
        # residual variance for each gene
        dm_nrow = dm_full.shape[0]
        dm_ncol = dm_full.shape[1]
        sigma2 = ((data_norm-dm_full.dot(coef))**2).sum(axis=0)/(dm_nrow-dm_ncol)

        q = dm_full.transpose().dot(dm_full)
        chol = np.linalg.cholesky(q)
        chol2inv = scipy.linalg.cho_solve((chol, False), np.eye(chol.shape[ 0 ]))
        std_dev = np.sqrt(np.diag(chol2inv))
        clusts = np.unique(leiden_cl)

        # compare clusters
        for _, k in enumerate(np.unique(leiden_cl)):
            ref_cl = clusts[k]
            ref_coef = coef.iloc[k,:]

            # recompute coefficients for contrasts
            # https://genomicsclass.github.io/book/pages/interactions_and_contrasts.html
            con = np.zeros((coef.shape[0], k))
            np.fill_diagonal(con, -1)
            con[k,] = 1

            std_new = np.sqrt((std_dev**2).dot(con**2))
            std_err = std_new**2*sigma2

            for i in np.arange(k-1):
                target_cl = clusts[i]

                # log2 fold change, reminder: log2(A/B)=log2(A)-log2(B)
                cur_lfc = ref_coef - coef.iloc[i,:]
                cur_lfc.index = std_err.index

                # compute p-values
                cur_t = cur_lfc/np.sqrt(std_err)
                t_dist = scipy.stats.t(resid_df)
                left = t_dist.cdf(cur_t)
                right = 1 - left

                # two sided test
                pv = np.minimum(left,right)

        log_debug('Exiting findMarkers()')
