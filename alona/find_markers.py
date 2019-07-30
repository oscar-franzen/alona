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
import joblib

import numpy as np
import statsmodels.api as sm
import scipy.linalg
import pandas as pd
import patsy

from .log import (log_info, log_debug, log_error, log_warning)
from .stats import p_adjust_bh
from .constants import (OUTPUT)
from .celltypes import AlonaCellTypePred

class AlonaFindmarkers(AlonaCellTypePred):
    """
    Find markers class.
    """

    def __init__(self):
        super().__init__()

    def findMarkers(self):
        """
        Finds differentially expressed (DE) genes between clusters by fitting a linear
        model (LM) to gene expression (response variables) and clusters (explanatory
        variables). Model coefficients are estimated via the ordinary least squares
        method and p-values are calculated using t-statistics. One benefit of using LM for
        DE is that computations are vectorized and therefore very fast.

        The ideas behind using LM to explore DE have been extensively covered in the
        limma R package.

        Useful references:
        https://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf
        https://newonlinecourses.science.psu.edu/stat555/node/12/
        """

        log_debug('Entering findMarkers()')
        data_norm = self.data_norm.transpose()
        leiden_cl = self.leiden_cl
        clusters_targets = self.clusters_targets
        
        # remove clusters with too few cells
        data_norm = data_norm[np.isin(leiden_cl, clusters_targets)]
        leiden_cl = np.array(leiden_cl)[np.isin(leiden_cl, clusters_targets)]

        #joblib.dump([data_norm, leiden_cl, clusters_targets], 'testing.joblib')
        #sys.exit()

        # full design matrix
        dm_full = patsy.dmatrix('~ 0 + C(cl)', pd.DataFrame({'cl' : leiden_cl}))
        resid_df = dm_full.shape[0] - dm_full.shape[1]

        # gene expression should be the response
        lm = sm.regression.linear_model.OLS(endog=data_norm, # response
                                            exog=dm_full, # design matrix of clusters
                                           )
        res = lm.fit()
        coef = res.params # coefficients

        # computing standard errors
        # https://stats.stackexchange.com/questions/44838/how-are-the-standard-errors-of-coefficients-calculated-in-a-regression
        # http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/lm.series.html
        # residual variance for each gene
        dm_nrow = dm_full.shape[0]
        dm_ncol = dm_full.shape[1]
        sigma2 = ((data_norm-dm_full.dot(coef))**2).sum(axis=0)/(dm_nrow-dm_ncol)

        q = dm_full.transpose().dot(dm_full)
        chol = np.linalg.cholesky(q)
        chol2inv = scipy.linalg.cho_solve((chol, False), np.eye(chol.shape[0]))
        std_dev = np.sqrt(np.diag(chol2inv))
        clusts = np.unique(leiden_cl)

        # compare clusters
        comparisons = []
        out_pv = []

        for _, k in enumerate(np.unique(leiden_cl)):
            ref_cl = clusts[k]
            ref_coef = coef.iloc[k, :]

            # recompute coefficients for contrasts
            # https://genomicsclass.github.io/book/pages/interactions_and_contrasts.html
            con = np.zeros((coef.shape[0], k))
            np.fill_diagonal(con, -1)
            con[k,] = 1

            std_new = np.sqrt((std_dev**2).dot(con**2))

            for i in np.arange(k):
                std_err = std_new[i]**2*sigma2
                target_cl = clusts[i]

                # log2 fold change, reminder: log2(A/B)=log2(A)-log2(B)
                cur_lfc = ref_coef - coef.iloc[i, :]
                cur_lfc.index = std_err.index

                # compute p-values
                cur_t = cur_lfc/np.sqrt(std_err)
                t_dist = scipy.stats.t(resid_df)
                left = t_dist.cdf(cur_t)
                right = 1 - left

                # two sided test
                pv = np.minimum(left, right)
                
                comparisons.append('%s_vs_%s' % (k,i))
                out_pv.append(pd.Series(pv))
                
        out_merged = pd.concat(out_pv,axis=1)
        out_merged.columns = comparisons
        out_merged.index = data_norm.columns
        
        fn = self.get_wd() + OUTPUT['FILENAME_ALL_T_TESTS']
        out_merged.to_csv(fn, sep=',')
        
        lab = []

        for i in out_pv:
            lab.append(pd.Series(data_norm.columns))

        pval = pd.concat(out_pv)
        ll = pd.DataFrame({'gene' : pd.concat(lab), 'pval' : pval,
                           'padj' : p_adjust_bh(pval)})
        
        fn = self.get_wd() + OUTPUT['FILENAME_ALL_T_TESTS_LONG']
        ll.to_csv(fn, sep=',', index=False)

        log_debug('Exiting findMarkers()')
