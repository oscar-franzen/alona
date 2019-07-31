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

        joblib.dump([data_norm, leiden_cl, clusters_targets], 'testing.joblib')
        #sys.exit()

        # full design matrix
        dm_full = patsy.dmatrix('~ 0 + C(cl)', pd.DataFrame({'cl' : leiden_cl}))
        resid_df = dm_full.shape[0] - dm_full.shape[1]

        # gene expression should be the response
        lm = sm.regression.linear_model.OLS(endog=data_norm, # response
                                            exog=dm_full # design matrix of clusters
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
        
        # mean gene expression for every gene in every cluster
        mge = [data_norm.iloc[leiden_cl == cl, :].mean() for cl in np.unique(leiden_cl)]

        # compare clusters
        comparisons = []
        out_t_stats = []
        out_pv = []
        out_lfc = []
        out_mge_g1 = []
        out_mge_g2 = []

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
                
                # times two since it is a two sided p-value
                pv = np.minimum(left, right)*2
                
                # cdf precision problem relating to floating point precision
                # https://stackoverflow.com/questions/6298105/precision-of-cdf-in-scipy-stats
                # https://github.com/scipy/scipy/issues/2238
                # sf = 1-cdf
                # check those with 0 with sf
                genes_recheck = pv==0
                pp = t_dist.sf(cur_t[genes_recheck])
                # put back
                pv[genes_recheck] = pp
                
                comparisons.append('%s_vs_%s' % (k,i))
                out_pv.append(pd.Series(pv))
                out_t_stats.append(pd.Series(cur_t))
                out_lfc.append(pd.Series(cur_lfc))
                
                out_mge_g1.append(mge[k])
                out_mge_g2.append(mge[i])
        
        out_merged = pd.concat(out_pv,axis=1)
        out_merged.columns = comparisons
        out_merged.index = data_norm.columns
        
        fn = self.get_wd() + OUTPUT['FILENAME_ALL_T_TESTS']
        out_merged.to_csv(fn, sep=',')
        
        lab1 = []
        lab2 = []

        for q in comparisons:
            lab1.append(pd.Series([q]*data_norm.columns.shape[0]))
            lab2.append(pd.Series(data_norm.columns))

        pval = pd.concat(out_pv, ignore_index=True)
        ll = pd.DataFrame({'comparison_A_vs_B' : pd.concat(lab1, ignore_index=True),
                   'gene' : pd.concat(lab2, ignore_index=True),
                   'p_val' : pval,
                   'FDR' : p_adjust_bh(pval),
                   't_stat' : pd.concat(out_t_stats, ignore_index=True),
                   'logFC' : pd.concat(out_lfc, ignore_index=True),
                   'mean.A' : pd.concat(out_mge_g1, ignore_index=True),
                   'mean.B' : pd.concat(out_mge_g2, ignore_index=True)})
        
        fn = self.get_wd() + OUTPUT['FILENAME_ALL_T_TESTS_LONG']
        ll.to_csv(fn, sep=',', index=False)

        log_debug('Exiting findMarkers()')
