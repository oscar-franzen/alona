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
        
    def _choose_leftright_pvalues(self, left, right, direction):
        if direction == 'up':
            pv = right
        elif direction == 'down':
            pv = left
        else:
            pv = np.minimum(left,right)*2
        return pv
        
    def discover_markers(self, lfc=0, direction='up'):
        pval_mat = self.fit_lm_tt(lfc, direction)
        self.combine_tests(pval_mat)

    def fit_lm_tt(self, lfc, direction):
        """
        Finds differentially expressed (DE) genes between clusters by fitting a linear
        model (LM) to gene expression (response variables) and clusters (explanatory
        variables) and performs pairwise t-tests for significance. Model coefficients are
        estimated via the ordinary least squares method and p-values are calculated using
        t-statistics. One benefit of using LM for DE is that computations are vectorized
        and therefore very fast.

        The ideas behind using LM to explore DE have been extensively covered in the
        limma R package.
        
        Arguments
        =========
        direction : Can be 'any' for any direction 'up' for up-regulated and 'down'
                    for down-regulated.

        Useful references:
        https://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf
        https://newonlinecourses.science.psu.edu/stat555/node/12/
        """

        log_debug('Entering fit_lm_tt()')
        data_norm = self.data_norm.transpose()
        
        leiden_cl = self.leiden_cl
        clusters_targets = self.clusters_targets
        
        # remove clusters with too few cells
        data_norm = data_norm[np.isin(leiden_cl, clusters_targets)]
        leiden_cl = np.array(leiden_cl)[np.isin(leiden_cl, clusters_targets)]

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

        # perform all pairwise comparisons of clusters (t-tests)
        comparisons = []
        out_t_stats = []
        out_pv = []
        out_lfc = []
        out_mge_g1 = []
        out_mge_g2 = []

        for k, _ in enumerate(clusts):
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
                if lfc == 0:
                    cur_t = cur_lfc/np.sqrt(std_err)
                    t_dist = scipy.stats.t(resid_df)
                    
                    left = t_dist.cdf(cur_t)
                    right = t_dist.sf(cur_t)
                    
                    pv1 = self._choose_leftright_pvalues(left, right, direction)
                    pv2 = self._choose_leftright_pvalues(right, left, direction)
                
                comparisons.append('%s_vs_%s' % (clusts[k], clusts[i]))
                comparisons.append('%s_vs_%s' % (clusts[i], clusts[k]))
                
                out_pv.append(pd.Series(pv1))
                out_pv.append(pd.Series(pv2))
            
                out_t_stats.append(pd.Series(cur_t))
                out_t_stats.append(pd.Series(cur_t))
                
                out_lfc.append(pd.Series(cur_lfc))
                out_lfc.append(pd.Series(cur_lfc*-1))
                
                out_mge_g1.append(mge[k])
                out_mge_g2.append(mge[i])
                
                out_mge_g1.append(mge[i])
                out_mge_g2.append(mge[k])
        
        out_merged = pd.concat(out_pv,axis=1)
        out_merged.columns = comparisons
        out_merged.index = data_norm.columns
        
        fn = self.get_wd() + OUTPUT['FILENAME_ALL_T_TESTS']
        out_merged.to_csv(fn, sep=',')
        
        lab1 = []
        lab2 = []
        anno = []

        for q in comparisons:
            lab1.append(pd.Series([q]*data_norm.columns.shape[0]))
            lab2.append(pd.Series(data_norm.columns))
            if type(self.anno)==pd.core.frame.DataFrame:
                anno.append(self.anno['desc'])

        pval = pd.concat(out_pv, ignore_index=True)
        ll = pd.DataFrame({'comparison_A_vs_B' : pd.concat(lab1, ignore_index=True),
                   'gene' : pd.concat(lab2, ignore_index=True),
                   'p_val' : pval,
                   'FDR' : p_adjust_bh(pval),
                   't_stat' : pd.concat(out_t_stats, ignore_index=True),
                   'logFC' : pd.concat(out_lfc, ignore_index=True),
                   'mean.A' : pd.concat(out_mge_g1, ignore_index=True),
                   'mean.B' : pd.concat(out_mge_g2, ignore_index=True)})
        if anno:
            ll['annotation'] = pd.concat(anno, ignore_index=True)

        fn = self.get_wd() + OUTPUT['FILENAME_ALL_T_TESTS_LONG']
        ll.to_csv(fn, sep='\t', index=False)

        log_debug('Exiting fit_lm_tt()')
        return out_merged

    def combine_tests(self, pval_mat):
        """ Uses Simes' method for combining p-values.
        
        Arguments
        =========
        pval_mat : A data frame with p-values.

        Simes, R. J. (1986). An improved Bonferroni procedure for multiple tests of
        significance. Biometrika, 73(3):751-754. 
        """
        
        test_mat = pval_mat
        clusters_targets = self.clusters_targets
        
        res = []
        for cl in clusters_targets:
            idx = test_mat.columns.str.match('^%s_vs'%cl)
            subset_mat = test_mat.iloc[:,idx]
            
            r = subset_mat.rank(axis=1)
            T = (subset_mat.shape[1]*subset_mat/r).min(axis=1).sort_values()
            T[T>1] = 1
            
            df = pd.DataFrame({'cluster' : [cl]*len(T), 'gene' : T.index, 'pvalue.Simes' : T})
            res.append(df)
        
        fn = self.get_wd() + OUTPUT['FILENAME_MARKERS']
        pd.concat(res).to_csv(fn, sep='\t', index=False)
