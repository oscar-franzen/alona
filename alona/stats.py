"""
 alona

 Description:
 Statistical routines used by alona.

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

from glm.glm import GLM
from glm.families import Gamma

def p_adjust_bh(p):
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(float(len(p)), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]

def hvg_brennecke(norm_data, norm_ERCC=None, fdr=0.1, minBiolDisp=0.5):
    """
    Implements the method of Brennecke et al. (2013) to identify highly variable genes.
    Largely follows the function Brennecke_getVariableGenes from the R package M3Drop.
    
    The below code uses regular Gamma regression instead of GLM by Fisher Scoring with
    Levenberg Damping.
    
    Brennecke et al. (2013) Accounting for technical noise in single-cell RNA-seq
    experiments. Nature Methods 10.1038/nmeth.2645
    """

    if norm_ERCC == None:
        norm_ERCC = norm_data
    
    norm_data = 2**norm_data-1
    norm_ERCC = 2**norm_ERCC-1
    
    norm_ERCC = norm_ERCC.dropna(axis=1, how='all')
    
    # technical gene (spikes)
    meansSp = norm_ERCC.mean(axis=1)
    varsSp = norm_ERCC.var(axis=1)
    cv2Sp = varsSp/meansSp**2
    
    # biological genes
    meansGenes = norm_data.mean(axis=1)
    varsGenes = norm_data.var(axis=1)
    cv2Genes = varsGenes/meansGenes**2

    minMeanForFit = np.quantile(meansSp[cv2Sp>0.3], 0.8)
    useForFit = meansSp >= minMeanForFit
    
    if np.sum(useForFit) < 20:
        meansAll = norm_data.mean(axis=1)
        cv2All = norm_data.var(axis=1)
        minMeanForFit = np.quantile(meansAll[cv2All>0.3], 0.8)
        useForFit = meansSp >= minMeanForFit
    
    gamma_model = GLM(family=Gamma())
    
    x = pd.DataFrame({ 'a0' : [1]*len(meansSp[useForFit]), 'a1tilde' : 1/meansSp[useForFit] })
    
    # modified to use the identity link function
    gamma_model.fit(np.array(x),y=np.array(cv2Sp[useForFit]))
    a0 = gamma_model.coef_[0]
    a1 = gamma_model.coef_[1]

    psia1theta = a1
    minBiolDisp = minBiolDisp**2

    m = norm_ERCC.shape[1]
    cv2th = a0+minBiolDisp+a0*minBiolDisp

    testDenom = (meansGenes*psia1theta+(meansGenes**2)*cv2th)/(1+cv2th/m)
    q = varsGenes * (m - 1)/testDenom
    
    p = 1-scipy.stats.chi2.cdf(q,m-1)
    padj = p_adjust_bh(p)
    res = pd.DataFrame( {'gene': meansGenes.index, 'pvalue' : p, 'padj' : padj} )
    
    return res[res['padj']<fdr]['gene']
