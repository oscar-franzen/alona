"""
 alona

 Description:
 Methods to identify highly variable genes.

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

from scipy.optimize import least_squares
from scipy.stats import gaussian_kde
from scipy.stats import norm
from scipy.optimize import minimize

from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression

from .glm.glm import GLM
from .glm.families import Gamma

from .log import (log_info, log_debug, log_error, log_warning)
from .stats import p_adjust_bh

class AlonaHighlyVariableGenes():
    """
    HVG class.
    """

    def __init__(self, hvg_method, hvg_n, data_norm, data_ERCC):
        self.hvg_method = hvg_method
        self.hvg_n = hvg_n
        self.data_norm = data_norm
        self.data_ERCC = data_ERCC

    def find(self):
        """ Finds HVG. Returns an array of HVG. """
        if self.hvg_method == 'seurat':
            hvg = self.hvg_seurat()
        elif self.hvg_method == 'Brennecke2013':
            hvg = self.hvg_brennecke()
        elif self.hvg_method == 'scran':
            hvg = self.hvg_scran()
        elif self.hvg_method == 'Chen2016':
            hvg = self.hvg_chen2016()
        elif self.hvg_method == 'M3Drop_smartseq2':
            hvg = self.hvg_M3Drop_smartseq2()
        elif self.hvg_method == 'M3Drop_UMI':
            hvg = self.hvg_M3Drop_UMI()
        else:
            log_error('Unknown hvg method specified.')
        return hvg

    @staticmethod
    def _exp_mean(mat):
        # Axis 0 will act on all the ROWS in each COLUMN
        # Axis 1 will act on all the COLUMNS in each ROW
        return mat.mean(axis=1)

    def hvg_brennecke(self, norm_ERCC=None, fdr=0.1, minBiolDisp=0.5):
        """
        Implements the method of Brennecke et al. (2013) to identify highly variable genes.
        Largely follows the function BrenneckeGetVariableGenes from the R package M3Drop.

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

    def hvg_scran(self):
        """
        This function implements the approach from scran to identify highly variable genes.

        Expression counts should be normalized and on a log scale.

        Outline of the steps:

        1. fits a polynomial regression model to mean and variance of the technical genes
        2. decomposes the total variance of the biological genes by subtracting the
           technical variance predicted by the fit
        3. sort based on biological variance

        Reference:
        Lun ATL, McCarthy DJ, Marioni JC (2016). “A step-by-step workflow for low-level
        analysis of single-cell RNA-seq data with Bioconductor.” F1000Research,
        doi.org/10.12688/f1000research.9501.2
        """

        if self.data_ERCC.empty:
            log_error(None, 'Running "--hvg scran" requires ERCC spikes in the dataset. \
these should begin with ERCC- followed by numbers.')

        norm_data = self.data_norm
        norm_ERCC = self.data_ERCC
        norm_ERCC = norm_ERCC.dropna(axis=1, how='all')

        means_tech = norm_ERCC.mean(axis=1)
        vars_tech = norm_ERCC.var(axis=1)

        to_fit = np.log(vars_tech)
        arr = [list(item) for item in zip(*sorted(zip(means_tech, to_fit)))]

        x = arr[0]
        y = arr[1]

        poly_reg = PolynomialFeatures(degree=4)
        x_poly = poly_reg.fit_transform(np.array(x).reshape(-1, 1))
        pol_reg = LinearRegression()
        pol_reg.fit(x_poly, y)

        #plt.scatter(x, y, color='red')
        #plt.plot(x, pol_reg.predict(poly_reg.fit_transform(np.array(x).reshape(-1,1))), color='blue')
        #plt.xlabel('mean')
        #plt.ylabel('var')
        #plt.show()

        # predict and remove technical variance
        bio_means = norm_data.mean(axis=1)
        vars_pred = pol_reg.predict(poly_reg.fit_transform(np.array(bio_means).reshape(-1, 1)))
        vars_bio_total = norm_data.var(axis=1)

        # biological variance component
        vars_bio_bio = vars_bio_total - vars_pred
        vars_bio_bio = vars_bio_bio.sort_values(ascending=False)
        return vars_bio_bio.head(self.hvg_n).index.values

    def hvg_chen2016(norm_data):
        """
        This function implements the approach from Chen (2016) to identify highly variable
        genes.
        https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-2897-6

        The code is essentially a direct translation of the R code from:
        https://github.com/hillas/scVEGs/blob/master/scVEGs.r

        Expression counts should be normalized and not on a log scale.
        """

        norm_data = 2**norm_data-1
        avg = norm_data.mean(axis=1)
        norm_data = norm_data[avg > 0]

        rows = norm_data.shape[0]
        avg = norm_data.mean(axis=1)
        std = norm_data.std(axis=1)
        cv = std / avg

        xdata = avg
        ydata = np.log10(cv)

        A = np.vstack([np.log10(xdata), np.ones(len(xdata))]).T
        res = np.linalg.lstsq(A, ydata, rcond=None)[0]

        def predict(k, m, x):
            y = k*x+m
            return y

        #plt.clf()
        #plt.scatter(np.log10(xdata), ydata, color='red')
        #plt.scatter(np.log10(xdata), predict(*res, np.log10(xdata)), color='blue')
        #plt.show()

        xSeq = np.arange(min(np.log10(xdata)), max(np.log10(xdata)), 0.005)

        def h(i):
            a = np.log10(xdata) >= (xSeq[i] - 0.05)
            b = np.log10(xdata) < (xSeq[i] + 0.05)
            return np.sum((a & b))

        gapNum = [h(i) for i in range(0, len(xSeq))]
        cdx = np.nonzero(np.array(gapNum) > rows*0.005)[0]
        xSeq = 10 ** xSeq

        ySeq = predict(*res, np.log10(xSeq))
        yDiff = np.diff(ySeq, 1)
        ix = np.nonzero((yDiff > 0) & (np.log10(xSeq[0:-1]) > 0))[0]

        if len(ix) == 0:
            ix = len(ySeq) - 1
        else:
            ix = ix[0]

        xSeq_all = 10**np.arange(min(np.log10(xdata)), max(np.log10(xdata)), 0.001)

        xSeq = xSeq[cdx[0]:ix]
        ySeq = ySeq[cdx[0]:ix]

        reg = LinearRegression().fit(np.log10(xSeq).reshape(-1, 1), ySeq)

        #plt.clf()
        #plt.scatter(np.log10(xSeq), ySeq, color='red')
        #plt.scatter(xSeq, reg.predict(np.array(xSeq).reshape(-1,1)), color='blue')
        #plt.show()

        ydataFit = reg.predict(np.log10(xSeq_all).reshape(-1, 1))

        logX = np.log10(xdata)
        logXseq = np.log10(xSeq_all)

        cvDist = []

        for i in range(0, len(logX)):
            cx = np.nonzero((logXseq >= (logX[i] - 0.2)) & (logXseq < (logX[i] + 0.2)))[0]
            tmp = np.sqrt((logXseq[cx] - logX[i])**2 + (ydataFit[cx] - ydata[i])**2)
            tx = np.argmin(tmp)

            if logXseq[cx[tx]] > logX[i]:
                if ydataFit[cx[tx]] > ydata[i]:
                    cvDist.append(-1*tmp[tx])
                else:
                    cvDist.append(tmp[tx])
            elif logXseq[cx[tx]] <= logX[i]:
                if ydataFit[cx[tx]] < ydata[i]:
                    cvDist.append(tmp[tx])
                else:
                    cvDist.append(-1*tmp[tx])

        cvDist = np.log2(10**np.array(cvDist))
        dor = gaussian_kde(cvDist)
        dor_y = dor(cvDist)
        distMid = cvDist[np.argmax(dor_y)]
        dist2 = cvDist - distMid

        a = dist2[dist2 <= 0]
        b = abs(dist2[dist2 < 0])
        c = distMid
        tmpDist = np.concatenate((a, b))
        tmpDist = np.append(tmpDist, c)

        # estimate mean and sd using maximum likelihood
        distFit = norm.fit(tmpDist)

        pRaw = 1-norm.cdf(cvDist, loc=distFit[0], scale=distFit[1])
        pAdj = p_adjust_bh(pRaw)

        res = pd.DataFrame({'gene': norm_data.index, 'pvalue' : pRaw, 'padj' : pAdj})
        res = res.sort_values(by='pvalue')

        filt = res[res['padj'] < 0.10]['gene']

        return np.array(filt.head(self.hvg_n))

    def hvg_M3Drop_smartseq2(self, norm_data):
        """
        This function implements the approach from M3Drop to identify highly variable genes
        and takes an alternative approach by using genes' dropout rates instead of 
        variance. In short, this method calculates dropout rates and mean expression for
        every gene, then models these with the Michaelis-Menten equation (parameters
        are estimated with maximum likelihood optimization). The basis for using MM
        is because most dropouts are caused by failure of the enzyme reverse transcriptase,
        thus the dropout rate can be modelled with theory developed for enzyme reactions.

        Use this method with SMART-seq2 data.

        The method is briefly described here: https://doi.org/10.1093/bioinformatics/bty1044
        R code: https://github.com/tallulandrews/M3Drop

        Expression counts should be normalized and not on a log scale.
        """

        norm_data = 2**norm_data-1
        ncells = norm_data.shape[1]

        gene_info_p = 1-np.sum(norm_data > 0, axis=1)/ncells
        gene_info_p_stderr = np.sqrt(gene_info_p*(1-gene_info_p)/ncells)

        gene_info_s = norm_data.mean(axis=1)
        gene_info_s_stderr = np.sqrt((np.mean(norm_data**2, axis=1) - gene_info_s**2)/ncells)

        xes = np.log(gene_info_s)/np.log(10)

        # maximum likelihood estimate of model parameters
        s = gene_info_s
        p = gene_info_p

        def neg_loglike(theta):
            krt = theta[0]
            sigma = theta[1]
            R = p-(1-(s/(krt+s)))
            R = norm.logpdf(R, 0, sigma)
            return -1*np.sum(R)

        theta_start = np.array([3, 0.25])
        res = minimize(neg_loglike, theta_start, method='Nelder-Mead', options={'disp': True})
        est = res.x

        krt = est[0]
        res_err = est[1]

        # testing step
        p_obs = gene_info_p
        always_detected = p_obs == 0
        p_obs[p_obs == 0] = min(p_obs[p_obs > 0])/2
        p_err = gene_info_p_stderr
        K_obs = krt
        S_mean = gene_info_s
        K_equiv = p_obs*S_mean/(1-p_obs)
        S_err = gene_info_s_stderr
        K_equiv_err = abs(K_equiv)*np.sqrt((S_err/S_mean)**2 + (p_err/p_obs)**2)
        K_equiv_log = np.log(K_equiv)
        thing = K_equiv-K_equiv_err
        thing[thing <= 0] = 10**-100
        K_equiv_err_log = abs(np.log(thing)-K_equiv_log)
        K_equiv_err_log[K_equiv-K_equiv_err <= 0] = 10**10
        K_obs_log = np.log(krt)
        K_err_log = np.std(K_equiv_log-K_obs_log)/np.sqrt(len(K_equiv_log))
        Z = (K_equiv_log - K_obs_log)/np.sqrt(K_equiv_err_log**2+K_err_log**2)

        pval = 1 - norm.cdf(Z)
        pval[always_detected] = 1
        padj = p_adjust_bh(pval)
        #pval[is.na(pval)] <- 1; # deal with never detected
        #effect_size <- K_equiv/fit$K;
        #effect_size[is.na(effect_size)] <- 1; # deal with never detected

        res = pd.DataFrame({'gene': norm_data.index, 'pvalue' : pval, 'padj' : padj})
        res = res.sort_values(by='pvalue')
        filt = res[res['padj'] < 0.10]['gene']

        return np.array(filt.head(self.hvg_n))

    def hvg_M3Drop_UMI(self, data):
        """
        This function implements the approach from M3Drop to identify highly variable genes
        and takes an alterntive approach to identify highly variable genes by using
        the dropout rate instead of the variance. Briefly, this method calculate the dropout
        rates and mean expression for every gene. Then it finds the parameters for the
        Michaelis-Menten equation using maximum likelihood optimization.

        Use this method with UMI data.

        The method is briefly described here: https://doi.org/10.1093/bioinformatics/bty1044
        R code: https://github.com/tallulandrews/M3Drop

        Expression counts should be raw read counts.
        """

        tjs = data.sum(axis=1) # no. mol/gene
        tis = data.sum(axis=0) # no. mol/cell
        
        djs = data.shape[1]-np.sum(data>0,axis=1) # dropouts no. per gene
        dis = data.shape[0]-np.sum(data>0,axis=0) # dropouts no. per cell

        nc = data.shape[1]
        ng = data.shape[0]
        
        # total sampled molecules
        total = sum(tis)
        
        min_size = 10**-10
        my_rowvar = []
        
        for i, row in data.iterrows():
            mu_is = tjs[i]*tis/total
            my_rowvar.append(np.var(row-mu_is))

        my_rowvar = np.array(my_rowvar)
        size = tjs**2*(sum(tis**2)/total**2)/((nc-1)*my_rowvar-tjs)

        max_size = 10*max(size)
        size[size < 0] = max_size
        size[size < min_size] = min_size
        
        size_g = size
        forfit = (size < max(size_g)) & (tjs > 0) & (size_g > 0)
        higher = (np.log(tjs/nc)/np.log(2)) > 4
        
        if sum(higher==True) > 2000:
            forfit = higher & forfit
        
        rg = LinearRegression()
        X = np.array(np.log((tjs/nc)[forfit])).reshape(-1,1)
        y = np.array(np.log(size_g[forfit]))
        rg.fit(X=X, y=y)
        coef_1 = rg.intercept_
        coef_2 = rg.coef_[0]
        
        exp_size = np.exp(coef_1 + coef_2 * np.log(tjs/nc))
        
        droprate_exp = []
        droprate_exp_err = []
        
        for i in range(0,ng):
            mu_is = tjs[i]*tis/total
            p_is = (1+mu_is/exp_size[i])**(-exp_size[i])
            p_var_is = p_is*(1-p_is)
            droprate_exp.append(sum(p_is)/nc)
            droprate_exp_err.append(np.sqrt(sum(p_var_is)/(nc**2)))
            
        droprate_exp = np.array(droprate_exp)
        droprate_exp[droprate_exp < (1/nc)] = 1/nc
        droprate_obs = djs/nc
        droprate_obs_err = np.sqrt(droprate_obs*(1-droprate_obs)/nc)

        diff = droprate_obs-droprate_exp
        
        droprate_exp_err = np.array(droprate_exp_err)
        combined_err = np.sqrt(droprate_exp_err**2+droprate_obs_err**2)

        Zed = diff/combined_err
        pvalues = 1-norm.cdf(Zed)
        pvalues = pd.Series(pvalues, index=data.index)
        
        padj = p_adjust_bh(pvalues)

        res = pd.DataFrame({'gene': data.index, 'pvalue' : pvalues, 'padj' : padj})
        res = res.sort_values(by='pvalue')
        filt = res[res['padj'] < 0.10]['gene']

        return np.array(filt.head(self.hvg_n))
