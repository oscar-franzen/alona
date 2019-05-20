"""
 alona

 Description:
 An analysis pipeline for scRNA-seq data.

 How to use:
 https://github.com/oscar-franzen/alona/

 Details:
 https://alona.panglaodb.se/

 Contact:
 Oscar Franzen <p.oscar.franzen@gmail.com>
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import sklearn.manifold

from sklearn.cluster import DBSCAN

import alona.irlbpy

from .log import (log_info, log_debug, log_error)
from .constants import OUTPUT

class AlonaAnalysis():
    """
    Analysis class.
    """

    def __init__(self, alonacell):
        self._alonacell = alonacell
        self.top_hvg = None
        self.pca_components = None
        self.embeddings = None

    @staticmethod
    def _exp_mean(mat):
        # Axis 0 will act on all the ROWS in each COLUMN
        # Axis 1 will act on all the COLUMNS in each ROW
        return mat.mean(axis=1)

    def find_variable_genes(self):
        """
        Retrieves a list of highly variable genes. Mimics Seurat's `FindVariableGenes`.
        The function bins the genes according to average expression, then calculates
        dispersion for each bin as variance to mean ratio. Within each bin, Z-scores are
        calculated and returned. Z-scores are ranked and the top 1000 are selected.
        """
        log_debug('Finding variable genes')

        # number of bins
        num_bin = 20
        top_genes = 1000

        gene_mean = self._exp_mean(self._alonacell.data_norm)
        # equal width (not size) of bins
        bins = pd.cut(gene_mean, num_bin)

        ret = []

        for _, sliced in self._alonacell.data_norm.groupby(bins):
            # Axis 0 will act on all the ROWS in each COLUMN
            # Axis 1 will act on all the COLUMNS in each ROW
            dispersion = sliced.var(axis=1)/sliced.mean(axis=1)
            zscores = (dispersion-dispersion.mean())/dispersion.std()
            ret.append(zscores)

        ret = pd.concat(ret)
        ret = ret.sort_values(ascending=False)
        self.top_hvg = ret.head(top_genes)

        wd = self._alonacell.alonabase.get_working_dir()
        self.top_hvg.to_csv(wd + OUTPUT['FILENAME_HVG'], header=False)

    def PCA(self):
        """
        Calculate principal components using irlba.

        The augmented implicitly restarted Lanczos bidiagonalization algorithm (IRLBA)
        finds a few approximate largest singular values and corresponding singular
        vectors using a method of Baglama and Reichel.

        A fast and memory-efficient way to compute a partial SVD, principal
        components, and some specialized partial eigenvalue decompositions.

        Reference:
        Baglama, James, and Lothar Reichel. “Augmented implicitly restarted Lanczos
        bidiagonalization methods.” SIAM Journal on Scientific Computing 27.1 (2005):
        19-42.

        Some useful notes about the R implementation:
        http://bwlewis.github.io/irlba/
        """

        log_debug('Running PCA...')

        index_v = self._alonacell.data_norm.index.isin(self.top_hvg.index)
        sliced = self._alonacell.data_norm[index_v]
        lanc = alona.irlbpy.lanczos(sliced, nval=75, maxit=1000)

        # weighing by var (Seurat-style)
        self.pca_components = np.dot(lanc.V, np.diag(lanc.s))
        log_debug('Finished PCA')

    def tSNE(self, tsne_out_path):
        """ Projects data to a two dimensional space using the tSNE algorithm. """
        log_debug('Running t-SNE...')

        tsne = sklearn.manifold.TSNE(n_components=2, n_iter=2000)
        self.embeddings = tsne.fit_transform(pd.DataFrame(self.pca_components))
        pd.DataFrame(self.embeddings).to_csv(path_or_buf=tsne_out_path, sep=',',
                                             header=None, index=False)

        log_debug('Finished t-SNE')

    def cluster(self):
        """ Cluster cells. """
        pass
        #db = DBSCAN(eps=0.3, min_samples=10)
        #db.fit(np.rot90(pca.components_)))
        #np.unique(db.labels_, return_counts=True)
    
    def cell_scatter_plot(self):
        plt.clf()
        figure(num=None, figsize=(5, 5))
        df = pd.DataFrame(self.embeddings)
        plt.scatter(df[0], df[1], s=1)
        plt.ylabel('tSNE1')
        plt.xlabel('tSNE2')

        plt.savefig(self._alonacell.alonabase.get_working_dir() + \
            OUTPUT['FILENAME_CELL_SCATTER_PLOT'], bbox_inches='tight')
        plt.close()
