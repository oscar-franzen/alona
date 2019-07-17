# Implementation of the Scrublet algorithm for scRNA-seq doublet detection.
# Supposedly to function as a doublet detecting routine in alona. After testing it seems
# like the proposed algorithm is very sensitive to the initial number of specified
# PCA components. I'm saving this code for putative future use; the code may be
# incorporated into alona at a later point.
#
# Attention: preprint is not identical to the published paper.
# Preprint:  https://www.biorxiv.org/content/biorxiv/early/2018/07/09/357368.full.pdf
# Published: https://www.cell.com/cell-systems/pdfExtended/S2405-4712(18)30474-5
#
# Oscar Franzen <p.oscar.franzen@gmail.com>, July 2019

import random
import joblib

from ctypes import cdll
import ctypes
import numpy.ctypeslib as npct
from scipy.sparse import coo_matrix
from sklearn import decomposition
import pandas as pd
import numpy as np
import alona.irlbpy
import alona.cell
from alona.hvg import AlonaHighlyVariableGenes

input_mat = 'input_matrix.csv'

# non-normalized read counts
data_counts = pd.read_csv(input_mat, sep='\t')
data_counts.index = data_counts['cellid']
data_counts = data_counts.drop('cellid', axis=1)

# normalize
data_norm = alona.cell.AlonaCell.normalization(None, data=data_counts, mrnafull = False,
                                               input_type = 'raw',
                                               remove_low_quality=False)

hvg_finder = AlonaHighlyVariableGenes(hvg_method='seurat',
                                              hvg_n=1000,
                                              data_norm=data_norm,
                                              data_ERCC=None)

hvg = hvg_finder.find()

data_norm_hvg = data_norm[data_norm.index.isin(hvg)]

# Z-score (z=(x-mean)/sd)
data_norm_zscore = data_norm_hvg.apply(lambda x : (x-x.mean())/x.std(), axis=1)
gene_mean = data_norm_hvg.mean(axis=1)
gene_std = data_norm_hvg.std(axis=1)

# PCA
pca = decomposition.PCA(n_components=30)
pca.fit(data_norm_zscore.transpose())
eigen1 = pca.transform(data_norm_zscore.transpose())
# same as:
# data_norm_zscore.transpose().dot(pca.components_.transpose())

# the number of doublets to simulate, relative to the number of observed transcriptomes
sim_doublet_ratio = 2
n_doublets = data_norm_hvg.shape[1]*sim_doublet_ratio

# simulate doublets (done at the count level in the published paper)
pairs = [(random.randint(1,data_norm_hvg.shape[1]),
          random.randint(1,data_norm_hvg.shape[1])) for i in range(n_doublets)]
pairs = pd.DataFrame(pairs)

cells1 = data_counts.iloc[:, pairs[0]-1]
cells2 = data_counts.iloc[:, pairs[1]-1]

l0 = ['doublet-']*n_doublets
l1 = list(range(0, n_doublets))

doublet_ids = [ i+str(j) for i,j in zip(l0,l1) ]

# add real counts together to create the simulated doublets
doublets = pd.DataFrame(cells1.values+cells2.values, index=cells1.index,
                        columns = doublet_ids)

# normalize simulated doublets
doublets_norm = alona.cell.AlonaCell.normalization(None, data=doublets, mrnafull = False,
                                                   input_type = 'raw',
                                                   remove_low_quality=False)

# take hvgs from observed cells
doublets_norm_hvg = doublets_norm[doublets_norm.index.isin(hvg)]

# apply zscore using the same means and sd as measured in the observed transcriptomes
doublets_norm_hvg = doublets_norm_hvg.sub(gene_mean, axis=0)
doublets_norm_hvg_zscore = doublets_norm_hvg.div(gene_std, axis=0)

# project doublets in the same PCA space as the observed cells
eigen2 = pca.transform(doublets_norm_hvg_zscore.transpose())

# union of observed cells and simulated doublets
merged_mat = np.concatenate((eigen1,eigen2))

# NN search
k = round(0.5*merged_mat.shape[0]**0.5)
k_adj = round(k*(1+sim_doublet_ratio))

libpath = 'annlib.so'
lib = cdll.LoadLibrary(libpath)

#pca_rotated = np.rot90(self.pca_components)
npa = np.ndarray.flatten(merged_mat)
npa2 = npa.astype(np.double).ctypes.data_as(ctypes.POINTER(ctypes.c_double))

lib.get_NN_2Set.restype = None
lib.get_NN_2Set.argtypes = [ctypes.POINTER(ctypes.c_double),
                            ctypes.POINTER(ctypes.c_double),
                            ctypes.POINTER(ctypes.c_int),
                            ctypes.POINTER(ctypes.c_int),
                            ctypes.POINTER(ctypes.c_int),
                            ctypes.POINTER(ctypes.c_int),
                            ctypes.POINTER(ctypes.c_double),
                            ctypes.POINTER(ctypes.c_int),
                            ctypes.POINTER(ctypes.c_int),
                            npct.ndpointer(dtype=np.double, ndim=1,
                                           flags='CONTIGUOUS'),
                            npct.ndpointer(ctypes.c_int),
                            npct.ndpointer(dtype=np.double, ndim=1,
                                           flags='CONTIGUOUS')]
no_cells = merged_mat.shape[0]
no_comps = merged_mat.shape[1]

out_index = np.zeros(no_cells*k_adj, 'int32')
out_dists = np.zeros(no_cells*k_adj)

lib.get_NN_2Set(npa2,
                npa2,
                ctypes.c_int(no_comps),
                ctypes.c_int(no_cells),
                ctypes.c_int(no_cells),
                ctypes.c_int(k_adj),
                ctypes.c_double(0),  # EPS
                ctypes.c_int(1),     # SEARCHTYPE
                ctypes.c_int(1),     # USEBDTREE
                np.ndarray(0),       # SQRAD
                out_index,
                out_dists)

out_index_mat = np.reshape(out_index, (no_cells, k_adj))
#out_dists_mat = np.reshape(out_dists, (no_cells, k_adj))

knn_graph = pd.DataFrame(out_index_mat)
melted = knn_graph.melt(id_vars=[0])[[0, 'value']]

col1 = []
col2 = []

# Calculate doublet scores
doublet_rate = 0.10 # 0.05-0.1
for idx in range(1, merged_mat.shape[0]):
    ss = melted[melted[melted.columns[0]]==idx]
    # count number of simulated neighbors
    no_sim = np.sum(ss['value']>data_norm_zscore.shape[1])
    q = (no_sim+1)/(k_adj+2)
    score = ((q*doublet_rate)/sim_doublet_ratio)/(1-doublet_rate-q*(1-doublet_rate-(doublet_rate/sim_doublet_ratio)))
    
    if idx>data_norm_zscore.shape[1]:
        status = 'sim'
    else:
        status = 'obs'
        
    col1.append(status)
    col2.append(score)

pd.DataFrame({ 'status' : col1, 'score' : col2}).to_csv('score.txt',sep='\t', index=False)
