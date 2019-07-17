# Implementation of the DoubletFinder algorithm for scRNA-seq doublet detection.
# Preprint: https://www.biorxiv.org/content/biorxiv/early/2018/07/19/352484.full.pdf
# Published: https://doi.org/10.1016/j.cels.2019.03.003

import random
from ctypes import cdll
import ctypes
import numpy.ctypeslib as npct

import pandas as pd
import numpy as np
from sklearn import decomposition

import alona.cell
from alona.hvg import AlonaHighlyVariableGenes

input_mat = '/home/rand/Bioinformatik/Proj/single_cell_db/alona/python-alona/q/input.mat.C'

# non-normalized read counts
data_counts = pd.read_csv(input_mat, sep='\t')
data_counts.index = data_counts['cellid']
data_counts = data_counts.drop('cellid', axis=1)

prop_doublets = 0.25
n_doublets = round(data_counts.shape[1]*prop_doublets)

# simulate doublets (done at the count level in the published paper)
pairs = [(random.randint(1,data_counts.shape[1]),
          random.randint(1,data_counts.shape[1])) for i in range(n_doublets)]
pairs = pd.DataFrame(pairs)

cells1 = data_counts.iloc[:, pairs[0]-1]
cells2 = data_counts.iloc[:, pairs[1]-1]

l0 = ['doublet-']*n_doublets
l1 = list(range(0, n_doublets))

doublet_ids = [ i+str(j) for i,j in zip(l0,l1) ]

# add real counts together to create the simulated doublets
doublets = pd.DataFrame(cells1.values+cells2.values, index=cells1.index,
                        columns = doublet_ids)
                        
# combine real and simulated data
merged_data_counts = pd.concat([data_counts, doublets], axis=1)

# normalize merged data
merged_data_norm = alona.cell.AlonaCell.normalization(None, data=merged_data_counts,
                                                      mrnafull = False, input_type = 'raw',
                                                      remove_low_quality=False)

pca = decomposition.PCA(n_components=10)
pca.fit(merged_data_norm.transpose())
comp = pca.transform(merged_data_norm.transpose())

libpath = '/home/rand/Bioinformatik/Proj/single_cell_db/alona/python-alona/alona/ANN/annlib.so'
lib = cdll.LoadLibrary(libpath)

npa = np.ndarray.flatten(comp)
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

no_cells = comp.shape[0]
no_comps = comp.shape[1]

k = round(merged_data_norm.shape[1]**0.5)

out_index = np.zeros(no_cells*k, 'int32')
out_dists = np.zeros(no_cells*k)

lib.get_NN_2Set(npa2,
                npa2,
                ctypes.c_int(no_comps),
                ctypes.c_int(no_cells),
                ctypes.c_int(no_cells),
                ctypes.c_int(k),
                ctypes.c_double(0),  # EPS
                ctypes.c_int(1),     # SEARCHTYPE
                ctypes.c_int(1),     # USEBDTREE
                np.ndarray(0),       # SQRAD
                out_index,
                out_dists)

out_index_mat = np.reshape(out_index, (no_cells, k))
knn_graph = pd.DataFrame(out_index_mat)
melted = knn_graph.melt(id_vars=[0])[[0, 'value']]

col1 = []
col2 = []

# Calculate doublet scores
doublet_rate = 0.10 # 0.05-0.1
for idx in range(1, data_counts.shape[1]+1):
    ss = melted[melted[melted.columns[0]]==idx]
    # count number of simulated neighbors
    no_sim = np.sum(ss['value']>data_counts.shape[1])
    q = no_sim/k
    
    col1.append(idx)
    col2.append(q)

pd.DataFrame({ 'status' : col1, 'score' : col2}).to_csv('score.txt',sep='\t', index=False)
