"""
 alona

 Description:
 An analysis pipeline for scRNA-seq data.

 How to use:
 https://github.com/oscar-franzen/alona/

 Contact:
 Oscar Franzen <p.oscar.franzen@gmail.com>
"""

ORANGE = '#E69F00'

# Output files
OUTPUT = {
    'FILENAME_BARPLOT_RAW_READ_COUNTS' : '/plots/barplot_rrc.pdf',
    'FILENAME_BARPLOT_GENES_EXPRESSED' : '/plots/barplot_ge.pdf',
    'FILENAME_CELL_SCATTER_PLOT_PREFIX' : '/plots/2d_plot_',
    'FILENAME_CELL_VIOLIN_GE_PLOT' : '/plots/ge_violin.pdf',
    'FILENAME_CELL_VIOLIN_TOP' : '/plots/ge_violin_top.pdf',
    'FILENAME_MARKER_HEATMAP' : '/plots/marker_heatmap.png',
    'FILENAME_PCA' : '/csvs/pca.csv',
    'FILENAME_EMBEDDING_PREFIX' : '/csvs/embeddings_',
    'FILENAME_HVG' : '/csvs/highly_variable_genes.tsv',
    'FILENAME_ALL_T_TESTS' : '/csvs/all_t_tests.csv',
    'FILENAME_ALL_T_TESTS_LONG' : '/csvs/all_t_tests_long.tsv',
    'FILENAME_SNN_GRAPH' : '/csvs/snn_graph.csv',
    'FILENAME_CLUSTERS_LEIDEN' : '/csvs/clusters_leiden.csv',
    'FILENAME_MEDIAN_EXP' : '/csvs/median_exp.tsv',
    'FILENAME_MEAN_EXP' : '/csvs/mean_exp.tsv',
    'FILENAME_CTA_RANK_F' : '/csvs/CTA_RANK_F/cell_type_pred_full_table.txt',
    'FILENAME_CTA_RANK_F_BEST' : '/csvs/CTA_RANK_F/cell_type_pred_best.txt',
    'FILENAME_SETTINGS' : '/settings.txt',
    'FILENAME_QC_SCORE' : '/csvs/Mahalanobis.csv',
    'FILENAME_KNN_map' : '/KNN.joblib'
}

# Reference data
GENOME = {
    'ENTREZ_GENE_IDS' : '/genome/MGI_Gene_Model_Coord.txt.C',
    'MOUSE_GENOME_ANNOTATIONS' : '/genome/Mus_musculus.GRCm38.gencode.vM17.primary_assembly.annotation.gene_level.ERCC.gtf',
    'HUMAN_GENE_SYMBOLS_TO_ENTREZ' : '/genome/hgnc_complete_set.txt',
    'MOUSE_VS_HUMAN_ORTHOLOGS' : '/genome/human_to_mouse_1_to_1_orthologs.tsv',
    'MOUSE_EXON_LENGTHS' : '/genome/Mus_musculus.GRCm38.gencode.vM17.primary_assembly.annotation.gene_level.gtf.exon_lengths',
    'MOUSE_GENE_SYMBOLS' : '/genome/mouse_gene_symbols.txt'
}

MARKERS = {
    'PANGLAODB' : '/markers/markers.tsv'
}

# For the terminal
# https://github.com/s0md3v/Photon/blob/master/core/colors.py
WHITE = '\033[97m'
GREEN = '\033[92m'
RED = '\033[91m'
YELLOW = '\033[93m'
END = '\033[0m'
BACK = '\033[7;91m'
INFO = '\033[93m[!]\033[0m'
QUE = '\033[94m[?]\033[0m'
BAD = '\033[91m[-]\033[0m'
GOOD = '\033[92m[+]\033[0m'
RUN = '\033[97m[~]\033[0m'
