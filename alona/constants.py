ORANGE = '#E69F00'

# Output files
OUTPUT = {
    'FILENAME_BARPLOT_RAW_READ_COUNTS' : '/plots/barplot_rrc.pdf',
    'FILENAME_BARPLOT_GENES_EXPRESSED' : '/plots/barplot_ge.pdf',
    'FILENAME_CELL_SCATTER_PLOT' : '/plots/scatter.pdf',
    
    'FILENAME_EMBEDDINGS' : '/csvs/embeddings.csv',
    'FILENAME_HVG' : '/csvs/highly_variable_genes.csv'
}

# Reference data
GENOME = {
    'ENTREZ_GENE_IDS' : './genome/MGI_Gene_Model_Coord.txt.C',
    'MOUSE_GENOME_ANNOTATIONS' : './genome/Mus_musculus.GRCm38.gencode.vM17.primary_assembly.annotation.gene_level.ERCC.gtf',
    'HUMAN_GENE_SYMBOLS_TO_ENTREZ' : './genome/hgnc_complete_set.txt',
    'MOUSE_VS_HUMAN_ORTHOLOGS' : './genome/human_to_mouse_1_to_1_orthologs.tsv',
    'MOUSE_EXON_LENGTHS' : './genome/Mus_musculus.GRCm38.gencode.vM17.primary_assembly.annotation.gene_level.gtf.exon_lengths'
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
