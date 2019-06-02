```
  ##   #       ####  #    #   ##   
 #  #  #      #    # ##   #  #  #    A pipeline for cell type prediction from
#    # #      #    # # #  # #    #        single cell RNA sequencing data.
###### #      #    # #  # # ###### 
#    # #      #    # #   ## #    #           "Life is better at the beach"
#    # ######  ####  #    # #    #          
```

# Description
`alona` is a Python-based software pipeline for analysis of single cell RNA sequencing data. `alona` performs normalization, quality control, clustering and cell type annotation of single cell RNA-seq data ([1][1]).

`alona` also exists as a parallel cloud-based service ([2][2]).

# Installation
### Requirements
* Linux (alona should work on MacOS too, but it is untested)
* Python >= 3.6

### From GitHub and pip3
```bash
# Clone the repository
git clone https://github.com/oscar-franzen/alona/

# Enter the directory
cd alona

# Install the package
pip3 install .
```

# Usage
```bash
python3 -m alona \
        --dark_bg \
        --hvg 2000 \
        --leiden_res 0.1 \
        --output test \
        --loglevel debug \
        --header yes \
        --minexpgenes 0.001 \
        --nomito input.mat
```

# All command line options
```bash
$[~/bio]> python3 -m alona --help
Usage: alona.py [OPTIONS] FILENAME

Options:
  -o, --output TEXT               Specify name of output directory
  -df, --dataformat [raw|rpkm|log2]
                                  Data format.
                                  (raw read counts, rpkm, log2
                                  normalized data). Default: raw
  -mr, --minreads INTEGER         Minimum number of reads per cell to keep the
                                  cell. Default: 1000
  -mg, --minexpgenes FLOAT        Minimum number of expressed genes as percent
                                  of all cells, i.e. genes expressed in fewer
                                  cells than this are removed. Default: 0.01
  --mrnafull                      Data come from a full-length protocol, such
                                  as SMART-seq2.
  -d, --delimiter [auto|tab|space]
                                  Data delimiter. The character used to
                                  separate data values. The default setting is
                                  to autodetect this character. Default: auto
  -h, --header [auto|yes|no]      Data has a header line. The default setting
                                  is to autodetect if a header is present or
                                  not. Default: auto
  -m, --nomito                    Exclude mitochondrial genes from analysis.
  --hvg [seurat|brennecke]        Method to use for identifying highly
                                  variable genes.Default: seurat
  --hvg_n INTEGER                 Number of top highly variable genes to use.
                                  Default: 1000
  --nn_k INTEGER                  k in the nearest neighbour search. Default:
                                  10
  --prune_snn FLOAT               Threshold for pruning the SNN graph, i.e.
                                  the edges with lower value (Jaccard index)
                                  than this will be removed. Set to 0 to
                                  disable pruning. Increasing this value will
                                  result in fewer edges in the graph. Default:
                                  0.067
  --leiden_partition [RBERVertexPartition|ModularityVertexPartition]
                                  Partitioning algorithm to use. Can be
                                  RBERVertexPartition or
                                  ModularityVertexPartition. Default:
                                  RBERVertexPartition
  --leiden_res FLOAT              Resolution parameter for the Leiden
                                  algorithm (0-1). Default: 0.8
  --ignore_small_clusters INTEGER
                                  Ignore clusters with fewer or equal to N
                                  cells. Default: 1
  --perplexity INTEGER            The perplexity parameter in the t-SNE
                                  algorithm. Default: 30
  -s, --species [human|mouse]     Species your data comes from. Default: mouse
  --dark_bg                       Use dark background in scatter plots.
                                  Default: False
  --color_labels TEXT             Plot cell type labels with the same color as
                                  the corresponding cell cluster cells.
                                  Default: True
  --cleanup                       Perform cleanup of temporary files.
  -lf, --logfile TEXT             Name of log file. Set to /dev/null if you
                                  want to disable logging to a file. Default:
                                  alona.log
  -ll, --loglevel [regular|debug]
                                  Set how much runtime information is written
                                  to the log file. Default: regular
  -n, --nologo                    Hide the logo.
  --version                       Display version number.
  --help                          Show this message and exit.
```

# Detailed help for all command line options
option | detailed description
--- | ---
`-out, --output [TEXT]` | Specify name of output directory
`-df, --dataformat [raw\|rpkm\|log2]` | test
`--mrnafull` | Data come from a full-length protocol, such as SMART-seq2. This option is important if data represent full mRNAs. Drop-seq/10X and similar protocols sequence the *ENDS* of an mRNA, it is therefore not necessary to normalize for gene *LENGTH*. However, if we sequence the complete mRNA then we must also normalize measurements for the length of the gene, since longer genes have more mapped reads. If this option is not set, then cell type prediction may give unexpected results when analyzing full-length mRNA data.
`--hvg [seurat|brennecke]` | Method to use for identifying highly variable genes. This option specifies the method to be used for identifying variable genes. Default: seurat.

## Contact
* Oscar Franzen <p.oscar.franzen@gmail.com>

## Cite
A manuscript is in preparation.

# License
GPLv3

[1]: https://en.wikipedia.org/wiki/Single-cell_transcriptomics
[2]: http://alona.panglaodb.se/
