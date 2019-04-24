# About this document
This is a brief tutorial on how to pre-process a single cell RNA-seq dataset from [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra/) for input into `alona` (https://alona.panglaodb.se). Your data will likely not come from the NCBI SRA, but it serves as an example for this tutorial and it is easy to change.

# Installation
## Aligner (HISAT2)
There are multiple aligners. Here we will use [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) - a not too memory-hungry aligner. Download the source code, unpack and compile it:

```bash
wget http://ccb.jhu.edu/software/hisat2/dl/hisat2-2.1.0-source.zip

unzip hisat2-2.1.0-source.zip

cd hisat2-2.1.0

make

# wait for a while, on my computer it took ~20 min.
```
