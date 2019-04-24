# About this document
This is a brief tutorial on how to pre-process a single cell RNA-seq dataset from [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra/) for input into `alona` (https://alona.panglaodb.se). Your data will likely not come from the NCBI SRA, but it serves as an example for this tutorial and it can easily be changed to fit your data. If nothing else is specified, all code are terminal commands.

# Contact
* Oscar Franzén p.oscar.franzen@gmail.com

# Installation
## Aligner (HISAT2)
There are multiple aligners. Here we will use [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) - a not too memory-hungry aligner. Download the source code, unpack and compile it:

```bash
mkdir ~/Temp

cd ~/Temp

wget http://ccb.jhu.edu/software/hisat2/dl/hisat2-2.1.0-source.zip

unzip hisat2-2.1.0-source.zip

cd hisat2-2.1.0

make

# wait for a while, on my computer it took ~5 min.

cd ..

rm -v hisat2-2.1.0-source.zip
```

# Prepare the reference genome
Download the reference genome of your species. Here I will download and build an index of the mouse `GRCm38` genome. It is important not to use the reference genome containing complete haplotype sequences, because otherwise some genes located in these blocks will get zero expression as the aligner flag the corresponding reads as multimappers.

```
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/GRCm38.primary_assembly.genome.fa.gz

gunzip GRCm38.primary_assembly.genome.fa.gz

# check number of processors
cat /proc/cpuinfo | grep processor | wc -l

# build the index with 2 CPUs, runtime might be an hour or more.
~/Temp/hisat2-2.1.0/hisat2-build -p 2 GRCm38.primary_assembly.genome.fa GRCm38.primary_assembly.genome.fa.hisat2
```
