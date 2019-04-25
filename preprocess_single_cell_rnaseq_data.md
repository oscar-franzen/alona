# About this document
This is a brief tutorial on how to pre-process a single cell RNA-seq dataset from [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra/) for input into `alona` (https://alona.panglaodb.se). Your data will likely not come from the NCBI SRA, but it serves as an example for this tutorial and it can easily be changed to fit your data. If nothing else is specified, all code are terminal commands.

# Contact
* Oscar Franzén p.oscar.franzen@gmail.com

# Installation
## sratoolkit
Only used for interacting with the NCBI SRA. Choose your system version from https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/. Here I will use the Linux version:

```bash
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz

# unpack
tar -zxvf sratoolkit.current-ubuntu64.tar.gz

rm -v sratoolkit.current-ubuntu64.tar.gz
```

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
Download the reference genome of your species. Here I will download and build an index of the mouse `GRCm38` genome. It is important not to use the reference genome containing complete haplotype sequences, because otherwise some genes located in these blocks will get zero expression as the aligner flag the corresponding reads as multimappers. Finally, to increase alignment sensitivity around splice junctions, you might instead want to consider using an aligner such as [STAR](https://github.com/alexdobin/STAR), which can use exisiting genome annotations when creating the index to improve alignment accuracy around splice junctions.

```bash
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/GRCm38.primary_assembly.genome.fa.gz

gunzip GRCm38.primary_assembly.genome.fa.gz

# check number of processors
cat /proc/cpuinfo | grep processor | wc -l

# build the index with 2 CPUs, runtime might be an hour or more.
~/Temp/hisat2-2.1.0/hisat2-build -p 2 GRCm38.primary_assembly.genome.fa GRCm38.primary_assembly.genome.fa.hisat2
```

# Preparing genome annotations
```bash
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/gencode.vM21.annotation.gtf.gz

gunzip gencode.vM21.annotation.gtf.gz

# collapse transcripts into "meta genes"
python3 gencode_meta_genes.py > gencode.vM21.annotation.meta_genes.gtf
```

# Download data for this tutorial (or use your own)
We will use the mouse lung dataset [SRS4031561](https://www.ncbi.nlm.nih.gov/sra/?term=SRR8176398), which was generated with Drop-seq [[1]](https://www.cell.com/abstract/S0092-8674(15)00549-8). The dataset consists of two sequencing runs, but for this example we will just use one of the runs:

```bash
./sratoolkit.2.9.6-ubuntu64/bin/prefetch SRR8176398

mv -v ~/ncbi/public/sra/SRR8176398.sra .

# dump it to a fastq file
./sratoolkit.2.9.6-ubuntu64/bin/fastq-dump SRR8176398.sra
```

# Prepare the fastq data
Identify and extract barcodes and unique molecular identifiers (UMI) and move them to the fastq header. In this experiment barcodes are 12 bp long and located at the start of the sequence, the next 8 bp is the UMI, and the remaining sequence is the transcript:

```
CATGAGTTCGTACGTGGATCTTTTTTTTTGTTGGGGGAGGTAATGATGAGGCTAGGTAAGTGAAGGTGGATTTGGCAACTG
^^^^^^^^^^
 barcode
          ^^^^^^^^
            UMI
                  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                                              transcript
```

## Count barcodes
First count the number of available barcodes. We expect some noise so we don't need to deal with barcodes with very few sequences. This step can be performed with any scripting language, here we will use `awk` - a swiss-knife language for text manipulation and analysis. Set `bc_len` to the barcode length.

```bash
awk -v bc_len=12 '$0~/^\@/ { getline; lines[substr($0,0,bc_len)]++; getline; getline; } END { for (i in lines) { print(i,lines[i]) } }' SRR8176398.fastq > SRR8176398.fastq.bc
```

The output will be a file with two columns (1. the barcode; 2. number of sequences with the barcode).
