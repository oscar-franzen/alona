# About this document
This is a brief tutorial on how to pre-process a single cell RNA-seq dataset from [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra/) for input into `alona` (https://alona.panglaodb.se). Your data will likely not come from the NCBI SRA, but it serves as an example for this tutorial and it can easily be changed to fit your data. If nothing else is specified, all code are terminal commands.

# Disclaimer
This guide serves to show the basics. There are several steps that can be optimized or changed.

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

## samtools
Samtools is used for SAM to BAM conversion and sorting.
```bash
wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2

tar -jxvf samtools-1.9.tar.bz2

cd samtools-1.9 && mkdir installed && cd installed && pwd && cd ../

./configure --prefix=/full/path/to/samtools-1.9/installed

make && make install

# run from terminal or modify your `~/.zshrc` or `~/.bashrc` and restart terminal
export PATH=${PATH}:/full/path/to/samtools-1.9/installed/bin
```

## subread
We will use `featureCounts` (part of `subread`) for counting alignments. [Htseq](https://htseq.readthedocs.io/en/release_0.11.1/) is another popular counter, but featureCounts is **much** faster. Download the source code of subread ([subread-1.6.4-source.tar.gz](https://github.com/oscar-franzen/alona/raw/master/preprocessing_tutorial/tools_src/subread-1.6.4-source.tar.gz)).

```bash
tar -zxvf subread-1.6.4-source.tar.gz && cd subread-1.6.4-source/src

# Choose your makefile flavor, I use Linux so I use:
make -f Makefile.Linux
```

## UMI-tools
Tools for collapsing unique molecular identifiers, see [this](http://genome.cshlp.org/content/early/2017/01/18/gr.209601.116.abstract) paper.
```bash
pip install umi_tools

# or (Python3)

pip3 install umi_tools
```

## R
Follow instructions from [https://www.r-project.org/](https://www.r-project.org/).

# Prepare the reference genome
Download the reference genome of your species. Here I will download and build an index of the mouse `GRCm38` genome. It is important not to use the reference genome containing complete haplotype sequences, because otherwise some genes located in these blocks will get zero expression as the aligner flag the corresponding reads as multimappers. Finally, to increase alignment sensitivity around splice junctions, you might instead want to consider using an aligner such as [STAR](https://github.com/alexdobin/STAR), which can use exisiting genome annotations when creating the index to improve alignment accuracy around splice junctions.

```bash
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/GRCm38.primary_assembly.genome.fa.gz

gunzip GRCm38.primary_assembly.genome.fa.gz

# check number of processors
cat /proc/cpuinfo | grep processor | wc -l

# build the index with 2 CPUs, runtime might be an hour or more.
time ~/Temp/hisat2-2.1.0/hisat2-build -p 2 \
      GRCm38.primary_assembly.genome.fa GRCm38.primary_assembly.genome.fa.hisat2

4817.61s user 28.27s system 202% cpu 39:49.91 total
```

# Preparing genome annotations
[`gencode_meta_genes.py`](https://github.com/oscar-franzen/alona/blob/master/preprocessing_tutorial/gencode_meta_genes.py) is available in this repository.

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
First count the number of available barcodes. We expect some noise so we don't need to deal with barcodes with very few sequences. This step can be performed with any scripting language, here we will use a one-liner with `awk` - a swiss-knife language for text manipulation. Set `bc_len` to the barcode length. **Note on awk's `substr(...)`, the first character of the  string has position 1.**

```bash
awk -v bc_len=12 \
  '$0 ~ /^\@/ {
    getline;
    lines[substr($0,1,bc_len)]++;
    getline;
    getline;
  }
  
  END {
    for (i in lines) {
      print(i,lines[i])
    }
  }' SRR8176398.fastq > SRR8176398.fastq.bc
```

The output will be a file with two columns (1. the barcode; 2. number of sequences with the barcode). If we take a look at the output file we see that there are numerous singleton barcodes as well as barcodes with just a few reads, these can be safely excluded.

Let's examine the top 10 barcodes:

```bash
sort -k 2,2nr SRR8176398.fastq.bc | head -n 10
AAAAAAAAAAAA 1930334
ACCCCGTTCCAT 841302
AGCACTCTGTTC 660082
GTGGCCGGTACA 589671
CTGTCTCTTATA 567723
ATAACATACAGG 463122
ATTGTTAGTAGG 424067
CCGACGTCTTAG 399235
TGGGTTCGTGAA 383446
CGATCCGGCACT 373765
```

The top ranking barcode "AAAAAAAAAAAA" seems suspicious, indicating that in the next step we should remove barcodes consisting of 100% one nucleotide.

## Organize the fastq file
We will move barcode and UMI to the FASTQ sequence header using the helper script [`organize_fastq.py`](https://github.com/oscar-franzen/alona/blob/master/preprocessing_tutorial/organize_fastq.py) (available in this repository).

```bash
python3 organize_fastq.py -i SRR8176398.fastq \
                            -b SRR8176398.fastq.bc \
                            -c 1000 > SRR8176398.clean.fastq
```

# Align the FASTQ file
Confirm the quality score encoding of your FASTQ data. This small script, [`determine_phred.pl`](https://raw.githubusercontent.com/oscar-franzen/alona/master/preprocessing_tutorial/determine_phred.pl), can be used to determine the quality score format (will require installation of the `List::MoreUtils` Perl module). For SRA data, this encoding is Phred33. The additional pipe in the script below is a filter, which keeps unique alignments. Filtering in this way prevents writing of uninformative (multimapping) sequence alignments to your disk, which may save time. **Note that quality scores in your FASTQ file are unrelated with the mapping qualities in the SAM file.**

```bash
time (hisat2 --phred33 -p 2 \
      -x GRCm38.primary_assembly.genome.fa.hisat2 \
      -U SRR8176398.clean.fastq | awk '$5 >= 60' > SRR8176398.sam)

10777.77s user 162.88s system 206% cpu 1:28:17.88 total
```

### A directory of FASTQ files
In many cases we don't have one single FASTQ file, rather we have a directory of multiple individual FASTQ files. Instead of writing multiple commands we can cut and paste a `for` loop in the terminal:

```bash
# create an output directory to prevent litter
mkdir aln

# data can be compressed with gzip to save space
for file in ./data/*.gz
do
    # truncates the path and first two file extensions (.fastq.gz)
    out=${file##*/}
    out=${out%.*}
    out=${out%.*}

    echo $file

    hisat2 --phred33 -p 20 -x /path/to/GRCm38.primary_assembly.genome.fa.hisat2 -U \
    ${file} | awk '$5 >= 60' > ./aln/${out}.sam
done
```

# Convert to `BAM`
`BAM` is the compressed version of `SAM`. The format is based on gzip and we use it to save space. It has become a standard format for saving sequence alignments from "NGS" technologies.

```bash
time samtools view -@ 2 \
                     -T GRCm38.primary_assembly.genome.fa \
                     -bS SRR8176398.sam > SRR8176398.bam
```

### A directory of `sam` files
Iterate through all files.

```bash
for file in ./aln/*.sam
do
samtools view -@ 2 \
                     -T /path/to/GRCm38.primary_assembly.genome.fa \
                     -bS ${file} > ${file%.*}.bam
done
```

# Sort the BAM file
Sorting the BAM file by alignment coordinates. This is a prerequisite in many downstream operations. `-T .` selects the present directory as the place to place temporary files. `-m 2G` restricts memory usage. `-@ 2` it will use two sorting threads.

```bash
time samtools sort -T . \
                     -m 2G \
                     -@ 2 SRR8176398.bam > \
                      SRR8176398.sorted.bam

630.24s user 28.51s system 229% cpu 4:47.47 total
```

### A directory of BAM files
```bash
for file in ./aln/*.bam
do
    samtools sort -T . -m 2G -@ 2 ${file} > ${file%.*}.sorted.bam
done
```

# Create an index of the BAM file
An index is used for faster lookups. It's needed if we want to open the BAM file in a browser such as the [IGV](https://software.broadinstitute.org/software/igv/).

```
time ./samtools-1.9/installed/bin/samtools index SRR8176398.sorted.bam

44.52s user 0.67s system 99% cpu 45.218 
```

### A directory of sorted BAM files
```
for file in ./aln/*.sorted.bam
do
    samtools index ${file}
done
```

# Add gene tag to SAM file
This step annotates each alignment with the gene most likely to be the transcribed source of the sequence. **Note, the output is another BAM file.**

```bash
time ./subread-1.6.4-source/bin/featureCounts -R BAM \
                                                --tmpDir . \
                                                -T 2 \
                                                -F GTF \
                                                -a gencode.vM21.annotation.meta_genes.gtf \
                                                -o SRR8176398.sorted.bam.featureCounts SRR8176398.sorted.bam

227.62s user 4.38s system 198% cpu 1:56.91 total
```

### A directory of sorted and indexed BAM files
```bash
# save output files in a separate directory
mkdir counts

for file in ./aln/*.sorted.bam
do
    echo $file
    
    featureCounts -R BAM \
                  --tmpDir . \
                  -T 2 \
                  -F GTF \
                  -a /path/to/PlasmoDB-25_Pfalciparum3D7.gtf \
                  -o ./counts/${file##*/}.featureCounts ${file} 2>/dev/null
done
```

# Sort and index again
Another round is sorting is needed.

```bash
time samtools sort -T . \
                     -m 2G \
                     -@ 2 SRR8176398.sorted.bam.featureCounts.bam > \
                      SRR8176398.sorted.bam.featureCounts.sorted.bam

631.08s user 31.53s system 230% cpu 4:47.10 total

samtools index SRR8176398.sorted.bam.featureCounts.sorted.bam
```

# umi-tools
Perform UMI deduplication.

```bash
time umi_tools dedup --no-sort-output \
                     --method unique \
                     --gene-tag=XT \
                     --per-gene \
                     --per-cell \
                     --log2stderr \
                     -I SRR8176398.sorted.bam.featureCounts.sorted.bam > \
                     SRR8176398.sorted.bam.featureCounts.sorted.umi_dedup.bam
                     
1345.55s user 4.82s system 100% cpu 22:30.08 total
```

# Sort and index final time
```bash
time samtools sort -T . \
                     -m 2G \
                     -@ 2 SRR8176398.sorted.bam.featureCounts.sorted.umi_dedup.bam > \
                      SRR8176398.sorted.bam.featureCounts.sorted.umi_dedup.sorted.bam

631.08s user 31.53s system 230% cpu 4:47.10 total

samtools index SRR8176398.sorted.bam.featureCounts.sorted.umi_dedup.sorted.bam
```

# Count reads
```bash
time umi_tools count --per-gene \
                       --gene-tag=XT \
                       --per-cell \
                       -I SRR8176398.sorted.bam.featureCounts.sorted.umi_dedup.sorted.bam \
                       -S SSRR8176398.sorted.bam.featureCounts.sorted.umi_dedup.sorted.bam.counts
                       
410.74s user 4.60s system 100% cpu 6:54.94 total

# confirm output
head SSRR8176398.sorted.bam.featureCounts.sorted.umi_dedup.sorted.bam.counts
gene	cell	count
0610005C13Rik_ENSMUSG00000109644	AATAGTGGTAGC	1
0610005C13Rik_ENSMUSG00000109644	AGGTGGGGGGAT	1
0610005C13Rik_ENSMUSG00000109644	AGTGCGTTCACT	1
0610005C13Rik_ENSMUSG00000109644	CATTGCCCGATC	1
0610005C13Rik_ENSMUSG00000109644	CTCATATACGGT	1
0610005C13Rik_ENSMUSG00000109644	TCAGAAGCTGGT	1
0610005C13Rik_ENSMUSG00000109644	TCTGTCGCGCTG	1
0610005C13Rik_ENSMUSG00000109644	TTCACGGAGCAA	1
0610005C13Rik_ENSMUSG00000109644	TTGAGTATATTG	1
```

# Reshape data into a matrix
We will borrow some R code to make this table into a matrix. [`linear_to_matrix.R`](https://github.com/oscar-franzen/alona/blob/master/preprocessing_tutorial/linear_to_matrix.R) is in the repository.

```bash
R --no-save \
   --no-restore \
   --slave \
   --args SSRR8176398.sorted.bam.featureCounts.sorted.umi_dedup.sorted.bam.counts \
   SSRR8176398.mat < linear_to_matrix.R
   
# ls SSRR8176398.mat.txt SSRR8176398.mat.RDS
 19M -rw-rw-r--  1 qwerty  19M 29 Apr 2019 10:46 SSRR8176398.mat.RDS
1.1G -rw-rw-r--  1 qwerty 1.1G 29 Apr 2019 10:46 SSRR8176398.mat.txt
```
As you can see the actual matrix eats up a lot of disk space owing to being zero inflated. We don't need to write zeros if we assume missing fields are zeros. The second file ending with RDS contains the sparse R representation of this matrix, which can be loaded into R with `readRDS(...)`.

A near similar trick is just to compress the plain text table:

```bash
time gzip SSRR8176398.mat.txt
7.60s user 0.30s system 99% cpu 7.904 total

ls -lh SSRR8176398.mat.txt.gz
-rw-rw-r-- 1 qwerty qwerty 16M Apr 29 10:46 SSRR8176398.mat.txt.gz

```

# Done!
We are done. `SSRR8176398.mat.txt.gz` can be uploaded into `alona`.
