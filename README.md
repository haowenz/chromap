## <a name="started"></a>Getting Started
```sh
git clone https://github.com/haowenz/chromap.git
cd chromap && make
# create an index first and then map
./chromap -i -r test/ref.fa -o ref.index
./chromap -x ref.index -r test/ref.fa -1 test/read1.fq -2 test/read2.fq -o test.bed
# use presets (no test data)
./chromap --preset atac -x index -r ref.fa -1 read1.fq -2 read2.fq -o aln.bed       # ATAC-seq reads
./chromap --preset atac -x index -r ref.fa -1 read1.fq -2 read2.fq -o aln.bed \
 -b barcode.fq.gz --barcode-whitelist whitelist.txt                                       # scATAC-seq reads
./chromap --preset chip -x index -r ref.fa -1 read1.fq -2 read2.fq -o aln.bed       # ChIP-seq reads
./chromap --preset hic -x index -r ref.fa -1 read1.fq -2 read2.fq -o aln.pairs      # Hi-C reads and pairs output
./chromap --preset hic -x index -r ref.fa -1 read1.fq -2 read2.fq --SAM -o aln.sam  # Hi-C reads and SAM output
```
## Table of Contents

- [Getting Started](#started)
- [User Guide](#uguide)
  - [Installation](#install)
  - [General usage](#general)
  - [Use cases](#cases)
    - [Map ChIP-seq short reads](#map-chip)
    - [Map ATAC-seq/scATAC-seq short reads](#map-atac)
    - [Map Hi-C short reads](#map-hic)
  - [Getting help](#help)

## <a name="uguide"></a>User Guide

Chromap is an ultrafast method for aligning and preprocessing high throughput
chromatin profiles. Typical use cases include: (1) trimming sequencing adapters,
mapping bulk ATAC-seq or ChIP-seq genomic reads to the human genome and removing
duplicates; (2) trimming sequencing adapters, mapping single cell ATAC-seq genomic
reads to the human genome, correcting barcodes, removing duplicates and performing
Tn5 shift; (3) split alignment of Hi-C reads against a reference genome. In all
these three cases, Chromap is 10-20 times faster while being accurate.

### <a name="install"></a>Installation

You can acquire precompiled binaries from
the [release page][release] with:

```sh
curl -L https://github.com/haowenz/chromap/releases/download/v0.1/chromap-0.1_x64-linux.tar.bz2 | tar -jxvf -
./chromap-0.1_x64-linux/chromap
```
If you want to compile from the source, you need to have the GCC compiler, GNU make
and zlib development files installed. Then type `make` in the source code
directory to compile. 

Chromap is also available form conda, including [bioconda](https://anaconda.org/bioconda/chromap). You can install Chromap with `conda install -c bioconda chromap` or `conda install -c liulab-dfci chromap`.

### <a name="general"></a>General usage
Before mapping, an index of the reference needs to be created and saved on the disk:

```sh
chromap -i -r ref.fa -o index
```
The users can input the min fragment length expected in their sequencing experiments, e.g. read length, by **--min-frag-length**. Then Chromap will choose proper k-mer length and window size to build the index. For human genome, it only takes a few minutes to build the index. 
Without any preset parameters, Chromap takes a reference database and a query sequence
file as input and produce approximate mapping, without base-level alignment in the [BED format][bed]:

```sh
chromap -x index -r ref.fa -1 query.fq -o approx-mapping.bed
```
You may ask Chromap to output alignments in the [SAM format][sam]:

```sh
chromap -x index -r ref.fa -1 query.fq --SAM -o alignment.sam
```
But note that the the processing of SAM files is not fully optimized and can be slow. Thus generating the output in SAM format is not preferred and should be avoided when possible. Chromap can take multiple input read files:

```sh
chromap -x index -r ref.fa -1 query1.fq,query2.fq,query3.fq --SAM -o alignment.sam
```
Chromap works with gzip'd FASTA and FASTQ formats as input. You don't need to convert between FASTA and FASTQ or decompress gzip'd files first. 

***Importantly***, it should be noted that once you build the index, indexing
parameters such as **-k**, **-w** and **--min-frag-length** can't be changed during
mapping. If you are running Chromap for different data types, you will
probably need to keep multiple indexes generated with different parameters.
This makes Chromap different from BWA which always uses the same index
regardless of query data types. Chromap can build the human genome index file in about 10 minutes.

### <a name="cases"></a>Use cases

To support different data types (e.g. ChIP-seq, Hi-C, ATAC-seq),
Chromap needs to be tuned for optimal performance and accuracy. It is usually
recommended to choose a preset with option **--preset**, which sets multiple
parameters at the same time.

#### <a name="map-chip"></a>Map ChIP-seq short reads

```sh
chromap --preset chip -x index -r ref.fa -1 read1.fq.gz -2 read2.fq.gz -o aln.bed      # ChIP-seq reads
```
This set of parameters is tuned for mapping ChIP-seq reads. Chromap will trim the
adapters on 3' end, map the paired-end reads with max insert size (**-l**) up to
2000 and then remove duplicates.

#### <a name="map-atac"></a>Map ATAC-seq/scATAC-seq short reads

```sh
chromap --preset atac -x index -r ref.fa -1 read1.fq.gz -2 read2.fq.gz -o aln.bed      # ATAC-seq reads
chromap --preset atac -x index -r ref.fa -1 read1.fq.gz -2 read2.fq.gz -o aln.bed\
 -b barcode.fq.gz --barcode-whitelist whitelist.txt                                    # scATAC-seq reads
```
When barcodes and a whitelist are given as input, by default Chromap will
estimate barcode abundance and use this information to perform barcode
correction with up to 1 Hamming distance from a whitelist barcode. By setting
**--bc-error-threshold** to 2, Chromap is able to correct barcodes with up to 2
Hamming distance from a whitelist barcode. User can also increase the probability
threshold to make a correction by setting **--bc-probability-threshold**
(set to 0.9 by default) to a large value (e.g., 0.975) to only make reliable
corrections. For scATAC-seq data with multiple read and barcode files, you can
use "," to concatenate multiple input files as the example [above](#general). 

Chromap also supports user-defined barcode format, including mixed barcode and genomic 
data case. User can specify the sequence structure through option **--read-format**. The value
is comma-separated string, each field is also semi-comma-splitted string: [r1|r2|bc]:start:end.
The start and end(inclusive, -1 means to the read end). For the example that the barcode is in read1's 
first 16bp, one can use the option 
`-1 read1.fq.gz -2 read2.fq.gz --barcode read1.fq.gz --read-format bc:0:15,r1:16:-1`

The BED format (fragment file) for bulk and single-cell is different except for the first
three columns. For bulk data, the columns are

    chrom chrom_start chrom_end N mapq strand
For single-cell data, the columns are 
    
    chrom chrom_start chrom_end barcode duplicate_count
as the definition of the fragment file in 
[CellRanger](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments). 
Note that chrom_end is open-end.

#### <a name="map-hic"></a>Map Hi-C short reads

```sh
chromap --preset hic -x index -r ref.fa -1 read1.fa -2 read2.fa -o aln.pairs           # Hi-C reads and pairs output
```
Chromap will perform split alignment on Hi-C reads and output mappings
in [pairs][pairs] format, which is used in [4DN Hi-C data processing pipeline][4DN]. 
Some Hi-C data analysis pipelines may require the reads are sorted in specific chromosome order
other than the one in the index. Therefore, Chromap provides the option **--chr-order** 
to specify the alignment order, and **--pairs-natural-chr-order** for flipping the pair 
in the pairs format. 

### <a name="help"></a>Getting help

Detailed description of Chromap command line options and optional tags 
can be displayed by running Chromap with **-h**. If you encounter bugs or have further questions or requests,
you can raise an issue at the [issue page][issue].


[bed]: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
[paf]: https://github.com/lh3/miniasm/blob/master/PAF.md
[sam]: https://samtools.github.io/hts-specs/SAMv1.pdf
[pairs]: https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md
[4DN]: https://data.4dnucleome.org/resources/data-analysis/hi_c-processing-pipeline
[minimap]: https://github.com/lh3/minimap
[release]: https://github.com/haowenz/chromap/releases
[issue]: https://github.com/haowenz/chromap/issues
