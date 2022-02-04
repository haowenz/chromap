[![CI](https://github.com/haowenz/chromap/actions/workflows/ci.yml/badge.svg)](https://github.com/haowenz/chromap/actions/workflows/ci.yml)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/chromap/badges/license.svg)](https://anaconda.org/bioconda/chromap)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/chromap/badges/version.svg)](https://anaconda.org/bioconda/chromap)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/chromap/badges/platforms.svg)](https://anaconda.org/bioconda/chromap)

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
  - [Citing Chromap](#cite)

## <a name="uguide"></a>User Guide

Chromap is an ultrafast method for aligning and preprocessing high throughput
chromatin profiles. Typical use cases include: (1) trimming sequencing adapters,
mapping bulk ATAC-seq or ChIP-seq genomic reads to the human genome and removing
duplicates; (2) trimming sequencing adapters, mapping single cell ATAC-seq genomic
reads to the human genome, correcting barcodes, removing duplicates and performing
Tn5 shift; (3) split alignment of Hi-C reads against a reference genome. In all
these three cases, Chromap is 10-20 times faster while being accurate.

### <a name="install"></a>Installation

To compile from the source, you need to have the GCC compiler with version>=7.3.0, 
GNU make and zlib development files installed. Then type `make` in the source code
directory to compile. 

Chromap is also available on [bioconda][bioconda]. 
Thus you can easily install Chromap with Conda
```sh
conda install -c bioconda -c conda-forge chromap
```

### <a name="general"></a>General usage
Before mapping, an index of the reference needs to be created and saved on the disk:

```sh
chromap -i -r ref.fa -o index
```
The users can input the min fragment length expected in their sequencing experiments,
e.g. read length, by **--min-frag-length**. Then Chromap will choose proper k-mer 
length and window size to build the index. For human genome, it only takes a few 
minutes to build the index. Without any preset parameters, Chromap takes a reference 
database and a query sequence file as input and produce approximate mapping, without
base-level alignment in the [BED format][bed]:

```sh
chromap -x index -r ref.fa -1 query.fq -o approx-mapping.bed
```
You may ask Chromap to output alignments in the [SAM format][sam]:

```sh
chromap -x index -r ref.fa -1 query.fq --SAM -o alignment.sam
```
But note that the the processing of SAM files is not fully optimized and can be slow. 
Thus generating the output in SAM format is not preferred and should be avoided when 
possible. Chromap can take multiple input read files:

```sh
chromap -x index -r ref.fa -1 query1.fq,query2.fq,query3.fq --SAM -o alignment.sam
```
Chromap works with gzip'd FASTA and FASTQ formats as input. You don't need to convert 
between FASTA and FASTQ or decompress gzip'd files first. 

***Importantly***, it should be noted that once you build the index, indexing
parameters such as **-k**, **-w** and **--min-frag-length** can't be changed during
mapping. If you are running Chromap for different data types, you will
probably need to keep multiple indexes generated with different parameters.
This makes Chromap different from BWA which always uses the same index
regardless of query data types. Chromap can build the human genome index file in a few minutes.

Detailed explanations for the options can be found at the [manpage][manpage].

### <a name="cases"></a>Use cases

To support different data types (e.g. ChIP-seq, Hi-C, ATAC-seq),
Chromap needs to be tuned for optimal performance and accuracy. It is usually
recommended to choose a preset with option **--preset**, which sets multiple
parameters at the same time.

#### <a name="map-chip"></a>Map ChIP-seq short reads

```sh
chromap --preset chip -x index -r ref.fa -1 read1.fq.gz -2 read2.fq.gz -o aln.bed      # ChIP-seq reads
```
This set of parameters is tuned for mapping ChIP-seq reads. Chromap will map the 
paired-end reads with max insert size up to 2000 (**-l 2000**) and then remove
duplicates (**--remove-pcr-duplicates**) using the low memory mode
(**--low-mem**). The output is in BED format (**--BED**). In the output BED file,
each row is a mapping of a fragment (i.e., a read pair) and the columns are

    chrom chrom_start chrom_end N mapq strand
The strand here is the strand of the first read in a read pair (specified by **-1**).
If the mapping start and end locations of each read in a read pair are desired,
**--TagAlign** should be used to overide **--BED** in the preset parameters as following
```sh
chromap --preset chip -x index -r ref.fa -1 read1.fq.gz -2 read2.fq.gz --TagAlign -o aln.tagAlign      # ChIP-seq reads
```
For each read pair, there will be two rows in the output file, one for each read in the pair
respectively. The meaning of the columns remains the same.

#### <a name="map-atac"></a>Map ATAC-seq/scATAC-seq short reads

```sh
chromap --preset atac -x index -r ref.fa -1 read1.fq.gz -2 read2.fq.gz -o aln.bed      # ATAC-seq reads
chromap --preset atac -x index -r ref.fa -1 read1.fq.gz -2 read2.fq.gz -o aln.bed\
 -b barcode.fq.gz --barcode-whitelist whitelist.txt                                    # scATAC-seq reads
```
This set of parameters is tuned for mapping ATAC-seq/scATAC-seq reads.
Chromap will trim the adapters on 3' end (**--trim-adapters**), map the 
paired-end reads with max insert size up to 2000 (**-l 2000**) and then
remove duplicates at cell level (**--remove-pcr-duplicates-at-cell-level**).
Tn5 shift will also be applied to the fragments (**--Tn5-shift**). The 
forward mapping start positions are increased by 4bp and the reverse
mapping end positions are decreased by 5bp. The processing is run in
the low memory mode (**--low-mem**).

If no barcode whitelist file is given, Chromap will skip barcode correction. 
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
is a comma-separated string, each field in the string is also a semi-comma-splitted string

    [r1|r2|bc]:start:end:strand
The start and end are inclusive and -1 means the end of the read. The strand is presented by '+' and '-' symbol, if '-' the barcode will be reverse-complemented after extraction. The strand symbol can be omitted if it is '+' and is ignored on r1 and r2. For example,
when the barcode is in the first 16bp of read1, one can use the option 
`-1 read1.fq.gz -2 read2.fq.gz --barcode read1.fq.gz --read-format bc:0:15,r1:16:-1`

The output file formats for bulk and single-cell data are different except for the first
three columns. For bulk data, the columns are

    chrom chrom_start chrom_end N mapq strand
For single-cell data, the columns are 
    
    chrom chrom_start chrom_end barcode duplicate_count
the same as the definition of the fragment file in [CellRanger][cellranger]. 
Note that chrom_end is open-end. This output fragment file can be used as input of
downstream analysis tools such as [MAESTRO][MAESTRO], [ArchR][ArchR], [signac][signac]
and etc.

Besides, Chromap can translate input cell barcodes to another set of barcodes. 
Users can specify the translation file through the option **--barcode-translate**. 
The translation file is a two-column tsv/csv file with the translated barcode on the first column
and the original barcode on the second column. This is useful for 10x Multiome data, where
scATAC-seq and scRNA-seq data use different sets of barcodes. This option also supports
combinatorial barcoding, such as SHARE-seq. Chromap can translate each barcode segment provided
in the second column to the ID in the first column and add "-" to concatenate the IDs in the
output.

#### <a name="map-hic"></a>Map Hi-C short reads

```sh
chromap --preset hic -x index -r ref.fa -1 read1.fa -2 read2.fa -o aln.pairs           # Hi-C reads and pairs output
```
Chromap will perform split alignment (**--split-alignment**) on Hi-C reads and output mappings
in [pairs][pairs] format (**--pairs**), which is used in [4DN Hi-C data processing pipeline][4DN]. 
Some Hi-C data analysis pipelines may require the reads are sorted in specific chromosome order
other than the one in the index. Therefore, Chromap provides the option **--chr-order** 
to specify the alignment order, and **--pairs-natural-chr-order** for flipping the pair 
in the pairs format. 

### <a name="help"></a>Getting help

Detailed description of Chromap command line options and optional tags 
can be displayed by running Chromap with **-h** or be found at the
[manpage][manpage]. If you encounter bugs or have further questions or 
requests, you can raise an issue at the [issue page][issue].

### <a name="cite"></a>Citing Chromap

If you use Chromap, please cite:

> Zhang, H., Song, L., Wang, X., Cheng, H., Wang, C., Meyer, C. A., ... & Li, H. (2021). Fast alignment and preprocessing of chromatin profiles with Chromap. Nature communications, 12(1), 1-6.
> https://doi.org/10.1038/s41467-021-26865-w

[bed]: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
[paf]: https://github.com/lh3/miniasm/blob/master/PAF.md
[sam]: https://samtools.github.io/hts-specs/SAMv1.pdf
[pairs]: https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md
[4DN]: https://data.4dnucleome.org/resources/data-analysis/hi_c-processing-pipeline
[minimap]: https://github.com/lh3/minimap
[release]: https://github.com/haowenz/chromap/releases
[issue]: https://github.com/haowenz/chromap/issues
[cellranger]: https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments
[manpage]: https://zhanghaowen.com/chromap/chromap.html
[bioconda]: https://anaconda.org/bioconda/chromap
[ArchR]: https://www.archrproject.com/index.html
[MAESTRO]: https://github.com/liulab-dfci/MAESTRO
[signac]: https://satijalab.org/signac/articles/pbmc_vignette.html
