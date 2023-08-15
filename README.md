# CIRI3

[![Latest Release](https://img.shields.io/github/release/Xinglab/espresso.svg?label=Latest%20Release)](https://github.com/Xinglab/espresso/releases/latest)
[![Total GitHub Downloads](https://img.shields.io/github/downloads/Xinglab/espresso/total.svg?label=Total%20GitHub%20Downloads)](https://github.com/Xinglab/espresso/releases)
[![Total Bioconda Installs](https://img.shields.io/conda/dn/bioconda/espresso.svg?label=Total%20Bioconda%20Installs)](https://anaconda.org/bioconda/espresso)

## About

CIRI3 is a comprehensive analysis package for the detection and quantification of circRNAs in RNA-Seq data, while providing differential expression analysis of circRNAs at multiple levels.

## Table of contents

* [Dependencies](#dependencies)
* [Usage](#usage)
  + [Snakemake](#snakemake)
  + [Basic Usage](#basic-usage)
  + [Example](#example)
  + [Preparing Input Files](#preparing-input-files)
  + [All Arguments](#all-arguments)
* [Output](#output)

## Dependencies

ESPRESSO requires the following to be installed and available on $PATH:

* [samtools](http://www.htslib.org/) >= 1.6

ESPRESSO requires a sorted SAM or BAM file as input. If the data is in another format then other tools can be used to create the sorted alignment file:

* [minimap2](https://github.com/lh3/minimap2)
* ONT guppy basecaller

The [visualization](visualization/) scripts require

* [UCSC](http://hgdownload.soe.ucsc.edu/admin/exe/) tools
  + bedGraphToBigWig
  + faToTwoBit
  + twoBitInfo
* Python 3
  + [NumPy](https://numpy.org/)

If using the [Snakemake](#snakemake) then its installation script can install the dependencies (except for guppy which requires a login).

## Usage

### Step1. Mapping of the fastq files with BWA MEM
Recommended protocols for running BWA-MEM:
```
bwa index -a bwtsw ref.fa
bwa mem –T 19 ref.fa reads.fq > aln-se.sam (single-end reads)
bwa mem –T 19 ref.fa read1.fq read2.fq > aln-pe.sam (paired-end reads)
```
### Step2. CircRNA detection and quantification

CIRI3 provides multiple input options for the identification and quantification of circRNAs, including single-sample input, multiple-sample input, multiple-sample files containing RNase R treated information. In addition, users have the option to input a collection of circRNAs of interest, allowing CIRI3 to quantify the BSJ and FSJ of these circRNAs within the samples.

A small test data set is provided in [test_data/test_data_espresso_sirv.tar.gz](test_data/test_data_espresso_sirv.tar.gz). The unpacked files are:

#### 1)Single-SAM/BAM files as input
CIRI3
```
## without the annotation gtf
Java -jar CIRI3.jar -I ./data/circRNA/SAM/0/sample.sam -O ./data/circRNA/SAM/0/result.txt -F ./data/circRNA/chr1.fa
## with the annotation gtf
Java -jar CIRI3.jar -I sample.sam -O result.txt -F chr1.fa -A chr1.gtf
```

#### 2)Multiple-SAM/BAM files as input

```
## without the annotation gtf
Java -jar CIRI3.jar -I ./data/circRNA/SAM/0/sample.sam -O ./data/circRNA/SAM/0/result.txt -F ./data/circRNA/chr1.fa
## with the annotation gtf
Java -jar CIRI3.jar -I sample.sam -O result.txt -F chr1.fa -A chr1.gtf
```
#### 3)Multiple-SAM/BAM files containing RNase R treated information as input

```
## without the annotation gtf
Java -jar CIRI3.jar -I ./data/circRNA/SAM/0/sample.sam -O ./data/circRNA/SAM/0/result.txt -F ./data/circRNA/chr1.fa
## with the annotation gtf
Java -jar CIRI3.jar -I sample.sam -O result.txt -F chr1.fa -A chr1.gtf
```
#### 4)CircRNA collections of interest and sam/bam files as inputs

```
## without the annotation gtf
Java -jar CIRI3.jar -I ./data/circRNA/SAM/0/sample.sam -O ./data/circRNA/SAM/0/result.txt -F ./data/circRNA/chr1.fa
## with the annotation gtf
Java -jar CIRI3.jar -I sample.sam -O result.txt -F chr1.fa -A chr1.gtf
```

### Step3. Differential Expression Analysis

CIRI3 provides multiple input options for the identification and quantification of circRNAs, including single-sample input, multiple-sample input, multiple-sample files containing RNase R treated information. In addition, users have the option to input a collection of circRNAs of interest, allowing CIRI3 to quantify the BSJ and FSJ of these circRNAs within the samples.

A small test data set is provided in [test_data/test_data_espresso_sirv.tar.gz](test_data/test_data_espresso_sirv.tar.gz). The unpacked files are:

#### 1)Single-SAM/BAM files as input
CIRI3
```
## without the annotation gtf
Java -jar CIRI3.jar -I ./data/circRNA/SAM/0/sample.sam -O ./data/circRNA/SAM/0/result.txt -F ./data/circRNA/chr1.fa
## with the annotation gtf
Java -jar CIRI3.jar -I sample.sam -O result.txt -F chr1.fa -A chr1.gtf
```
#### 1)Single-SAM/BAM files as input
CIRI3
```
## without the annotation gtf
Java -jar CIRI3.jar -I ./data/circRNA/SAM/0/sample.sam -O ./data/circRNA/SAM/0/result.txt -F ./data/circRNA/chr1.fa
## with the annotation gtf
Java -jar CIRI3.jar -I sample.sam -O result.txt -F chr1.fa -A chr1.gtf
```
#### 1)Single-SAM/BAM files as input
CIRI3
```
## without the annotation gtf
Java -jar CIRI3.jar -I ./data/circRNA/SAM/0/sample.sam -O ./data/circRNA/SAM/0/result.txt -F ./data/circRNA/chr1.fa
## with the annotation gtf
Java -jar CIRI3.jar -I sample.sam -O result.txt -F chr1.fa -A chr1.gtf
```



### Preparing Input Files

#### Base Calling

The following is one possible procedure for base calling. Other base calling procedures may be preferable for different data sets. If using Nanopore data then the base calling could be done using guppy with a command similar to:
```
guppy_basecaller --input_path /path/to/fast5 --save_path /path/to/fastq/ --cpu_threads_per_caller 20 --config dna_r9.4.1_450bps_fast.cfg
```

That command assumes that the data is from a Nanopore r9.4.1 cDNA library and it uses the fast mode of guppy with 20 threads. For a Nanopore r9.4.1 direct RNA library `--config rna_r9.4.1_70bps_fast.cfg` could be used.


All the separate .fastq files can be combined into a single file with:
```
cat /path/to/fastq/*.fastq > combined.fastq
```

### All Arguments

```
perl ESPRESSO_S.pl --help

Arguments:

    -L, --list_samples
          tsv list of sample(s) (each file in a line with 1st column as sorted
          BAM/SAM file and 2nd column as sample name; required)
    -F, --fa
          FASTA file of all reference sequences. Please make sure this file is
          the same one provided to mapper. (required)
    -A, --anno
          input annotation file in GTF format (optional)
    -B, --SJ_bed
          input custom reliable splice junctions in BED format (optional; each
          reliable SJ in one line, with the 1st column as chromosome, the 2nd
          column as upstream splice site 0-base coordinate, the 3rd column as
          downstream splice site and 6th column as strand)
    -O, --out
          work directory (existing files in this directory may be OVERWRITTEN;
          default: ./)

    -H, --help
          show this help information

    -N, --read_num_cutoff
          min perfect read count for denovo detected candidate splice junctions
          (default: 2)
    -R, --read_ratio_cutoff
          min perfect read ratio for denovo detected candidate splice junctions:
          Set this as 1 for completely GTF-dependent processing (default: 0)

    -C, --cont_del_max
          max continuous deletion allowed; intron will be identified if longer
          (default: 50)
    -M, --chrM
          tell ESPRESSO the ID of mitochondrion in reference file (default:
          chrM)

    -T, --num_thread
          thread number (default: minimum of 5 and sam file number)
    -Q, --mapq_cutoff
          min mapping quality for processing (default: 1)
```

```
perl ESPRESSO_C.pl --help

Arguments:

    -I, --in
          work directory (generated by ESPRESSO_S)
    -F, --fa
          FASTA file of all reference sequences. Please make sure this file is
          the same one provided to mapper. (required)
    -X, --target_ID
          ID of sample to process (required)

    -H, --help
          show this help information

    -T, --num_thread
          thread number (default: 5)
```

```
perl ESPRESSO_Q.pl --help

Arguments:

    -L, --list_samples
          tsv list of multiple samples (each bam in a line with 1st column as
          sorted bam file, 2nd column as sample name in output, 3rd column as
          directory of ESPRESSO_C results; this list can be generated by
          ESPRESSO_S according to the initially provided tsv list; required)
    -A, --anno
          input annotation file in GTF format (optional)
    -O, --out_dir
          output directory (default: directory of -L)
    -V, --tsv_compt
          output tsv for compatible isoform(s) of each read (optional)
    -T --num_thread
          how many threads to use (default: 5)

    -H, --help
          show this help information

    -N, --read_num_cutoff
          min perfect read count for all splice junctions of novel isoform
          (default: 2)
    -R, --read_ratio_cutoff
          min perfect read ratio for all splice junctions of novel isoform
          (default: 0)
```

## Output

The three main output files are:
* `sample_N2_R0_updated.gtf` is an updated GTF annotation for detected isoforms.
  + Each detected isoform is reported as a transcript.
  + The source column indicates for each isoform whether it is a `"novel_isoform"` or an `"annotated_isoform"`.
* `sample_N2_R0_abundance.esp` is a tsv file for expression of detected isoforms.
  + Each detected isoform is reported on a separate line.
  + The first columns are `transcript_ID`, `transcript_name`, `gene_ID`.
  + There is an additional column for each sample name provided in samples.tsv. Those columns show the number of reads from that sample which were counted toward this isoform.
  + The counts are assigned by expectation maximization. Each input read contributes at most 1 count, either to a single isoform or distributed as fractional counts to multiple isoforms.
* `sample_N2_R0_compatible_isoform.tsv` is a tsv file for compatible isoforms of each read.
  + This file is only produced if the -V parameter of ESPRESSO_Q is used.
  + The columns are `read_id`, `sample_name`, `read_classification`, `compatible_isoforms`.
  + The possible classifications are NIC/NNC, NCD, ISM, FSM, single-exon.

## Visualization

Visualizations can be created for the ESPRESSO output files similar to [test_data/visualization_cd44.png](test_data/visualization_cd44.png) and [test_data/visualization_sirv.png](test_data/visualization_sirv.png). The script [visualization/visualize.py](visualization/visualize.py) can generate files that can then be viewed in a genome browser

### Visualization Arguments

```
python3 visualization/visualize.py --help

usage: visualize.py [-h] --genome-fasta GENOME_FASTA --updated-gtf UPDATED_GTF
                    --abundance-esp ABUNDANCE_ESP --target-gene TARGET_GENE
                    --minimum-count MINIMUM_COUNT [--normalize-counts-to-cpm]
                    --descriptive-name DESCRIPTIVE_NAME --output-dir
                    OUTPUT_DIR

Generate files for visualizing ESPRESSO output

optional arguments:
  -h, --help            show this help message and exit
  --genome-fasta GENOME_FASTA
                        the fasta input to ESPRESSO
  --updated-gtf UPDATED_GTF
                        the *_updated.gtf file output by ESPRESSO
  --abundance-esp ABUNDANCE_ESP
                        the *_abundance.esp file output by ESPRESSO
  --target-gene TARGET_GENE
                        the name of the gene to visualize. transcripts with
                        name like {target-gene}-{number} or gene_id like
                        {target-gene}.* will have output generated. Use the
                        gene_id to match novel isoforms output by ESPRESSO
  --minimum-count MINIMUM_COUNT
                        only isoforms where the (normalized) count for a
                        sample meets the minimum count are considered
  --normalize-counts-to-cpm
                        convert raw counts to counts per million
  --descriptive-name DESCRIPTIVE_NAME
                        used as a label in the visualization and for filenames
  --output-dir OUTPUT_DIR
                        where to write visualization files
```

## References

1. Li H. Minimap2: pairwise alignment for nucleotide sequences[J]. Bioinformatics, 2018, 34(18): 3094-3100.
