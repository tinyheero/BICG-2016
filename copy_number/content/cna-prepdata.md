# Module 5 Copy Number Alterations - Data Preparation

## Introduction

This lab uses a number of publicly available datasets.  Below we give details about how to obtain these datasets.  We also include instructions for the time intensive processing has been applied to these datasets prior to the lab.  The instructions below will not be needed during the lab, they are included for reference purposes in case you wish to reproduce the lab in your own time.

## Environment

The following environment variables are required.

`PICARD_DIR=/usr/local/picard` picard tools binaries directory
`HMMCOPY_DIR=/usr/local/HMMcopy` hmm copy binaries directory

Bcftools install:

```
SAMTOOLS_DIR=/usr/local/opt/samtools
```

SnpEff install:

```
SNPEFF_DIR=/usr/local/snpEff
```

Summary of all the above environment commands (for copy and pasting convenience):

```
SAMTOOLS_DIR=/usr/local/opt/samtools
SNPEFF_DIR=/usr/local/snpEff
```

## Geting/Preparing Data for Sequencing Analysis

## Affymetrix SNP 6.0 HCC1395 cell line data

The cel files containing the probe intensity data for cell line HC1395 can be downloaded as follows.

    mkdir -p data/cel
    cd data/cel
    wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM337nnn/GSM337662/suppl/GSM337662%2ECEL%2Egz
    wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM337nnn/GSM337641/suppl/GSM337641%2ECEL%2Egz
    gunzip *.gz
    cd ../..

## Whole genome sequencing reference genome

Whole genome reference data will be stored in the `genome` directory.

    mkdir genome
    cd genome

The reference genome fasta file *must* match the reference genome which was used to create the HCC1143 bam files.  Conventiently the url is stored in the bam header.

    samtools view -H ../HCC1143/G15511.HCC1143.1.chr21.bam

Download the genome fasta.  Also create a fasta index using samtools and 

    wget http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta

Create a samtools fasta index for constant time lookup into the sequence.

    samtools faidx Homo_sapiens_assembly19.fasta

Create a dict file for use with the GATK.

    java -jar $PICARD_DIR/CreateSequenceDictionary.jar R=Homo_sapiens_assembly19.fasta O=Homo_sapiens_assembly19.dict

## HMMCopy genome setup

Two major systematic biases in whole genome sequence data are due to sequence mappability and GC content.  For HMMCopy we need to precalculate the GC content and mappability accross each chromosome.  

To calculate GC for chromosome 21, run HMMCopy's `gcCounter` tool.

    $HMMCOPY_DIR/bin/gcCounter -c 21 Homo_sapiens_assembly19.fasta > hg19.21.gc.wig

Mappability information can be downloaded from ucsc in the form of a mappability bigWig file.  Run HMMCopy's `mapCounter` to convert to HMMCopy format.

    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig
    $HMMCOPY_DIR/bin/mapCounter -c chr21 wgEncodeCrgMapabilityAlign100mer.bigWig > hg19.21.map.wig

Finished working in the genome directory.

    cd ..

## Extract chromosome 21 alignments

The next few steps will manipulate the HCC1143 bam files.

    cd HCC1143

We will focus on a single chromosome in this lab as the whole genome analysis is too time consuming.  Extract chromosome 21 using samtools.

    samtools view -b G15511.HCC1143.1.bam 21 > G15511.HCC1143.1.chr21.bam
    samtools view -b G15511.HCC1143_BL.1.bam 21 > G15511.HCC1143_BL.1.chr21.bam

Use samtools to index the bam files.

    samtools index G15511.HCC1143.1.chr21.bam
    samtools index G15511.HCC1143_BL.1.chr21.bam

To get a preliminary look at the data, create a `.tdf` file using `igvtools` and view the bam file in igv.

    igvtools count G15511.HCC1143.1.chr21.bam G15511.HCC1143.1.chr21.bam.tdf hg19
    igv G15511.HCC1143.1.chr21.bam

The `.tdf` calculates read depth that igv displays as a histogram.

## Calculate read depth for HMMCopy

Copy number prediction is based on read depth.  Calculate read depth for the tumour and normal bam files and store in `.wig` format.  Specify chromosome 21.

    $HMMCOPY_DIR/bin/readCounter -c 21 G15511.HCC1143.1.chr21.bam > G15511.HCC1143.1.chr21.wig
    $HMMCOPY_DIR/bin/readCounter -c 21 G15511.HCC1143_BL.1.chr21.bam > G15511.HCC1143_BL.1.chr21.wig

## Convert to pileup

For OncoSNP-SEQ, we require pileup files.  These can take some time to prepare, so they have been created for you.

    samtools mpileup G15511.HCC1143_BL.1.chr21.bam > G15511.HCC1143_BL.1.chr21.pileup
    samtools mpileup G15511.HCC1143.1.chr21.bam > G15511.HCC1143.1.chr21.pileup

## OncoSNP-Seq setup

    mkdir -p genome/oncoseq
    cd genome/oncoseq

Download dbSNP for chromosome 21 and extract.

    wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b141_GRCh37p13/BED/bed_chr_21.bed.gz
    gunzip bed_chr_21.bed.gz

Download the GC content files from the [OncoSNP-SEQ google sites page](https://sites.google.com/site/oncosnpseq/downloads).  If the following wget does not work, download using your browser.

    wget --no-check-certificate "https://doc-08-7c-docs.googleusercontent.com/docs/securesc/ha0ro937gcuc7l7deffksulhg5h7mbp1/b3vbuuo4e743kg87p1qehbamqsft7lpo/1401026400000/03938588792210094527/*/0B_XFGmx3odi4MzlPTjhqTlktS3c?h=16653014193614665626&e=download" -O b37.tar.gz
    tar -xzvf b37.tar.gz

# Preparing Data of Sequencing 

As mentioned in the lab, we need the following ...

## Calculate tumour and normal read depth using HMMCopy

Copy number prediction is based on read depth.  Calculate read depth for the tumour and normal bam files and store in `.wig` format.  We use the -c parameter to specify all autosome chromosomes:

```
mkdir -p hmmCopy/wig 
/usr/bin/readCounter -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 HCC1395/exome/HCC1395_exome_tumour.bam \ 
	> hmmCopy/wig/HC1395_exome_tumour.wig
/usr/bin/readCounter -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 HCC1395/exome/HCC1395_exome_normal.bam \
	> hmmCopy/wig/HC1395_exome_normal.wig 
```

## Retrieve allele counts of heterozygous positions from the tumour and normal

The next step we need to do is retrieve the allele count data from heterozygous positions in both the tumour and normal. The first we will do is get a list of all heterozygous positions in the normal:

```
mkdir -p titan/bcftools/vcf 
$SAMTOOLS_DIR/samtools mpileup -u -I -f ref_data/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa HCC1395/exome/HCC1395_exome_normal.bam | \
	$SAMTOOLS_DIR/bcftools/bcftools view -vcg - | \
	java -jar $SNPEFF_DIR/SnpSift.jar filter "isHet(GEN[0]) & (QUAL >= 20)" > \
	titan/bcftools/vcf/HCC1395_exome_normal.var.het.vcf
```

Here we used samtools and bcftools (1.18) to retrieve these positions. We impose a `QUAL>=20` filter that only retrieves positions with a quality of > 20 to enhance quality. You can use other programs (e.g. GATK) to retrieve these heterozygous positions. 

Now that we have the 