# Module 5 Lab - Copy Number Analysis
> # Data Preparation

## Introduction

This lab uses a number of publically available datasets.  Below we give details about how to obtain these datasets.  We also include instructions for the time intensive processing has been applied to these datasets prior to the lab.  The instructions below will not be needed during the lab, they are included for reference purposes in case you wish to reproduce the lab in your own time.

## Environment

The following environment variables are required.

`PICARD_DIR=/usr/local/picard` picard tools binaries directory
`HMMCOPY_DIR=/usr/local/HMMcopy` hmm copy binaries directory

## SNP6.0 HCC1143 cell line data

The cel files containing the probe intensity data for cell line HCC1143 can be downloaded as follows.

    mkdir -p data/cel
    cd data/cel
    wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM337nnn/GSM337662/suppl/GSM337662%2ECEL%2Egz
    wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM337nnn/GSM337641/suppl/GSM337641%2ECEL%2Egz
    gunzip *.gz
    cd ../..

## Whole genome sequencing HCC1143 cell line data

The WGS HCC1143 cell line data can be downloaded from the [TCGA benchmarking website](https://cghub.ucsc.edu/datasets/benchmark_download.html).  The downloading instructions on the website are reproduced here.

You require the [GeneTorrent](https://cghub.ucsc.edu/software/downloads.html) software to download the data.  If you are running Ubuntu 13.10 download the binaries as follows.  Note that GeneTorrent depends on [boost](http://www.boost.org/).

    wget https://cghub.ucsc.edu/software/downloads/GeneTorrent/3.8.5a/GeneTorrent-download-3.8.5a-94-Ubuntu13.10.x86_64.tar.gz
    tar -xzvf GeneTorrent-download-3.8.5a-94-Ubuntu13.10.x86_64.tar.gz

The tumour and normal HCC1143 data is the first on the list on the TCGA benchmarking website.  Download the xml uuid files as follows.

    wget --no-check-certificate https://cghub.ucsc.edu/cghub/metadata/analysisAttributes?analysis_id=ad3d4757-f358-40a3-9d92-742463a95e88 -O ad3d4757-f358-40a3-9d92-742463a95e88.xml
    wget --no-check-certificate https://cghub.ucsc.edu/cghub/metadata/analysisAttributes?analysis_id=f0eaa94b-f622-49b9-8eac-e4eac6762598 -O f0eaa94b-f622-49b9-8eac-e4eac6762598.xml

Download the bam files using GeneTorrent.

    gtdownload -c https://cghub.ucsc.edu/software/downloads/cghub_public.key -vv -d ad3d4757-f358-40a3-9d92-742463a95e88.xml
    gtdownload -c https://cghub.ucsc.edu/software/downloads/cghub_public.key -vv -d f0eaa94b-f622-49b9-8eac-e4eac6762598.xml

The files will appear in two subdirectories named after the uuid for each dataset.  For convencience, move the files to a subdirectory named HCC1143.

    mkdir HCC1143
    mv ad3d4757-f358-40a3-9d92-742463a95e88/* HCC1143
    mv f0eaa94b-f622-49b9-8eac-e4eac6762598/* HCC1143

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




 momac23:copy_number amcphers$ ls
GenomeWideSNP_6,Full.CDF    apt-1.15.1-src.zip      data                old             output              packages
momac23:copy_number amcphers$ wget http://bioinfo-out.curie.fr/projects/freec/data/testChr19.zip
--2014-05-23 17:56:19--  http://bioinfo-out.curie.fr/projects/freec/data/testChr19.zip
Resolving bioinfo-out.curie.fr (bioinfo-out.curie.fr)... 193.49.205.28
Connecting to bioinfo-out.curie.fr (bioinfo-out.curie.fr)|193.49.205.28|:80... connected.
HTTP request sent, awaiting response... 200 OK
Length: 1334106500 (1.2G) [application/zip]
Saving to: ‘testChr19.zip’

100%[==================================================================================================================================================================================================================================================================================>] 1,334,106,500  331KB/s   in 37m 12s

2014-05-23 18:33:32 (584 KB/s) - ‘testChr19.zip’ saved [1334106500/1334106500]

momac23:copy_number amcphers$ unzip testChr19.zip
Archive:  testChr19.zip
   creating: testChr19/
  inflating: testChr19/GC_profile.cnp
  inflating: testChr19/run.sh
  inflating: testChr19/hg19_snp131.SingleDiNucl.1based.txt
  inflating: testChr19/config_chr19.txt
  inflating: testChr19/chr_19.noDup0.pileup.gz
  inflating: testChr19/hs19_chr19.len
  inflating: testChr19/README.txt
  inflating: testChr19/makeGraph_Chromosome.R
  inflating: testChr19/chr_19.noDup0.pileup.gz_sample.cpn
momac23:copy_number amcphers$
momac23:copy_number amcphers$ ls
GenomeWideSNP_6,Full.CDF    apt-1.15.1-src.zip      data                old             output              packages            testChr19           testChr19.zip
momac23:copy_number amcphers$ pwd
/Users/amcphers/Analysis/cbw_tutorial/copy_number
momac