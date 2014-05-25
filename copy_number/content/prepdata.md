# Module 5 Lab - Copy Number Analysis
> # Data Preparation

## Introduction

This lab uses a number of publically available datasets.  Below we give details about how to obtain these datasets.  We also include instructions for the time intensive processing has been applied to these datasets prior to the lab.  The instructions below will not be needed during the lab, they are included for reference purposes in case you wish to reproduce the lab in your own time.

## SNP6.0 HCC1143 cell line data

The cel files containing the probe intensity data for cell line HCC1143 can be downloaded as follows.

    mkdir -p data/cel
    cd data/cel
    wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM337nnn/GSM337662/suppl/GSM337662%2ECEL%2Egz
    wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM337nnn/GSM337641/suppl/GSM337641%2ECEL%2Egz
    gunzip *.gz
    cd ../..

## Obtaining whole genome sequencing HCC1143 cell line data

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

The next few steps will manipulate the HCC1143 bam files.

    cd HCC1143

### Extract chromosome 21 

We will focus on a single chromosome in this lab as the whole genome analysis is too time consuming.  Extract chromosome 21 using samtools and index the bam files.

    samtools view -b G15511.HCC1143.1.bam 20 > G15511.HCC1143.1.chr21.bam
    samtools view -b G15511.HCC1143_BL.1.bam 20 > G15511.HCC1143_BL.1.chr21.bam
    
    samtools index G15511.HCC1143.1.chr21.bam
    samtools index G15511.HCC1143_BL.1.chr21.bam

To get a preliminary look at the data, create a `.tdf` file using `igvtools` and view the bam file in igv.

    igvtools count G15511.HCC1143.1.chr21.bam G15511.HCC1143.1.chr21.bam.tdf hg19
    igv G15511.HCC1143.1.chr21.bam

The `.tdf` calculates read depth that igv displays as a histogram.

## Calculate read depth for HMMCopy

Copy number prediction is based on read depth.  Calculate read depth for the tumour and normal bam files and store in `.wig` format.

    $HMMCOPY_DIR/readCounter G15511.HCC1143.1.chr21.bam > G15511.HCC1143.1.chr21.wig
    $HMMCOPY_DIR/readCounter G15511.HCC1143_BL.1.chr21.bam > G15511.HCC1143_BL.1.chr21.wig

## Convert to pileup





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