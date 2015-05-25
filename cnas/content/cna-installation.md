# Lab Module 5 - Software Installation Instructions

This page describes where to get and how to install the software used for this lab module. 

## Installing Array Normalization Packages

We will use the procedure described on the [penncnv site](http://www.openbioinformatics.org/penncnv/penncnv_tutorial_affy_gw6.html). This page also contains links to most of the software we need to download.

> Note: You will need to register with Affymetrix to download the library files, specifically the .cdf file, from the [affimetrix site](http://www.affymetrix.com/support/technical/byproduct.affx?product=genomewidesnp_6). The [aroma-project](http://www.aroma-project.org/docs) has a copy available which we will download since it can be done from the cloud.

```
cd $INSTALL_DIR
wget http://www.openbioinformatics.org/penncnv/download/penncnv.latest.tar.gz
tar -xzvf penncnv.latest.tar.gz

wget http://www.openbioinformatics.org/penncnv/download/gw6.tar.gz
tar -xzvf gw6.tar.gz
    
wget http://media.affymetrix.com/Download/updates/apt-1.17.0-x86_64-intel-linux.zip
unzip apt-1.17.0-x86_64-intel-linux.zip
    
wget http://www.aroma-project.org/data/annotationData/chipTypes/GenomeWideSNP_6/GenomeWideSNP_6.cdf.gz
gunzip GenomeWideSNP_6.cdf.gz
```

## Installing OncoSNP

For predicting CNVs from the array data we will use OncoSNP in this lab. The files for OncoSNP can be downloaded from <https://sites.google.com/site/oncosnp/user-guide/downloads>. Before you can use OncoSNP, you will need to register with the author and he will supply a password which you can use to unlock the downloaded files.

Here we will download OncoSNP version 1.4:

```
wget http://www.well.ox.ac.uk/~cyau/oncosnp/oncosnp_v1.4.run
bash oncosnp_v1.4.run          # enter your oncosnp password
```

Then follow the on-screen instructions. 

You'll also need to get the matching MATLAB MCR file to be used with OncoSNP. This will be described on the website. 

```
wget ftp://ftp.stats.ox.ac.uk/pub/yau/oncosnp/mcr/MCRinstaller.run.zip
unzip MCRinstaller.run.zip
bash MCRinstaller.run          # enter your oncosnp password, and select installation directory
```

The `MCRinstaller.run` installer will ask you where to install the MATLAB MCR files. You should `export MCR_DIR=` to create an environment variable pointing to the location to which you installed MATLAB MCR.  The MCR will be located in a subdirectory `MATLAB_Compiler_Runtime/vXXX` for version XXX of the compiler.

We will also need to download GC content files for the relevant build of the genome. For this tutorial we will use the hg19 (build 37) files.

```
wget ftp://ftp.stats.ox.ac.uk/pub/yau/quantisnp2/download/b36.tar.gz
tar -zxvf b36.tar.gz
export GC_DIR=$INSTALL_DIR/b36
```

## Installing HMMCopy

For this tutorial we will install HMMCopy. It is a Bioconductor package http://www.bioconductor.org/packages/release/bioc/html/HMMcopy.html.

We will start R.

    R

and inside the R environment execute

    source("http://bioconductor.org/biocLite.R")
    biocLite("HMMcopy")

where the last line was copied from the HMMcopy Bioconductor page. Follow along and accept any suggested updates.

HMMcopy also has C package associated with it. We will need to download this from <http://compbio.bccrc.ca/software/hmmcopy>.

    cd $INSTALL_DIR/src
    wget http://compbio.bccrc.ca/files/2013/04/HMMcopy_0.1.1_14Dec2012.tar_.gz
    mv HMMcopy_0.1.1_14Dec2012.tar_.gz HMMcopy_0.1.1_14Dec2012.tar.gz
    tar -zxvf HMMcopy_0.1.1_14Dec2012.tar.gz
    cd HMMcopy_0.1.1_14Dec2012
    make

## Installing TITAN

TITAN is available as a Bioconductor R package. We can install it in R. First open R:

```
R
```

In R, enter the following commands:

```
...
```

