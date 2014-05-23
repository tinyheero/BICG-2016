% Module 5 Lab - Software Installation Instructions
% Andrew Roth

## Setup
Before we start please read the general guide on installing software in Linux.

First we will setup a directory to work in. For this workshop we have been using the `~/workspace` directory. For your own use you will need to modify this path.

We will create an installation directory for Module5 and enter it. To save some typing we will first create a *shell variable* to store this directory.

```
export INSTALL_DIR=~/workspace/Module5/install

mkdir -p $INSTALL_DIR

cd $INSTALL_DIR
```

We will also create a place to put the *source* files we download.

```
mkdir src
```

Finally we want to update the system variables to look in our personal install directory.

```
export PATH=$INSTALL_DIR/bin:$PATH

export LD_LIBRARY_PATH=$INSTALL_DIR/lib:$LD_LIBRARY_PATH
```

## Installing Python
The first piece of software we will need is Python. This is a programming language like R, that is used by some of the helper scripts in the lab.

The latest version of Python can be found at <http://www.python.org/download>. For this tutorial we will use the 2.7.5 version of Python. We will download the software and unpack it.

> Note: There are two major versions of Python, versions 2 and 3. Version 3 is a major revision of the language which breaks backwards compatibility. As a result many packages still use version 2, which is still being updated.

```
cd $INSTALL_DIR/src

wget http://www.python.org/ftp/python/2.7.5/Python-2.7.5.tar.bz2

tar -jxvf Python-2.7.5.tar.bz2
```

Now we can enter the directory and perform the build process. We want to install the files to our personal install folder so we use the `--prefix` flag for configure.

```
cd Python-2.7.5

./configure --prefix=$INSTALL_DIR

make

make install
```

Finally check which version of Python the system is using.

```
which python
```

it should be

```
/home/ubuntu/workspace/Module5/install/bin/python
```

### Installing Python Packages
For this tutorial we will need the PyVCF package. It can be obtained from <https://pypi.python.org>. We will go through the usual Python package installation procedure. PyVCF requires a package called setuptools which we will install first.

```
cd $INSTALL_DIR/src

wget https://pypi.python.org/packages/source/s/setuptools/setuptools-0.6c11.tar.gz

tar -zxvf setuptools-0.6c11.tar.gz

cd setuptools-0.6c11

python setup.py install
```

Now to get PyVCF.

```
cd $INSTALL_DIR/src

wget https://pypi.python.org/packages/source/P/PyVCF/PyVCF-0.6.3.tar.gz

tar -zxvf PyVCF-0.6.3.tar.gz

cd PyVCF-0.6.3

python setup.py install
```

## Installing R
We will follow a similar procedure for installing R.

Fetch the software from http://www.r-project.org/ and unpack it.

```
cd $INSTALL_DIR/src

wget http://mirror.its.dal.ca/cran/src/base/R-3/R-3.0.1.tar.gz

tar -zxvf R-3.0.1.tar.gz

cd R-3.0.1

./configure --prefix=$INSTALL_DIR --with-x=no

make

make install
```

> Note: We pass the `--with-x=no flag` to `./configure` because R is being built (compiled) on a machine without a graphical interface. If you were building R on your own computer or a server which allowed you to run graphical tools you could omit this flag.

And check which version of R the system is using.

```
which R
```

Which should hopefully print.

```
/home/ubuntu/workspace/Module5/install/bin/R
```

## Installing Array Normalisation Packages
We will use the procedure described at <http://www.openbioinformatics.org/penncnv/penncnv_tutorial_affy_gw6.html>. This page also contains links to most of the software we need to download.

> Note: You will need to register with Affymetrix to download the library files, specifically the .cdf file, from <http://www.affymetrix.com/support/technical/byproduct.affx?product=genomewidesnp_6>. The [aroma-project](http://www.aroma-project.org/docs) has a copy available which we will download since it can be done from the cloud.

You will need to change some variables in the `run_PennCNV.sh` to point to the folders you are creating below. I will point these out as we go.

```
cd $INSTALL_DIR/src

wget http://www.openbioinformatics.org/penncnv/download/penncnv.latest.tar.gz

 # Change penncnvDir variable to point to the directory created.
tar -zxvf penncnv.latest.tar.gz

wget http://www.openbioinformatics.org/penncnv/download/gw6.tar.gz

 # Change the gw6Dir variable to point to folder created 
tar -zxvf gw6.tar.gz

wget http://media.affymetrix.com/Download/updates/apt-1.15.1-src.zip

 # Change the aptDir variable to point to the folder created
unzip apt-1.15.1-src.zip

wget http://www.aroma-project.org/data/annotationData/chipTypes/GenomeWideSNP_6/GenomeWideSNP_6,Full.CDF.gz

 # Change the cdfFile variable to point to the file created
gunzip GenomeWideSNP_6,Full.CDF.gz
```

## Installing OncoSNP
For predicting CNVs from the array data we will use OncoSNP in this tutorial. The files for OncoSNP can be downloaded from <https://sites.google.com/site/oncosnp/user-guide/downloads>. Before you can use OncoSNP you will need to register with the author and he will supply a password which you can use to unlock the downloaded files.

We will also need to download GC content files for the relevant build of the genome. For this tutorial we will use the hg18 (build 36) files.

```
cd $INSTALL_DIR/src

wget ftp://ftp.stats.ox.ac.uk/pub/yau/oncosnp/executables/oncosnp_v2.28.run

chmod +x oncosnp_v2.28.run

 # This will require that you have registered and received a password.
./oncosnp_v2.28.run

wget ftp://ftp.stats.ox.ac.uk/pub/yau/oncosnp/mcr/MCRinstaller.run.zip

unzip MCRinstaller.run.zip

chmod +x MCRinstaller.run

 # You will need the same password from the previous step here.
./MCRinstaller.run

wget ftp://ftp.stats.ox.ac.uk/pub/yau/quantisnp2/download/b36.tar.gz

tar -zxvf b36.tar.gz
```

The `MCRinstaller.run` program will ask you where to install the MATLAB install files. You should edit the `mcrDir` in the `run_oncosnp.sh` script to point here. You will also need to edit the `installDir` variable to point to where OncoSNP installs and the `gcDir` variable to point to the folder created by decompressing `b36.tar.gz`. 

## Installing HMMcopy
For this tutorial we will install HMMcopy. It is a Bioconductor package http://www.bioconductor.org/packages/release/bioc/html/HMMcopy.html.

We will start R.

```
R
```

and inside the R environment execute

```
source("http://bioconductor.org/biocLite.R")
biocLite("HMMcopy")
```

where the last line was copied from the HMMcopy Bioconductor page. Follow along and accept any suggested updates.

HMMcopy also has C package associated with it. We will need to download this from <http://compbio.bccrc.ca/software/hmmcopy>.

```
cd $INSTALL_DIR/src

wget http://compbio.bccrc.ca/files/2013/04/HMMcopy_0.1.1_14Dec2012.tar_.gz

mv HMMcopy_0.1.1_14Dec2012.tar_.gz HMMcopy_0.1.1_14Dec2012.tar.gz

tar -zxvf HMMcopy_0.1.1_14Dec2012.tar.gz

cd HMMcopy_0.1.1_14Dec2012

make
```

## Installing APOLLOH

We will download APOLLOH from <http://compbio.bccrc.ca/software/apolloh>. APOLLOH is built using MATLAB which requires us to install some additional files.

In the first step we will download the MATLAB library files and install them. The last line below will start the interactive installation process. You will be required to specify a path where the files should be put. Choose `/home/ubuntu/workspace/Module5/install/opt/matlab/apolloh`.

```
cd $INSTALL_DIR/src

wget ftp://ftp.bcgsc.ca/public/shahlab/Apolloh/MCRInstaller.bin

chmod +x MCRInstaller.bin

./MCRInstaller.bin
```

When you run APOLLOH you need to edit the line

```
mcrDir=/opt/MATLAB/MATLAB_Component_Runtime/v77
``` 

in the `run_APOLLOH.sh` to be

```
mcrDir=/home/ubuntu/workspace/Module5/install/opt/matlab/apolloh/v77
```


Now we will install APOLLOH.

```
wget http://compbio.bccrc.ca/files/2011/12/APOLLOH_0.1.1_01Oct2012.tar_.gz

mv APOLLOH_0.1.1_01Oct2012.tar_.gz APOLLOH_0.1.1_01Oct2012.tar.gz

tar -zxvf APOLLOH_0.1.1_01Oct2012.tar.gz
```

If you wanted to use this version of APOLLOH you would need to modify the `run_APOLLOH.sh` script so that

```
installDir=/usr/local/APOLLOH_0.1.1/
```

becomes

```
installDir=/home/ubuntu/workspace/Module5/install/APOLLOH_0.1.1/
```

## Installing GATK
Installing [GATK](http://www.broadinstitute.org/gatk) is straightforward, you need to register and download the package. Once you have downloaded and decompressed the file you need to point any calls to GATK here.

For example to the `runGATK.sh` script used for this lab has the line

```
java -jar ~/Documents/teaching/cbw/2013/install/src/GATK/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper -baq RECALCULATE -L 21 -I $normalBam -I $tumourBam -o $outFile
```

This would need to be changed to 

```
 java -jar /path/to/GATK/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper -baq RECALCULATE -L 21 -I $normalBam -I $tumourBam -o $outFile
```

assuming `/path/to/GATK` is the folder you decompressed GATK to.