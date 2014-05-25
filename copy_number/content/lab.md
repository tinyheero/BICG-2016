# Module 5 Lab - Copy Number Analysis

## Setup

First login into the cloud.

Now enter the ~/workspace directory

    cd ~/workspace

Now we will download the scripts we need for this module from the wiki page.  I have copied the link from the wiki page below. If the command below does not work try copying the link from the wiki and pasting yourself.

wget http://bioinformatics.ca/workshop_wiki/images/0/0c/Module5.tar.gz

Now we will "unzip" the files.

    tar -zxvf Module5.tar.gz

This should create a folder called Module5. Check this is true.

    ls -lh

We can now remove the compressed ".tar.gz" file to keep our workspace clean.

    rm Module5.tar.gz

Enter the Module5 directory

    cd Module5/


Link to the scripts




## SNP6.0 Analysis

### Fetch Array Data

> *Note* the array data has already been downloaded for you.  The following instructions are provided for reference purposes only

To link to the existing data:

    mkdir data
    ln -s /home/ubuntu/CourseData/CG_data/Module5/data/cel data/cel

To download the cel data:

    mkdir -p data/cel
    cd data/cel
    wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM337nnn/GSM337662/suppl/GSM337662%2ECEL%2Egz
    wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM337nnn/GSM337641/suppl/GSM337641%2ECEL%2Egz
    gunzip *.gz
    cd ../..

Create a list of the cel files to be used by downstream tools.  In practice we would normalize many arrays in a batch.  For demonstration purposes we use just a single tumour and normal.

    echo cel_files > file_list1
    echo `pwd`/GSM337641.CEL >> file_list1
    echo `pwd`/GSM337662.CEL >> file_list1

### Array Normalisation and LRR/BAF Extraction

The first step in array analysis is to normalise the data and extract the log R and BAF (B Allele Frequencies). The following steps will create normalised data from Affymetrix SNP6 data.

We require a number of data files that define the SNP6.0 arrays.

The sketch file gives a subset of the probes that work well for normalization.

    SKETCH_FILE=$GW6_DIR/lib/hapmap.quant-norm.normalization-target.txt

The cluster file defines genotype clusters from hapmap, and is used for small batches.

    CLUSTER_FILE=$GW6_DIR/lib/hapmap.genocluster

Chromosome positions for each probe.

    LOC_FILE=$GW6_DIR/lib/affygw6.hg18.pfb

#### Step 1: Probe set summarization

    $APT_DIR/bin/apt-probeset-summarize --cdf-file $SNP6_CDF --analysis quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true --target-sketch $SKETCH_FILE --out-dir results/apt --cel-files data/cel/file_list.txt --chip-type GenomeWideEx_6 --chip-type GenomeWideSNP_6

#### Step 2: B-allele and log ratios

    $GW6_DIR/bin/normalize_affy_geno_cluster.pl $CLUSTER_FILE results/apt/quant-norm.pm-only.med-polish.expr.summary.txt -locfile $LOC_FILE -out output/gw6.lrr_baf.txt

#### Step 3: Split results into a single file per sample

    $PENN_CNV_DIR/kcolumn.pl results/gw6.lrr_baf.txt split 2 -tab -head 3 -name --output results/gw6

The normalised files will be placed in `results/gw6*`.

    ls -lh results/gw6*

The file structure is one probe per line, giving the position, normalized log R and BAF for each probe.

    less -S results/gw6.GSM337641

### CNV Calling and Visualisation 

Now that we have the BAF and LRR data we will use OncoSNP to analyse this data.  Create a working directory for OncoSNP.

    mkdir oncosnp

OncoSNP requires an input file describing the dataset to be analyzed.  A template is available with the oncoSNP install.

    cp /usr/local/oncosnp/demo/example-batch-file.txt oncosnp/batch-file.txt

Next we edit this file to point to our data.  Note the tumour data is in gw6.GSM337641 and the normal data in gw6.GSM337662.  Set the name of the sample to `HCC113`.

    nano -w oncosnp/batch-file.txt

OncoSNP has many command line parameters, and most will not change between runs of different datasets.  For convencience we have provided the run_oncosnp.sh.  This will be easier to execute than a large command line.

    less -S scripts/run_oncosnp.sh

Create an output directory for oncosnp.  Run OncoSNP using the script provided.  

We will run OncoSNP in the background since it is slow.  While it goes we can take a few more minutes to study the script for launching it. Please ask questions.

    mkdir oncosnp/results
    bash scripts/run_oncosnp.sh oncosnp/batch-file.txt oncosnp/results &

Rather then print to screen, the script sends the output of OncoSNP to a log file. We can monitor the progress of the program by examining this file.

    less -S oncosnp/results/run.log

Similarly the errors are also sent to a file which we can explore.

    less -S oncosnp/results/run.err

We can see if the script is still running by looking at our background jobs

    jobs

When the program finishes we can go to the output folder and browser the results.

    ls -lh oncosnp/results

The first key file is the .qc file which outputs some basic quality control values and some parameters. Probably the most interesting value is the stromal contamination i.e. fraction of normal cells. Two values are reported by default because OncoSNP does multiple analysis runs. The first value is the most probable.

    less -S oncosnp/results/HCC1143.qc

Next is .cnvs file which contains the smoothed segments with there copy number prediction.

    less -S oncosnp/results/HCC1143.cnvs

The last file we will look at is the .cnv file. This is essentially a more informative version of the .cnvs file. One column of particular interest is the "Tumour State" column. This is an integer >= 1 which represents the most likely state of the HMM for that segment. 

    less -S oncosnp/results/HCC1143.cnv

One downside of OncoSNP is that it does not give parental copy number directly. We can recover this information by taking the tumour state for a segment and looking it up in the file passed to OncoSNP after the --tumourstatesfile flag.  The states correspond to rows in this file. The minor copy number is given by the TumourBalleles2 column, the major copy number is given by the tumourBalleles3 and the total copy number by TumourBalleles4.

    less -S /usr/local/oncosnp/configuration/tumourStates.dat

We've included a script to parse this information out.

    python scripts/parse_segments.py oncosnp/results/HCC1143.cnv oncosnp/results/HCC1143.pcn.seg /usr/local/oncosnp/configuration/tumourStates.dat
    less -S oncosnp/results/HCC1143.pcn.seg

The final intersting file that OncoSNP produces is the plots HCC1143.ps.gz.  Since I would like to avoid dealing with transfering from the cloud in this tutorial I have included a copy of the figure in scripts for this lab.  Download this file to your own computer decompress it. The plot is under figures/oncosnp.

## Whole Genome Sequencing (WGS) Analysis

The workflow for WGSS data is not dramatically different. We still need to do some normalisation and B allele extraction.

The dataset we will be using is a breast cancer cell line sequenced by the TCGA.  The data has been downloaded and partially processed for you (See data preparation).

    ln -s /home/ubuntu/CourseData/CG_data/TCGA/HCC1143/

It is very important that we use the same genome fasta as was used to generate the bams.  The genome fasta has been downloaded for you.

    ln -s /home/ubuntu/CourseData/CG_data/Module5/genome/

### HMMCopy for total copy number

HMMcopy requires several input files all in .wig format. For most analyses you can use the GC content and mappability files on the HMMcopy website.  For this lab we have created versions for chr21.

    less -S /home/ubuntu/CourseData/CG_data/Module5/genome/hg19.21.gc.wig
    less -S /home/ubuntu/CourseData/CG_data/Module5/genome/hg19.21.map.wig

HMMCopy also requires read depth for tumour and normal data.  These are obtained by counting reads overlapping 1000bp segments in the human genome. 

    less -S /home/ubuntu/CourseData/CG_data/Module5/HCC1143/G15511.HCC1143.1.chr21.wig
    less -S /home/ubuntu/CourseData/CG_data/Module5/HCC1143/G15511.HCC1143_BL.1.chr21.wig

An R script is provided to run HMMCopy.  

    less -S scripts/run.R

Start R,

    R

and use the source function to execute the script.  Then quit.

    source('scripts/hmmcopy/run.R')
    quit()

The `run.R` script will create two directories: `results/hmmcopy/results` containts the segment files and `results/hmmcopy/plots` contains the various plots generated.  The `.igv.seg` segment file can be loaded into IGV.

    less -S results/hmmcopy/results/tumour.igv.seg

<!-- The plots file are in the scripts you downloaded earlier under figures/hmmcopy. -->

### APOLLOH for allele specific copy number

To obtain BAF data from bam files, we first need to run a genotyper to call heterozygous SNP positions in the normal.  GATK has been run for you (see data preparation).  The output of GATK is a VCF file.

    less -S data/vcf/HCC1143.GATK.vcf

Given heterozygous SNPs, we now need to extract the reference and alternate allele counts from the bam file.  The script `build_apolloh_allelic_counts_file.py` does this for you.

    mkdir -p data/apolloh
    python scripts/build_apolloh_allelic_counts_file.py data/vcf/HCC1143.GATK.vcf data/apolloh/baf.txt --normal_column 1

The BAF file is a tab delimited file with one row per SNP.

    less -S data/apolloh/baf.txt

We are now in a position to run APOLLOH.

    mkdir -p results/apolloh
    bash scripts/run_apolloh.sh data/apolloh/baf.txt results/hmmcopy/results/tumour.apolloh.seg results/apolloh/

There are three files. params.txt contains the APOLLOH model paramters. Most interesting Normal cell contamination

The `params.txt` file contains the APOLLOH model parameters.

    less -S results/apolloh/params.txt

The `loh.txt` file contains information about the state of each SNP.

    less -S loh.txt

The `segs.txt` contains information about the segments.

    less -S segs.txt


### OncoSNP-SEQ 



    perl /usr/local/oncosnpseq/scripts/process_pileup.pl --infile HCC1143/G15511.HCC1143.1.chr21.pileup  --outfile test --snpfile genome/oncoseq/bed_chr_21.bed
    perl /usr/local/oncosnpseq/scripts/process_pileup.pl --infile HCC1143/G15511.HCC1143.1.chr21.pileup  --outfile test --snpfile genome/oncoseq/bed_chr_21.bed




## Visualizing Datasets In IGV

For this part please download the METABRIC dataset from the wiki and open it in IGV.

