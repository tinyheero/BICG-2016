# Lab Module 5 - Copy Number Analysis

## Setup

First login into the cloud.

Now enter the ~/workspace directory

```
cd ~/workspace
```

Download the module package we need for this module from the wiki page onto your amazon instance.

```
wget http://bioinformatics.ca/workshop_wiki/images/0/0c/Module5.tar.gz
```

Now we will "unzip" the files.

```
tar -zxvf Module5.tar.gz
```

This should create a folder called Module5. Check this is true.

```
ls -lh
```

We can now remove the compressed ".tar.gz" file to keep our workspace clean.

```
rm Module5.tar.gz
```

Enter the Module5 directory

```
cd Module5/
```

Link to the scripts

## Environment

```
INSTALL_DIR=/home/ubuntu/CourseData/CG_data/Module5/install/
```

Normalization files:

```
GW6_DIR=$INSTALL_DIR/gw6
```

Affimetrix power tools:

```
APT_DIR=$INSTALL_DIR/apt-1.17.0-x86_64-intel-linux
```

Cell definition file for SNP6.0:

```
SNP6_CDF=$INSTALL_DIR/GenomeWideSNP_6.cdf
```

Oncosnp and oncosnp-seq locations:

```
export ONCOSNP_DIR=/usr/local/oncosnp
export ONCOSNPSEQ_DIR=/usr/local/oncosnpseq/
export MCR_DIR=/home/ubuntu/CourseData/software/MATLAB/MCR/v82
```

GC content files for oncosnp and oncosnpseq:

```
export GC_DIR=/home/ubuntu/CourseData/CG_data/Module5/install/b37
```

Summary of all the above environment commands (for copy and pasting convenience):

```
INSTALL_DIR=/home/ubuntu/CourseData/CG_data/Module5/install/
GW6_DIR=$INSTALL_DIR/gw6
APT_DIR=$INSTALL_DIR/apt-1.17.0-x86_64-intel-linux
SNP6_CDF=$INSTALL_DIR/GenomeWideSNP_6.cdf
export ONCOSNP_DIR=/usr/local/oncosnp
export ONCOSNPSEQ_DIR=/usr/local/oncosnpseq/
export MCR_DIR=/home/ubuntu/CourseData/software/MATLAB/MCR/v82
export GC_DIR=/home/ubuntu/CourseData/CG_data/Module5/install/b37
```

## Affymetrix SNP 6.0 Analysis

### Fetch Array Data

For calling copy number variants from Affymetrix SNP 6.0 data, we will be using breast cancer cell-line (HC1395). The array data for HCC1395 has already been downloaded for you. Create a directory and copy the array data to your workspace.

```
mkdir -p data/cel
cp /home/ubuntu/CourseData/CG_data/HCC1395/cel/* data/cel
```

Create a list of the cel files to be used by downstream tools.  In practice we would normalize many arrays in a batch.  For demonstration purposes we use just a single tumour.

```
echo cel_files > data/cel/file_list.txt
echo `pwd`/data/cel/GSM888107.CEL >> data/cel/file_list.txt
```

### Step 1 - Array Normalisation and LRR/BAF Extraction

The first step in array analysis is to normalize the data and extract the log R and BAF (B Allele Frequencies). The following steps will create normalised data from Affymetrix SNP6 data.

We require a number of data files that define the SNP6.0 arrays.

The sketch file gives a subset of the probes that work well for normalization.

```
SKETCH_FILE=$GW6_DIR/lib/hapmap.quant-norm.normalization-target.txt
```

The cluster file defines genotype clusters from hapmap, and is used for small batches.

```
CLUSTER_FILE=$GW6_DIR/lib/hapmap.genocluster
```

Chromosome positions for each probe.

```
LOC_FILE=$GW6_DIR/lib/affygw6.hg19.pfb
```

#### Step 1: Probe set summarization

```
$APT_DIR/bin/apt-probeset-summarize --cdf-file $SNP6_CDF \
    --analysis quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true \
    --target-sketch $SKETCH_FILE --out-dir results/apt \
    --cel-files data/cel/file_list.txt --chip-type GenomeWideEx_6 \
    --chip-type GenomeWideSNP_6
```

#### Step 2: B-allele and log ratios

```
mkdir -p results/array
$GW6_DIR/bin/normalize_affy_geno_cluster.pl $CLUSTER_FILE \
    results/apt/quant-norm.pm-only.med-polish.expr.summary.txt \
    -locfile $LOC_FILE -out results/array/gw6.lrr_baf.txt
```

#### Step 3: Split results into a single file per sample

```
perl scripts/penncnv/kcolumn.pl results/array/gw6.lrr_baf.txt split 2 -tab -head 3 \
    -name --output results/array/gw6
```

The normalised files will be placed in `results/array/gw6*`.

```
ls -lh results/array/gw6*
```

The file structure is one probe per line, giving the position, normalized log R and BAF for each probe.

```
less -S results/array/gw6.GSM337641
```


### CNV Calling and Visualisation 

Now that we have the BAF and LRR data we will use OncoSNP to analyse this data.  Create a working directory for OncoSNP.

```
mkdir -p results/oncosnp
```

OncoSNP requires an input file describing the dataset to be analyzed.  A template is available with the OncoSNP install.

```
cp /usr/local/oncosnp/demo/example-batch-file.txt results/oncosnp/batch-file.txt
```

Next we edit this file to point to our data.  Note the tumour data is in gw6.GSM337641 and the normal data in gw6.GSM337662.  Set the name of the sample to `HCC1143`.

```
nano -w results/oncosnp/batch-file.txt
```

OncoSNP has many command line parameters, and most will not change between runs of different datasets.  For convenience we have provided the run_oncosnp.sh.  This will be easier to execute than a large command line.

```
less -S scripts/run_oncosnp.sh
```

Create an output directory for oncosnp.  Run OncoSNP using the script provided.  

We will run OncoSNP only on chromosome 21 since it can take awhile for the whole genome. While it goes we can take a few more minutes to study the script for launching it. Please feel free to ask questions.

```
bash scripts/run_oncosnp.sh results/oncosnp/batch-file.txt results/oncosnp &
```

Rather then print to screen, the script sends the output of OncoSNP to a log file. We can monitor the progress of the program by examining this file.

```
less -S results/oncosnp/run.log
```

Similarly the errors are also sent to a file which we can explore.

```
less -S results/oncosnp/run.err
```

We can see if the script is still running by looking at our background jobs

```
jobs
```

To bring the job back into the foreground, type

```
fg
```

To put it in the background again, suspend it using conrol-z, then type

```
bg
```

When the program finishes we can go to the output folder and browse the results.

```
ls -lh results/oncosnp
```

The first key file is the .qc file which outputs some basic quality control values and some parameters. Probably the most interesting value is the stromal contamination i.e. fraction of normal cells. Two values are reported by default because OncoSNP does multiple analysis runs. The first value is the most probable.

```
less -S results/oncosnp/HCC1143.qc
```

Next is .cnvs file which contains the smoothed segments with there copy number prediction.

```
less -S results/oncosnp/HCC1143.cnvs
```

The last file we will look at is the .cnv file. This is essentially a more informative version of the .cnvs file. One column of particular interest is the "Tumour State" column. This is an integer >= 1 which represents the most likely state of the HMM for that segment. 

```
less -S results/oncosnp/HCC1143.cnvs
```

The final interesting file that OncoSNP produces is the plots HCC1143.*.ps.gz.  This file can be found in the module package under `content/figures/oncosnp` in case you have trouble copying the file from your amazon instance.

## Whole Genome Sequencing (WGS) Analysis

The workflow for WGS data is not dramatically different. We still need to do some normalisation and B allele extraction.

The dataset we will be using is a breast cancer cell line sequenced by the TCGA.  The data has been downloaded and partially processed for you (See data preparation).

```
ln -s /home/ubuntu/CourseData/CG_data/TCGA/HCC1143/
```

It is very important that we use the same genome fasta as was used to generate the bams.  The genome fasta has been downloaded for you.

```
ln -s /home/ubuntu/CourseData/CG_data/Module5/genome/
```

### HMMCopy for total copy number

HMMcopy requires several input files all in .wig format. For most analyses you can use the GC content and mappability files on the HMMcopy website.  For this lab we have created versions for chr21.

```
less -S /home/ubuntu/CourseData/CG_data/Module5/genome/hg19.21.gc.wig
less -S /home/ubuntu/CourseData/CG_data/Module5/genome/hg19.21.map.wig
```

HMMCopy also requires read depth for tumour and normal data.  These are obtained by counting reads overlapping 1000bp segments in the human genome. 

```
less -S /home/ubuntu/CourseData/CG_data/Module5/HCC1143/G15511.HCC1143.1.chr21.wig
less -S /home/ubuntu/CourseData/CG_data/Module5/HCC1143/G15511.HCC1143_BL.1.chr21.wig
```

An R script is provided to run HMMCopy.  

```
less -S scripts/hmmcopy/run.R
```

Run R.

```
R
```

and use the source function to execute the script.  Then quit.

```
source('scripts/hmmcopyrun.R')
quit()
```

The `run.R` script will create two directories: `results/hmmcopy/results` containts the segment files and `results/hmmcopy/plots` contains the various plots generated.  The `.igv.seg` segment file can be loaded into IGV.

```
less -S results/hmmcopy/results/tumour.igv.seg
```

### APOLLOH for allele specific copy number

To obtain BAF data from bam files, we first need to run a genotyper to call heterozygous SNP positions in the normal.  GATK has been run for you (see data preparation).  The output of GATK is a VCF file.

```
mkdir data/vcf
cp /home/ubuntu/CourseData/CG_data/Module5/data/vcf/HCC1143.GATK.vcf data/vcf
less -S data/vcf/HCC1143.GATK.vcf
```

Given heterozygous SNPs, we now need to extract the reference and alternate allele counts from the bam file.  The script `build_apolloh_allelic_counts_file.py` does this for you.

```
mkdir -p data/apolloh
python scripts/build_apolloh_allelic_counts_file.py data/vcf/HCC1143.GATK.vcf \
    data/apolloh/baf.txt --normal_column 1
```

The BAF file is a tab delimited file with one row per SNP.

```
less -S data/apolloh/baf.txt
```

We are now in a position to run APOLLOH.

```
mkdir -p results/apolloh
bash scripts/run_apolloh.sh data/apolloh/baf.txt results/hmmcopy/results/tumour.apolloh.seg \
    results/apolloh/
```

There are three files. params.txt contains the APOLLOH model paramters. Most interesting Normal cell contamination

The `params.txt` file contains the APOLLOH model parameters.

```
less -S results/apolloh/params.txt
```

The `loh.txt` file contains information about the state of each SNP.

```
less -S results/apolloh/loh.txt
```

The `segs.txt` contains information about the segments.

```
less -S results/apolloh/segs.txt
```

### OncoSNP-SEQ 

OncoSNP-SEQ is another tool that predicts allele specific copy number of WGS.  The input is a set of SNP positions a pileup file from the tumour and normal genomes.  Pileup files take a considerable amount of time to generate, and the pileups for chromosome 21 have been generated for you.  They provide information of the read bases that match each reference base with one line per position in the genome.

```
less -S HCC1143/G15511.HCC1143.1.chr21.pileup
```

The pileups are processed by oncosnpseq using the process_pileup.pl script.  Run these in parallel using the `&`.

```
mkdir results/oncosnpseq
perl /usr/local/oncosnpseq/scripts/process_pileup.pl \
    --infile HCC1143/G15511.HCC1143.1.chr21.pileup \
    --outfile results/oncosnpseq/tumour.txt \
    --snpfile genome/oncoseq/bed_chr_21.bed \
    > results/oncosnpseq/tumour.log &
perl /usr/local/oncosnpseq/scripts/process_pileup.pl \
    --infile HCC1143/G15511.HCC1143_BL.1.chr21.pileup \
    --outfile results/oncosnpseq/normal.txt \
    --snpfile genome/oncoseq/bed_chr_21.bed \
    > results/oncosnpseq/normal.log &
```

OncoSNP-SEQ will be run in a similar way to oncosnp.  For convencience we have provided the run_oncosnpseq.sh.  This will be easier to execute than a large command line.

```
less -S scripts/run_oncosnpseq.sh
```

Create an output directory for oncosnpseq.  Run OncoSNP using the script provided.  

We will run OncoSNP in the background since it is slow.  While it goes we can take a few more minutes to study the script for launching it. Please ask questions.

```
bash scripts/run_oncosnpseq.sh results/oncosnpseq/tumour.txt \
    results/oncosnpseq/normal.txt HCC1143 results/oncosnpseq/ &
```

The results are similar to oncosnp and can be found in the module package under `results/oncosnpseq`.

## Visualizing Datasets In IGV

For this part please download the METABRIC dataset from the wiki and open it in IGV.

