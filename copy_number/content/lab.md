# Lab Module 5 - Copy Number Analysis

## Setup

First login into the cloud.

Now enter the ~/workspace directory

```
cd ~/workspace
```

Create a directory for this module and enter this directory:

```
mkdir Module5
cd Module5
```

Now let's create a link to some helper scripts we will need for this module:

```
ln -s /home/ubuntu/CourseData/CG_data/Module5/scripts
```


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

## Analysis Of CNAs using Arrays

For calling copy number variants from Affymetrix SNP 6.0 data, we will be using breast cancer cell-line (HC1395). The array data for HCC1395 has already been downloaded for you. Create a directory and copy the array data to your workspace.

### Fetch Array Data

```
ln -s /home/ubuntu/CourseData/CG_data/HCC1395
```

Create a list of the cel files to be used by downstream tools.  In practice we would normalize many arrays in a batch.  For demonstration purposes we use just a single tumour.

```
echo cel_files > cel_file_list.txt
echo `pwd`/HCC1395/cel/GSM888107.CEL >> cel_file_list.txt
```

### Step 1 - Array Normalization

The first step in array analysis is to normalize the data and extract the log R and BAF (B Allele Frequencies). The following steps will create normalized data from Affymetrix SNP 6.0 data. We require a number of data files that define the Affymetrix SNP 6.0 arrays.

The sketch file gives a subset of the probes that work well for normalization.

```
SKETCH_FILE=$GW6_DIR/lib/hapmap.quant-norm.normalization-target.txt
```

The cluster file defines genotype clusters from HapMap, and is used for small batches.

```
CLUSTER_FILE=$GW6_DIR/lib/hapmap.genocluster
```

Chromosome positions for each probe.

```
LOC_FILE=$GW6_DIR/lib/affygw6.hg19.pfb
```

Once these reference files have been defined, we can now perform probeset summarization:

```
$APT_DIR/bin/apt-probeset-summarize --cdf-file $SNP6_CDF \
    --analysis quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true \
    --target-sketch $SKETCH_FILE --out-dir results/apt \
    --cel-files cel_file_list.txt --chip-type GenomeWideEx_6 \
    --chip-type GenomeWideSNP_6
```

### Step 2 -  Extract BAF & LRR

Now that normalization is complete, we can extract the B-allele frequencies (BAF) and log R ratios (LRR).

```
mkdir -p results/array
$GW6_DIR/bin/normalize_affy_geno_cluster.pl $CLUSTER_FILE \
    results/apt/quant-norm.pm-only.med-polish.expr.summary.txt \
    -locfile $LOC_FILE -out results/array/gw6.lrr_baf.txt
```

The BAF and LRR values for every sample in the batch will be placed into a single file. The next step will be split them into sample specific BAF and LRR files for downstream analyses (even though we only have one sample, we will still do this to follow a consistent workflow):

```
perl scripts/penncnv/kcolumn.pl results/array/gw6.lrr_baf.txt split 2 -tab -head 3 \
    -name --output results/array/gw6
```

The sample-specific BAF and LRR files will be placed in `results/array/gw6*`. The file structure is one probe per line, giving the position, normalized log R and BAF for each probe.

```
less -S results/array/gw6.GSM888107
```

| Name          | Chr | Position | GSM888107.CEL Log R Ratio | GSM888107.CEL.B Allele Freq |
|---------------|-----|----------|---------------------------|-----------------------------|
| SNP_A-2131660 | 1   | 1156131  | 0.3040                    | 0.0501                      |
| SNP_A-1967418 | 1   | 2234251  | -0.0355                   | 1.0000                      |
| SNP_A-1969580 | 1   | 2329564  | -0.2625                   | 0.9403                      |
| SNP_A-4263484 | 1   | 2553624  | -0.3366                   | 0.9780                      |
| SNP_A-1978185 | 1   | 2936870  | -0.0276                   | 0.0528                      |
| SNP_A-4264431 | 1   | 2951834  | -0.1812                   | 0.0272                      |
| SNP_A-1980898 | 1   | 3095126  | 0.0830                    | 0.9793                      |

### Step 3 - Call CNV

Now that we have the BAF and LRR data we will use OncoSNP to analyze this data.  Create a working directory for OncoSNP.

```
mkdir -p results/oncosnp
```

OncoSNP has many command line parameters, and most will not change between runs of different datasets. Below is an example of how you could run it:

```
$ONCOSNP_DIR/run_oncosnp.sh $MCR_DIR \
	--sampleid HCC1395 \
	--tumour-file results/array/gw6.GSM888107 \
	--output-dir results/oncosnp \
	--fulloutput --plot \
	--gcdir $GC_DIR \
	--paramsfile $ONCOSNP_DIR/data/cyau/temp/oncosnp/configuration/hyperparameters-affy.dat \
	--levelsfile $ONCOSNP_DIR/configuration/levels-affy.dat \
	--subsample 30 \
	--emiters 1 \
	--female \
	--trainingstatesfile $ONCOSNP_DIR/configuration/trainingStates.dat \
	--tumourstatesfile $ONCOSNP_DIR/configuration/tumourStates.dat \
	--chr 21 \
	--hgtables $ONCOSNP_DIR/configuration/hgTables_b37.txt \
	> results/oncosnp/run.log 2> results/oncosnp/run.err &
```

Some important parameters to consider:

* --tumour-file: Specify the location of where the BAF and LRR values are for the sample
* --chr: Specify the chromosome you want to run on. In this example, we run only on chromosome 21 since it can take awhile for the whole genome. Don't specify this parameter for whole genome analysis.
* --stroma: This parameter can be specified for normal content adjustment. As this is a cell-line, we did not set this. 
* --intratumor: This parameter can be specified for correcting intratumor heterogeneity. As this is a cell-line, we did not set this. 

The `&` character at the end of the above command sends the job to run in the background. Rather then print the progress of the job to screen, this command will send output of OncoSNP to a log file. We can monitor the progress of the program by examining this file.

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
less -S results/oncosnp/HCC1395.qc
```

| LogRRatioShift | NormalContent | Copy Number (Average) | Log-likelihood | OutlierRate | LogRRatioStd | BAlleleFreqStd | PloidyNo |
|----------------|---------------|-----------------------|----------------|-------------|--------------|----------------|----------|
| -0.1710        | 0.0           | 2.0                   | 1156.17639     | 0.010       | 0.237        | 0.041          | 1        |
| -0.1099        | 0.0           | 1.9                   | 1098.29187     | 0.010       | 0.237        | 0.041          | 2        |

Next is .cnvs file which contains the smoothed segments with there copy number prediction.

```
less -S results/oncosnp/HCC1395.cnvs
```

| Chromosome | StartPosition | EndPosition | CopyNumber | LOH | Rank | Loglik       | nProbes | NormalFraction | TumourState | PloidyNo | MajorCopyNumber | MinorCopyNumber |
|------------|---------------|-------------|------------|-----|------|--------------|---------|----------------|-------------|----------|-----------------|-----------------|
| 21         | 10913441      | 11039570    | 2          | 0   | 1    | 11.732701    | 8       | 0.0            | 3           | 1        | 1               | 1               |
| 21         | 14369207      | 48084747    | 2          | 0   | 1    | 31023.688404 | 12520   | 0.0            | 3           | 1        | 1               | 1               |
| 21         | 38323528      | 40232808    | 1          | 1   | 3    | 634.355920   | 763     | 0.0            | 2           | 1        | 1               | 0               |
| 21         | 47133549      | 48084747    | 3          | 0   | 3    | 314.794518   | 265     | 0.0            | 4           | 1        | 2               | 1               |
| 21         | 43993615      | 44503173    | 2          | 2   | 4    | 56.394753    | 200     | 0.0            | 21          | 1        | 2               | 0               |
| 21         | 14369207      | 14775085    | 3          | 0   | 5    | 7.201958     | 24      | 0.0            | 4           | 1        | 2               | 1               |

The last file we will look at is the .cnv file. This is essentially a more informative version of the .cnvs file. One column of particular interest is the "Tumour State" column. This is an integer >= 1 which represents the most likely state of the HMM for that segment. 

```
less -S results/oncosnp/HCC1395.cnvs
```

The final interesting file that OncoSNP produces is the plots HCC1395.*.ps.gz.  This file can be found in the module package under `content/figures/oncosnp` in case you have trouble copying the file from your Amazon instance.

## Analysis Of CNAs using Sequencing Data

The workflow for analyzing CNAs using sequencing data is not dramatically different from microarrays. The major differences is starting with aligned sequencing data (bam files) rather than raw microarray data (e.g. CEL).

The sample we will be using is the same breast cancer cell line we used for the microarray analysis section. This time we will be using the matching exome sequencing data.

### Fetch Sequencing Data

The raw sequencing data has been downloaded and aligned for you (See data preparation). Create a hardlink to the folder containing these data (if not already completed from the "Analysis of CNA using Microarray Data" section already)

```
ln -s /home/ubuntu/CourseData/CG_data/HCC1395
```

### Get Input Data

We will be using TITAN, available as a R bioconductor package (TitanCNA) for the copy number analysis. The program has the ability to perform the normalization, extraction of LRR/BAF, and calling of CNV. But before we can use TITAN, we need to retrieve a few input files:

1. Tumour read count data
2. Normal read count data
3. Tumour allele counts for normal heterozygous positions
4. Genome reference mappability file
5. Genome reference GC content file

Generating these files can take a bit of time. So for this lab, they have been already been generated for you and can be copied for running (Please see the data preparation page for details on how these files were generated

Copy the tumour and normal read count data:

```
mkdir -p hmmCopy/wig
cp /home/ubuntu/CourseData/CG_data/Module5/hmmCopy/wig/* hmmCopy/wig
```

Copy the tumour allele count data:

```
mkdir -p titan/bcftools/tables
cp /home/ubuntu/CourseData/CG_data/Module5/titan/bcftools/tables/* titan/bcftools/tables/
```

Create a link for the folder containing the genome reference along with its mappability and GC content file:

```
ln -s /home/ubuntu/CourseData/CG_data/ref_data
```

### Running TITAN

Once these input files have been retrieved/generated, we can now run TITAN. A R script `scripts/run_titan.R` is provided to run TITAN.

```
Rscript scripts/run_titan.R &> run_titan.log &
```

Just like the oncosnp run, this will run in the background. You can check the progress of the job by going:

```
less -S run_titan.log
```

This will take a few minutes to run. Take this time to review the script itself. Please any questions regarding the content of the script:

```
less -S scripts/run_titan.R
```

The `run_titan.R` script will create the directory `results/titan` which contains the TITAN results and plots. 

## Visualizing Datasets In IGV

For this part please download the METABRIC dataset from the wiki and open it in IGV.

