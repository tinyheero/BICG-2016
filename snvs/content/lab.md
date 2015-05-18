# Module 6 Somatic Mutations - Lab

## Setup

First login into the cloud.

Now enter the ~/workspace directory

```
cd ~/workspace
```

Now we will download the scripts we need for this module from the wiki page.  I have copied the link from the wiki page below. If the command below does not work try copying the link from the wiki and pasting yourself.

```
wget http://bioinformatics.ca/workshop_wiki/images/7/77/Module6.tar.gz
```

Now we will "unzip" the files.

```
tar -zxvf Module6.tar.gz
```

This should create a folder called Module6. Check this is true.

```
ls -lh
```

We can now remove the compressed ".tar.gz" file to keep our workspace clean.

```
rm Module6.tar.gz
```

Enter the Module6 directory

```
cd Module6/
```

## Environment

Module directory:

```
MODULE_DIR=/home/ubuntu/CourseData/CG_data/Module6/
```

MutationSeq install:

```
MUTATIONSEQ_DIR=/usr/local/mutationSeq/
```

SnpEff install:

```
SNPEFF_DIR=/usr/local/snpEff
```

## Prepare Data

We will be using the exome data from the HCC1395 breast cancer cell line.  The tumour and normal bam files have been placed on the server. Create a link to the location of the bam files to allow us to access them quickly.

```
ln -s /home/ubuntu/CourseData/CG_data/HCC1395
```

Additionally, we will need the reference genome for the mutation calling. This has also been placed on the server. We shall create a link to this file for easy access:

```
ln -s /home/ubuntu/CourseData/CG_data/refdata
```

We will restrict our analysis to a 1 Mb region (7Mb and 8Mb) within chromosome 17. Use `samtools` to create bam files containing only alignments within this region.  The argument `17:7000000-8000000` specifies the region, and `-b` specifies the output is bam (default output is sam format).

```
mkdir data
samtools view -b HCC1395/exome/HCC1395_exome_tumour.bam 17:7000000-8000000 \
    > data/HCC1395_exome_tumour.17.7MB-8MB.bam
samtools view -b HCC1395/exome/HCC1395_exome_normal.bam 17:7000000-8000000 \
    > data/HCC1395_exome_normal.17.7MB-8MB.bam
```

Create an index for each bam file.

```
samtools index data/HCC1395_exome_tumour.17.7MB-8MB.bam
samtools index data/HCC1395_exome_normal.17.7MB-8MB.bam
```

If you are unfamiliar with the sam/bam format, take some time to look at the alignments and metadata.

The header information contains information about the reference genome and read groups (per read information about sample/flow cell/lane etc).

```
samtools view -H data/HCC1395_exome_tumour.17.7MB-8MB.bam | less -S
```

The main contents of the file contain read alignments, tab separated, one line per alignment.

```
samtools view data/HCC1395_exome_tumour.17.7MB-8MB.bam | less -S
```

Samtools will also calculate statistics about the reads and alignments.  Unfortunately this information is not cached, and thus this command will take considerable time on a regular sized bam.

```
samtools flagstat data/HCC1395_exome_tumour.17.7MB-8MB.bam
```

## Predicting SNVs

### Strelka

Create a local copy of the strelka config.  Strelka provides aligner specific config files for bwa, eland, and isaac.  Each file contains default configuration parameters that work well with the aligner.  The bam files we are working with were created using bwa, so we select that config file and make a local copy to make changes.

```
mkdir config
cp /usr/local/etc/strelka_config_bwa_default.ini config/strelka_config_bwa.ini
```

Since we will be using Exome data for this, we need to change the `isSkipDepthFilters` parameter in the strelka_config_bwa.ini file. Let's create a new config file for exome analysis:

```
cp config/strelka_config_bwa.ini config/strelka_config_bwa_exome.ini
```

Now let's edit the `config/strelka_config_bwa_exome.ini` and change the `isSkipDepthFilters = 0` to `isSkipDepthFilters = 1`. The reason why we do this is described on the [Strelka FAQ page](https://sites.google.com/site/strelkasomaticvariantcaller/home/faq):

> The depth filter is designed to filter out all variants which are called above a multiple of the mean chromosome depth, the default configuration is set to filter variants with a depth greater than 3x the chromosomal mean. If you are using exome/targeted sequencing data, the depth filter should be turned off...
> 
> However in whole exome sequencing data, the mean chromosomal depth will be extremely small, resulting in nearly all variants being (improperly) filtered out.

A strelka analysis is performed in 2 steps.  In the first step we provide Strelka with all the information it requires to run the analysis, including the tumour and normal bam filenames, the config, and the reference genome.  Strelka will create an output directory with the setup required to run the analysis.

```
perl $STRELKA_DIR/bin/configureStrelkaWorkflow.pl \
    --tumor HCC1395/exome/HCC1395_exome_tumour.17.7MB-8MB.bam \
    --normal HCC1395/exome/HCC1395_normal_tumour.17.7MB-8MB.bam \
    --ref refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
    --config config/strelka_config_bwa_exome.ini \
    --output-dir results/strelka/
```

The output directory will contain a _makefile_ that can be used with the tool _make_.  One benefit of the makefile style workflow is that it can be easily parallelized using _qmake_ on a grid engine cluster.  

To run the strelka analysis, use make and specify the directory constructed by `configureStrelkaWorkflow.pl` with make's `'-C'` option.

```
make -C results/strelka/
```

> If you access to a grid engine cluster, you can replace the command `make` with `qmake` to launch Strelka on the cluster.

Strelka has the benefit of calling SNVs and small indels.  Additionally, Strelka calculates variant quality and filters the data in two tiers.  The filenames starting with `passed` contain high quality candidates, and filenames starting with `all` contain high quality and marginal quality candidates.

The Strelka results are in VCF format, with additional fields explained on the [strelka website](https://sites.google.com/site/strelkasomaticvariantcaller/home/somatic-variant-output).

```
less -S results/strelka/results/passed.somatic.snvs.vcf
less -S results/strelka/results/passed.somatic.indels.vcf
```

### MutationSeq

Running MutationSeq is a one step process.  The tumour and normal bam files and reference genome are provided on the command line. MutationSeq uses supervised learning (random forest) to classify each mutation as true or artifact.  The training set consists of over 1000 known true and positive mutations validated by deep SNV sequencing.  The trained model is provided with the mutationseq package, but must be provided on the command line using the `model:` argument.

Additional parameters include

- `-i` to specify the region
- `-q 1` to remove ambiguously mapping reads
- `--all` to print variants that are classified as low probability.

```
mkdir -p results/mutationseq
python $MUTATIONSEQ_DIR/classify.py \
    tumour:data/HCC1395_exome_tumour.17.7MB-8MB.bam \
    normal:data/HCC1395_normal_tumour.17.7MB-8MB.bam \
    reference:refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
    model:$MUTATIONSEQ_DIR/model_v4.1.1.npz \
    -i 17:7000000-8000000 \
    -c $MUTATIONSEQ_DIR/metadata.config -q 1 -o results/mutationseq/HCC1395.vcf \
    -l results/mutationseq/HCC1395.log --all 
```
    
The results are provided in VCF format.

```
less -S results/mutationseq/HCC1395.vcf
```

Filtering on these results are left to the end-user. A threshold of 0.85 on PR (probability of being a true somatic mutation) will be used in this lab.

## Converting the VCF format into a tabular format

The VCF format is sometimes not useful for visualization and data exploration purposes which often requires the data to be in tabular format. We can convert from VCF format to tabular format using the extractField() function from SnpSift/SnpEff. Since we mutation caller has a different set of output values in the vcf file, the command needs be adjusted for the mutation caller. 

For example, to convert the Strelka vcf file into a tabular format:

```
java -jar $SNPEFF_DIR/SnpSift.jar extractFields -e "."  results/strelka/results/passed.somatic.snvs.vcf CHROM POS ID REF ALT QUAL FILER QSS TQSS NT QSS_NT TQSS_NT SGT SOMATIC GEN[0].DP GEN[1].DP GEN[0].FDP GEN[1].FDP GEN[0].SDP GEN[1].SDP GEN[0].SUBDP GEN[1].SUBDP GEN[0].AU GEN[1].AU GEN[0].CU GEN[1].CU GEN[0].GU GEN[1].GU GEN[0].TU GEN[1].TU > results/strelka/results/passed.somatic.snvs.txt 
```

To convert the MutationSeq vcf file into a tabular format:

```
java -jar $(SNPEFF_DIR)/SnpSift.jar extractFields -e "." results/museq/HCC1395.museq.vcf CHROM POS ID REF ALT QUAL QSS FILTER PR TR TA NR NA TC NI ND > results/museq/HCC1395.museq.txt
```

The -e parameter specifies how to represent empty fields. In this case, the "." character is placed for any empty fields. This facilities loading and completeness of data. For more details on the extractField() function see the [SnpSift documentation](http://snpeff.sourceforge.net/SnpSift.html#Extract).

## Data Exploration

### Integrative Genomics Viewer (IGV) 

Download the chromosome 17 bam files from the Module 6 wiki section.  We can investigate a few predicted positions in IGV:

* 17:7491818
* 17:7578406
* 17:7482929

### Exploration in R

We can further explore the results using R.

Further exploration of the data can be done using R. Open up the rmarkdown file 

  
The predictions for the whole HCC1143 dataset are contained in the Module 6 package at `Module6/content/data/HCC1143.vcf.gz`.  These can also be loaded into IGV.  Additionally, a plot provided by mutationseq for the HCC1143 tumour is located at `Module6/content/data/HCC1143_0.5_0.9_PASS.pdf`.



