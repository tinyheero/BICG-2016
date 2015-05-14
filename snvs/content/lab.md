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

Strelka install:

```
STRELKA_DIR=/usr/local/strelka/
```

MutationSeq install:

```
MUTATIONSEQ_DIR=/usr/local/mutationSeq/
```

GATK install:

```
GATK_DIR=/usr/local/GATK/
```

SnpEff install:

```
SNPEFF_DIR=/usr/local/snpEff
```

## Prepare Data

We will be using data from the HCC1143 breast cancer cell line.  The bam files for chromosome 21 are stored on the server.  Create a link to the location of the bam files to allow us to access them quickly.

```
ln -s /home/ubuntu/CourseData/CG_data/TCGA/HCC1143/
```

We will restrict our analysis to chromosome 21, between 19Mb and 20Mb.  Use samtools to create a bam file containing only alignments within this region.  The argument `21:19000000-20000000` specifies the region, and `-b` specifies the output is bam, not sam.

```
mkdir data
samtools view -b HCC1143/G15511.HCC1143.1.chr21.bam 21:19000000-20000000 \
    > data/G15511.HCC1143.1.chr21.19M-20M.bam
samtools view -b HCC1143/G15511.HCC1143_BL.1.chr21.bam 21:19000000-20000000 \
    > data/G15511.HCC1143_BL.1.chr21.19M-20M.bam
```

Create an index for each bam file.

```
samtools index data/G15511.HCC1143.1.chr21.19M-20M.bam
samtools index data/G15511.HCC1143_BL.1.chr21.19M-20M.bam
```

If you are unfamiliar with the sam/bam format, take some time to look at the alignments and metadata.

The header information contains information about the reference genome and read groups (per read information about sample/flow cell/lane etc).

```
samtools view -H data/G15511.HCC1143.1.chr21.19M-20M.bam | less -S
```

The main contents of the file contain read alignments, tab separated, one line per alignment.

```
samtools view data/G15511.HCC1143.1.chr21.19M-20M.bam | less -S
```

Samtools will also calculate statistics about the reads and alignments.  Unfortunately this information is not cached, and thus this command will take considerable time on a regular sized bam.

```
samtools flagstat data/G15511.HCC1143.1.chr21.19M-20M.bam
```

## Predicting SNVs

### Strelka

Create a local copy of the strelka config.  Strelka provides config files for the bwa, eland, and isaac.  Each file contains default configuration parameters that work well with the aligner.  The bam files we are working with were created using bwa, so we select that config file and make a local copy to make changes.

```
mkdir config
cp $STRELKA_DIR/etc/strelka_config_bwa_default.ini config/strelka_config_bwa.ini
```

The `binSize` config option allows a user to parallelize by genomic regions of a certain size.  Since we are not going to run jobs in parallel using a cluster, we will set this to a reasonably high value.  Edit `config/strelka_config_bwa.ini` and set `binSize` to `250000000`.  The file should then contain the line

```
binSize = 250000000
```

near the bottom of the file.

A strelka analysis is performed in 2 steps.  In the first step we provide strelka with all the information it requires to run the analysis, including the tumour and normal bam filenames, the config, and the reference genome.  Strelka will create an output directory with the setup required to run the analysis.

```
perl $STRELKA_DIR/bin/configureStrelkaWorkflow.pl \
    --tumor data/G15511.HCC1143.1.chr21.19M-20M.bam \
    --normal data/G15511.HCC1143_BL.1.chr21.19M-20M.bam \
    --ref HCC1143/Homo_sapiens_assembly19.fasta \
    --config config/strelka_config_bwa.ini \
    --output-dir results/strelka/
```

The output directory will contain a _makefile_ that can be used with the tool _make_.  One benefit of the makefile style workflow is that it can be easily parallelized using _qmake_ on a grid engine cluster.  

To run the strelka analysis, use make and specify the directory constructed by `configureStrelkaWorkflow.pl` with make's `'-C'` option.

```
make -C results/strelka/
```

Strelka has the benefit of calling SNVs and small indels.  Additionally, strelka calculates variant quality and filters the data in two tiers.  The filenames starting with `passed` contain high quality candidates, and filenames starting with `all` contain high quality and marginal quality candidates.

The strelka results are in VCF format, with additional fields explained on the [strelka website](https://sites.google.com/site/strelkasomaticvariantcaller/home/somatic-variant-output).

```
less -S results/strelka/results/passed.somatic.snvs.vcf
less -S results/strelka/results/passed.somatic.indels.vcf
```


### MutationSeq

Running mutationseq is a one step process.  The tumour and normal bam files and reference genome are provided on the command line.  MutationSeq uses supervised learning (random forest) to classify each mutation as true or artifact.  The training set consists of over 1000 known true and positive mutations validated by deep SNV sequencing.  The trained model is provided with the mutationseq package, but must be provided on the command line using the `model:` argument.

Additional parameters include

- `-i` to specify the region
- `-q 1` to remove ambiguously mapping reads
- `--all` to print variant's that are classified as low probability.

```
mkdir -p results/mutationseq
python $MUTATIONSEQ_DIR/classify.py \
    tumour:data/G15511.HCC1143.1.chr21.19M-20M.bam \
    normal:data/G15511.HCC1143_BL.1.chr21.19M-20M.bam \
    reference:HCC1143/Homo_sapiens_assembly19.fasta \
    model:$MUTATIONSEQ_DIR/model_v4.1.1.npz \
    -i 21:19000000-20000000 \
    -c $MUTATIONSEQ_DIR/metadata.config -q 1 -o results/mutationseq/HCC1143.vcf \
    -l results/mutationseq/HCC1143.log --all 
```
    
The results are provided in VCF format.

```
less -S results/mutationseq/HCC1143.vcf
```

## Converting the VCF format into a tabular format

The VCF format is sometimes not useful for visualization and data exploration purposes which often requires the data to be in tabular format. We can convert from VCF format to tabular format using the extractField() function from SnpSift/SnpEff.

```
SNPEFF_DIR=~/share/usr/snpEff/snpEff-4.0
java -jar $SNPEFF_DIR/SnpSift.jar extractFields -e "."  passed.somatic.snvs.vcf CHROM POS ID REF ALT QUAL QSS TQSS NT QSS_NT TQSS_NT SGT SOMATIC GEN[0].DP GEN[1].DP GEN[0].FDP GEN[1].FDP GEN[0].SDP GEN[1].SDP GEN[0].SUBDP GEN[1].SUBDP GEN[0].AU GEN[1].AU GEN[0].CU GEN[1].CU GEN[0].GU GEN[1].GU GEN[0].TU GEN[1].TU > passed.somatic.snvs.txt
```

The -e parameter specifies how to represent empty fields. In this case, the "." character is placed for any empty fields. This facilities loading and completeness of data. For more details on the extractField() function see the [SnpSift documentation](http://snpeff.sourceforge.net/SnpSift.html#Extract).

## Data Exploration

Download the chromosome 21 bam files from the Module 6 wiki section.  View these bam files in IGV.

Further exploration of the data can be done using R. Open up the rmarkdown file 

  
The predictions for the whole HCC1143 dataset are contained in the Module 6 package at `Module6/content/data/HCC1143.vcf.gz`.  These can also be loaded into IGV.  Additionally, a plot provided by mutationseq for the HCC1143 tumour is located at `Module6/content/data/HCC1143_0.5_0.9_PASS.pdf`.



