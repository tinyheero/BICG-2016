# Module 6 - Somatic Mutations Lab

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

## Predicting SNVs

GATK is a well used tool kit providing many bioinformatic functions, including SNV calling using the _UnifiedGenotyper_.  GATK is quite strict about its inputs, bam's and reference genomes must be constructed exactly as specified otherwise GATK will report an error and exit.  For example: the bam files and the genome must match exactly, and the reads must have properly defined read groups (per read information about sample/flow cell/lane etc).

Run GATK using java as follows.

```
mkdir -p results/gatk
java -jar $GATK_DIR/GenomeAnalysisTK.jar \
    -R HCC1143/Homo_sapiens_assembly19.fasta \
    -T UnifiedGenotyper \
    -baq RECALCULATE \
    -L 21:19000000-20000000 \
    -I data/G15511.HCC1143_BL.1.chr21.19M-20M.bam \
    -I data/G15511.HCC1143.1.chr21.19M-20M.bam \
    -o results/gatk/HCC1143.vcf
```

Command line options: 

- -R specifies the reference genome fasta file  
- -T specifies the tool we want to execute; we want UnifiedGenotyper 
- -L specifies the region we want to analyse
- -I specifies input BAM filename
- -baq RECALCULATE refers to Base Alignment Quality calculation, see the [samtools docs](http://samtools.sourceforge.net/mpileup.shtml) for details

The raw GATK results are provided in VCF format.

```
less -S results/gatk/HCC1143.vcf
```

GATK does not automatically call SNVs as somatic or germline.  Instead it calls variants in the tumour sample and normal sample independently.  To select somatic variant's, we use a simple script to select all variant's called in the tumour but not the normal.

```
python scripts/call_somatic_mutations.py results/gatk/HCC1143.vcf \
    --normal_column 1 --min_genotype_quality 30 > results/gatk/HCC1143.somatics.txt
```

The results is a simple text file containing chromosome and position of somatic variants.

```
less -S results/gatk/HCC1143.somatics.txt
```

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


