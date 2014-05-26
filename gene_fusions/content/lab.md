# Module 4 - Gene Fusions
># Lab

## Setup

First login into the cloud.

Now enter the ~/workspace directory

    cd ~/workspace

Now we will download the scripts we need for this module from the wiki page.  I have copied the link from the wiki page below. If the command below does not work try copying the link from the wiki and pasting yourself.

wget http://bioinformatics.ca/workshop_wiki/images/0/0c/Module4.tar.gz

Now we will "unzip" the files.

    tar -zxvf Module4.tar.gz

This should create a folder called Module4. Check this is true.

    ls -lh

We can now remove the compressed ".tar.gz" file to keep our workspace clean.

    rm Module4.tar.gz

Enter the Module4 directory

    cd Module4/

## Environment

Module directory:

    MODULE_DIR=/home/ubuntu/CourseData/CG_data/Module4/

ChimeraScan install:

    CHIMERASCAN_DIR=$MODULE_DIR/tutorial/install/chimerascan

deFuse install:

    DEFUSE_DIR=$MODULE_DIR/tutorial/install/defuse

TopHat-Fusion install:

    TOPHATFUSION_DIR=$MODULE_DIR/tutorial/install/tophat_fusion

Trinity install:

    TRINITY_DIR=$MODULE_DIR/tutorial/install/trinity

RSEM install:

    RSEM_DIR=$MODULE_DIR/tutorial/install/rsem/packages/rsem-1.2.12

## Fusion analysis on a simulated dataset

We will run 4 different tools on a simulated RNA-Seq dataset consisting of fused and wild type genes from chromosome 20.  

To make things easier we will store the paths to the fastq files in environment variables.

    FASTQ1=$MODULE_DIR/tutorial/data/SIM001/SIM001_1.fastq
    FASTQ2=$MODULE_DIR/tutorial/data/SIM001/SIM001_2.fastq

The information about the simulated fusions is stored in a tab separated file and can be viewed using the following command.

    column -n -t $MODULE_DIR/tutorial/data/SIM001/simulation/fusions.tsv | less -S

    cut -f8,16,17 $MODULE_DIR/tutorial/data/SIM001/simulation/fusions.tsv | uniq

### ChimeraScan

ChimeraScan requires that a number of its python libraries are in the python path.  Update the `PYTHONPATH` environment variable.

    export PYTHONPATH=$PYTHONPATH:$CHIMERASCAN_DIR/lib/python2.7/site-packages

Running ChimeraScan is a one step process.  Configuration can be in an xml configuration file, in addition to on the command line.  Options relate to parallelization and filtering.  Type the following to see available options.

    $CHIMERASCAN_DIR/bin/chimerascan_run.py --help

Provide the index directory containing the reference genome and transcriptome on the command line, in addition to the pair of fastq files.

    mkdir -p results/chimerascan
    python $CHIMERASCAN_DIR/bin/chimerascan_run.py $CHIMERASCAN_DIR/index $FASTQ1 $FASTQ2 results/chimerascan

### deFuse

Running deFuse is a one step process.  Configuration can be in a text file of key value pairs separated by `=`, with a default configuration provided in the package.

    mkdir -p results/defuse
    $DEFUSE_DIR/packages/defuse/scripts/defuse.pl -c $DEFUSE_DIR/data/config.txt \
        -1 $FASTQ1 -2 $FASTQ2 -o results/defuse

### TopHat-Fusion

Tophat-Fusion is run in two steps.  First we run tophat2 to obtain spliced alignments to the genome.  Providing a gene models file is always advised (`-G` option).

    mkdir -p results/tophat_fusion/tophat_SIM001
    tophat2 --mate-inner-dist 100 --mate-std-dev 30 --fusion-search \
        -G $TOPHATFUSION_DIR/data/ucsc_hg19.gtf \
        -o results/tophat_fusion/tophat_SIM001 \
        $TOPHATFUSION_DIR/data/ucsc_hg19 $FASTQ1 $FASTQ2

The second step is to post-process the `tophat2` results using `tophat-fusion-post`.  The `tophat-fusion-post` tool requires a specific working directory structure.  

    cd results/tophat_fusion

The `tophat-fusion-post` tool will search the working directory for directories prefixed with `tohat_` and will assume these directories contain the results of tophat for samples of interest.

Additionaly, `tophat-fusion-post` requires 3 datasets in specific subdirectories of the working directory, a blast database for filtering and transcript information from ensmbl and refseq.

    ln -s $TOPHATFUSION_DIR/data/blast_human/ blast
    ln -s $TOPHATFUSION_DIR/data/refGene.txt
    ln -s $TOPHATFUSION_DIR/data/ensGene.txt

We now run the `tophat-fusion-post` tool.

    tophat-fusion-post --num-fusion-reads 1 \
        --num-fusion-pairs 2 \
        --num-fusion-both 5 \
        $TOPHATFUSION_DIR/data/ucsc_hg19

### Trinity and gmap

    mkdir results/trinity
    cp $FASTQ1 results/reads_1.fq
    cp $FASTQ2 results/reads_2.fq
    $TRINITY_DIR/bin/trinity/Trinity --seqType fq --JM 12G --CPU 4 \
        --left results/reads_1.fq --right results/reads_2.fq \
        --output results/trinity

    gmap -f gff3_gene -D $TRINITY_DIR/data/gmap -d hg19 results/trinity/Trinity.fasta \
        > results/trinity/Trinity.gff3

    python ./scripts/gmap_extract_fusions.py results/trinity/Trinity.gff3 \
        results/trinity/Trinity.fasta results/trinity/Fusions.fasta

### RSEM for fusion expression

    mkdir -p results/rsem
    $RSEM_DIR/rsem-prepare-reference results/trinity/Trinity.fasta \
        results/rsem/trinity

    $RSEM_DIR/rsem-calculate-expression -p 8 --paired-end \
        $FASTQ1 $FASTQ2 results/rsem/trinity results/rsem/expression



