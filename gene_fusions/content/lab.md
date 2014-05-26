# Module 4 Lab - Gene Fusions

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

## Fusion analysis on a simulated dataset

We will run 4 different tools on a simulated RNA-Seq dataset consisting of fused and wild type genes from chromosome 20.  

To make things easier we will store the paths to the fastq files in environment variables.

    FASTQ1=/home/ubuntu/CourseData/CG_data/Module4/tutorial/data/SIM001/SIM001_1.fastq
    FASTQ2=/home/ubuntu/CourseData/CG_data/Module4/tutorial/data/SIM001/SIM001_2.fastq


# ChimeraScan install
export CHIMERASCAN_DIR=/home/ubuntu/CourseData/CG_data/Module4/tutorial/install/chimerascan

# deFuse install
export DEFUSE_DIR=/home/ubuntu/CourseData/CG_data/Module4/tutorial/install/defuse

# Tophat-fusion install
export TOPHATFUSION_DIR=/home/ubuntu/CourseData/CG_data/Module4/tutorial/install/tophat_fusion

# Trinity install
export TRINITY_DIR=/home/ubuntu/CourseData/CG_data/Module4/tutorial/install/trinity


### deFuse


    mkdir -p results/defuse
    $DEFUSE_DIR/bin/scripts/defuse.pl -c $DEFUSE_DIR/data/config.txt -1 $FASTQ1 -2 $FASTQ2 -o results/defuse

### ChimeraScan

    mkdir -p results/chimerascan
    python $CHIMERASCAN_DIR/bin/chimerascan_run.py $CHIMERASCAN_DIR/index $FASTQ1 $FASTQ2 results/chimerascan

### TopHat-Fusion

Tophat-Fusion is run in two steps.  First we run tophat2 to obtain spliced alignments to the genome.  Providing a gene models file is always advised.

******-G

    mkdir -p results/tophat_fusion/tophat_SIM001
    tophat2 --mate-inner-dist 100 --mate-std-dev 30 --fusion-search -o results/tophat_fusion/tophat_SIM001 $HG19_INDEX_BASE $FASTQ1 $FASTQ2

Tophat


    with Sentinal('run_tophat_'+sample_id) as sentinal:

        if sentinal.unfinished:

            utils.makedirs(tophat_fusion_info.sample_tophat_output_directory(sample_id))

            subprocess.check_call(['tophat2',
                                   '--mate-inner-dist', '100',
                                   '--mate-std-dev', '30',
                                   '--fusion-search',
                                   '-o', tophat_fusion_info.sample_tophat_output_directory(sample_id),
                                   info.hg19_index_base,
                                   info.fastq_filename(sample_id, '1'),
                                   info.fastq_filename(sample_id, '2')])

    with Sentinal('run_tophat_fusion_'+sample_id) as sentinal:

        if sentinal.unfinished:

            with utils.CurrentDirectory(tophat_fusion_info.sample_results_directory(sample_id)):

                utils.symlink(tophat_fusion_info.blast_human_directory)
                utils.symlink(tophat_fusion_info.ref_gene_filename)
                utils.symlink(tophat_fusion_info.ens_gene_filename)

                subprocess.check_call(['tophat-fusion-post',
                                       '--num-fusion-reads', '1',
                                       '--num-fusion-pairs', '2',
                                       '--num-fusion-both', '5',
                                       info.hg19_index_base])

### Trinity 


