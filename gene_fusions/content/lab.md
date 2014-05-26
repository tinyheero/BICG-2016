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


## Data Exploration

In the following section we will use the ucsc genome browser's [online blat](http://genome.ucsc.edu/cgi-bin/hgBlat?command=start) to explore a number of example positive and negative fusion transcripts.

### 313B RNA-Seq

#### Example 1

    >19455
    CAGGAAATCTTCAGCAAGCTGTCTTACTTCTTTTGGCCAAGTCGCACTCCACATTAGA
    GTTTGCCTATCAGGTCTTATTTGATCCACAATCTTCCTTATTTGGGGT|TCAAATCCC
    ATATCCAGCATCCTATCAGCTTCATCCAACACTAAGTACTTGCAGAAGTCTAATCC

Notice that in the alignment summary, the best matches overlap in the query sequence.  In deFuse this prediction is given a lower probability according to the breakpoint homology feature.

       ACTIONS      QUERY           SCORE START  END QSIZE IDENTITY CHRO STRAND  START    END      SPAN
    ---------------------------------------------------------------------------------------------------
    browser details 19455            132     1   141   171  97.2%    17   +   62499145  62499374    230
    browser details 19455             96     2   154   171  88.8%    22   +   38890737  38890976    240
    browser details 19455             67   104   171   171 100.0%     Y   -   15027123  15027592    470
    browser details 19455             62   108   171   171  98.5%     X   -   73350813  73350876     64
    browser details 19455             40    88   141   171  87.1%     6   -   74118977  74119030     54
    browser details 19455             32    84   116   171 100.0%    22   -   49607171  49613016   5846
    browser details 19455             22    92   114   171 100.0%     4   +   72709432  72709456     25
    browser details 19455             21    95   116   171 100.0%     6   +  131014581 131014603     23
    browser details 19455             20    93   112   171 100.0%     2   -  157047083 157047102     20
    browser details 19455             20   101   122   171  95.5%    14   -   29461851  29461872     22
    browser details 19455             20     1    20   171 100.0%     3   +  175393751 175393770     20

#### Example 2

    >29877
    CTGTTTCCTCTTTTACCAAGGACCCGCCAACATGGGCCG|GCTATCTTGTTGCGGAGCTT
    CTTGCTGGGGATAATGGCGATCTCCTTGCACACGCGCTTGTTCGTGTG

Blat the sequence.  Notice in the alignment summary, one part of the query sequence maps to many locations in the genome.  Although one of the alignments may represent a true fusion, the prediction is more likely a mapping artifact.

       ACTIONS      QUERY           SCORE START  END QSIZE IDENTITY CHRO STRAND  START    END      SPAN
    ---------------------------------------------------------------------------------------------------
    browser details YourSeq           68    40   107   107 100.0%    22   -   32435560  32435627     68
    browser details YourSeq           64    40   107   107  97.1%    15   +   82824392  82824459     68
    browser details YourSeq           64    40   107   107  97.1%    15   +   83208735  83208802     68
    browser details YourSeq           62    40   107   107  95.6%     5   -  116052023 116052090     68
    browser details YourSeq           60    40   107   107  94.2%    17   -   29158008  29158075     68
    browser details YourSeq           58    40   107   107  92.7%     1   +  167131923 167131990     68
    browser details YourSeq           57    44   102   107  98.4%     6   -   63257401  63257459     59
    browser details YourSeq           54    44   107   107  92.2%     6   +   50825228  50825291     64
    browser details YourSeq           54    41   106   107  85.8%     3   +  143574838 143574900     63
    browser details YourSeq           37    44    88   107  91.2%     7   -   37156486  37156530     45
    browser details YourSeq           35     1    35   107 100.0%    15   -   82824833  82824867     35
    browser details YourSeq           35     1    35   107 100.0%    15   -   83209176  83209210     35
    browser details YourSeq           34    40    75   107  97.3%    11   +  110976470 110976505     36
    browser details YourSeq           30     8    39   107  96.9%    22   +   32435452  32435483     32
    browser details YourSeq           20     1    20   107 100.0%     5   +  144765194 144765213     20

#### Example 3

    >40571
    TAGAATTAGAATTGTGAAGATGATAAGTGTAGAGGGAAGGTTAATAGTTGATATTGCTAG
    TGTGGCGCTTCCAATTAGGTGCATGAGTAGGTGGCCTGCAGTAAT|GTTAGCGACAGGGA
    GGGATGCGCGCCTGGGTGTAGTTGTGGGGGAGGAAGTGGCTAGCTCAGGGCTTCAGGGGA
    CAGACAGGGAGAGATGACTGAG

Blat the sequence and select the first alignment result.  Ensure you have the NUMT track turned on in UCSC.  This fusion prediction is actually a NUMT insertion in the patient's genome.

#### Example 4

    >33864
    CCAGGGCGCCATTGAGCGGCGAGGGGGTGAGGGGGTTGACGGTGGCGGTGGTCCTGGTCG
    CGGTGGAAAGCATCCCTAGCGAAGGGGACTTGGGCTCATGGCTCATGCCTG|CACCAGTA
    AGGTCTGGTCCGTCCTCCTCCCGGCTGCTCTGCAGACACTGTGCTGGCCTCAGCTCCTGG
    GCCATCCTGGGGCCTCTGGGCAG
    >17735
    CATGGGCACGCGCTTGGGTGTGCTGGCGGGGGAGCTGTGGTTGGTGGCCGGAGAGGACAC
    GGGGGACGACTCGCTGCTCAGTGAGGACC|CTGCACCAGTAAGGTCTGGTCCGTCCTCCT
    CCCGGCTGCTCTGCAGACACTGTGCTGGCCTCAGCTCCTGGGCCATCCTGGGGCCTCTGG
    GCAGGGTCTCCGTGGGGGCGCGTGGCCGGGTCTCGGACT

Blat both sequences and select the first alignment results.  These are examples of likely read through chimeras.

#### Example 5

    >5655
    TGATCAAGCAACTTCCCTGAGGATCCTCAACAATGGTCATGCTTTCAACGTGGAGTTTGA
    TGACTCTCAGGACAAAGC|AGAACGTAAGCTCCATGAGGACCAGGAAGTCTGTCTGCTTT
    GTTCACTGCTGGATCCCGTGACTCGGAACAGTGCACGTAACAGGTGTTCAATAAACCTTT
    GTTGAATGAATAAGTGAA
    >11908
    TCTGTTTCCTATGATCAAGCAACTTCCCTGAGGATCCTCAACAATGGTCATGCTTTCAAC
    GTGGAGTTTGATGACTCTCAGGACAAAGC|AGGGGCTCTTTCCAGGATTCCTGGGTGATG
    GTGCATGATTCTAACAAGCAACAACAGAGGATGAACCCCCGGCCAGATTCAGAAAACCCC
    ACGCCCCTTCCAGGCA

Blat both sequences.  Select browser for any alignment result and then browse to `chr1:110,721,056-110,731,655` and `chr8:86,373,054-86,382,253`.

#### Example 6

    >37434
    TTGGCATCAAATAGATGAACAGGAGAAAAGCTGTTTTAATGTATGTACTCACAGATGGGA
    ATCCCACAAGAATATGAGACTTAAAGAACAGGCCAGGT|TATTCCAGGATCTTTGGAGAC
    CCGAGGAAAGCCGTGTTGACCAAAAGCAAGACAAATGACTCACAGAGAAAAAAGATGGCA
    GAACCAAGG
    >26643
    CCGGGACAGTCTGAATCATGTCCTTCAGTAAGCCAGCCCATCTACCAGCTGTTCAGAACC
    TGACGGCTTTAGTTGCCCTTGGTTCTGCCATCTTTTTTCTCTGTGAGTCATTTGTCTTGC
    TTTTGGTCAACACGGCTTTCCTCGGGTCTCCAAAGATCCTGGAATAACCT|TCCTGGTGG
    AGTAGAAGTAGTCTATAGCTTCTCCTTGGTAGTCCAGATGGGTCTCCCCAGCCAATGCAT
    AACTCTCTCTTTGCCTTTTGATTCAGAGGCATGTGGAGCTCAGCGTGGCCAGGT
    >29118
    AGTAAGCCAGCCCATCTACCAGCTGTTCAGAACCT|AGAGGTCTTAGTTCCGGAGGGAGG
    AATGCTGCCACCAGGAGACACAACAATGATTCAATTAAACTAGAATTTACGACTGC
    >29539
    AGGCACACTCAAACAACGACTGGTCCTCACTCACAACTGATAAGGCTTCCTTGATATGAG
    CTGCTGGGTCCGGGACAGTCTGAATCATGTCCTTCAGTAAGCCAGCCCATCTACCAGCTG
    TTCAGAACCTGACGGCTTTAGTTGCCCTTGGTTCTGCCATCTTTTTTCTCTGTGAGTCAT
    TTGTCTTGCTTTTGGTCAACACGGCTTTCCTCGGGTCTC|CAAAGCCATCTTGCTGTTAT
    CAACAGCATCGAGTAATGATAGGTATCTGGAATGTTCAATATGACCTGCCGCGCTCCAGG
    CGGCGCTCCCCGCCCCTCGCCCTCCGCCTCCGCCTCCGCCTCCTGCTTAGCTCGCGCCTA
    CTCG
    >29538
    TTCAGAACCTGACGGCTTTAGTTGCCCTTGGTTCTGCCATCTTTTTTCTCTGTGAGTCAT
    TTGTCTTGCTTTTGGTCAACACGGCTTTCCTCGGGTCTCCAAAGATCCTGGAATA|ACCT
    GTCCAGTAGTTCTGTAGCGGAGCAGGGCAGGTCCTACTTCTTCAAAAGCACTCAGTAAAG
    GTGGGGAAGTCCTGAGCAACCT
    >29535
    AGGCTTCCTTGATATGAGCTGCTGGGTCCGGGACAGTCTGAATCATGTCCTTCAGTAAGC
    CAGCCCATCTACCAGCTGTTCAGAACCTGACGGCTTTAGTTGCCCTTGGTTCTGCCATCT
    TTTTTCTCTGTGAGTCATTTGTCTTGCTTTTGGTCAACACGGCTTTCCTCGGGTCTCCAA
    AGATCCTGGAATA|ACCTGCCGCGCCGCGCTCCTCACACCCGCTTTCACCTCCGGGCGGG
    GCAGGGGGCATCGGCGGGTCCCAGGCGCCCAGGTTCCCCTCCCCAGCCCGGACCCCGAGC
    CGGGACCCTGGTACCGGCGCCGCTCACCTGCCGCGCTCCAGGCGGCGCTCCCCGCCCCTC
    GCCCTCCGCCTCCGCCTCCGCCTCCTGCTTAGCTCGCGCCTACTCGGC
    >37354
    ACCCTCCAAAGCAACATGAAATGAAACCAAACCACAATAACAACCAAATGAAATAAGACT
    GACAAGAAGTATGCGGTCATGGCCAATACATGGCT|CGATTTTTTTTTCTTTAACATGCA
    CCTTCCTGAGCAAATAAAGGGCTTTTTTCCACCCCTTCCCGCTTGGCTTTAAATGACCAA
    AGAATATT
    >20051
    TGGAATGTTCAATATGACCTGCCGCGCTCCAGGCGGCGCTCCCCGCCCCTCGCCCTCCGC
    CTCCGCCTCCGCCTCCTGCTTAGCTCGCGCCTA|CTCCAGCGACTATGGACAGACTTCCA
    AGATGAGCCCACGCGTCCCTCAGCAGGATTGGCTGTCTCAACCCCCAGCCAGGGTCACCA
    TCAAAATGGAATGTAACCCTAGCCAGGTGAATGGCTCAAGGAACTCTCCTGATGAATGCA
    GTGTG

Blat all the sequences and select browser for any alignment result.  Then navigate to `TMPRSS2` and `ERG`.

### 313B 


