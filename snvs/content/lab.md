

Environment

    STRELKA_DIR=/usr/local/strelka/
    MUTATIONSEQ_DIR=/usr/local/mutationSeq/
    GATK_DIR=/usr/local/GATK/

Prepare Data

    samtools view -b HCC1143/G15511.HCC1143.1.chr21.bam 21:19000000-20000000 > data/G15511.HCC1143.1.chr21.19M-20M.bam
    samtools view -b HCC1143/G15511.HCC1143_BL.1.chr21.bam 21:19000000-20000000 > data/G15511.HCC1143_BL.1.chr21.19M-20M.bam
    samtools index data/G15511.HCC1143.1.chr21.19M-20M.bam
    samtools index data/G15511.HCC1143_BL.1.chr21.19M-20M.bam


Create a local copy of the strelka config.

    mkdir config
    cp $STRELKA_DIR/etc/strelka_config_bwa_default.ini config/strelka_config_bwa.ini

Edit `binSize` to 250000000 in `config/strelka_config_bwa.ini`.  The file should contain the line

    binSize = 250000000

near the bottom of the file.

Configure strelka

    perl $STRELKA_DIR/bin/configureStrelkaWorkflow.pl --tumor data/G15511.HCC1143.1.chr21.19M-20M.bam \
        --normal data/G15511.HCC1143_BL.1.chr21.19M-20M.bam --ref HCC1143/Homo_sapiens_assembly19.fasta \
        --config $STRELKA_DIR/etc/strelka_config_bwa_default.ini --output-dir results/strelka/


Run strelka

    make -C results/strelka/


Run MutationSeq

    python $MUTATIONSEQ_DIR/classify.py \
        tumour:data/G15511.HCC1143.1.chr21.19M-20M.bam \
        normal:data/G15511.HCC1143_BL.1.chr21.19M-20M.bam \
        reference:HCC1143/Homo_sapiens_assembly19.fasta \
        model:$MUTATIONSEQ_DIR/model_v4.1.1.npz \
        -i 21:19000000-20000000 \
        -c $MUTATIONSEQ_DIR/metadata.config -q 1 -o results/mutationseq/HCC1143.vcf \
        -l results/mutationseq/HCC1143.log --no_filter --all --coverage 0
    


Run GATK

requires dict file.


    mkdir -p results/gatk
    java -jar $GATK_DIR/GenomeAnalysisTK.jar -R HCC1143/Homo_sapiens_assembly19.fasta \
        -T UnifiedGenotyper -baq RECALCULATE -L 21:19000000-20000000 \
        -I data/G15511.HCC1143_BL.1.chr21.19M-20M.bam \
        -I data/G15511.HCC1143.1.chr21.19M-20M.bam \
        -o results/gatk/HCC1143.vcf

Command line options: 

* -R specifies the reference genome fasta file  
* -T specifies the tool we want to execute; we want UnifiedGenotyper 
* -L specifies the region we want to analyse
* -I specifies input BAM filename



