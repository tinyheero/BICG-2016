

Environment

    STRELKA_DIR=/usr/local/strelka/

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

