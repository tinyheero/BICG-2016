#!/bin/bash

ref=$1  # Reference genome

normalBam=$2 # BAM File

tumourBam=$3 # BAM File

outDir=$4

jsm.py train $ref $normalBam $tumourBam $outDir/params.txt --skip_size 100 >$outDir/train.log 2>$outDir/train.err

jsm.py classify $ref $normalBam $tumourBam --parameters_file $outDir/params.txt --somatic_threshold 0.01 --post_process --out_file $outDir/results.tsv >$outDir/classify.log 2>$outDir/classify.err
