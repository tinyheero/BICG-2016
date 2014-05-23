#!/bin/bash

ref=$1  # Reference genome

normalBam=$2 # BAM File

tumourBam=$3 # BAM File

outFile=$4

gatkJar=/usr/local/GATK/GenomeAnalysisTK.jar
#gatkJar=~/Documents/teaching/cbw/2013/install/src/GATK/GenomeAnalysisTK.jar

######### EXECUTE GATK UNIFIED GENOTYPER ###################################################
######### Requires Java to be installed ####################################################
######### -R specifies the reference genome fasta file #####################################
######### -T specifies the tool we want to execute; we want UnifiedGenotyper ###############
######### -L specifies the region we want to analyse
######### -I specifies input BAM filename
java -jar $gatkJar -R $ref -T UnifiedGenotyper -baq RECALCULATE -L 21:19000000-20000000 -I $normalBam -I $tumourBam -o $outFile
