# HCC1395 dataset location
# https://github.com/genome/gms/wiki/HCC1395-WGS-Exome-RNA-Seq-Data
#
.PHONY : ref hmmcopy.gc.wig hmmcopy.map hmmcopy.ref

# Create a subset of the genome to be used
ref : refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.17.fa 
hmmcopy.gc.wig : refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.gc.wig refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.17.gc.wig 
hmmcopy.map : refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.map.ws_1000.wig 
#	refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.17.fa.map.ws_1000.bw 
hmmcopy.ref : hmmcopy.gc.wig hmmcopy.map

GENERATEMAP_OPTS ?= 

# create index of genome fa for quick indexing
refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.fai : refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
	samtools faidx $<

refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.%.fa : refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.fai
	samtools faidx refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa $* > $@

###
# HMM-Copy references

## GC content wig file
refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.gc.wig : refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
	~/share/usr/HMMcopy-0.1.1/bin/gcCounter $< > $@

# GC content wig file (for specific chromosome)
refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.%.gc.wig : refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.%.fa
	~/share/usr/HMMcopy-0.1.1/bin/gcCounter $< > $@

# Mappability BigWig file
# This will fail if the bowtie index is not built yet. 
# To build, specify -b in the GENERATEMAP_OPTS. Once built need to re-run with
# 1) get BigWig mappability file, then 2) generate mappbility wig file
# the parameter
refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.map.bw : refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
	~/share/usr/HMMcopy-0.1.1/util/mappability/generateMap.pl $(GENERATEMAP_OPTS) $< -o $@

refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.%.fa.map.bw : refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.%.fa
	~/share/usr/HMMcopy-0.1.1/util/mappability/generateMap.pl $(GENERATEMAP_OPTS) $< -o $@

# Convert the mappability BigWig File into a wig file
refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.map.ws_1000.wig : refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.map.bw
	~/share/usr/HMMcopy-0.1.1/bin/mapCounter -w 1000 $< > $@
