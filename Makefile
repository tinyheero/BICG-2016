.PHONY : ref

# Create a subset of the genome to be used
ref : refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.17.fa 
hmmcopy.ref : refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.17.gc.wig refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.17.map.bw 

GENERATEMAP_OPTS ?= 

# create index of genome fa for quick indexing
refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.fai : refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
	samtools faidx $<

refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.%.fa : refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.fai
	samtools faidx refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa $* > $@

###
# HMM-Copy references

# GC content wig file
refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.%.gc.wig : refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.%.fa
	~/share/usr/HMMcopy-0.1.1/bin/gcCounter $< > $@

# Mappability BigWig file
# This will fail if the bowtie index is not built yet. 
# To build, specify -b in the GENERATEMAP_OPTS. Once built need to re-run with
# the parameter
refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.%.map.bw : refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.%.fa
	~/share/usr/HMMcopy-0.1.1/util/mappability/generateMap.pl $(GENERATEMAP_OPTS) $<
