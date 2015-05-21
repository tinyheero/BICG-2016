# The Makefile for the copy number module in CBW
.PHONY : oncosnp 

SAMTOOLS ?= ~/share/usr/samtools-0.1.18/samtools
BCFTOOLS ?= ~/share/usr/samtools-0.1.18/bcftools/bcftools
JAVA ?= ~/share/usr/jre1.7.0_25/bin/java
REF_FASTA ?= ref_data/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa

SNPEFF_DIR ?= ~/share/usr/snpEff/snpEff-4.0

VCF_FIELDS ?= CHROM POS REF ALT DP4
DBSNP ?= 

CHROMOSOMES = 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22

HMMCOPY_DIR ?= ~/share/usr/HMMcopy-0.1.1

, := ,
space :=
space +=
$(space) := 
$(space) +=
comma := ,

INIT = mkdir -p $(@D); set -o pipefail;

### 
## Lab Wiki
###

# convert markdown into a mediawiki format to upload onto the CBW wiki
lab.wiki : 
	pandoc -t mediawiki -o lab.wiki lab.md

###
## Oncosnp
## HCC1395 SNP 6.0 Microarray (http://165.112.7.20/geo/query/acc.cgi?acc=GSM888107) from the CCLE
###

oncosnp :
	$(MAKE) -f ~/share/modules/oncosnp.mk ONCOSNP_VERSION=1.3 STROMA=false INTRATUMOUR=false oncosnp

###
## HMMCopy
###
hmmcopy.readcount : hmmCopy/wig/HCC1395_exome_tumour.wig hmmCopy/wig/HCC1395_exome_normal.wig

READCOUNTER_OPTS ?= 

# generate tumour/normal read counts
hmmCopy/wig/%.wig :
	mkdir -p $(@D); \
	$(HMMCOPY_DIR)/bin/readCounter -c $(subst $( ),$(,),$(strip $(CHROMOSOMES))) HCC1395/exome/$*.bam > $@

###
## TITAN
###
SPLIT_CHR ?= true

.SECONDARY : $(foreach chr,$(CHROMOSOMES),titan/bcftools/vcf/HCC1395_exome_normal.var.het.$(chr).vcf) $(foreach chr,$(CHROMOSOMES),titan/bcftools/vcf/HCC1395_exome_tumour.var.het.$(chr).vcf)

.PHONY += titan/bcftools/tables/HCC1395_exome_tumour.var.het.table.tmp

titan.normal.hetpos : titan/bcftools/vcf/HCC1395_exome_normal.var.het.vcf
titan.tumour.hetpos : titan/bcftools/vcf/HCC1395_exome_tumour.var.het.vcf

# --- retrieve all heterozygous positions in the normal ---
# split by chromosome to speed up
# won't actually present this for the workshop
ifeq ($(SPLIT_CHR),true)

# add a QUAL >= 20 filter to improve the het position calls
titan/bcftools/vcf/HCC1395_exome_normal.var.het.%.vcf : HCC1395/exome/HCC1395_exome_normal.bam
	$(INIT) $(SAMTOOLS) view -uF 1024 $< $* | $(SAMTOOLS) mpileup -u -I -f $(REF_FASTA) - | $(BCFTOOLS) view -vcg - | $(JAVA) -jar ~/share/usr/snpEff/snpEff-4.0/SnpSift.jar filter "isHet(GEN[0]) & (QUAL >= 20)" > $@

# merge vcf files from different chromosome positions
titan/bcftools/vcf/HCC1395_exome_normal.var.het.vcf : $(foreach chr,$(CHROMOSOMES),titan/bcftools/vcf/HCC1395_exome_normal.var.het.$(chr).vcf)
	grep '#' $< > $@;\
	for file in $^; do\
		grep -v '#' $${file} >> $@;\
	done && rm $^

# build a bed file containing normal het positions to be 
# passed into the tumour
titan/bcftools/bed/HCC1395_exome_normal.var.het.%.bed : titan/bcftools/vcf/HCC1395_exome_normal.var.het.vcf
	$(INIT) awk '{OFS="\t"; if (!/^#/ && $$1 == "$*"){print $$1,$$2-1,$$2,$$4"/"$$5,"+"}}' $< > $@

titan/bcftools/vcf/HCC1395_exome_tumour.var.het.%.vcf : titan/bcftools/bed/HCC1395_exome_normal.var.het.%.bed HCC1395/exome/HCC1395_exome_tumour.bam
	$(SAMTOOLS) view -uF 1024 $(word 2,$^) $* | $(SAMTOOLS) mpileup -u -I -f $(REF_FASTA) -l $< - | $(BCFTOOLS) view -cg - > $@

# merge vcf files from different chromosome positions
titan/bcftools/vcf/HCC1395_exome_tumour.var.het.vcf : $(foreach chr,$(CHROMOSOMES),titan/bcftools/vcf/HCC1395_exome_tumour.var.het.$(chr).vcf)
	grep '#' $< > $@;\
	for file in $^; do\
		grep -v '#' $${file} >> $@;\
	done && rm $^

else
# this will be presented for the worshop
titan/bcftools/vcf/HCC1395_exome_normal.var.het.vcf : HCC1395/exome/HCC1395_exome_normal.bam
	$(INIT) $(SAMTOOLS) mpileup -u -I -f $(REF_FASTA) $< | $(BCFTOOLS) view -vcg - | $(JAVA) -jar ~/share/usr/snpEff/snpEff-4.0/SnpSift.jar filter "isHet(GEN[0]) & (QUAL >= 20)" > $@

# retrieve information for normal heterozygous positions in the tumour
titan/bcftools/bed/HCC1395_exome_normal.var.het.bed : titan/bcftools/vcf/HCC1395_exome_normal.var.het.vcf
	mkdir -p $(@D); \
	grep -v "^#" $< | awk '{OFS="\t"; print $$1,$$2-1,$$2,$$4"/"$$5,"+"}' > $@

titan/bcftools/vcf/HCC1395_exome_tumour.bcftools.vcf : titan/bcftools/vcf/het_pos.txt HCC1395/exome/HCC1395_exome_tumour.bam
	mkdir -p $(@D); \
	$(SAMTOOLS) view -uF 1024 $(word 2,$^) 17 | $(SAMTOOLS) mpileup -u -I -f $(REF_FASTA) -l $(word 1,$^) - | $(BCFTOOLS) view -cg - > $@
endif

# construct a tumour and normal table
# normal just for interest
titan/bcftools/tables/HCC1395_exome_normal.var.het.table.tmp : titan/bcftools/vcf/HCC1395_exome_normal.var.het.vcf
	mkdir -p $(@D); \
	java -jar $(SNPEFF_DIR)/SnpSift.jar extractFields $< $(VCF_FIELDS) > $@

titan/bcftools/tables/HCC1395_exome_normal.var.het.table.txt : titan/bcftools/tables/HCC1395_exome_normal.var.het.table.tmp
	perl getTitanAlleleReadCounts.pl $< > $@ && rm $<

# construct a table
titan/bcftools/tables/HCC1395_exome_tumour.var.het.table.tmp : titan/bcftools/vcf/HCC1395_exome_tumour.var.het.vcf
	mkdir -p $(@D); \
	java -jar $(SNPEFF_DIR)/SnpSift.jar extractFields $< $(VCF_FIELDS) > $@

titan/bcftools/tables/HCC1395_exome_tumour.var.het.table.txt : titan/bcftools/tables/HCC1395_exome_tumour.var.het.table.tmp
	perl getTitanAlleleReadCounts.pl $< > $@ && rm $<

# run titan
titan.timestamp :
	Rscript run_titan.R && touch titan.timestamp

# Merges the different plots per-chromosomes into a single pdf document for easy viewing
# Individual per-chromosome png files are removed once merged 
titan/results/HCC1395_exome_tumour.all.pdf : $(foreach chr,$(CHROMOSOMES),titan/results/HCC1395_exome_tumour.$(chr).png)
	convert $^ $@

