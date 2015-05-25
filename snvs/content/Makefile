include ~/share/modules/Makefile.inc 

LOGDIR ?= log

STRELKA_VERSION ?= 1.0.14
.PHONY : strelka museq

MUTATIONSEQ_DIR = /share/lustre/fong/usr/museq/4.2.1
STRELKA_DIR = /share/lustre/fong/usr/strelka/1.0.14
SNPEFF_DIR = ~/share/usr/snpEff/snpEff-4.0
REF_FASTA = ref_data/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa

CHROMOSOMES = 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 MT X Y

museq : results/museq/HCC1395.museq.vcf
museq.full : $(foreach chr,$(CHROMOSOMES),results/museq/HCC1395.full.museq.$(chr).vcf)
museq.table : results/museq/HCC1395.museq.txt results/museq/HCC1395.full.museq.txt

.SECONDARY : $(foreach chr,$(CHROMOSOMES),results/museq/HCC1395.full.museq.$(chr).vcf)


# convert markdown into a mediawiki format to upload onto the CBW wiki
snv-lab.wiki : 
	pandoc -t mediawiki -o snv-lab.wiki snv-lab.md

snv-data-analysis-package.tar.gz : 
	tar -cvzf snv-data-analysis-package.tar.gz --exclude=Makefile --exclude=*.swp snv-data-analysis-package

####
## Mutationseq
####

## Chromosome 17 (7-8 Mb)
results/museq/HCC1395.museq.vcf : 
	mkdir -p $(@D); \
	~/share/usr/anaconda/2.1.0/bin/python $(MUTATIONSEQ_DIR)/classify.py \
	tumour:bam/HCC1395_exome_tumour.17.7MB-8MB.bam \
	normal:bam/HCC1395_exome_normal.17.7MB-8MB.bam \
	reference:$(REF_FASTA) \
	model:${MUTATIONSEQ_DIR}/model_v4.1.1.npz \
	-i 17:7000000-8000000 \
	-c ${MUTATIONSEQ_DIR}/metadata.config -q 1 -o $@ \
	-l $(@D)/$(@F).log

## Whole genome
results/museq/HCC1395.full.museq.%.vcf : 
	$(call LSCRIPT_MEM_DRMAA,2G,4G,"~/share/usr/anaconda/2.1.0/bin/python $(MUTATIONSEQ_DIR)/classify.py -b 2G tumour:HCC1395/exome/HCC1395_exome_tumour.bam normal:HCC1395/exome/HCC1395_exome_normal.bam reference:$(REF_FASTA) model:$(MUTATIONSEQ_DIR)/model_v4.1.1.npz -i $* -c $(MUTATIONSEQ_DIR)/metadata.config -q 1 -o $@ -l $(@D)/$(@F).log --all")

# merge all the position_file chunks into one file
results/museq/HCC1395.full.museq.vcf : $(foreach chr,$(CHROMOSOMES),results/museq/HCC1395.full.museq.$(chr).vcf)
	grep '#' $< > $@;\
	for file in $^; do\
		grep -v '#' $${file} >> $@;\
	done

%.museq.txt : %.museq.vcf
	java -jar $(SNPEFF_DIR)/SnpSift.jar extractFields -e "." $< CHROM POS ID REF ALT QUAL FILTER PR TR TA NR NA TC NI ND > $@

####
## Strelka
####

strelka : results/strelka/task.complete results/strelka.full/task.complete
strelka.table : results/strelka/results/passed.somatic.snvs.txt results/strelka.full/results/passed.somatic.snvs.txt

## Chromosome 17 (7-8 Mb)
results/strelka/Makefile :
	perl $(STRELKA_DIR)/bin/configureStrelkaWorkflow.pl \
    --tumor HCC1395/exome/HCC1395_exome_tumour.17.7MB-8MB.bam \
    --normal HCC1395/exome/HCC1395_exome_normal.17.7MB-8MB.bam \
    --ref refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
    --config config/strelka_config_bwa_exome.ini \
    --output-dir $(@D)

results/strelka/task.complete : results/strelka/Makefile
	make -C $(@D) -j 4

results/strelka.full/Makefile :
	perl $(STRELKA_DIR)/bin/configureStrelkaWorkflow.pl \
    --tumor HCC1395/exome/HCC1395_exome_tumour.bam \
    --normal HCC1395/exome/HCC1395_exome_normal.bam \
    --ref refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
    --config config/strelka_config_bwa_exome.ini \
    --output-dir $(@D)

results/strelka.full/task.complete : results/strelka.full/Makefile
	make -C $(@D) -j 4

%.somatic.snvs.txt : %.somatic.snvs.vcf
	java -jar $(SNPEFF_DIR)/SnpSift.jar extractFields -e "." $< CHROM POS ID REF ALT QUAL QSS FILTER TQSS NT QSS_NT TQSS_NT SGT SOMATIC GEN[0].DP GEN[1].DP GEN[0].FDP GEN[1].FDP GEN[0].SDP GEN[1].SDP GEN[0].SUBDP GEN[1].SUBDP GEN[0].AU GEN[1].AU GEN[0].CU GEN[1].CU GEN[0].GU GEN[1].GU GEN[0].TU GEN[1].TU > $@

%.html : %.Rmd
	Rscript --slave -e "rmarkdown::render('$<')"
