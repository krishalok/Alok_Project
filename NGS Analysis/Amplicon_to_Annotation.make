#----------- Project Description ----------------#
# AMPLICON TO ANNOTATION MAKESCRIPT
#
# This script will perform the following:
#	1. Align paired-end reads and map them to the human genome
#	2. Sort the Alignments and create a sorted BAM file
#	3. Index the BAM file
#	4. Call Variants using Varscan 2.3.5
#   5. Filter Variant calls (SNP & Indels)
#		a. Calls must have >10x Coverage
#		b. Call must have p-value <= 0.01
#	6. Annotates the Variants with Annovar for the following:
#		a. refGene Information
#		b. Segmental Duplications
#		c. Gwas catalog
#		d. whole-exome LJBSIFT
#		e. PolyPhen
#		f. PhyloP
#		g. LRT
#		h. MutationTaster
#		i. GERP++ scores
#		j. COSMIC database version 65
#		k. Alternate allele frequencies esp6500 European Americans
#		l. Alternate allele frequencies esp6500 European Americans including indels and chrY
#		m. 1000g alternative allele frequency data in 1000 Genomes Project for ALL, AMR (admixed american), EUR (european), ASN (asian), AFR (african) populations
#		n. dbSNP137
#		o. dbSNP137NonFlagged (dbSNP with ANNOVAR index files, after removing those flagged SNPs (SNPs < 1% minor allele frequency (MAF) (or unknown), mapping only once to reference assembly, flagged in dbSnp as "clinically associated"))
#		p. nci60 (NCI-60 human tumor cell line panel exome sequencing allele frequency data)
#
#	Final results will be 2 annotation tables
#		<Sample name>.variant.annov.hg19_multianno.csv
#		<Sample name>.indel.annov.hg19_multianno.csv




#----------- PROJECT FASTQ FOLDER ----------------#
# Example FASTQ names:
#	CRC_SAMPLE_1_R1_001.fastq
#   CRC_SAMPLE_1_R2_002.fastq
#	CRC_SAMPLE_2_R1_001.fastq
#   CRC_SAMPLE_2_R2_002.fastq
# Place the uncompressed FASTQ files in folder specific to their project
# Then all that you'll need to do is enter the directory where the files reside
# For example, if you put the FASTQ files in the FQ folder
# You would fill out the:
# FASTQ_FOLDER = FQ

FASTQ_FOLDER=


#----------- PROJECT BED FILE ----------------#
# You will also need to provide a sorted bed file for the amplicons
# Where each line consists of (Chr	Start_Position	Stop_Position	Amplicon_Name)
# Ex:	chr1    242013622       242013947       EXO1-P1E1-FO
#
# If the BED file has not been sorted, please use BEDTools to sort the bed file.
# Example Commmand:  sortBed -i my_unsorted_bed.bed > my_sorted_bed.bed
#
# An example of a bed filename would be
# BED_FILE = my_sorted_bed.bed

BED_FILE=


#----------- DON'T EDIT BELOW THIS LINE ----------------#
######################## Program PATHS and Variables ######################
#Reference Genome
GENOME = Genome/ucsc.hg19.fasta

#NGS QC Toolkit FASTQ Quality Trimming 
TRIM = /local/scratch/PACKAGES/src/ngsqctoolkit/TrimmingReads.pl

#Default Java Installation
JAVA = /usr/bin/java

#FastQC Sequence Metrics
FASTQC	= /Users/plott/Genomic/Packages/FastQC.app/Contents/Resources/Java/fastqc

#Default VarScan Program
VARSCAN = /local/scratch/PACKAGES/src/VarScan.v2.3.6.jar


#PICARD CollectMultipleMetrics
CMM = /Users/plott/Genomic/Packages/picard/CollectMultipleMetrics.jar

#PICARD CollectTargetedPcrMetrics
CTPM = /Users/plott/Genomic/Packages/picard/CollectTargetedPcrMetrics.jar

#Default Annovar Table Annotation Program
ANNOVAR = /local/scratch/PACKAGES/src/annovar

SHELL:=/bin/bash

##########################
#Variant Calling Variables
##########################
PVALUE=0.01

COVERAGE=10


#####################  Functions and Commands ###########
FASTQ_FILES_PREFIX = $(shell ls -1 $(FASTQ_FOLDER)/*.fastq | sed 's/_L001_R[12]_001.fastq//' | sed 's/^$(FASTQ_FOLDER)\///' )

GERMLINE_SUFFIX =  $(addsuffix .germline.annov.hg19_multianno.csv, $(FASTQ_FILES_PREFIX))
GERMLINE = $(addprefix VARIANTS/, $(GERMLINE_SUFFIX))

FQ_QC = $(addprefix FASTQ_QC/, $(FASTQ_FILES_PREFIX))

BAM_QC_SUFFIX = $(addsuffix .alignment_summary_metrics, $(FASTQ_FILES_PREFIX))
BAM_QC = $(addprefix MAP_QC/, $(BAM_QC_SUFFIX))

TARGET_QC_SUFFIX = $(addsuffix .target_metrics, $(FASTQ_FILES_PREFIX))
TARGET_QC = $(addprefix MAP_QC/, $(TARGET_QC_SUFFIX))

all: $(GERMLINE) $(FQ_QC) $(BAM_QC) $(TARGET_QC)

.SECONDARY:
##########################
# Need to finish testing and implementing Trimming
##########################
#Perform Quality Trimming and Filtering of sequences
#$(FASTQ_FOLDER)/%_L001_R1_001_trimmed.fastq: $(FASTQ_FOLDER)/%_L001_R1_001.fastq
#	perl $(TRIM) -i $(FASTQ_FOLDER)/$(*)_L001_R1_001.fastq -q 30 -n 200 -o $(FASTQ_FOLDER)/$(*)_L001_R1_001_trimmed.fastq

##########################
#Perform BWA Alignment for the Normal FASTQ files (Paired end)
##########################
BAM/%_sorted.bam: $(FASTQ_FOLDER)/%_L001_R1_001.fastq $(FASTQ_FOLDER)/%_L001_R2_001.fastq
	mkdir -p BAM
	@echo "Aligning Paired-end Amplicon reads to the genome using BWA-MEM"
#	bwa mem -M -t 1 $(GENOME) $(FASTQ_FOLDER)/$(*)_L001_R1_001.fastq $(FASTQ_FOLDER)/$(*)_L001_R2_001.fastq | samtools view -bS - | samtools sort - BAM/$(*)_sorted


##########################	
#Index the Normal BAM file
##########################
BAM/%_sorted.bai: BAM/%_sorted.bam
	@echo "Index Normal BAM file for viewing with IGV"
#	samtools index BAM/$(*)_sorted.bam BAM/$(*)_sorted.bai


##########################
#Variant Calling using Varscan 2.3.6 (Germline Indel/SNP)
##########################
VARIANTS/%.combined.variants: BAM/%_sorted.bam  BAM/%_sorted.bai
	mkdir -p VARIANTS
	@echo "Generate Pileup and Calling SNP/INDEL Variants with Varscan v2.3.5"
#	samtools mpileup -B -f $(GENOME) BAM/$(*)_sorted.bam | $(JAVA) -Xmx2g -jar $(VARSCAN) mpileup2cns --variants --p-value $(PVALUE) --min-coverage $(COVERAGE) --strand-filter 1 --min-avg-qual 20 > VARIANTS/$(*).combined.variants


##########################
#Annotate the High Confidence Calls using Annovar
##########################
VARIANTS/%.germline.annov.hg19_multianno.csv: VARIANTS/%.combined.variants
	echo "Target: $@"
	echo "Prerequisites: $?"
	@echo "Convert VarScan output to Annovar"
#	Varscan2Annovar.pl -type germline -varscan VARIANTS/$(*).combined.variants -output VARIANTS/$(*).germline.tmp
	@echo "Annotating variants and indels with Annovar\n"
#	perl $(ANNOVAR)/table_annovar.pl VARIANTS/$(*).germline.tmp $(ANNOVAR)/humandb/ --buildver hg19 -out VARIANTS/$(*).germline.annov -otherinfo -remove -protocol refGene,segdup,gwascatalog,ljb23_all,cosmic68wgs,esp6500_ea,esp6500si_ea,1000g2012apr_all,snp138,snp138NonFlagged,nci60 -operation  g,r,r,f,f,f,f,f,f,f,f -nastring NA -csvout

##########################
#Get Quality Metrics for FASTQ files
##########################
FASTQ_QC/%:  $(FASTQ_FOLDER)/%_L001_R1_001.fastq  $(FASTQ_FOLDER)/%_L001_R2_001.fastq
	@echo "Generating Fastq QA metrics for $(*)"
	mkdir -p FASTQ_QC
	mkdir -p FASTQ_QC/$(*)_R1
	mkdir -p FASTQ_QC/$(*)_R2
	$(FASTQC) --outdir=FASTQ_QC/$(*)_R1  $(FASTQ_FOLDER)/$(*)_L001_R1_001.fastq
	$(FASTQC) --outdir=FASTQ_QC/$(*)_R2  $(FASTQ_FOLDER)/$(*)_L001_R2_001.fastq

##########################
#Get Multiple Alignment Summary Statistics for BAM Files 
##########################
MAP_QC/%.alignment_summary_metrics: BAM/%_sorted.bam
	@echo "Generating Alignment Summary Metrics for $(*)"
	mkdir -p MAP_QC
	mkdir -p MAP_QC/$(*)
	$(JAVA) -jar $(CMM) INPUT=BAM/$(*)_sorted.bam OUTPUT=MAP_QC/$(*)/$(*)

##########################
#Get Target PCR Metrics for Amplicons
##########################
MAP_QC/%.target_metrics: BAM/%_sorted.bam MAP_QC/coord.txt
	@echo "Generating Alignment Target Metrics for $(*)"
	mkdir -p MAP_QC/$(*)
	samtools view -H BAM/$(*)_sorted.bam > MAP_QC/header.txt
	cat MAP_QC/header.txt MAP_QC/coord.txt > MAP_QC/Amplicon_interval.list
	$(JAVA) -jar $(CTPM) INPUT=BAM/$(*)_sorted.bam OUTPUT=MAP_QC/$(*)/$(*).pcr_metrics AMPLICON_INTERVALS=MAP_QC/Amplicon_interval.list TARGET_INTERVALS=MAP_QC/Amplicon_interval.list PER_TARGET_COVERAGE=MAP_QC/$(*)/$(*).target_metrics REFERENCE_SEQUENCE=$(GENOME)

MAP_QC/coord.txt: $(BED_FILE)
	mkdir -p MAP_QC
	awk '{print $$1"\t"$$2"\t"$$3"\t+\t"$$4}' $(BED_FILE) > MAP_QC/coord.txt

##########################
# Combine QC Metrics(Fastq, Alignment, and Amplicon) for Each Sample into a single worksheet
##########################
