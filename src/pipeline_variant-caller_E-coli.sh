#!/bin/bash

# ONLY RUNS ON BRIDGES SUPERCOMPUTER (USES BRIDGES-SPECIFIC COMMANDS) #

# Description #
# This script runs the following software:
    # SRA toolkit: using SRA #, creates twp zipped FASTQ files for paired end reads
    # fastqc: quick check of raw reads, e.g. per base seq content, adapter content
    # trimmomatic: trims adapters, removes low quality seqs, keeps order of reads
    # fastqc: quick check of trimmed reads, e.g. per base quality, adapter content
    # bowtie2: FM index of reference, aligns reads end-to-end, creates SAM file
    # samtools: sorts SAM file, creates BAM file to align Illumina reads to an
      # E. coli reference genome, creates mpileup file
    # varscan2: calls short variants (SNPs, indels)
    # bcftools: left-aligns indels for accurate annotation & Sx variant calling
    # snpEff: annotates short variants
    # delly: calls structural variants such as large deletions, duplications, inversions

# E. coli and Sequencing Information #
# E. coli - Strain: K12 / Substrain: MG1655 / 216.5 million bases
# SRA Accession #: SRR8224895 (104.3 Mb file size)
# Specs: Illumina HiSeq 2500 (WGS) / Paired Ends / Adapters: Nextera XT


# removes pre-loaded modules and activates virtual environment for pipeline
# echo "module purge"
# module purge

# echo "activating anaconda 3"
# conda activate /home/monorato/anaconda3

echo "module load sra-toolkit/2.8.1-2"
module load sra-toolkit/2.8.1-2

# converts SRA file of E. coli reads into FASTQ file
# --split-files: creates separate files for paired end reads (i.e. 2+ reads/sample)
echo "fastq-dump --split-files --gzip SRR8224895"
fastq-dump --split-files --gzip SRR8224895


# cannot load both java7/8 @ same time -> java7 dependency for fastqc/0.11.3
echo "module load java7/jdk7u80"
module load java7/jdk7u80

# FastQC: initial check on raw sequence data, compatible with fastq files
# (input files), option to change Kmer size or use a contaminant filter,
# can be performed pre- and post-trim (if adapters still present after trim,
# they will appear as over-represented seqs in summary)
echo "module load fastqc/0.11.3"
module load fastqc/0.11.3

echo "running fastqc on SRR8224895_1.fastq.gz and SRR8224895_2.fastq.gz"
fastqc SRR8224895_1.fastq.gz SRR8224895_2.fastq.gz

echo "unzipping SRR8224895_1_fastqc.zip and SRR8224895_2_fastqc.zip"
unzip SRR8224895_1_fastqc.zip
unzip SRR8224895_2_fastqc.zip

echo "module load java/jdk8u73"
module load java/jdk8u73


# Trimmomatic: trims adapters (simple or palindromic), supports Nextera adapter
# library preps, removes low quality seqs, keeps order of reads (can work with
# all or subset of reads, e.g. only FWD AND REV reads survive trimming)
echo "module load trimmomatic/0.35"
module load trimmomatic/0.36

# calculates Phred scores (i.e. probability of incorrect base call)
# outputs fq.gz files (compressed) b/c FASTA files take up a lot of space
echo "creating .fa file for Nextera adapters; file to be used in trimmomatic"
echo $'>PrefixNX/1\nAGATGTGTATAAGAGACAG\n>PrefixNX/2\nAGATGTGTATAAGAGACAG\n>Trans1\nTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG\n>Trans1_rc\nCTGTCTCTTATACACATCTGACGCTGCCGACGA\n>Trans2\nGTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG\n>Trans2_rc\nCTGTCTCTTATACACATCTCCGAGCCCACGAGAC' > NexteraPE-PE.fa
echo "trimming reads with trimmomatic-0.36"
java -jar $TRIMMOMATIC_HOME/trimmomatic-0.36.jar PE -phred33 -trimlog matt_trim_log.fastq SRR8224895_1.fastq.gz SRR8224895_2.fastq.gz SRR8224895_forward_paired.fq.gz SRR8224895_forward_unpaired.fq.gz SRR8224895_reverse_paired.fq.gz SRR8224895_reverse_unpaired.fq.gz ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


# purges modules, mainly java8 so java7 can be reloaded properly
echo "module purge"
module purge

echo "module load java7/jdk7u80"
module load java7/jdk7u80

echo "module load fastqc/0.11.3"
module load fastqc/0.11.3

echo "running fastqc on trimmed reads, i.e. 4 fastq files"
fastqc SRR8224895_forward_paired.fq.gz SRR8224895_forward_unpaired.fq.gz SRR8224895_reverse_paired.fq.gz SRR8224895_reverse_unpaired.fq.gz

echo "unzipping 4 fastqc.zip files"
unzip SRR8224895_forward_paired_fastqc.zip
unzip SRR8224895_forward_unpaired_fastqc.zip
unzip SRR8224895_reverse_paired_fastqc.zip
unzip SRR8224895_reverse_unpaired_fastqc.zip


# Bowtie2: optimized for short reads, ex. Illumina reads of ~100 bps
# trimmed reads here are 90-100 bps long (PacBio: no upper limit on read length)
echo "module load bowtie2"
module load bowtie2

echo "FM indexing E. coli K12 reference genome"
bowtie2-build ecolik12.fasta ecolik12

echo "aligning paired and unpaired reads to reference genome based on FM index (end-to-end, not local)"
bowtie2 -p 8 -x ecolik12 -1 SRR8224895_forward_paired.fq.gz -2 SRR8224895_reverse_paired.fq.gz -U SRR8224895_forward_unpaired.fq.gz,SRR8224895_reverse_unpaired.fq.gz -S SRR8224895.sam
# 65.98% overall alignment rate

echo "module load samtools/1.9"
module load samtools/1.9

echo "sorting alignments of SAM file by leftmost coordinates and creating a BAM file"
samtools sort --reference ecolik12.fasta SRR8224895.sam > SRR8224895.bam

echo "creating index of BAM file (for fast access, like FM index)"
samtools index SRR8224895.bam
# SRR8224895.bam = 125,337,680 bytes vs. SRR8224895.bam.bai = 13,472 bytes
# 125,338 KB vs. 13 KB -> BAM file 10,000x bigger than BAM.bai file

# samtools stats SRR8224895.bam
# SN	raw total sequences:	2077660
# SN	reads mapped:	1370857
# SN	reads unmapped:	706803
# SN	mismatches:	2457394	# from NM fields
# SN	error rate:	1.904334e-02	# mismatches / bases mapped (cigar)
# SN	average length:	94
# RL read length range: 36-100 -> trimmomatic set to 36
# SN	average quality:	34.8

echo "generating mpileup, a summary file containing base calls of aligned reads to a reference sequence"
samtools mpileup -f ecolik12.fasta SRR8224895.bam --min-MQ 15 > SRR8224895.mpileup
# Why MapQV > 15? for structural variants (reads often partially mapped)
# mpileup file is largest of all: 320,285,791 (~300 MB)


echo "module load java7/jdk7u80"
module load java7/jdk7u80

echo "module load varscan/2.4.2"
module load varscan/2.4.2

echo "running mpileup2cns to call variants for SNPs and indels only"
java -jar $VARSCAN_HOME/VarScan.jar mpileup2cns SRR8224895.mpileup --variants --min-reads2 4 --p-value 99e-04 --output-vcf 1 > SRR8224895.vcf
# 57715 variant positions reported: 57313 SNP, 402 indel

echo "zipping VCF file to reduce size"
bgzip -c SRR8224895.vcf > SRR8224895.vcf.gz

echo "indexing zipped VCF file similar to bam.bai file"
tabix -p vcf SRR8224895.vcf.gz
# VCF file (609,682,538 KB) -> vcf.gz (9,411,593 KB) -> vcf.gz.tbi (2,674 KB)
# Total Size Reduction: ~66x * ~3,300x = ~225,000x

echo "module load bcftools/1.3.1"
module load bcftools/1.3.1

echo "left-aligning indels to improve annotation accuracy (and possibly SV calling)"
bcftools norm -cw -f ecolik12.fasta SRR8224895.vcf.gz


# snpEff: annotates variants and determines impact on known genes, requires an
# annotated E. coli database and a VCF file (both with matching CHROM names)

echo "downloading E. coli databases for all strains in snpEff"
snpEff download Escherichia_coli_str_k_12_substr_mg1655

echo "updating chromosome name from NC_000913.3 to Chromosome in SRR8224895.vcf (to match E. coli database)"
cat SRR8224895.vcf | sed "s/^NC_000913.3/Chromosome/" > SRR8224895_updated.vcf

echo "annotating variants with snpEff relative to reference strain"
snpEff eff -csvStats snpEff_summary.csv Escherichia_coli_str_k_12_substr_mg1655 SRR8224895_updated.ann.vcf
# Variants by type
# Type , Count , Percent
# DEL , 196 , 0.3396%
# INS , 206 , 0.356926%
# SNP , 57313 , 99.303474%

# Effects by functional class
# Type , Count , Percent
# MISSENSE , 7526 , 14.838036%
# NONSENSE , 50 , 0.098578%
# SILENT , 43145 , 85.063386%


# Delly: uses paired-ends, split-reads and read-depth to identify genomic
# rearrangements; requires sorted, indexed and duplicate marked bam file and an
# indexed reference genome
delly call -g ecolik12.fasta SRR8224895.bam
bcftools view sv.bcf > SRR8224895_sv.vcf
bgzip -c SRR8224895_sv.vcf > SRR8224895_sv.vcf.gz
tabix -p vcf SRR8224895_sv.vcf.gz

echo ">>> Genome assembly and variant calling pipeline complete. <<<"
