#!/bin/bash

# ONLY RUNS ON BRIDGES SUPERCOMPUTER (USES BRIDGES-SPECIFIC COMMANDS) #

# Description #
# This script runs the following software:
    # Prokka (rapid prokaryotic genome annotation): de novo annotation transfer between strains (no reference)
        # i) finds genes with prodigal: long ORFs, transcription initiation sites, and poly-A sites
        # ii) annotates gene function: protein sequence similarity & database hierarchy (curated > UniProt > RefSeq > HMM)
    # RATT (rapid annotation transfer tool): reference-based annotation transfer between strains
        # i-ii) finds and annotates genes: conserved synteny (i.e. genes with same genetic location)
    # biodiff: identifies variants between NC and FS genomes
    # bcftools: identifies variant types and plots variant statistics
    # snpEff: assigns functional information to DNA variants (e.g. non-synonymous substitution)
        # crucial step in correlating DNA changes to changes in phenotype

# S. aureus and Sequencing Information #
# S. aureus: strain: ATCC BAA-39 / 2.8 million bases
# GenBank Accession #: CP033505 (control), CP033506 (FS treatment)
# Specs: PacBio (long read WGS) / SMRT Link v.5 pipeline + HGAP v.4 (assembly)


echo "module load prokka/1.11"
module load prokka/1.11

echo "running Prokka: annotating NC genome (NC.fasta)"
prokka --outdir pratt_NC --prefix pratt_NC --addgenes --genus Staphylococcus --species aureus --strain ATCC-BAA-39 --gram pos $SCRATCH/group/input_files/NC.fasta
  # bases: 2791218
  # CDS: 2733
  # gene: 2795
  # tRNA: 61 // tmRNA: 1 // sig_peptide: 174 // repeat_region: 2


echo "installing RATT in directory ratt-code/"
svn co "https://svn.code.sf.net/p/ratt/code/" ratt-code

echo "setting RATT_HOME and configuration file in RATT directory"
RATT_HOME=~/programs/ratt-code/; export RATT_HOME
RATT_CONFIG=$RATT_HOME/RATT.config_bac; export RATT_CONFIG

# RATT needs Mummer to generate the sequence comparison
# Nucmer: genome alignment component of MUMmer -> aligns query and reference genomes
echo "module load mummer"
module load mummer

# features from reference are transferred to query if in same syntenic block
echo "running RATT: using synteny to annotate FS genome from NC genome"
$HOME/programs/ratt-code/start.ratt.sh ../../embl $SCRATCH/group/input_files/FS.fasta FS_pratt Strain
# no synteny: 1.45% vs. 1.73%
# 5720/5767	Elements transferred (47 not)
# 2711/2733	Gene models/CDSs transferred (22 not)

# Biological Interpretation #
# regions with no synteny: less than 200 base pairs = probable deletion
# if the region is bigger: it might be a gap or a real insertion in the query


# biodiff: exact comparison of biological sequences --> input: fasta, output: vcf
# bcftools norm: left-align and normalize indels to reference sequence (using --fasta-ref)
echo "running biodiff: annotating FS variants (vs. NC) then normalizing indels to minimize sequencing error (e.g. AA vs. AAA)"
biodiff ../genomes/NC/NC.fasta ../genomes/FS/FS.fasta | bcftools norm --fasta-ref ../genomes/NC/NC.fasta - > ../data/biodiff/vcf/NC_FS_biodiff.vcf

echo "zipping .vcf file then running bcftools stats to calculate variant counts"
bgzip -c ../data/biodiff/vcf/NC_FS_biodiff.vcf > ../data/biodiff/vcf/NC_FS_biodiff.vcf.gz
tabix -p vcf ../data/biodiff/vcf/NC_FS_biodiff.vcf.gz
bcftools stats ../data/biodiff/vcf/NC_FS_biodiff.vcf.gz > ../data/biodiff/vcf/NC_FS_biodiff.vchk
# bcftools stats: Parses VCF or BCF and produces text file stats (can be plotted using plot-vcfstats)
# number of records:	5434
# number of SNPs:	565
# number of MNPs:	223
# number of indels:	3325
# number of others:	1320

echo "creating viualizations of numerous variant statistics"
plot-vcfstats --prefix ../data/biodiff/variants_bcf --main-title NC_FS_variants ../data/biodiff/vcf/NC_FS_biodiff.vchk


# snpEff: annotates and predicts the effects of genetic variants
# Are variants in a gene? In an exon? Do they change protein coding? Do they cause
# premature stop codons? SnpEff can help answer these questions.

# user-defined databases are required by snpEff to produce annotations
# chromosome name in VCF file must match database (snpeff.config: 'Chromosome')
echo "renaming CP033505.1 to Chromosome in VCF file to match database"
cat $SCRATCH/group/data/NC_FS_biodiff.vcf | sed "s/^CP033505.1/Chromosome/" > $SCRATCH/group/data/NC_FS_biodiff.vcf

echo "creating a custom S. auerus databased of NC genes to annotate FS variants"
cd /home/monorato/anaconda3/share/snpeff-4.3.1t-2
mkdir data/CP033505

cd data/CP033505
cp $SCRATCH/group/prokka_NC/pratt_NC/pratt_NC.gbf /home/monorato/anaconda3/share/snpeff-4.3.1t-2/data/CP033505/
mv pratt_NC.gbf genes.gbk

cd ../..
echo 'CP033505.genome : CP033505' >> snpEff.config
snpEff build -genbank CP033505
# Protein check: CP033505, OK: 2733, Not found: 0, Errors: 0, Error percentage: 0.0%

echo "running snpEff to annotate FS variants"
snpEff eff -csvStats FS_variants.csv CP033505 $SCRATCH/group/data/NC_FS_biodiff.vcf > NC_FS_biodiff.ann.vcf

echo ">>> Genome annotation and variant calling pipeline complete. <<<"
