#!/bin/bash
# Intersect Snippy output with focal virulence genes

# Load software
m activate /fs/ess/PAS0471/jelmer/conda/bedtools-env

# Variables
genome_id=SM51-19

# Input files
snippy_tab=results/snippy/"$genome_id"/snps.tab
snippy_vcf=results/snippy/"$genome_id"/snps.filt.vcf

biofilm_gff=results/pseudofinder/virgenes/biofilm_"$genome_id"_all.gff
motility_gff=results/pseudofinder/virgenes/motility_"$genome_id"_all.gff
secretion_gff=results/pseudofinder/virgenes/secretion_"$genome_id"_all.gff

# Output files
biofilm_vcf=results/snippy/virgenes/biofilm.vcf
motility_vcf=results/snippy/virgenes/motility.vcf
secretion_vcf=results/snippy/virgenes/secretion.vcf

# Intersect to get VCF files
bedtools intersect -a "$snippy_vcf" -b "$biofilm_gff" > "$biofilm_vcf"
bedtools intersect -a "$snippy_vcf" -b "$motility_gff" > "$motility_vcf"
bedtools intersect -a "$snippy_vcf" -b "$secretion_gff" > "$secretion_vcf"

# Merge with TAB file to get TAB output
join <(awk -v OFS="\t" '{print $1 "_" $2, $0}' "$snippy_tab" | sort)

