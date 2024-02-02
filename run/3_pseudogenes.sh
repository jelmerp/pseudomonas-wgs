#!/bin/bash
# Runner script to identify pseudogenes and SNPs

# Variables
nr_fa=/fs/scratch/PAS0471/jelmer/dbs/nr/nr_database.fasta
nr_db=/fs/scratch/PAS0471/jelmer/dbs/nr/nr
asm_dir=results/spades/decontam
ref_pseudo=results/refgenomes/ASM14582v2.fasta
ref_snps_id="SM914-13"
ref_snps=results/rast_fixed_gbk/"$ref_snps_id".gbk

# ==============================================================================
#                           IDENTIFY PSEUDOGENES
# ==============================================================================
# Make a Diamond db for the NCBI NR database for Pseudofinder
sbatch mcic-scripts/misc/diamond_db.sh -i $nr_fa -o "$nr_db"

# Run RAST #NOTE: I did this on my laptop, where I was able to install Rast-tk
for assembly in "$asm_dir"/*fasta; do
    outdir=results/rast/$(basename "$assembly" .fasta)
    bash mcic-scripts/bact/rast_STUB.sh "$assembly" "$outdir" "Pseudomonas syringae"
done

# Fix RAST's Genbank file to work with Pseudofinder
for genbank in $(find results/rast -name "*gbk"); do
    sbatch scripts/fix_rast_gbk.sh "$genbank" results/rast_fixed_gbk
done

# Run Pseudofinder
for asm_gbk in results/rast_fixed_gbk/*gbk; do
    outdir=results/pseudofinder/"$(basename "$asm_gbk" .gbk)"
    sbatch mcic-scripts/bact/pseudofinder.sh \
        -i "$asm_gbk" --db "$nr_db" --ref "$ref_pseudo" -o "$outdir"
done
#> Pseudogene GFFs are in: results/pseudofinder/SM51-19/SM51-19_fixed_pseudos.gff

# Pseudogenes among virulence genes
sbatch scripts/intersect_pseudogenes.sh

# Summarize
cat results/pseudofinder/virgenes/*_pseudo.txt > results/pseudofinder/virgenes/all_pseudogenes.tsv

# ==============================================================================
#                           IDENTIFY SNPS
# ==============================================================================
asm_tsv=results/snippy/input_assemblies.tsv && mkdir -p results/snippy

# Make an input TSV for Snippy-multi
find data/fastq -name "*fastq.gz" | sort | \
    paste <(find data/fastq -name "*_R1_001.fastq.gz" -print0 | xargs -0 -n1 basename | sed -E 's/_S[0-9]+_L001_R1_001.fastq.gz//' | sort) - - |
    awk -v id="$ref_snps_id" '$1 != id' | sort > "$asm_tsv"

# Run snippy-multi
sbatch mcic-scripts/bact/snippy-multi.sh -i "$asm_tsv" -r "$ref_snps" -o results/snippy

# Intersect Snippy SNPs with focal virulence genes
bash scripts/intersect_snps.sh

# Summarize results across genomes
Rscript scripts/process_snps.R

# Get GFFs of focal genes for the reference genome
#biofilm_gff_in=results/pseudofinder/virgenes/biofilm_"$ref_snps_id"_all.gff
#biofilm_gff_ed=results/virgenes/biofilm_"$ref_snps_id".gff
#sed -e 's/%2C/,/g' -e 's/%3B/;/g' "$biofilm_gff_in" > "$biofilm_gff_ed"
