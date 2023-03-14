#!/bin/bash

# Variables
nr_fa=/fs/scratch/PAS0471/jelmer/dbs/nr/nr_database.fasta
nr_db=/fs/scratch/PAS0471/jelmer/dbs/nr/nr

asm_dir=results/spades/decontam

ref_pseudo=results/refgenomes/ASM14582v2.fasta
ref_snps_id="SM914-13"
#ref_snps="$asm_dir"/"$ref_snps_id".fasta
ref_snps=results/rast_fixed_gbk/"$ref_snps_id".gbk

# ==============================================================================
#                  IDENTIFY PLASMIDS WITH MOB-SUITE
# ==============================================================================
# Mob-suite - https://github.com/phac-nml/mob-suite
for asm in "$asm_dir"/*fasta; do
    outdir=results/mobsuite/$(basename "$asm" .fasta)
    sbatch mcic-scripts/bact/mob-suite.sh -i "$asm" -o "$outdir"
done

# Combine Mob-suite results for all assemblies
res_smr=results/mobsuite/mobtyper_results_all.tsv
contig_rep=results/mobsuite/contig_report_all.tsv

find results/mobsuite -name "mobtyper_results.txt" -print0 |
    xargs -0 awk 'FNR==1 && NR!=1{next} 1' > "$res_smr"

find results/mobsuite -name "contig_report.txt" | head -n1 | xargs head -n1 > "$contig_rep"
find results/mobsuite -name "contig_report.txt" -print0 | xargs -0 awk '$2 == "plasmid"' >> "$contig_rep"


# ==============================================================================
#                           IDENTIFY PSEUDOGENES
# ==============================================================================
# Make a Diamond db for the NCBI NR database for Pseudofinder
sbatch mcic-scripts/misc/diamond_db.sh -i $nr_fa -o "$nr_db"

# Run RAST #NOTE: I did this on my laptop, where I was able to install Rast-tk
for assembly in "$asm_dir"/*fasta; do
    outdir=results/rast/$(basename "$asm" .fasta)
    mcic-scripts/bact/rast_STUB.sh "$assembly" "$outdir" "Pseudomonas syringae"
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
    paste <(find data/fastq -name "*_R1_001.fastq.gz" -print0 | xargs -0 -n1 basename | sed -E 's/_S[0-9]+_L001_R1_001.fastq.gz//') - - |
    awk -v id="$ref_snps_id" '$1 != id' | sort > "$asm_tsv"

# Run snippy-multi
sbatch mcic-scripts/bact/snippy-multi.sh -i "$asm_tsv" -r "$ref_snps" -o results/snippy

# Explore Snippy results for a single assembly by means of example
restab=results/snippy/SM51-19/snps.tab
#restab=results/old_snippy/SM51-19/snps.tab

column -t "$restab" | less -S
tail -n+2 "$restab" | cut -f 3 | sort | uniq -c                         # Check mutation type - SNP vs indel etc
tail -n+2 "$restab" | cut -f 11  | cut -d " " -f 1 | sort | uniq -c     # Check SnpEff effects
tail -n+2 "$restab" | awk 'NF > 7 && $3 == "snp"' | head

results/pseudofinder/virgenes/biofilm_SM1031-4_all.gff