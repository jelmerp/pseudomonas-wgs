#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=intersect_snps
#SBATCH --output=slurm-intersect_snps-%j.out

# Intersect Snippy output with focal virulence genes

# ==============================================================================
#                                   SETUP
# ==============================================================================
# Load software
module load miniconda3/4.12.0-py39
source activate /fs/ess/PAS0471/jelmer/conda/bedtools-env

# Strict bash settings
#set -euo pipefail
LC_ALL=C ; LANG=C ; export LC_ALL LANG

# Variables
virgenes_dir=results/pseudofinder/virgenes
refgenome_id=SM914-13
snippy_dir=results/snippy
outdir_part=results/virgenes_snps

# constants
VARIANTS_REMOVE_REGEX="intragenic_variant|intergenic_region|synonymous_variant|intron_variant"

# Input files
biofilm_withpseudo_gff="$virgenes_dir"/biofilm_"$refgenome_id"_all.gff
motility_withpseudo_gff="$virgenes_dir"/motility_"$refgenome_id"_all.gff
secretion_withpseudo_gff="$virgenes_dir"/secretion_"$refgenome_id"_all.gff

biofilm_pseudo_gff="$virgenes_dir"/biofilm_"$refgenome_id"_pseudo.gff
motility_pseudo_gff="$virgenes_dir"/motility_"$refgenome_id"_pseudo.gff
secretion_pseudo_gff="$virgenes_dir"/secretion_"$refgenome_id"_pseudo.gff

# Output files - GFFs with pseudogenes excluded
biofilm_gff="$outdir_part"/biofilm_"$refgenome_id"_nopseudo.gff
motility_gff="$outdir_part"/motility_"$refgenome_id"_nopseudo.gff
secretion_gff="$outdir_part"/secretion_"$refgenome_id"_nopseudo.gff

# Report
echo -e "\n# Starting script intersect_snps.sh"
date
echo
echo "Snippy results dir:                   $snippy_dir"
echo "Virulence genes dir:                  $virgenes_dir"
echo "Base output dir:                      $outdir_part"
echo "=========================================================================="
echo

# ==============================================================================
#                                   RUN
# ==============================================================================
# Remove pseudogenes
bedtools intersect -v -a "$biofilm_withpseudo_gff" -b "$biofilm_pseudo_gff" > "$biofilm_gff"
bedtools intersect -v -a "$motility_withpseudo_gff" -b "$motility_pseudo_gff" > "$motility_gff"
bedtools intersect -v -a "$secretion_withpseudo_gff" -b "$secretion_pseudo_gff" > "$secretion_gff"

ngenes_biofilm1=$(grep -cv "^#" "$biofilm_withpseudo_gff")
ngenes_biofilm2=$(grep -cv "^#" "$biofilm_gff")
echo "Nr biofilm genes before/after pseudogene removal: $ngenes_biofilm1 / $ngenes_biofilm2"
ngenes_motility1=$(grep -cv "^#" "$motility_withpseudo_gff")
ngenes_motility2=$(grep -cv "^#" "$motility_gff")
echo "Nr motility genes before/after pseudogene removal: $ngenes_motility1 / $ngenes_motility2"
ngenes_secretion1=$(grep -cv "^#" "$secretion_withpseudo_gff")
ngenes_secretion2=$(grep -cv "^#" "$secretion_gff")
echo "Nr secretion genes before/after pseudogene removal: $ngenes_secretion1 / $ngenes_secretion2"

# Loop over the genomes
echo -e "\n# Looping over the genomes..."
for genome_dir in "$snippy_dir"/SM*; do
    
    genome_id=$(basename "$genome_dir")
    echo -e "\nFocal genome: $genome_id"

    # Input files
    snippy_tab="$genome_dir"/snps.tab
    snippy_vcf="$genome_dir"/snps.filt.vcf

    # Output files
    outdir="$outdir_part"/"$genome_id"
    snippy_tab_ed="$outdir"/snps_nonsyn.tab

    # Output files - VCFs with snippy SNPs only for focal genes
    biofilm_vcf="$outdir"/biofilm_all.vcf
    motility_vcf="$outdir"/motility_all.vcf
    secretion_vcf="$outdir"/secretion_all.vcf
    
    # Output files - tab files with nonsyn mutations only
    biofilm_tab="$outdir"/biofilm_nonsyn.tab
    motility_tab="$outdir"/motility_nonsyn.tab
    secretion_tab="$outdir"/secretion_nonsyn.tab

    # Output files - summary mutation-count-per-genes files
    biofilm_countpergene="$outdir"/biofilm_nonsyn_countspergene.txt
    motility_countpergene="$outdir"/motility_nonsyn_countspergene.txt
    secretion_countpergene="$outdir"/secretion_nonsyn_countspergene.txt

    # Create the output dir
    mkdir -p "$outdir"

    # Intersect to get VCF files
    bedtools intersect -a "$snippy_vcf" -b "$biofilm_gff" > "$biofilm_vcf"
    bedtools intersect -a "$snippy_vcf" -b "$motility_gff" > "$motility_vcf"
    bedtools intersect -a "$snippy_vcf" -b "$secretion_gff" > "$secretion_vcf"

    n_total=$(grep -v "^#" "$snippy_vcf" | tail -n +2 | wc -l)
    n_biofilm=$(grep -v "^#" "$biofilm_vcf" | tail -n +2 | wc -l)
    n_motility=$(grep -v "^#" "$motility_vcf" | tail -n +2 | wc -l)
    n_secretion=$(grep -v "^#" "$secretion_vcf" | tail -n +2 | wc -l)
    echo "Nr of total SNPs for this genome: $n_total"
    echo "Nr of total SNPs in biofilm/motility/secretion: $n_biofilm / $n_motility / $n_secretion"

    # Remove synonymous (etc) SNPs
    tail -n+2 "$snippy_tab" | grep -Ev "$VARIANTS_REMOVE_REGEX" > "$snippy_tab_ed"
    n_nonsyn=$(wc -l < "$snippy_tab_ed")
    echo "Nr of nonsyn SNPs for this genome: $n_nonsyn"

    # Merge with TAB file to get TAB output - biofilm
    head -n 1 "$snippy_tab" > "$biofilm_tab"
    join -1 1 -2 1 -t $'\t' \
        <(awk -v OFS="\t" '{print $1 "_" $2, $0}' "$snippy_tab_ed" | sort -k1b,1) \
        <(awk -v OFS="\t" '{print $1 "_" $2, $0}' "$biofilm_vcf" | sort -k1b,1 | cut -f 1) |
        cut --complement -f 1 >> "$biofilm_tab"

    ngenes_biofilm=$(tail -n +2 "$biofilm_tab" | cut -f 14 | sort | uniq | wc -l)
    tail -n +2 "$biofilm_tab" | cut -f 14 | sort | uniq -c > "$biofilm_countpergene"

    # Merge with TAB file to get TAB output - motility
    head -n 1 "$snippy_tab" > "$motility_tab"
    join -1 1 -2 1 -t $'\t' \
        <(awk -v OFS="\t" '{print $1 "_" $2, $0}' "$snippy_tab_ed" | sort -k1b,1) \
        <(awk -v OFS="\t" '{print $1 "_" $2, $0}' "$motility_vcf" | sort -k1b,1 | cut -f 1) |
        cut --complement -f 1 >> "$motility_tab"

    ngenes_motility=$(tail -n +2 "$motility_tab" | cut -f 14 | sort | uniq | wc -l)
    tail -n +2 "$motility_tab" | cut -f 14 | sort | uniq -c > "$motility_countpergene"

    # Merge with TAB file to get TAB output - serectiom
    head -n 1 "$snippy_tab" > "$secretion_tab"
    join -1 1 -2 1 -t $'\t' \
        <(awk -v OFS="\t" '{print $1 "_" $2, $0}' "$snippy_tab_ed" | sort -k1b,1) \
        <(awk -v OFS="\t" '{print $1 "_" $2, $0}' "$secretion_vcf" | sort -k1b,1 | cut -f 1) |
        cut --complement -f 1 >> "$secretion_tab"

    ngenes_secretion=$(tail -n +2 "$secretion_tab" | cut -f 14 | sort | uniq | wc -l)
    tail -n +2 "$secretion_tab" | cut -f 14 | sort | uniq -c > "$secretion_countpergene"

    # Count and report
    n_biofilm_nonsyn=$(tail -n+2 "$biofilm_tab" | wc -l)
    n_motility_nonsyn=$(tail -n+2 "$motility_tab" | wc -l)
    n_secretion_nonsyn=$(tail -n+2 "$secretion_tab" | wc -l)
    echo "Nr of nonsyn. SNPs in biofilm/motility/secretion: $n_biofilm_nonsyn / $n_motility_nonsyn / $n_secretion_nonsyn"
    echo "Nr of genes with nonsyn. SNPs in biofilm/motility/secretion: $ngenes_biofilm / $ngenes_motility / $ngenes_secretion"
done

# Report
echo -e "\nDone with script"
date

# Explore Snippy results for a single assembly by means of example
#restab=results/snippy/SM51-19/snps.tab
#column -t "$restab" | less -S
#tail -n+2 "$restab" | cut -f 3 | sort | uniq -c                         # Check mutation type - SNP vs indel etc
#tail -n+2 "$restab" | cut -f 11  | cut -d " " -f 1 | sort | uniq -c     # Check SnpEff effects
#tail -n+2 "$restab" | awk 'NF > 7 && $3 == "snp"' | head
