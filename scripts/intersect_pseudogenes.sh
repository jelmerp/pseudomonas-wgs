#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=intersect_pseudogenes
#SBATCH --output=slurm-intersect_pseudogenes-%j.out

# Load software
module load miniconda3/4.12.0-py39
source activate /fs/ess/PAS0471/jelmer/conda/bedtools-env

# Strict bash settings, and locale for sorting
set -euo pipefail
LC_ALL=C ; LANG=C ; export LC_ALL LANG   # Needed for proper sorting for 'join'

# Constants
MIN_INTERSECT_PROP=0.1                  # Min. reciprocal proportion of intersection among genes for Bedtools
assembly_dir=results/spades/decontam
pseudofinder_dir=results/pseudofinder
virgenes_dir=results/virgenes
rast_dir=results/rast
outdir="$pseudofinder_dir"/virgenes

# Report
echo -e "\n# Starting script intersect_pseudogenes.sh"
date
echo
echo "Assembly dir:                         $assembly_dir"
echo "Pseudofinder results dir:             $pseudofinder_dir"
echo "Virulence genes dir:                  $virgenes_dir"
echo "RAST results dir:                     $rast_dir"
echo "Output dir:                           $outdir"
echo "=========================================================================="

# Create the output dir
echo -e "\n# Creating the output dirs..."
mkdir -pv "$outdir"/logs

# Run
echo -e "\n# Looping over the genomes and virulence gene types..."
for genome in "$assembly_dir"/*fasta; do

    genome_id=$(basename "$genome" .fasta | sed 's/Spades//')
    pseudo="$pseudofinder_dir"/"$genome_id"/"$genome_id"_pseudos.gff
    [[ ! -f "$pseudo" ]] && echo "Can't find file $pseudo" && exit 1
    echo
    echo "Input genome: $genome_id"

    # Loop over virulence gene types
    for genetype in biofilm motility secretion; do

        # Define input files
        genetype_list="$virgenes_dir"/"$genetype".txt
        gff_in="$rast_dir"/"$genome_id"/"$genome_id".gff
        [[ ! -f "$genetype_list" ]] && echo "Can't find file $genetype_list" && exit 1
        [[ ! -f "$gff_in" ]] && echo "Can't find file $gff_in" && exit 1

        # Define output files
        gff_focal="$outdir"/"$genetype"_"$genome_id"_all.gff
        genes_focal="$outdir"/"$genetype"_"$genome_id"_all.txt
        gff_pseudo="$outdir"/"$genetype"_"$genome_id"_pseudo.gff
        genes_pseudo="$outdir"/"$genetype"_"$genome_id"_pseudo.txt
        
        # Create a GFF file with an added 1st column with 'sanitized' gene names for merging
        grep -v "^#" "$gff_in" > tmp_gff_nohead
        sed -e 's/%2C/,/g' -e 's/%3B/;/g' -e 's/.*Name=//' -e 's/;.*//' tmp_gff_nohead > tmp_gff_list
        paste tmp_gff_list tmp_gff_nohead | sort -k 1b,1 > tmp_gff

        # Merge the GFF with the list of genes of interest
        join -1 1 -2 1 -t $'\t' tmp_gff <(sort -k 1b,1 $genetype_list) > tmp_gff_merged
        cut --complement -f 1 tmp_gff_merged > "$gff_focal"
        cut -f 1 tmp_gff_merged > "$genes_focal"

        # Run bedtools to intersect pseudo- and virulence genes
        bedtools intersect -f "$MIN_INTERSECT_PROP" -r \
            -a "$gff_focal" -b "$pseudo" -wo > "$gff_pseudo"
        
        # Get a simple list of pseudogenes and their counts
        sed -e 's/%2C/,/g' -e 's/%3B/;/g' -e 's/.*Name=//' -e 's/;.*//' "$gff_pseudo" |
            cut -f 1 |
            sed 's/ /__/g' |
            sort | uniq -c |
            awk -v OFS="\t" -v type="$genetype" -v id="$genome_id" \
                '{print type, id, $1, $2}' |
            sed 's/__/ /g' \
            > "$genes_pseudo"

        # Report
        n_list=$(wc -l < "$genetype_list")
        n_found=$(wc -l < "$gff_focal")
        n_found_uniq=$(sort "$genes_focal" | uniq | wc -l)
        n_pseudo=$(wc -l < "$gff_pseudo")
        n_pseudo_uniq=$(wc -l < "$genes_pseudo")
        echo "# input: $n_list   # found (uniq): $n_found ($n_found_uniq)   # pseudo (uniq): $n_pseudo ($n_pseudo_uniq)   Gene type: $genetype"

        # Remove temporary files
        rm tmp_gff tmp_gff_list tmp_gff_nohead tmp_gff_merged
    done
    
    echo "=================================================================="
done

# Report
echo -e "\nDone with script"
date
