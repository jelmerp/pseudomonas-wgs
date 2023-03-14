#!/bin/bash

# Software
module load miniconda3/4.12.0-py39
source activate /fs/ess/PAS0471/jelmer/conda/seqkit

# Strict bash settings
set -euo pipefail

# Parse command-line args
genome_id=$1
outfile=$2

# Input files
spades=results/spades/decontam/Spades"$genome_id".fasta
prokka=results/prokka/Spades"$genome_id"/Spades"$genome_id".fna

# Output files
outdir=$(dirname "$outfile")
spades_list="$outdir"/tmp/spades_"$genome_id".tsv
prokka_list="$outdir"/tmp/prokka_"$genome_id".tsv
spades_exclusive="$outdir"/tmp/spades_exclusive_"$genome_id".tsv
prokka_exclusive="$outdir"/tmp/prokka_exclusive_"$genome_id".tsv

# Report
echo
echo "Starting script contig_ids.sh"
date
echo
echo "Spades FASTA:                         $spades"
echo "Prokka FASTA:                         $prokka"
echo "Output file:                          $outfile"
echo
echo "# Listing the input files:"
ls -lh "$spades" "$prokka"
echo

# Make output dirs
mkdir -p "$outdir"/logs "$outdir"/tmp

# Sequence lengths
seqkit fx2tab --length --name --gc "$spades" --base-count C | \
    awk '{print $1 "\t" $2 "_" $3 "_" $4}' \
    > "$spades_list"

seqkit fx2tab --length --name --gc "$prokka" --base-count C | \
    awk '{print $1 "\t" $2 "_" $3 "_" $4}' |
    sed 's/gnl|Prokka|//' \
    > "$prokka_list"

# Join contig lists
join -t $'\t' -1 2 -2 2 "$spades_list" "$prokka_list" | cut -f 2,3 | sort > "$outfile"

# Also get non-matching lists
join -t $'\t' -1 2 -2 2 -v 1 "$spades_list" "$prokka_list" | cut -f 2,3 > "$spades_exclusive"
join -t $'\t' -1 2 -2 2 -v 2 "$spades_list" "$prokka_list" | cut -f 2,3 > "$prokka_exclusive"

# Report
echo "Statistics:"
echo "Nr contigs in Spades input:           $(wc -l < "$spades_list")"
echo "Nr contigs in Prokka input:           $(wc -l < "$prokka_list")"
echo "Nr contigs in outfile:                $(wc -l < "$outfile")"
echo "Nr unmatched contigs from Spades:     $(wc -l < "$spades_exclusive")"
echo "Nr unmatched contigs from Prokka:     $(wc -l < "$prokka_exclusive")"

if [[ -s "$spades_exclusive" ]]; then
    echo -e "\n# Showing the unmatched Spades contigs:"
    cat -n "$spades_exclusive"
fi
if [[ -s "$prokka_exclusive" ]]; then
    echo -e "\n# Showing the unmatched Spades contigs:"
    cat -n "$prokka_exclusive"
fi
echo -e "\n# Showing the first few lines of the output file:"
head "$outfile"

echo -e "\n# Listing the output file:"
ls -lh "$outfile"
echo -e "\nDone with script"
