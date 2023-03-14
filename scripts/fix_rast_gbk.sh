#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=fix_rast_gbk
#SBATCH --output=slurm-fix_rast_gbk-%j.out

# Strict base settings
set -euo pipefail

# Parse command-line args
infile=$1                   #infile=results/rast/SpadesSM51-19/SpadesSM51-19.gbk
outdir=$2

# Get the assembly ID and define the outdir
outfile="$outdir"/$(basename "$infile")

# Report
echo "# Starting script rast.sh"
echo "This script will Fix RAST's Genbank file to work with Pseudofinder"
echo "We need a slightly different locus header, 'locus_tag' instead of 'db_xref',"
echo "and most importantly, a 'gene' feature, which corresponds to each 'CDS' feature,"
echo "so this script is copying those coordinates and it's 'locus_tag'."
date
echo
echo "Input Genbank file:   $infile"
echo "Outdir:               $outdir"
echo "Output Genbank file:  $outfile"

# Create the output dir
echo -e "\n# Creating the output dirs..."
mkdir -pv "$outdir"

# Prepare for annotation by creating a Rast-specific file
echo -e "\n# Fixing the GFF file..."
grep -v "EC_number=" "$infile" |
    sed -e 's/db_xref/locus_tag/' \
        -e 's/dna     linear/DNA     linear/' |
    sed 's/linear   UNK/linear       01-FEB-2023/' |
    awk '/     CDS             / {x=$0; x2=$0; sub("CDS ", "gene", x); print x; getline; y = $0; print y; print x2} {print}' \
    > "$outfile"

# Report
echo -e "\n# Listing the output file:"
ls -lh "$(realpath "$outfile")"
echo
echo "Done with script"
date
