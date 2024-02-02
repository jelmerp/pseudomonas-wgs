#!/bin/bash
# Runner script to submit the genome assemblies to NCBI

# Settings
minlen=200  # Min contig length

# Input dirs and files
asm_dir=results/spades/decontam
plasmid_results=results/mobsuite/contig_report_all.tsv
gff_dir=results/rast
virgene_dir=results/pseudofinder/virgenes

# Output dirs and files
subdir=results/ncbi_submission/fasta"$minlen" && mkdir -p "$subdir"
plasmid_contigs=results/ncbi_submission/plasmid_contigs.txt
lenfiltdir=results/spades/decontam_lenfilt"$minlen" && mkdir -p "$lenfiltdir"
genecheck_dir=results/ncbi_submission/genecheck && mkdir -p "$genecheck_dir"

# Create a file with plasmid contigs that has assembly ID and contig concatenated
tail -n+2 "$plasmid_results" | awk '{print $1 "_" $5}' > "$plasmid_contigs"

# Remove contigs shorten than the min. length
for asm in "$asm_dir"/*fasta; do
    asm_out="$lenfiltdir"/"$(basename "$asm")"
    seqkit seq --min-len "$minlen" "$asm" > "$asm_out" 2>/dev/null
    n_in=$(grep -c ">" "$asm"); n_out=$(grep -c ">" "$asm_out")
    n_removed=$(( n_in - n_out ))
    echo "$asm: Removed $n_removed out of $n_out contigs"
done

# Check how many genes were removed
for asm in "$asm_dir"/*fasta; do
    # Inputs
    asm_id=$(basename "$asm" .fasta)
    gff="$gff_dir"/"$asm_id"/"$asm_id".gff
    virgenes="$virgene_dir"/biofilm_"$asm_id"_all.gff

    # Outputs
    short_contigs="$genecheck_dir"/"$asm_id"_short.txt
    short_genes="$genecheck_dir"/"$asm_id"_shortgenes.gff
    numgenes="$genecheck_dir"/"$asm_id"_numgenes.txt

    # Get IDs of too-short contigs
    seqkit seq --max-len $(( minlen - 1 )) "$asm" 2>/dev/null | grep ">" |
        sed -E 's/>(NODE_[0-9]+)_.*/\1/' > "$short_contigs"
    
    # Get genes on too-short contigs
    grep -f "$short_contigs" "$gff" > "$short_genes"
    n_removed=$(wc -l < "$short_genes")
    n_genes=$(grep -v -f "$short_contigs" "$gff" | awk '$3 == "CDS"' | wc -l)
    echo -e "${asm_id}\t${n_genes}" > "$numgenes"
    # Report
    echo -e "========\n$asm_id -- genes removed: $n_removed / genes remaining: $n_genes"

    # Cross-check with virulence genes
    bedtools intersect -a "$short_genes" -b "$virgenes"
done
cat "$genecheck_dir"/*_numgenes.txt > "$genecheck_dir"/all_numgenes.tsv

# Create FASTA files with the proper extension & definition lines, e.g.:
# ncbi_submission/fasta/SM914-13.fsa: >contig001 [organism=Pseudomonas syringae pv. syringae] [isolate=SM914-13]
for asm in "$lenfiltdir"/*fasta; do
    # Define variables and output files
    outfile="$subdir"/$(basename "$asm" | sed -e 's/Spades//' -e 's/.fasta$/.fsa/')
    asm_id=$(basename "$outfile" .fsa)
    echo -e "===========\n$asm_id"

    # Copy the original assembly file to the output file, getting rid of any header
    # line content after a space, and adding organism and isolate info
    awk '{print $1}' "$asm" |
        sed "/>/s/$/ [organism=Pseudomonas syringa pv. syringae] [isolate=$asm_id]/" \
        > "$outfile"

    # Plasmids - first create lookup with assembly ID and contig concatenated,
    asm_id2=Spades$(basename "$asm" .fasta)
    awk -v asm_id="$asm_id2" '{print asm_id "_" $0, $0}' \
        <(grep ">" "$asm" | sed 's/>//') > tmp_"${asm_id}".contigs
    # then search for plasmids
    grep -f "$plasmid_contigs" tmp_"${asm_id}".contigs | while read -r _ contig; do
        (( counter++))
        echo "Plasmid contig detected: $contig (nr $counter)"
        sed -i "/>$contig/s/$/ [plasmid-name=unnamed$counter]/" "$outfile"
    done

    # Contig renaming
    sed -i -E -e 's/>NODE_([0-9])_[^[:space:]]+/>contig00\1/' \
        -e 's/>NODE_([0-9]{2})_[^[:space:]]+/>contig0\1/' \
        -e 's/>NODE_([0-9]{3})_[^[:space:]]+/>contig\1/' \
        "$outfile"

    # Finalize: remove temp files, print header lines to check
    rm tmp_"${asm_id}".contigs
    grep ">" "$outfile" | head -n2
    grep -E "plasmid" "$outfile"
done

# Remove excluded genome
rm -v "$subdir"/*SM225*

# Get general stats for Sochina
seqkit stats -T "$subdir"/*fsa > ncbi_submission/asm_stats.tsv

# See if the checkm results changed after removing small contigs
sbatch mcic-scripts/assembly/checkm.sh -i "$subdir" -o results/checkm200 -x fsa
