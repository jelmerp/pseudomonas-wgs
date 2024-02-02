#!/bin/bash

# ==============================================================================
#                  REMOVE CONTAMINATION FROM ASSEMBLIES
# ==============================================================================
# Run Kraken to detect contamination
asm_dir=asm_dir=results/spades/all_assemblies
kraken_db=/fs/project/PAS0471/jelmer/refdata/kraken/std-plus-fungi
for asm in "$asm_dir"/*fasta; do
    sbatch mcic-scripts/meta/kraken.sh -i "$asm" -o results/kraken --db-dir "$kraken_db"
done
#grep "Homo sapiens" results/kraken/*main.txt
#grep "Homo sapiens" results/kraken/*report.txt

# Remove Kraken-identified contaminants
conda activate /fs/project/PAS0471/jelmer/conda/krakentools-1.2
mkdir -p results/spades/decontam results/spades/contam
taxid=9606 # => Human
for krak_out in results/kraken/*main.txt; do 
    smp=$(basename "$krak_out" _main.txt)
    asm=results/spades/all_assemblies/"$smp".fasta
    fa_decontam=results/spades/decontam/"$smp".fasta
    fa_contam=results/spades/contam/"$smp".fasta
    echo "Sample: $smp"
    # Create a decontaminated assembly
    extract_kraken_reads.py -k "$krak_out" -s "$asm" -o "$fa_decontam" -t "$taxid" --exclude
    # Get a FASTA with only human-contaminant contigs
    extract_kraken_reads.py -k "$krak_out" -s "$asm" -o "$fa_contam" -t "$taxid"
    echo -e "--------------------\n\n"
done

# Check how many contigs were removed
grep -c ">" results/spades/contam/SpadesSM*fasta

# BLAST contaminant contigs
db=/fs/project/PAS0471/jelmer/refdata/blast/nt-db/nt
# Run BLAST search
for fa in results/spades/contam/*fasta; do
    ID=$(basename "$fa" .fasta)
    sbatch mcic-scripts/ncbi/blast.sh -i "$fa" -o results/kraken/blast/"$ID".out -l -d "$db"
done
# Add species identity to BLAST hits
for blast_out in results/kraken/blast/*.out; do
    sbatch mcic-scripts/ncbi/blast-process.sh -i "$blast_out" -o "$blast_out".proc -t 1
    sleep 30s
done

# Check whether any genes were removed (none!)
outdir=results/spades/contam_genecheck && mkdir -p "$outdir"
for asm in results/spades/contam/*fasta; do
    asm_id=$(basename "$asm" .fasta | sed 's/Spades//')
    gff=results/rast/"$asm_id"/"$asm_id".gff
    virgenes=results/pseudofinder/virgenes/biofilm_"$asm_id"_all.gff

    grep ">" "$asm" | sed -E 's/>(NODE_[0-9]+)_.*/\1/' > "$outdir"/"$asm_id"_contam.txt
    grep -f "$outdir"/"$asm_id"_contam.txt "$gff" > "$outdir"/"$asm_id"_contamgenes.gff
    n_removed=$(wc -l < "$outdir"/"$asm_id"_contamgenes.gff)
    n_genes=$(grep -v -f "$outdir"/"$asm_id"_contam.txt "$gff" | awk '$3 == "CDS"' | wc -l)
    
    echo -e "========\n$asm_id -- genes removed: $n_removed / genes remaining: $n_genes" 
    bedtools intersect -a "$outdir"/"$asm_id"_contamgenes.gff -b "$virgenes"
done
