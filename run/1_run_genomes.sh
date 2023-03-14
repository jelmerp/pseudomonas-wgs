#!/bin/bash

asm_dir=results/spades/all_assemblies

# Copy & rename the assembly FASTAs for easier access
mkdir -p results/spades/all_assemblies
for dir in results/spades/Spades*; do
    cp -v "$dir"/contigs.fasta "$asm_dir"/"$(basename "$dir")".fasta
done

# ==============================================================================
#                               ASSEMBLY QC
# ==============================================================================
# Run Quast to check genome quality -- for all assemblies at once
sbatch mcic-scripts/assembly/quast.sh -d results/spades/all_assemblies -o results/quast

# Run checkM to check genome quality -- for all assemblies at once
sbatch mcic-scripts/assembly/checkm.sh -i results/spades/all_assemblies -o results/checkm

# ==============================================================================
#                               KSNP3
# ==============================================================================
# Run kSNP3 to get an alignment of core SNPs -- own genomes only
mkdir -p results/ksnp3/own_asm
paste <(ls "$PWD"/results/spades/all_assemblies/*) \
    <(ls results/spades/all_assemblies/ | sed -E 's/Spades(.*).fasta/\1/') \
    >results/ksnp3/own_asm/assembly_list.txt
sbatch mcic-scripts/trees/ksnp3.sh -i results/ksnp3/own_asm/assembly_list.txt -o results/ksnp3/own_asm
Rscript mcic-scripts/trees/ggtree.R -i results/ksnp3/own_asm/tree.core.tre -o results/ksnp3/own_asm/tree.core.png


# ==============================================================================
#                  REMOVE CONTAMINATION FROM ASSEMBLIES
# ==============================================================================
# Run Kraken to detect contamination
kraken_db=/fs/project/PAS0471/jelmer/refdata/kraken/std-plus-fungi
for asm in results/spades/all_assemblies/*fasta; do
    sbatch mcic-scripts/metagenomics/kraken-run.sh -i "$asm" -o results/kraken -d "$kraken_db" -n
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


# ==============================================================================
#                  IDENTIFY CLOSEST REFERENCE GENOMES
# ==============================================================================
#> Outlier sample: SM225-1

# 1. Taxonomic classfication (LCA) with sourmash, using GTDB database
db=/fs/project/PAS0471/jelmer/refdata/sourmash/gtdb-rs207.genomic.k31.lca.json.gz
for asm in results/spades/decontam/*fasta; do
    sbatch mcic-scripts/misc/sourmash_classify.sh -i "$asm" -d "$db" -o results/sourmash/classify
done
#> All except SM225-1 are classified as Ps; and SM225-1 as P. moraviensis.
#> For SpadesSM225-1.fasta (all others ID'ed as syringae):
#> results/spades/decontam/SpadesSM225-1.fasta,found,d__Bacteria,p__Proteobacteria,c__Gammaproteobacteria,o__Pseudomonadales,f__Pseudomonadaceae,g__Pseudomonas_E,s__Pseudomonas_E moraviensis_A,

# 2. Use sourmash-search to identify the closest Pseudomonas genome to each assembly
# 2a. Download all Pseudomonas genomes (Downloaded 23,721 out of 23,723 genomes)
sbatch scripts/dl-all-pseudo.sh results/ref_all_pseudomonas

# 2b. Make sourmash database for all Pseudomonas genomes
sbatch mcic-scripts/misc/sourmash_db.sh -i results/ref_all_pseudomonas \
    -o results/sourmash/search/all_pseudo/db -d ref_all_pseudomonas

# 2c. Look for closest match of each assembly in the database with all Pseudomonas genomes
db=results/sourmash/ref_all_pseudomonas/db/ref_all_pseudomonas.sbt.zip
outdir=results/sourmash/search/all_pseudo/output
for asm in results/spades/decontam/*fasta; do
    sbatch mcic-scripts/misc/sourmash_search.sh -i "$asm" -d "$db" -o "$outdir"
done

# 3. Sourmash ANI analysis with own assemblies and initial set of reference genomes
# Note: this could be done with any set, e.g. also with closest Pseudomonas genomes found by sourmash in step 3
indir=results/sourmash/ani/fasta_in && mkdir -p "$indir"
cp results/spades/decontam/*fasta results/refgenomes/*fasta "$indir" # Copy all genomes to be used into an input dir
sbatch mcic-scripts/misc/sourmash_ani.sh -i "$indir" -o results/sourmash/ani

# 4. Check previously generated Kraken results, especially for SM225
# less results/kraken/SpadesSM225-1_report.txt
