#!/bin/bash
# Runner script with main analyses using the genome assemblies

# Input dirs and settings
ref_genomes=(ASM1669475v2 ASM780v1 ASM1224v1 ASM1220v1 ASM14582v2 ASM14584v2 ASM2327812v1 ASM15699v2 ASM276365v1 ASM45244v3 ASM1860349v1 ASM14594v2 ASM290581v2 ASM2327794v1 CFBP4215 CFBP2118 ASM98848v1 ASM98839v1 ASM51227v1)

# Output dirs
asm_dir=results/spades/all_assemblies && mkdir -p "$asm_dir"   # Dir with Spades assemblies
ref_dir=data/ref && mkdir -p "$ref_dir"                        # Dir with reference (downloaded) genomes
refasm_list=metadata/refgenomes_GCA.txt                        # List of reference (downloaded) genomes

# ==============================================================================
#                                   SETUP
# ==============================================================================
# Copy & rename the assembly FASTAs for easier access
for dir in results/spades/Spades*; do
    cp -v "$dir"/contigs.fasta "$asm_dir"/"$(basename "$dir")".fasta
done

# ==============================================================================
#                   DOWNLOAD REFERENCE GENOMES FROM NCBI
# ==============================================================================
# Convert ASM ids to GCA ids
>"$refasm_list"
for id in "${ref_genomes[@]}"; do
    echo "$id"
    GCA=$(esearch -db assembly -query "$id" | esummary | grep "Genbank.*GCA" |
            sed -e 's/<Genbank>//' -e 's@</Genbank>@@' | awk '{print $1}')
    echo "$GCA" | tee -a "$refasm_list"
    echo "----"
done

# 2023-10-14: Download outgroup used in https://www.frontiersin.org/articles/10.3389/fmicb.2017.02422/full
# P. aeruginosa DSM 50071T = NZ_CP012001.1 = https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001045685.1/
sbatch mcic-scripts/dbs/dl-genomes.sh --accession GCA_001045685.1 -o data/ref --include genome
echo "GCA_001045685.1" >> "$refasm_list"

# Download the genomes
sbatch mcic-scripts/dbs/dl-genomes.sh --accession_file "$refasm_list" -o data/ref --include genome

# ==============================================================================
#                               ASSEMBLY QC
# ==============================================================================
# Run Quast to check genome quality -- for all assemblies at once
sbatch mcic-scripts/assembly/quast.sh -d results/spades/all_assemblies -o results/quast

# Run checkM to check genome quality -- for all assemblies at once
sbatch mcic-scripts/assembly/checkm.sh -i results/spades/all_assemblies -o results/checkm

# ==============================================================================
#                  IDENTIFY PLASMIDS WITH MOB-SUITE
# ==============================================================================
# Mob-suite
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
sbatch mcic-scripts/misc/sourmash_compare.sh -i "$indir" -o results/sourmash/ani --ani

# 4. Check previously generated Kraken results, especially for SM225
# less results/kraken/SpadesSM225-1_report.txt

# ==============================================================================
#                               Roary
# ==============================================================================
sbatch mcic-scripts/bact/roary.sh -i results/prokka/gff_withseqs_noSM225 -o results/roary

# ==============================================================================
#                         kSNP3 - own assemblies only
# ==============================================================================
# Run kSNP3 to get an alignment of core SNPs -- own genomes only
mkdir -p results/ksnp3/own_asm
paste <(ls "$PWD"/"$asm_dir"/*) <(ls "$asm_dir" | sed -E 's/Spades(.*).fasta/\1/') \
    >results/ksnp3/own_asm/assembly_list.txt
sbatch mcic-scripts/bact/ksnp3.sh -i results/ksnp3/own_asm/assembly_list.txt -o results/ksnp3/own_asm
Rscript mcic-scripts/trees/ggtree.R -i results/ksnp3/own_asm/tree.core.tre -o results/ksnp3/own_asm/tree.core.png

# ==============================================================================
#                         kSNP3 - with ref genomes
# ==============================================================================
# Renamed FASTAs for kSNP3
mkdir -p data/ref/renamed_links
for fna in "$PWD"/"$ref_dir"/*fna; do
    ln -sv "$fna" data/ref/renamed_links/"$(basename "${fna%%.*}")".fasta
done

# Prep kSNP3 input file
ksnp_dir=results/ksnp3/withrefs && mkdir -p "$ksnp_dir"
paste <(ls "$PWD"/"$asm_dir"/*fasta) <(ls "$asm_dir" | sed -E 's/Spades(.*).fasta/\1/') |
    grep -v "SM225-1" > "$ksnp_dir"/assembly_list.txt
paste <(ls "$PWD"/data/ref/renamed_links/*) <(ls data/ref/renamed_links | sed -E 's/.fasta//') \
    >> "$ksnp_dir"/assembly_list.txt

# Run kSNP3
sbatch mcic-scripts/bact/ksnp3.sh -i "$ksnp_dir"/assembly_list.txt -o "$ksnp_dir"

# Plot the tree
ls "$asm_dir" | sed -E 's/Spades(.*).fasta/\1/' | grep -v "SM225" > metadata/genomes.txt
micromamba activate /fs/ess/PAS0471/jelmer/conda/r_tree
Rscript mcic-scripts/trees/ggtree.R -i "$ksnp_dir"/tree.core.tre
