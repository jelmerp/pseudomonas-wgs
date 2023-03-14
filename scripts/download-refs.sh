#!/bin/bash

## Load software
module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/miniconda3/envs/blast-env
export NCBI_API_KEY=34618c91021ccd7f17429b650a087b585f08

## Constants
outdir=results/refgenomes
ncbi_ids=(ASM1669475v2 ASM780v1 ASM1224v1 ASM1220v1 ASM14582v2 ASM14584v2 ASM2327812v1 ASM15699v2 ASM276365v1 ASM45244v3 ASM1860349v1 ASM14594v2 ASM290581v2 ASM2327794v1 CFBP4215 CFBP2118 ASM98848v1 ASM98839v1)

## Create output dir
mkdir -p "$outdir"

## Download ref FASTA files
for ncbi_id in "${ncbi_ids[@]}"; do 
    echo -e "\n$ncbi_id"
    esearch -db assembly -query "$ncbi_id" |
        elink -target nucleotide -name assembly_nuccore_insdc |
        efetch -format fasta > "$outdir"/"$ncbi_id".fasta
    ls -lh "$outdir"/"$ncbi_id".fasta
done

echo "Done."
