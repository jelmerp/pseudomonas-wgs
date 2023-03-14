#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=24:00:00
#SBATCH --job-name=dl-refs
#SBATCH --output=slurm-dl-all-pseudo-%j.out

## Args
outdir=$1

## Make output dir
mkdir -p "$outdir"

## Report
echo "## Starting script dl-all-pseudo.sh"
date
echo
echo "## Output dir:            $outdir"
echo -e "----------------------\n"

# https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=286&lvl=3&lin=f&keep=1&srchmode=1&unlock

## Get list of refs
esearch -db assembly -query 'txid286[Organism:exp]' |
    esummary |
    xtract -pattern DocumentSummary -element FtpPath_GenBank \
    >"$outdir"/ftp_list.txt

## Download refs
cat "$outdir"/ftp_list.txt | while read -r line; do
    fname=$(echo "$line" | grep -o 'GCA_.*' | sed 's/$/_genomic.fna.gz/')
    echo "$fname"
    curl -L -o "$outdir"/"$fname" "$line/$fname"
    echo -e "-------------\n"
done

## Report
echo -e "\n--------------------------"
echo "## Number of files in output dir: $(ls "$outdir" | wc -l)"
echo "## Done with script dl-all-pseudo.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
