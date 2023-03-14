
# SET-UP -----------------------------------------------------------------------
# Command-line args
genome_id <- "SM51-19"
virgenes_gff_dir <- "results/pseudofinder/virgenes"
snippy_dir <- "results/snippy"
outdir <- "results/snippy/virgenes"

# Packages
library(tidyverse)

# Input files
biofilm_infile <- file.path(virgenes_gff_dir,
                            paste0("biofilm_", genome_id, "_all.gff"))
motility_infile <- file.path(virgenes_gff_dir,
                            paste0("motility_", genome_id, "_all.gff"))
secretion_infile <- file.path(virgenes_gff_dir,
                            paste0("secretion_", genome_id, "_all.gff"))

# Output files
biofilm_outfile <- file.path(outdir, paste0("biofilm_", genome_id, "_snps.tsv"))

# Make outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


# PREP INPUT FILES -------------------------------------------------------------
biofilm_all <- read_tsv(biofilm_infile, col_names = FALSE, show_col_types = FALSE)

# INTERSECT --------------------------------------------------------------------


