# SET-UP -----------------------------------------------------------------------
# Packages
library(tidyverse)
library(janitor)

# Input files
indir <- "results/virgenes_snps"

# Output files
outdir <- "results/virgenes_snps"
outfile_recoded <- file.path(outdir, "allgenomes_nonsyn_recoded.tab")
outfile_genotypes <- file.path(outdir, "allgenomes_nonsyn_genotypes.tab")

# Make outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# PROCESS ----------------------------------------------------------------------
infiles <- list.files(path = indir, pattern = "_nonsyn.tab",
                      recursive = TRUE, full.names = TRUE)
infiles <- infiles[!grepl("snps_", infiles)]

df_long <- read_tsv(infiles, show_col_types = FALSE, id = "file_id") |>
  janitor::clean_names() |> 
  mutate(assembly = gsub(".*virgenes_snps/(.*)/\\w+_nonsyn.tab", "\\1", file_id),
         gene_type = gsub(".*virgenes_snps/.*/(\\w+)_nonsyn.tab", "\\1", file_id)) |>
  select(assembly, gene_type, chrom, pos, type, strand, locus_tag, gene, product,
         ref, alt)

df_wide <- df_long |>
  pivot_wider(names_from = assembly, values_from = alt) |>
  mutate(alt = do.call(coalesce, across(starts_with("SM")))) |>
  relocate(alt, .after = ref)
  
df_recoded <- df_long |>
  mutate(alt2 = ifelse(is.na(alt), 0, 1)) |> 
  pivot_wider(names_from = assembly, values_from = alt2, values_fill = 0)

# Write the output files
write_tsv(df_recoded, outfile_recoded)
write_tsv(df_wide, outfile_genotypes)
