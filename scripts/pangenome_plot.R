# Load packages
library(pagoo)

# Define input files
roary_outfile <- "results/roary/gene_presence_absence.csv"

# Define output files
outdir <- "results/plots"
plotfile <-  file.path(outdir, "pangenome_rarefaction.png")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Load the input
pg <- roary_2_pagoo(roary_outfile)

# Create the pangenome rarefaction graph (https://iferres.github.io/pagoo/articles/Methods_Plots.html)
pg$gg_curves() + 
  geom_point() + 
  facet_wrap(~Category, scales = 'free_y') +
  scale_color_brewer(palette = "Accent", guide = "none") +
  scale_y_continuous(labels = scales::comma) +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank())

ggsave(plotfile, width = 7.5, height = 6, dpi = "retina")
