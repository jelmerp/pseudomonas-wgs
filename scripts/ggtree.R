# SETUP ------------------------------------------------------------------------
# Working dir should be /fs/ess/PAS0471/jelmer/assist/01_archive/2022-06_sochina

# Load packages
library(tidyverse)
library(ggtree)
library(ggtreeExtra)   # geom_fruit()
library(ggnewscale)    # Allow for multiple color scales
library(ggforce)
library(RColorBrewer)

# Define the input files
indir <- "results/ksnp3/w_moraviensis/"  # This is the final choice of tree/kSNP3 run
#indir <- "results/ksnp3/no_moraviensis" # Alternative tree, without GCA_000512275 (moraviensis)
#indir <- "results/ksnp3/withrefs"       # Alternative tree, with aeruginosa
tree_file <- file.path(indir, "tree.core.tre")
ncbi_meta_file <- "data/ref/metadata/meta_sel.tsv"
ncbi_host_file <- "metadata/refgenomes_hostdata.tsv"
ncbi_phylogroup_file <- "data/ref/metadata/meta_phylogroup.tsv"
own_meta_file <- "metadata/owngenome_metadata.tsv"

# Define the output file
figfile <- file.path(indir, "tree.core_withmeta.png")

# Settings
root_id <- "GCA_000512275"
species_levels <- c("syringae", "aeruginosa", "coronafaciens",
                    "moraviensis", "tremae", "amygdali", "savastanoi")


# PREP THE INPUT DATA ----------------------------------------------------------
# Read the tree, and reroot it
tree <- read.tree(tree_file)
tree <- ape::root(tree, outgroup = root_id)

# Read and process the metadata on the NCBI genomes
ncbi_host <- read_tsv(ncbi_host_file, show_col_types = FALSE) |>
  mutate(ID = sub("\\.\\d+$", "", assembly)) |>
  select(ID, host)
ncbi_phylo <- read_tsv(ncbi_phylogroup_file, show_col_types = FALSE) |>
  mutate(ID = sub("\\.\\d+$", "", ID),
         phylogroup = sub("Phylogroup ", "", phylogroup))
ncbi_meta <- read_tsv(ncbi_meta_file, show_col_types = FALSE) |>
  janitor::clean_names() |>
  mutate(ID = sub("\\.\\d+$", "", assembly_accession)) |> 
  select(ID, organism_name) |>
  mutate(species = sub("Pseudomonas (\\w+).*", "\\1", organism_name),
         pathovar = ifelse(grepl("pv.", organism_name), organism_name, NA),
         pathovar = sub("Pseudomonas \\w+ pv. ", "", pathovar),
         pathovar = ifelse(species != "syringae", NA, pathovar),
         pathovar = ifelse(organism_name == "Pseudomonas syringae Cit 7", "Cit 7", pathovar),
         pathovar_short = sub(" .*", "", pathovar),
         pathovar_short = ifelse(organism_name == "Pseudomonas syringae Cit 7",
                                 "Cit 7", pathovar_short),
         pathovar_short = ifelse(is.na(pathovar_short) & species == "syringae",
                                 "(unknown)", pathovar_short)) |>
  left_join(ncbi_host, by = "ID") |>
  left_join(ncbi_phylo, by = "ID")

# Read and process the metadata on the own genomes
own_meta <- read_tsv(own_meta_file, show_col_types = FALSE) |> 
  mutate(species = "syringae",
         pathovar_short = "syringae")

# Combine metadata
meta <- bind_rows(ncbi_meta, own_meta) |>
  mutate(species = factor(species, level = species_levels),
         pathovar_short = relevel(factor(pathovar_short), ref = "syringae"))

# Colors
species_cols <- c(
  "grey40",
  brewer.pal(n = length(species_levels) - 1, name = "Dark2")
  )
pathovar_cols <- c(
  "grey40",
  as.vector(paletteer::paletteer_d("basetheme::ink"))[1:length(levels(meta$pathovar_short)) - 1]
  )

# Make very long branch shorter, and give the shortenend branch a different line-type
## (https://yulab-smu.top/treedata-book/faq.html#change-colors-or-line-types-of-arbitrarily-selected-branches)
## ggtree(tree) + geom_text(aes(label = node), hjust = -.3) # Get tip and node labels
tree$edge.length[which(tree$edge.length == max(tree$edge.length))] <- 0.25
focal_branch <- 33
linetype_df <- data.frame(node=1:treeio::Nnode2(tree), lty = 1)
linetype_df[focal_branch, 2] <- 2


# PLOT THE TREE ----------------------------------------------------------------
ggtree(tree,
       aes(linetype = I(lty)),
       layout = "roundrect") %<+%
  # Add main metadata df
  meta %<+%
  # Add linetype df
  linetype_df +
  # Add an edge at the root
  geom_rootedge(rootedge = sum(tree$edge.length) / 100) +
  # Add bootstrap support
  geom_nodepoint(aes(subset = !is.na(as.numeric(label)),
                     fill = as.numeric(label)),
                 shape = 21, size = 2, color = "grey30") +
  # Fill scale for bootstrap support
  scale_fill_viridis_c(option = "magma", name = "Bootstrap\nsupport",
                       breaks = c(0, 0.5, 1)) +
  # Main tip-labels
  geom_tiplab(align = TRUE, offset = 0.02, size = 3.5, linesize = 0.1, color = "grey40") +
  # Phylogroup labels
  geom_tiplab(aes(label = phylogroup, color = phylogroup), fontface = "bold",
              align = TRUE, linetype = "blank", offset = 0.13, size = 3.5) +
  ggsci::scale_color_jama(guide = "none") +
  # Species labels, colored by species
  new_scale_color() +
  geom_tiplab(aes(label = species, color = species), fontface = "italic",
              align = TRUE, linetype = "blank", offset = 0.16, size = 3.5) +
  scale_color_manual(name = "Species", values = species_cols, guide = "none") +
  # Pathovar labels (only for P.s.), colored by pathovar
  new_scale_color() +
  geom_tiplab(aes(label = pathovar_short, color = pathovar_short),
              align = TRUE, linetype = "blank", offset = 0.24, size = 3.5,
              fontface = "italic") +
  scale_color_manual(values = pathovar_cols, guide = "none") +
  # Host labels, colored by host
  new_scale_color() +
  geom_tiplab(aes(label = host, color = host),
              align = TRUE, linetype = "blank", offset = 0.325, size = 3.5) +
  scale_color_discrete(guide = "none") +
  # Disease severity heatmap
  new_scale_fill() +
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = severity, y = ID),
             offset = 0.01, pwidth = 0.02, width = 0.03, color = "grey20") +
  scale_fill_viridis_c(name = "Disease\nSeverity", na.value = "white") +
  # Manually add the 'column' labels for species, pathovar, and host data
  annotate(geom = "text",
           x = 0.277, y = length(tree$tip.label) + 1,
           label = "ID", fontface = "bold", size = 3.5) +
  annotate(geom = "text",
           x = 0.37, y = length(tree$tip.label) + 1,
           label = "Phylogroup", fontface = "bold", size = 3.5) +
  annotate(geom = "text",
           x = 0.433, y = length(tree$tip.label) + 1,
           label = "Species", fontface = "bold", size = 3.5) +
  annotate(geom = "text",
           x = 0.530, y = length(tree$tip.label) + 1,
           label = "P. s. pathovar", fontface = "bold", size = 3.5) + #hjust = 0, 
  annotate(geom = "text",
           x = 0.589, y = length(tree$tip.label) + 1,
           label = "Host", fontface = "bold", size = 3.5) + #hjust = 0
  # This avoids the labels being cut off despite wide right margin
  coord_cartesian(clip = "off") +
  # Tree formatting - e.g. wide right-hand margin needed
  theme(plot.margin = margin(0.2, 3, 0.2, 0.2, "cm"),
        legend.position = "top",
        legend.title = element_text(vjust = 0.8, face = "bold"),
        legend.box.margin = margin(1, 5, 1, 1),
        legend.background = element_rect(fill = "grey95"),
        legend.box.background = element_rect(fill = "grey95"))

ggsave(figfile, width = 10, height = 7, dpi = "retina")


# ------------------------------------------------------------------------------
# Alternative ways of showing bootstrap support - with text labels:

#geom_text2(aes(label = label,
#               subset = !is.na(as.numeric(label)) & as.numeric(label) < 0.80),
#           color = "grey50", nudge_y = 0.1, nudge_x = -0.03, size = 2.5)

#geom_point2(aes(subset = !is.na(as.numeric(label)),
#                fill = cut(as.numeric(label), c(0, 0.5, 0.8, 1))), 
#            shape = 21, size = 3)
