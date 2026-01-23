### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(gggenes)))
suppressWarnings(suppressPackageStartupMessages(library(ggtree)))
suppressWarnings(suppressPackageStartupMessages(library(tibble)))
suppressWarnings(suppressPackageStartupMessages(library(ggtreeExtra)))
suppressWarnings(suppressPackageStartupMessages(library(patchwork)))
suppressWarnings(suppressPackageStartupMessages(library(ape)))
suppressWarnings(suppressPackageStartupMessages(library(phangorn)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

source("functions.R")

# Define output folder
table.dir <- "tables"
output.dir <- file.path(visual.dir, "part5_rem", "rem-equiv-ecfs")
ecf.subseqs.dir <- file.path(output.dir, "ecf-subseqs")
ecf.trees.dir <- file.path(output.dir, "ecf-trees")
make.dir(ecf.trees.dir)


all.ecfs <- list.files(ecf.subseqs.dir)
all.ecfs <- unique(substr(all.ecfs, 1, 9))


for(index in 1:length(all.ecfs)){
  ecf.name <- all.ecfs[index]
  ecf.tree.filename <- file.path(ecf.trees.dir, sprintf("%s-tree.png", ecf.name))
  if(!file.exists(ecf.tree.filename)){
    # Load phylogenetic tree
    tree_file <- file.path(ecf.subseqs.dir, sprintf("%s.aln.treefile", ecf.name))
    tree.unrooted <- read.tree(tree_file)
    tree <- midpoint(tree.unrooted)
    
    # Load ECOD domain data
    ecod_file <- file.path(ecf.subseqs.dir, sprintf("%s-ecod.txt", ecf.name))
    ecod_data <- read.csv(ecod_file, stringsAsFactors = FALSE)
    
    # Load PFAM domain data
    pfam_file <- file.path(ecf.subseqs.dir, sprintf("%s-pfam.txt", ecf.name))
    pfam_data <- read.csv(pfam_file, stringsAsFactors = FALSE)
    
    # Load conserved fragments and protein lengths
    frag_file <- file.path(ecf.subseqs.dir, sprintf("%s-frag.txt", ecf.name))
    fragment_data <- read.csv(frag_file, stringsAsFactors = FALSE)
    
    # Convert tree into a tibble for alignment
    tree_data <- tibble(label = tree$tip.label)
    
    # Ensure domain & fragment data align with tree
    pfam_data$y <- match(pfam_data$name, tree$tip.label)
    ecod_data$y <- match(ecod_data$name, tree$tip.label)
    fragment_data$y <- match(fragment_data$name, tree$tip.label)
    
    # Generate distinct colors for all unique ECOD and PFAM IDs
    set.seed(42)  # Ensures consistent color assignment
    unique_ecod_ids <- unique(ecod_data$id)
    unique_pfam_ids <- unique(pfam_data$id)
    
    # Assign colors dynamically for ECOD and PFAM **IDs**
    ecod_colors <- setNames(colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(unique_ecod_ids)), unique_ecod_ids)
    pfam_colors <- setNames(colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(unique_pfam_ids)), unique_pfam_ids)
    
    # Create Phylogenetic Tree Plot
    tree_plot <- ggtree(tree) +
      theme_tree2()
    
    # **PFAM domain plot with protein length as a grey line**
    pfam_plot <- ggplot(pfam_data, aes(xmin = start, xmax = end, y = name, fill = id)) +
      # Grey line for full protein length
      geom_segment(data = fragment_data, aes(x = 0, xend = length, y = name, yend = name), 
                   color = "grey40", size = 1, inherit.aes = FALSE) +
      # PFAM domains as arrows
      geom_gene_arrow() +  
      # Conserved fragment as a red segment
      geom_segment(data = fragment_data, 
                   aes(x = start, xend = end, y = name, yend = name), 
                   color = "red", size = 2, alpha = 0.8, inherit.aes = FALSE) +
      scale_fill_manual(values = pfam_colors) + 
      theme_minimal() +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
      labs(x = "PFAM Domains", fill = "PFAM ID")
    
    # **ECOD domain plot with protein length as a grey line**
    ecod_plot <- ggplot(ecod_data, aes(xmin = start, xmax = end, y = name, fill = id)) +
      # Grey line for full protein length
      geom_segment(data = fragment_data, aes(x = 0, xend = length, y = name, yend = name), 
                   color = "grey40", size = 1, inherit.aes = FALSE) +
      # ECOD domains as arrows
      geom_gene_arrow() +  
      # Conserved fragment as a red segment
      geom_segment(data = fragment_data, 
                   aes(x = start, xend = end, y = name, yend = name), 
                   color = "red", size = 2, alpha = 0.8, inherit.aes = FALSE) +
      scale_fill_manual(values = ecod_colors) + 
      theme_minimal() +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
      labs(x = "ECOD Domains", fill = "ECOD ID")
    
    # Combine tree, PFAM, and ECOD plots
    final_plot <- tree_plot + pfam_plot + ecod_plot + plot_layout(widths = c(1, 1.5, 1.5))
    ggsave(file = ecf.tree.filename, final_plot, width = 15, height = 9)  
  }
  
}
