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
suppressWarnings(suppressPackageStartupMessages(library(phytools)))
suppressWarnings(suppressPackageStartupMessages(library(Biostrings)))
suppressWarnings(suppressPackageStartupMessages(library(DECIPHER)))


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
    p <- ggtree(tree)
    tree.taxa.order <- get_taxa_name(p)
    
    # Load alignment
    aln_file <- file.path(ecf.subseqs.dir, sprintf("%s.aln", ecf.name))
    aln <- readAAStringSet(aln_file)  # must be aligned
    
    # Cluster closely related sequences
    D <- DistanceMatrix(aln, minCoverage = RE.MIN.OVERLAP, verbose = F)
    D[is.na(D)] <- 1
    hc <- hclust(as.dist(D), method = "complete")
    cl <- cutree(hc, h = RE.MAX.DISTANCE)
    cl.sizes <- table(cl)
    cluster.names.above.thr <- which(cl.sizes >= RE.MIN.CLUSTER.SIZE)
    closely.related.taxa <- names(cl[which(cl %in% cluster.names.above.thr)])
    tip_metadata <- data.frame(
      label = tree$tip.label,
      highlight = tree$tip.label %in% closely.related.taxa
    )

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
    tree_data <- tibble(label = tree.taxa.order)
    
    # Ensure domain & fragment data align with tree
    pfam_data$y <- match(pfam_data$name, tree.taxa.order)
    ecod_data$y <- match(ecod_data$name, tree.taxa.order)
    fragment_data$y <- match(fragment_data$name, tree.taxa.order)
    
    # Ensure Y axis names align with the tree
    pfam_data$name <- factor(pfam_data$name, levels = rev(tree.taxa.order))
    ecod_data$name <- factor(ecod_data$name, levels = rev(tree.taxa.order))
    fragment_data$name <- factor(fragment_data$name, levels = rev(tree.taxa.order))

    # Generate distinct colors for all unique ECOD and PFAM IDs
    set.seed(42)  # Ensures consistent color assignment
    unique_ecod_ids <- unique(ecod_data$id)
    unique_pfam_ids <- unique(pfam_data$id)
    
    # Assign colors dynamically for ECOD and PFAM **IDs**
    ecod_colors <- setNames(colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(unique_ecod_ids)), unique_ecod_ids)
    pfam_colors <- setNames(colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(unique_pfam_ids)), unique_pfam_ids)
    
    # Create Phylogenetic Tree Plot
    tree_depth <- max(nodeHeights(tree))
    tree_plot <- ggtree(tree) %<+% tip_metadata +
      theme_tree2() + geom_tiplab(
        aes(color = highlight),
        size = 3
      ) + scale_color_manual(
        values = c("FALSE" = "black", "TRUE" = "red"),
        guide = "none"
      ) + 
      xlim(0, tree_depth*1.35)

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
