# Load required libraries
library(ggtree)
library(gggenes)
library(ggplot2)
library(dplyr)
library(tibble)
library(ggtreeExtra)
library(patchwork)

# Simulated phylogenetic tree
set.seed(42)
tree <- rtree(5)  # Tree with 5 proteins
tree$tip.label <- paste0("protein", 1:5)  # Rename tips to protein names

# Simulated PFAM domain data
pfam_data <- data.frame(
  protein = rep(tree$tip.label, each = 3),
  domain = rep(c("PF00069", "PF00520", "PF12345"), times = 5),
  start = rep(c(5, 40, 80), times = 5),
  end = rep(c(35, 75, 110), times = 5)
)

# Simulated ECOD domain data
ecod_data <- data.frame(
  protein = rep(tree$tip.label, each = 2),
  domain = rep(c("ECOD_001", "ECOD_002"), times = 5),
  start = rep(c(10, 50), times = 5),
  end = rep(c(30, 90), times = 5)
)

# Simulated Conserved Fragment Data
fragment_data <- data.frame(
  protein = tree$tip.label,  # Same proteins as tree
  frag_start = c(15, 25, 60, 10, 45),  # Example conserved fragment positions
  frag_end = c(50, 55, 85, 40, 80)
)

# Convert tree into a tibble with positions for alignment
tree_data <- tibble(label = tree$tip.label)
pfam_data$y <- match(pfam_data$protein, tree$tip.label)
ecod_data$y <- match(ecod_data$protein, tree$tip.label)
fragment_data$y <- match(fragment_data$protein, tree$tip.label)  # Align fragments too

# Create the tree plot
tree_plot <- ggtree(tree) +
  theme_tree2()


# PFAM domain plot with gene arrows + conserved fragment as a segment
pfam_plot <- ggplot(pfam_data, aes(xmin = start, xmax = end, y = protein, fill = domain)) +
  geom_gene_arrow() +  # Use gene arrows for PFAM domains
  geom_segment(data = fragment_data, 
               aes(x = frag_start, xend = frag_end, y = protein, yend = protein), 
               color = "black", size = 2, alpha = 0.6, inherit.aes = FALSE) +  # Mark conserved fragment
  scale_fill_manual(values = c("PF00069" = "blue", "PF00520" = "red", "PF12345" = "green")) + 
  theme_minimal() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  labs(x = "PFAM Domains", fill = "PFAM")

# ECOD domain plot with gene arrows + conserved fragment as a segment
ecod_plot <- ggplot(ecod_data, aes(xmin = start, xmax = end, y = protein, fill = domain)) +
  geom_gene_arrow() +  # Use gene arrows for ECOD domains
  geom_segment(data = fragment_data, 
               aes(x = frag_start, xend = frag_end, y = protein, yend = protein), 
               color = "black", size = 2, alpha = 0.6, inherit.aes = FALSE) +  # Highlight conserved fragment
  scale_fill_manual(values = c("ECOD_001" = "purple", "ECOD_002" = "orange")) + 
  theme_minimal() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  labs(x = "ECOD Domains", fill = "ECOD")

# Combine tree, PFAM, and ECOD plots
tree_plot + pfam_plot + ecod_plot + plot_layout(widths = c(1, 1.5, 1.5))