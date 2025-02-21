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
domain_data <- data.frame(
  protein = rep(tree$tip.label, each = 3),  # Each protein has 3 domains
  domain = rep(c("PF00069", "PF00520", "PF12345"), times = 5),  # Example PFAMs
  start = rep(c(5, 40, 80), times = 5),  # Domain start positions
  end = rep(c(35, 75, 110), times = 5)  # Domain end positions
)

# Convert tree into a tibble with positions for alignment
tree_data <- tibble(label = tree$tip.label)
domain_data$y <- match(domain_data$protein, tree$tip.label)  # Align to tree order

# Create the tree plot
tree_plot <- ggtree(tree) +
  theme_tree2()

# Create the domain plot without arrows
domain_plot <- ggplot(domain_data, aes(xmin = start, xmax = end, y = protein, fill = domain)) +
  geom_rect(aes(ymin = as.numeric(y) - 0.4, ymax = as.numeric(y) + 0.4)) +  # Simple rectangles
# geom_gene_arrow() +  # Or use arrows instead
  scale_fill_manual(values = c("PF00069" = "blue", "PF00520" = "red", "PF12345" = "green")) + 
  theme_minimal() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  labs(x = "Protein Position", fill = "PFAM Domain")

# Combine tree and domain plots
tree_plot + domain_plot + plot_layout(widths = c(1, 1.5))  # Adjust widths for readability
