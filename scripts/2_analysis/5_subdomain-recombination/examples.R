library(ggplot2)

# Define parameters
protein1_length <- 100       # Length of Protein 1
protein2_length <- 200       # Length of Protein 2

aln1_start <- 50         # Start of alignment in Protein 1
aln1_end <- 100          # End of alignment in Protein 1
aln2_start <- 150        # Start of alignment in Protein 2
aln2_end <- 200          # End of alignment in Protein 2

ecf1_start <- 70
ecf1_end <- 90
ecf2_start <- 170
ecf2_end <- 190

aln_identity <- 0.2              # Percentage identity (e.g., 0.8 for 80%)
ecf_identity <- 0.7

# Customization options for the alignment polygon
aln_colour <- "lightblue" # Color of the polygon
ecf_colour <- "darkblue" # Color of the polygon

# Create a base plot with horizontal lines for proteins
p <- ggplot() +
  geom_segment(aes(x = 0, xend = protein1_length, y = 1, yend = 1), size = 2, color = "black") +
  geom_segment(aes(x = 0, xend = protein2_length, y = 2, yend = 2), size = 2, color = "black") +
  geom_text(aes(x = protein1_length + 10, y = 1, label = "Protein1"), size = 5, vjust = 0.5, hjust = 0) +
  geom_text(aes(x = protein2_length + 10, y = 2, label = "Protein2"), size = 5, vjust = 0.5, hjust = 0)

# Add the connecting rectangle for the alignment
p <- p + geom_polygon(data = data.frame(
  x = c(aln1_start, aln1_end, aln2_end, aln2_start),
  y = c(1, 1, 2, 2)
), aes(x = x, y = y), fill = aln_colour, alpha = aln_identity)

# Add the connecting rectangle for the alignment
p <- p + geom_polygon(data = data.frame(
  x = c(ecf1_start, ecf1_end, ecf2_end, ecf2_start),
  y = c(1, 1, 2, 2)
), aes(x = x, y = y), fill = ecf_colour, alpha = ecf_identity)

# Customize plot
p <- p + theme_minimal() +
  labs(title = paste("Protein Alignment Visualization -", round(identity * 100), "% Identity"), 
       x = "Protein Sequence Position", y = "") +
  scale_y_continuous(breaks = c(1, 2), labels = c("Protein1", "Protein2")) +
  theme(axis.text.y = element_text(size = 12))

print(p)