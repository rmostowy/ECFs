# Funciton to visulaise the protein alignment
plot_protein_alignment <- function(protein1_name, protein2_name,
                                   protein1_length, protein2_length, 
                                   aln1_start, aln1_end, aln2_start, aln2_end, 
                                   ecf1_start, ecf1_end, ecf2_start, ecf2_end, 
                                   aln_identity, protein1_description, protein2_description) {
  
  # Customization options for the alignment polygon) {
  
  # Customization options for the alignment polygon
  aln_colour <- "lightblue"  # Color of the alignment polygon
  domain_colour <- "red"     # Color of the domain regions
  
  # Create a base plot with horizontal lines for proteins
  p <- ggplot() +
    # Draw Protein 1 and Protein 2 as black horizontal lines
    geom_segment(aes(x = 0, xend = protein1_length, y = 2, yend = 2), linewidth = 2, color = "black") +
    geom_segment(aes(x = 0, xend = protein2_length, y = 1, yend = 1), linewidth = 2, color = "black")
  
  # Add the connecting rectangle for the alignment
  p <- p + geom_polygon(data = data.frame(
    x = c(aln1_start, aln1_end, aln2_end, aln2_start),
    y = c(2, 2, 1, 1)
  ), aes(x = x, y = y), fill = aln_colour, alpha = aln_identity) +
    # Add text annotations for protein descriptions
    annotate("text", x = 0, y = 1.03, label = protein2_description, hjust = 0, size = 5, vjust = 0) +
    annotate("text", x = 0, y = 1.97, label = protein1_description, hjust = 0, size = 5, vjust = 1)
  
  # Customize plot
  p <- p +
    labs(x = "Position", y = "") +
    scale_y_continuous(breaks = c(2, 1), labels = c(protein1_name, protein2_name)) +
    theme(axis.text.y = element_text(size = 12),
          plot.margin = margin(20, 30, 20, 30))  # Adjust margins to fit text
  
  return(p)
}

select.best.ecod.hits <- function(ecod.hits){
  ecod.hits.out <- NULL
  all.ecod.names <- ecod.hits$sname
  all.ecod.names <- all.ecod.names[!is.na(all.ecod.names)]
  n <- length(all.ecod.names)
  if(n > 0){
    setorder(ecod.hits, qname, t.id, -bitscore)
    ecod.hits <- ecod.hits[, .SD[1], by = .(qname, t.id)]
    ecod.hits[, name := sprintf("(%.2f) %s", scov, t.id)]
    ecod.hits.out <- ecod.hits %>% select(qstart, qend, name)
  }
  return(ecod.hits.out)
}

select.all.ecf.hits <- function(ecf.hits){
  ecf.hits.out <- NULL
  all.ecf.names <- ecf.hits$sname
  all.ecf.names <- all.ecf.names[!is.na(all.ecf.names)]
  n <- length(all.ecf.names)
  if(n > 0) ecf.hits.out <- ecf.hits %>% select(qstart, qend, name = sname, family = EF)
  return(ecf.hits.out)
}

select.best.ecf.hits <- function(ecf.hits){
  ecf.hits.out <- NULL
  all.ecf.names <- ecf.hits$sname
  all.ecf.names <- all.ecf.names[!is.na(all.ecf.names)]
  n <- length(all.ecf.names)
  if(n > 0){
    setorder(ecf.hits, qname, EF, -bitscore)
    ecf.hits <- ecf.hits[, .SD[1], by = .(qname, EF)]  
    ecf.hits[, ecf.name := sprintf("%s:%s", EF, sname)]
    ecf.hits.out <- ecf.hits %>% select(qstart, qend, name = ecf.name)
  }
  return(ecf.hits.out)
}

# Function to combine domain hits
combine_domain_hits <- function(..., types) {
  # List of input tables
  input_tables <- list(...)
  
  # Check for NULL input
  if (all(sapply(input_tables, is.null))) {
    return(NULL)  # Return NULL if all input tables are NULL
  }
  
  # Combine non-NULL tables with their types
  rbindlist(
    lapply(seq_along(input_tables), function(i) {
      if (!is.null(input_tables[[i]])) {
        data.table(type = types[i], input_tables[[i]])
      } else {
        NULL
      }
    }),
    use.names = TRUE, fill = TRUE
  )
}

# Function to visualise ECF domains and alignment on the query
plot_ecf_domains_with_alignment <- function(protein_name, protein_length, aln1_start, aln1_end, protein_ecf_hits_all) {
  
  # Check if domain_data is NULL or empty
  if (is.null(protein_ecf_hits_all) || nrow(protein_ecf_hits_all) == 0) {
    message("No domains detected for this protein.")
    return(NULL)
  }
  
  # Add a y-axis level for each ECF domain, reducing vertical spread for clarity
  protein_ecf_hits_all$y_level <- seq_len(nrow(protein_ecf_hits_all)) * 0.5
  
  # Define colors for families
  family_colors <- unique(protein_ecf_hits_all$family)
  color_palette <- setNames(RColorBrewer::brewer.pal(length(family_colors), "Set3"), family_colors)
  
  # Create the plot
  p <- ggplot() +
    # Plot the query protein as a black line with the aligned segment in red
    geom_segment(aes(x = 0, xend = protein_length, y = 0, yend = 0), linewidth = 2, color = "black") +
    geom_segment(aes(x = aln1_start, xend = aln1_end, y = 0, yend = 0), linewidth = 2, color = "red") +
    # Plot each ECF domain as a coloured segment, positioned on its y_level
    geom_segment(data = protein_ecf_hits_all, aes(x = qstart, xend = qend, y = y_level, yend = y_level, color = family),
                 linewidth = 4) +
    # Add labels for each ECF domain
    geom_text(data = protein_ecf_hits_all, aes(x = (qstart + qend) / 2, y = y_level + 0.1, label = name),
              size = 3, hjust = 0.5, color = "black") +
    # Define custom colours for families
    scale_color_manual(values = color_palette) +
    # Customise the plot
    labs(
      x = sprintf("Position along protein: %s", protein_name),
      y = "",
      title = "ECF Domains and Query Alignment",
      color = "Family"
    ) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.margin = margin(20, 30, 20, 30)
    ) +
    coord_cartesian(xlim = c(0, protein_length))
  
  return(p)
}

# Function to visualise protein domains with protein name in x-axis label
plot_protein_domains_compact <- function(protein_length, aln1_start, aln1_end, domain_data, protein_name) {
  
  # Check if domain_data is NULL or empty
  if (is.null(domain_data) || nrow(domain_data) == 0) {
    message("No domains detected for this protein.")
    return(NULL)
  }
  
  # Add a y-axis level for each domain, with reduced spacing for compactness
  domain_data$y_level <- seq_len(nrow(domain_data)) * 0.5  # Adjust spacing
  
  # Define colours for domain types
  domain_colours <- c("ECOD" = "red", "ECF" = "green")
  
  # Create the plot
  p <- ggplot(data = domain_data) +
    # Plot the full protein as a horizontal black line at the bottom
    geom_segment(aes(x = 0, xend = protein_length, y = 0, yend = 0),
                 linewidth = 2, color = "black") +
    geom_segment(aes(x = aln1_start, xend = aln1_end, y = 0, yend = 0), linewidth = 2, color = "red") +
    # Plot each domain as a segment coloured by its type
    geom_segment(aes(x = qstart, xend = qend, y = y_level, yend = y_level, color = type),
                 linewidth = 4) +
    # Add text labels for each domain's name above the segment
    geom_text(aes(x = qstart, y = y_level + 0.1, label = name),
              size = 4, hjust = 0, vjust = 0, color = "blue") +
    # Define custom colours for domain types
    scale_color_manual(values = domain_colours) +
    # Customise the plot
    labs(
      x = sprintf("%s", protein_name),  # Add protein name to x-axis
      y = "",
      color = "Domain Type"
    ) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.margin = margin(20, 30, 20, 30)
    ) +
    coord_cartesian(xlim = c(0, protein_length))  # Ensure x-axis matches protein length
  
  return(p)
}
