# Function to calculate overlap details
calculate_overlap_parms <- function(start1, end1, start2, end2) {
  overlap.start <- pmax(start1, start2)
  overlap.end <- pmin(end1, end2)
  overlap.length <- pmax(0, overlap.end - overlap.start + 1)
  return(list(overlap.start = overlap.start, overlap.end = overlap.end, overlap.length = overlap.length))
}

# Function to calculate overlap between two intervals
calculate_overlap_jaccard <- function(start1, end1, start2, end2) {
  # Calculate the intersection (overlap length)
  overlap <- pmax(0, pmin(end1, end2) - pmax(start1, start2) + 1)
  
  # Calculate the union (total combined region covered by both CDSs)
  total_region <- (end1 - start1 + 1) + (end2 - start2 + 1) - overlap
  
  # Compute Jaccard similarity (Intersection over Union)
  jaccard <- overlap / total_region
  
  # Handle cases where there's no valid overlap
  jaccard[is.na(jaccard)] <- 0  # If NA, set to 0
  
  return(jaccard)
}

get.fragment.domain.overlap <- function(fragment_start, fragment_end, domain_start, domain_end){
  intersection_start <- max(fragment_start, domain_start)
  intersection_end <- min(fragment_end, domain_end)
  intersect_length <- max(0, intersection_end - intersection_start + 1)
  union_start <- min(fragment_start, domain_start)
  union_end <- max(fragment_end, domain_end)
  union_length <- union_end - union_start + 1
  overlap <- intersect_length/union_length
  return(overlap)
}

dt.round.digits <- function(DT, digits = 4, byref = TRUE, cols = NULL, include_integer64 = FALSE) {
  if (!is.data.table(DT)) stop("DT must be a data.table")
  
  # Pick candidate columns
  all_cols <- if (is.null(cols)) names(DT) else intersect(cols, names(DT))
  
  is_int64 <- function(x) inherits(x, "integer64")
  is_double <- function(x) is.double(x) && !is_int64(x)
  
  # Start with plain doubles
  target <- all_cols[vapply(DT[, ..all_cols], is_double, logical(1))]
  
  # Optionally include integer64 (with precision warning!)
  if (include_integer64) {
    int64_cols <- all_cols[vapply(DT[, ..all_cols], is_int64, logical(1))]
    target <- c(target, int64_cols)
  }
  
  if (length(target) == 0L) {
    if (byref) return(invisible(DT))
    else return(copy(DT))
  }
  
  if (byref) {
    DT[, (target) := lapply(.SD, signif, digits = digits), .SDcols = target]
    invisible(DT)
  } else {
    out <- copy(DT)
    out[, (target) := lapply(.SD, signif, digits = digits), .SDcols = target]
    out
  }
}

map.ecf.coords.to.genome <- function(aa_start, aa_end,
                             gene_start, gene_end,
                             strand, protlen) {

  # convert AA coords to nucleotide offsets inside CDS
  nt_start <- (aa_start - 1) * 3 + 1
  nt_end   <- (aa_end   - 1) * 3 + 3

  # quick sanity check vs protein length
  if (aa_end > protlen)
    warning("ECF exceeds protein length")

  # gene bounds on genome
  gmin <- min(gene_start, gene_end)
  gmax <- max(gene_start, gene_end)

  # map to genome coordinates
  if (strand == "+") {
    ecf_start <- gmin + nt_start - 1
    ecf_end   <- gmin + nt_end   - 1
  } else {
    ecf_start <- gmax - nt_start + 1
    ecf_end   <- gmax - nt_end   + 1
  }

  # return genome interval (start < end)
  c(start = min(ecf_start, ecf_end),
    end   = max(ecf_start, ecf_end))
}

get.min.parms <- function(hhtable){
  col.names2check <- c("qname", "sname", "qcov", "scov", "eval", "prob", "hit.length")
  hh.table.cols <- colnames(hhtable)
  
  missing.colnames <- col.names2check[!col.names2check %in% hh.table.cols]
  if(length(missing.colnames) > 0){
    cat("Missing columns:\n", paste(missing.colnames, collapse = ", "), "\n")
  } else{
    cat("All columns in place!\n")
  }
  
  # helper: safe min/max with NA handling
  safe_min <- function(x) if(length(x) == 0 || all(is.na(x))) NA else min(x, na.rm = TRUE)
  safe_max <- function(x) if(length(x) == 0 || all(is.na(x))) NA else max(x, na.rm = TRUE)
  
  # first query
  if("qname" %in% hh.table.cols){
    cat("first query:", hhtable[["qname"]][1], "\n")
  }
  
  # first subject
  if("sname" %in% hh.table.cols){
    cat("first subject:", hhtable[["sname"]][1], "\n")
  }
  
  # min prob
  if("prob" %in% hh.table.cols){
    cat("min prob:", safe_min(hhtable[["prob"]]), "\n")
  }
  
  # max eval
  if("eval" %in% hh.table.cols){
    cat("max eval:", safe_max(hhtable[["eval"]]), "\n")
  }
  
  # min qcov
  if("qcov" %in% hh.table.cols){
    cat("min qcov:", safe_min(hhtable[["qcov"]]), "\n")
  }
  
  # min scov
  if("scov" %in% hh.table.cols){
    cat("min scov:", safe_min(hhtable[["scov"]]), "\n")
  }
  
    # min hit length
  if("hit.length" %in% hh.table.cols){
    cat("min hit.length:", safe_min(hhtable[["hit.length"]]), "\n")
  }
  
  invisible(NULL)
}

plot.ecf.distr.over.query <- function(protein_name, ecf.hits.dt, title = "ECF Domains and Query Alignment") {
  
  # Check if domain_data is NULL or empty
  if (is.null(ecf.hits.dt) || all(is.na(ecf.hits.dt$qstart))) {
    message("No domains detected for this protein.")
    return(NULL)
  }
  
  protein_length <- ecf.hits.dt$qlength[1]
  
  # sort by ECF family
  setorder(ecf.hits.dt, family, qstart)
  
  # Add a y-axis level for each ECF domain, reducing vertical spread for clarity
  ecf.hits.dt$y_level <- seq_len(nrow(ecf.hits.dt)) * 0.5
  
  # Define colors for families
  family_colors <- unique(ecf.hits.dt$family)
  no.fam <- length(family_colors)
  # if(no.fam < 3) no.fam <- 3
  # color_palette <- setNames(RColorBrewer::brewer.pal(no.fam, "Set3"), family_colors)
  base_cols <- RColorBrewer::brewer.pal(12, "Set3")
  color_palette <- setNames(
    colorRampPalette(base_cols)(length(family_colors)),
    family_colors
  )
  
  # Create the plot
  p <- ggplot() +
    # Plot the query protein as a black line with the aligned segment in red
    geom_segment(aes(x = 0, xend = protein_length, y = 0, yend = 0), linewidth = 2, color = "black") +
    # geom_segment(aes(x = aln1_start, xend = aln1_end, y = 0, yend = 0), linewidth = 2, color = "red") +
    # Plot each ECF domain as a coloured segment, positioned on its y_level
    geom_segment(data = ecf.hits.dt, aes(x = qstart, xend = qend, y = y_level, yend = y_level, color = family),
                 linewidth = 4) +
    # Add labels for each ECF domain
    geom_text(data = ecf.hits.dt, aes(x = (qstart + qend) / 2, y = y_level + 0.1, label = sname),
              size = 3, hjust = 0.5, color = "black") +
    # Define custom colours for families
    scale_color_manual(values = color_palette) +
    # Customise the plot
    labs(
      x = sprintf("rHMM position: %s", protein_name),
      y = "",
      title = title,
      color = "Family"
    ) +
    # theme_classic(base_size = 12) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    ) +
    coord_cartesian(xlim = c(0, protein_length))
  
  
  return(p)
}
