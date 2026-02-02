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

map_ecf_to_genome <- function(aa_start, aa_end,
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