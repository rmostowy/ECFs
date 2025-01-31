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