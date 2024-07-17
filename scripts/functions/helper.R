# Function to calculate overlap details
calculate_overlap_parms <- function(start1, end1, start2, end2) {
  overlap.start <- pmax(start1, start2)
  overlap.end <- pmin(end1, end2)
  overlap.length <- pmax(0, overlap.end - overlap.start)
  return(list(overlap.start = overlap.start, overlap.end = overlap.end, overlap.length = overlap.length))
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