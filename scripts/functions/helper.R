# Function to calculate overlap details
calculate_overlap_parms <- function(start1, end1, start2, end2) {
  overlap.start <- pmax(start1, start2)
  overlap.end <- pmin(end1, end2)
  overlap.length <- pmax(0, overlap.end - overlap.start)
  return(list(overlap.start = overlap.start, overlap.end = overlap.end, overlap.length = overlap.length))
}
