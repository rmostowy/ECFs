# Function to calculate overlap and coverage
calculate_overlap_and_coverage <- function(start1, end1, start2, end2) {
  overlap_length <- pmax(0, pmin(end1, end2) - pmax(start1, start2))
  coverage <- overlap_length / (end1 - start1)
  return(list(overlap_length = overlap_length, coverage = coverage))
}
