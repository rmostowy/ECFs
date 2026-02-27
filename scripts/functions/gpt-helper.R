reduce_ranges <- function(dt) {
  dt <- as.data.table(dt)               # ensure input is a data.table (copies if needed)
  
  # no rows or all NA → empty reduced set
  if (nrow(dt) == 0L || all(is.na(dt$qstart)) || all(is.na(dt$qend))) {
    return(data.table(qstart = integer(), qend = integer()))  # standardized empty return
  }
  
  # keep only valid ranges
  dt <- dt[!is.na(qstart) & !is.na(qend)]  # drop incomplete intervals
  if (nrow(dt) == 0L) return(data.table(qstart = integer(), qend = integer()))  # nothing left
  
  # enforce start <= end
  dt[qstart > qend, c("qstart","qend") := .(qend, qstart)]  # swap endpoints if reversed
  
  setorder(dt, qstart, qend)            # sort by start then end for linear merge pass
  
  out_start <- integer()                # accumulator for merged interval starts
  out_end   <- integer()                # accumulator for merged interval ends
  
  cur_start <- dt$qstart[1]             # current merged interval start
  cur_end   <- dt$qend[1]               # current merged interval end
  
  if (nrow(dt) > 1L) {                  # only loop if there are additional intervals
    for (i in 2:nrow(dt)) {             # scan intervals in sorted order
      if (dt$qstart[i] <= cur_end + 1L) {          # overlaps or is directly adjacent (merge)
        cur_end <- max(cur_end, dt$qend[i])        # extend current interval end if needed
      } else {                                     # disjoint interval (flush current)
        out_start <- c(out_start, cur_start)       # save finished merged interval
        out_end   <- c(out_end, cur_end)
        cur_start <- dt$qstart[i]                  # start new merged interval
        cur_end   <- dt$qend[i]
      }
    }
  }
  
  out_start <- c(out_start, cur_start)  # flush last interval after loop
  out_end   <- c(out_end, cur_end)
  
  data.table(qstart = out_start, qend = out_end)    # return merged, non-overlapping intervals
}

covered_length <- function(dt) {
  dt <- as.data.table(dt)               # ensure data.table
  if (nrow(dt) == 0L) return(0L)        # empty input covers 0 bp
  dt <- dt[!is.na(qstart) & !is.na(qend)]  # drop incomplete intervals
  if (nrow(dt) == 0L) return(0L)        # nothing valid left
  sum(dt$qend - dt$qstart + 1L)         # inclusive length per interval, summed
}

intersection_length <- function(a, b) {
  a <- as.data.table(a)                 # ensure data.table
  b <- as.data.table(b)
  
  a <- a[!is.na(qstart) & !is.na(qend)] # drop incomplete intervals
  b <- b[!is.na(qstart) & !is.na(qend)]
  
  if (nrow(a) == 0L || nrow(b) == 0L) return(0L)  # no overlap possible
  
  total <- 0L                         # accumulator for total intersecting length (bp)
  for (i in seq_len(nrow(a))) {       # iterate intervals in a
    for (j in seq_len(nrow(b))) {     # iterate intervals in b (O(|a|*|b|))
      start <- max(a$qstart[i], b$qstart[j])  # overlap start
      end   <- min(a$qend[i],   b$qend[j])    # overlap end
      if (start <= end) total <- total + (end - start + 1L)  # add inclusive overlap
    }
  }
  total                               # total intersecting length across all interval pairs
}