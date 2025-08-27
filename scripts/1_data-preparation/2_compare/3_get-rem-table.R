### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
# suppressWarnings(suppressPackageStartupMessages(library(Biostrings)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

### define functions
compute_overlap_pairs <- function(DT, start_col, end_col, by_cols = c("qname","sname"), thr = 0.10) {
  # Prepare interval table
  tmp <- DT[, .(
    qname, sname, hsp_id = .I,
    start = get(start_col),
    end   = get(end_col)
  )]
  # HSP length (inclusive)
  tmp[, len := pmax(0L, end - start + 1L)]
  
  # Keys needed for foverlaps
  setkeyv(tmp, c(by_cols, "start", "end"))
  
  # Self-join to find overlaps within each (qname,sname)
  # This function returns every pair of rows whose overlaps overlap at all
  # (given the same qname/sname), including itself, but it does not compute
  # percentage overlaps; this is done below
  ov <- foverlaps(tmp, tmp, type = "any", nomatch = 0L)
  
  # Drop self-matches and keep unique unordered pairs (i.hsp_id < hsp_id)
  ov <- ov[i.hsp_id > hsp_id]
  
  # Compute overlap length
  ov[, ov_len := pmin(i.end, end) - pmax(i.start, start) + 1L]
  ov <- ov[ov_len > 0L] # redundant
  
  # Overlap fraction relative to shorter HSP
  ov[, ov.frac := ov_len / pmin(i.len, len)]
  
  # Return the (qname,sname) pairs that have any overlap >= threshold
  unique(ov[ov.frac >= thr, .(qname, sname, ov.frac)])
}

# Calculate total number of unique q/s pairs in a table
get.unique.pairs <- function(DT, self = FALSE){
  all.comb <- rbind(DT[, .(qname,sname)], DT[, .(qname = sname,sname = qname)])
  if(self){
    pairs.final <- unique(all.comb[qname <= sname])
  } else{
    pairs.final <- unique(all.comb[qname < sname])
  }
  return(pairs.final)
}


# -------- helper: swap q<->s columns for rows where we chose B->A --------
swap_qs_cols <- function(DT) {
  # Pairs of columns to swap if present
  pairs <- list(
    c("qname","sname"),
    c("qstart","sstart"),
    c("qend","send"),
    c("qlength","slength"),
    c("q.alnlen","s.alnlen"),
    c("qcov","scov"),
    c("qcov.full","scov.full")
  )
  for (pq in pairs) {
    if (all(pq %chin% names(DT))) {
      tmp <- DT[[pq[1]]]
      DT[, (pq[1]) := DT[[pq[2]]]]
      DT[, (pq[2]) := tmp]
    }
  }
  invisible(DT)
}

####################################################################
############################### MAIN ###############################
####################################################################

##########################
### 1. LOAD & PROCESS DATA
##########################

# Load refseq table filtered by NCBI
reprseq.table <- fread(reprseq.table.output)
all.reprseq <- unique(reprseq.table$reprseq)

# Load all-by-all blast & ensure integer coords
rvr.blast <- fread(reprseq.vs.reprseq.blast.output)
coord.cols <- intersect(c("qstart","qend","sstart","send"), names(rvr.blast))
rvr.blast[, (coord.cols) := lapply(.SD, as.integer), .SDcols = coord.cols]

# Filter the full blast table by pident of 30%
# rvr.cov <- fread(reprseq.vs.reprseq.cov.output)
# setnames(rvr.cov, c("query", "subject"), c("qname", "sname"))
# setkey(rvr.cov, qname, sname)

##############################################
### 2. FILTER BLAST TABLE BY 
### REMOVING PAIRS WITH OVERLAPPING HITS
### AND MORE THAN 2 ALIGNMENTS PER PAIR
##############################################

# Identify all duplicated pairs
duplicated.pairs <- rvr.blast[, .N, by = .(qname, sname)][N > 1L, .(qname, sname)] # dups_pairs
rvr.blast.duplicated <- rvr.blast[duplicated.pairs, on = .(qname, sname), nomatch = 0L]

# Overlap checks on QUERY and SUBJECT coordinates separately ---
overlapping.q.pairs <- compute_overlap_pairs(rvr.blast.duplicated, "qstart", "qend", thr = 0)
overlapping.s.pairs <- compute_overlap_pairs(rvr.blast.duplicated, "sstart", "send", thr = 0)

# Union of pairs to flag/remove if either side shows overlap
overlapping.pairs <- unique(rbind(overlapping.q.pairs[, .(qname, sname)], overlapping.s.pairs[, .(qname, sname)], use.names = TRUE, fill = TRUE))

# divide the rvr.blast.duplicated into:
# a) overlapping and N = 2 (problematic)
# b) overlapping and N > 2 (problematic)
# c) non-overlapping and N = 2 (non-problematic)
# d) non-overlapping and N > 2 (problematic)
rvr.blast.duplicated.overlapping <- rvr.blast.duplicated[overlapping.pairs, on = .(qname, sname)]
rvr.blast.duplicated.nonoverlapping <- rvr.blast.duplicated[!overlapping.pairs, on = .(qname, sname)]

pairs.duplicated.overlapping.ne2 <- rvr.blast.duplicated.overlapping[, .N, by = .(qname, sname)][N == 2L, .(qname, sname)]
pairs.duplicated.overlapping.ng2 <- rvr.blast.duplicated.overlapping[, .N, by = .(qname, sname)][N > 2L, .(qname, sname)]
pairs.duplicated.nonoverlapping.ne2 <- rvr.blast.duplicated.nonoverlapping[, .N, by = .(qname, sname)][N == 2L, .(qname, sname)]
pairs.duplicated.nonoverlapping.ng2 <- rvr.blast.duplicated.nonoverlapping[, .N, by = .(qname, sname)][N > 2L, .(qname, sname)]

pairs.problematic <- rbind(pairs.duplicated.overlapping.ne2, pairs.duplicated.overlapping.ng2, pairs.duplicated.nonoverlapping.ng2)
pairs.nonproblematic <- rbind(pairs.duplicated.nonoverlapping.ne2)

# filter blast hits by removing all problematic sequence pairs
rvr.blast.filtered <- rvr.blast[!pairs.problematic, on = .(qname, sname)]

# Report number of pairs for the methods
no.all.pairs <- nrow(get.all.pairs(rvr.blast))
no.duplicated.pairs <- nrow(get.all.pairs(rvr.blast.duplicated))
no.pairs.duplicated.overlapping.ne2 <- nrow(get.all.pairs(pairs.duplicated.overlapping.ne2))
no.pairs.duplicated.overlapping.ng2 <- nrow(get.all.pairs(pairs.duplicated.overlapping.ng2))
no.pairs.duplicated.nonoverlapping.ne2 <- nrow(get.all.pairs(pairs.duplicated.nonoverlapping.ne2))
no.pairs.duplicated.nonoverlapping.ng2 <- nrow(get.all.pairs(pairs.duplicated.nonoverlapping.ng2))
no.pairs.duplicated.overlapping <- no.pairs.duplicated.overlapping.ne2 + no.pairs.duplicated.overlapping.ng2
no.problematic.pairs <- no.pairs.duplicated.overlapping.ne2 + no.pairs.duplicated.overlapping.ng2 + no.pairs.duplicated.nonoverlapping.ng2

text4methods <- sprintf("MMseqs2 aligned %d pairs of which %d were found to be duplicated,
        and of those we removed %d pairs which had overlapping alignments
        within the sequence, and %d non-overlapping pairs with more than 2
        alignments per pair, altogether %d problematic pairs, yielding %d aligned
        pairs after filtering", 
        no.all.pairs, no.duplicated.pairs, 
        no.pairs.duplicated.overlapping, no.pairs.duplicated.nonoverlapping.ng2,
        no.problematic.pairs, no.all.pairs - no.problematic.pairs)
cat(text4methods)

# fwrite(rvr.blast.filtered, file = reprseq.vs.reprseq.blast.filtered.output)

##############################################
### 3. TRANSFORM BLAST (A/B AND B/A) TABLE INTO
### UNIQUE, ORDERED PAIRS (A/B)
### (written by GPT-5)
##############################################

# 0) Start from your already-filtered BLAST hits
BL <- copy(rvr.blast.filtered)

# 1) Drop self-hits and add unordered pair labels
BL <- BL[qname != sname]
BL[, `:=`(
  a   = pmin(qname, sname),          # unordered left
  b   = pmax(qname, sname),          # unordered right
  dir = fifelse(qname <= sname, "qs", "sq")  # which way this HSP was found
)]

# 2) Summarize per unordered pair & direction
dir_stats <- BL[, .(
  n_hsps       = .N,
  sum_bitscore = sum(bitscore, na.rm = TRUE),
  max_bitscore = max(bitscore, na.rm = TRUE)
), by = .(a, b, dir)]

# 3) Choose ONE direction per unordered pair:
#    - prefer larger sum_bitscore
#    - tie-break by larger max_bitscore
#    - final tie-break: keep "qs" for determinism
choice <- dir_stats[
  order(a, b, -sum_bitscore, -max_bitscore, dir)
][
  , .SD[1L], by = .(a, b)
][
  , .(a, b, dir_chosen = dir)
]

# 4) Keep ALL HSPs from the chosen direction
BL_kept <- BL[choice, on = .(a, b, dir = dir_chosen), nomatch = 0L]

# 5) Normalize so final rows always read q=a, s=b
#    (Swap columns for rows chosen from "sq")
sq_rows <- BL_kept[dir == "sq"]
if (nrow(sq_rows)) {
  swap_qs_cols(sq_rows)
  sq_rows[, dir := "qs"]  # after swap, they read as qs
}
qs_rows <- BL_kept[dir == "qs"]

rvr.unordered <- rbind(qs_rows, sq_rows, use.names = TRUE)
rvr.unordered[, `:=`(qname = a, sname = b)]  # enforce canonical labels
rvr.unordered[, c("a","b","dir") := NULL][]
rvr.unordered[, q.hit.length := qend - qstart + 1]
rvr.unordered[, s.hit.length := send - sstart + 1]
rvr.unordered[, pair.index.global := .GRP, by = .(qname, sname)]

fwrite(rvr.unordered, file = reprseq.vs.reprseq.blast.filtered.unordered.output)

##############################################
### 3. Filter BLAST table to identify REMs
##############################################

# Filter to create a table of putative recent recombinations
rem.table <- rvr.unordered %>%
  filter(q.hit.length >= RECENT.REC.MIN.Q.HIT.LENGTH &
           s.hit.length >= RECENT.REC.MIN.S.HIT.LENGTH &
           pident >= RECENT.REC.MIN.PIDENT &
           qcov.full <= RECENT.REC.MAX.QCOV &
           scov.full <= RECENT.REC.MAX.SCOV)

rem.table.final <- rem.table %>%
  select(qname, sname, pident, qstart, qend, qlength, sstart, send, slength,
         hit.length.query = q.hit.length, hit.length.subject = s.hit.length,
         qcov, scov, q.alnlen, s.alnlen, qcov.full, scov.full,
         pair.index.global)
setorder(rem.table.final, qname, sname)

rem.table.final <- copy(rem.table.final)  # or force a clean selfref
rem.table.final[, pair.index.rem := .GRP, by = .(qname, sname)]

# write to a file
fwrite(rem.table.final, file = rem.table.output)
