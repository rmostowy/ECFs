### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(Biostrings)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

### functions (written by GPT-5)

# Calculate the total length covered by a set of intervals,
# making sure that overlapping parts are not counted twice.
#
# Example: intervals [1–4], [3–6], [10–12]
#   - The union is [1–6] and [10–12]
#   - Total length = 6 (from 1–6) + 3 (from 10–12) = 9
#
# Important:
# - Intervals are inclusive (so [1–4] has length 4, not 3).
# - The input vectors `s` (starts) and `e` (ends) should be sorted
#   by start position before calling this function:
#   e.g., setorder(DT, qname, sname, qstart, qend)
# - If an interval has start > end, it will be swapped automatically.
# - If there are NAs, they are ignored. If all are NA, result is 0.
#
# Arguments:
#   s : vector of start coordinates
#   e : vector of end coordinates
# Returns:
#   An integer giving the total length covered by all intervals
#   without double-counting overlaps.

union_len_sorted <- function(s, e) {
  # Remove intervals with missing coordinates
  keep <- !(is.na(s) | is.na(e))
  if (!any(keep)) return(0L)
  
  # Work with numeric to avoid integer overflow/NA-coercion edge cases.
  s <- as.numeric(s[keep])
  e <- as.numeric(e[keep])
  
  # Fix intervals where start > end (swap them)
  swap <- e < s
  if (any(swap)) {
    tmp <- s[swap]
    s[swap] <- e[swap]
    e[swap] <- tmp
  }
  if (!length(s)) return(0L)
  
  # We assume intervals are already sorted by start
  # (done outside this function with setorder)
  
  # Keep track of the maximum end seen so far
  cm <- cummax(e)
  
  # For each interval, look at the maximum end of all previous intervals
  # (first one gets -Inf because there is nothing before it)
  prev_end <- c(-Inf, head(cm, -1L))  # lag(cummax(e)) without using shift()
  
  # For each interval, work out the new part that was not already covered:
  # - If the start is inside an earlier interval, shift it forward
  # - If the interval is fully inside earlier coverage, length = 0
  new_start <- pmax(s, prev_end + 1)
  seg_len   <- pmax(0, e - new_start + 1)
  
  # Add up all new contributions to get the total union length
  as.integer(sum(seg_len))
}

####################################################################
############################### MAIN ###############################
####################################################################

# Load refseq table filtered by NCBI
reprseq.table <- fread(reprseq.table.output)
all.reprseq <- unique(reprseq.table$reprseq)

# Load reprseq sequences
reprseq.seqs <- readAAStringSet(reprseq.seqs.filename)
mmseqs.tmp.dir <- "/Users/rmostowy/Data/tmp-ecfs/recent-rec/mmseqs_tmp"
mmseqs.db <- file.path(mmseqs.tmp.dir, "reprseqDB")
mmseqs.res <- file.path(mmseqs.tmp.dir, "reprseq_vs_reprseq")
mmseqs.res.m8 <- file.path(mmseqs.tmp.dir, "reprseq_vs_reprseq.m8")

# Ensure the tmp directory exists
dir.create(mmseqs.tmp.dir, showWarnings = FALSE, recursive = TRUE)

# Step 1: Create MMseqs2 Database
mmseqs.createdb.cmd <- paste("mmseqs createdb", reprseq.seqs.filename, mmseqs.db)
system(mmseqs.createdb.cmd)

# Step 2: Run MMseqs2 All-vs-All Search
mmseqs.search.cmd <- paste(
  "mmseqs search", mmseqs.db, mmseqs.db, mmseqs.res, mmseqs.tmp.dir,
  "-s 7.5 --alt-ali 10 --cov-mode 0 -c 0.0 -e 1e-3 --max-seqs 2000 --threads 20"
)
system(mmseqs.search.cmd)

# Step 3: Convert MMseqs2 Alignments to M8 Format (BLAST-like)
mmseqs.convert.cmd <- paste(
  "mmseqs convertalis", mmseqs.db, mmseqs.db, mmseqs.res, mmseqs.res.m8,
  "--format-output 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits,qcov,tcov'"
)
system(mmseqs.convert.cmd)

# Step 4: Load MMseqs2 results into a data.table
reprseq.vs.reprseq <- fread(mmseqs.res.m8, header = FALSE, sep = "\t")
# setorder(reprseq.vs.reprseq, qname, sname)

# Assign meaningful column names
setnames(reprseq.vs.reprseq, c(
  "qname", "sname", "pident", "alnlen", "mismatch", "gapopen",
  "qstart", "qend", "qlength", "sstart", "send", "slength",
  "evalue", "bitscore", "qcov", "scov"
))
# fwrite(reprseq.vs.reprseq, file = reprseq.vs.reprseq.blast.output)

# Step 5: Calculate full query and subject coverage values based on the MMseqs2 output

# 1) QUERY-side union: sort once by (qname, sname, qstart, qend)
setorder(reprseq.vs.reprseq, qname, sname, qstart, qend)
cov_q <- reprseq.vs.reprseq[
  , .(
    q_union = union_len_sorted(qstart, qend),
    qlen    = suppressWarnings(max(qlength, na.rm = TRUE))
  ),
  by = .(qname, sname)
]

# 2) SUBJECT-side union: sort once by (qname, sname, sstart, send)
setorder(reprseq.vs.reprseq, qname, sname, sstart, send)
cov_s <- reprseq.vs.reprseq[
  , .(
    s_union = union_len_sorted(sstart, send),
    slen    = suppressWarnings(max(slength, na.rm = TRUE))
  ),
  by = .(qname, sname)
]

# 3) Merge and compute fractions; attach to every row (ordered pair)
cov_full <- cov_q[cov_s, on = .(qname, sname)]
cov_full[, `:=`(
  q.alnlen = as.integer(q_union),
  s.alnlen = as.integer(s_union),
  qcov.full = fifelse(is.finite(qlen) & qlen > 0, q_union / qlen, NA_real_),
  scov.full = fifelse(is.finite(slen) & slen > 0, s_union / slen, NA_real_)
)]

# 4) Update the final blast table
reprseq.vs.reprseq.with.full.cov <- cov_full[
  reprseq.vs.reprseq, on = .(qname, sname)
][
  , .(qname, sname, pident, alnlen, mismatch, gapopen,
      qstart, qend, qlength, sstart, send, slength,
      evalue, bitscore, q.alnlen, s.alnlen, qcov, scov, qcov.full, scov.full)
]

# Step 6: Calculate full query and subject coverage values based on the MMseqs2 output

# How many decimals MMseqs2 shows for qcov/scov
digits_cov <- 3L
# Tolerance: differences smaller than this are just rounding noise
eps <- 0.5 * 10^(-digits_cov)   # e.g. 0.0005 for 3 decimals

DT <- reprseq.vs.reprseq.with.full.cov

# Round full coverages to MMseqs2 precision
DT[, `:=`(qcov.full = round(qcov.full, digits_cov),
          scov.full = round(scov.full, digits_cov))]

# If the rounded full < MMseqs2 value by more than epsilon,
# bump the full up to MMseqs2's value so (full - local) won't be negative
DT[, `:=`(
  qcov.full = fifelse(qcov.full + eps < qcov, qcov, qcov.full),
  scov.full = fifelse(scov.full + eps < scov, scov, scov.full)
)]

# Clamp to [0, 1] just in case
DT[, `:=`(
  qcov.full = pmin(1, pmax(0, qcov.full)),
  scov.full = pmin(1, pmax(0, scov.full))
)]

# Optional sanity check: after this, deltas shouldn't be negative beyond epsilon
# stopifnot(DT[, all(qcov.full - qcov >= -eps, na.rm = TRUE)])
# stopifnot(DT[, all(scov.full - scov >= -eps, na.rm = TRUE)])

reprseq.vs.reprseq.with.full.cov <- DT

# Step 7: Save results
fwrite(reprseq.vs.reprseq.with.full.cov, file = reprseq.vs.reprseq.blast.output)
