### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

#################################
### Load all data:
### 1. ECF tables
### 2. reprseq lengths
### 3. reprseq annotations
### 4. reprseq vs ECF table
#################################

# ECF tables
ecf.table <- fread(ecf.table.output)
ecf.family.table <- fread(ecf.family.table.output)

# reprseq lengths
reprseq.lengths <- fread(reprseq.lengths.filename)

# Load functional annotations
reprseq.functions.unique <- fread(reprseq.functions.unique.output)
reprseq.functions <- fread(reprseq.functions.all.output)

# Load reprseq vs ECFs
reprseq.vs.ecfs.all <- fread(reprseq.ecfs.output)

#################################
### Process loaded tables
#################################

setkey(ecf.family.table, ecf)

setnames(reprseq.lengths, "name", "qname")
setkey(reprseq.lengths, qname)

# setkey(reprseq.functions, qname)
setkey(reprseq.functions.unique, reprseq)

reprseq.vs.ecfs <- reprseq.vs.ecfs.all %>% 
  filter(prob >= MIN.PROB.FOR.DOMAIN.DETECT.IN.RHMM,
         eval <= MAX.EVAL.FOR.DOMAIN.DETECT.IN.RHMM, 
         qcov >= MIN.RHMM.COV.FOR.DOMAIN.DETECT.IN.RHMM,
         scov >= MIN.DOMAIN.COV.FOR.DOMAIN.DETECT.IN.RHMM,
         hit.length >= MIN.HIT.LENGTH.FOR.DOMAIN.DETECT.IN.RHMM
  )
reprseq.vs.ecfs <- cbind(reprseq.vs.ecfs, family = ecf.family.table[.(reprseq.vs.ecfs$sname)]$family)
setkey(reprseq.vs.ecfs, qname)

###########################################
### Calculate protein architecture
### metrics from ECF/EF distributions:
###  1) coverage_fraction
###  2) n_families_total
###  3) mean_eff_families (Simpson effective #, per position)
###  4) segregation_index (1 - mean within-family overlap of position ranges)
###  5) n_dom_segments (dominant-family segmentation, optional purity filter)
###########################################

reprseq.ecf.hits <- copy(reprseq.vs.ecfs) #prev hits

# Expand to per-position table (heavy but intuitive)
# Each hit contributes 1 to each covered position (multiplicity kept).
reprseq.ecf.hits.by.pos <- reprseq.ecf.hits[, .(pos = qstart:qend, family), by = .(qname, qstart, qend, family)] # hit_pos -> reprseq.ecf.hits.by.pos
reprseq.ecf.hits.by.pos <- reprseq.ecf.hits.by.pos[, .(qname, pos, family)]  # drop interval cols early

# 1) covered.length  = (# positions covered by ≥1 hit)
covered.positions <- unique(reprseq.ecf.hits.by.pos[, .(qname, pos)]) # prev covered_pos
reprseq.covered.lengths <- covered.positions[, .(covered.length = .N), by = qname] # covered_len -> covered.length, cov -> reprseq.covered.lengths
setkey(reprseq.covered.lengths, qname)

# 2) n.fams = distinct families in protein (over all hits)
reprseq.families.total <- unique(reprseq.ecf.hits.by.pos[, .(qname, family)])[, .(n.fams = .N), by = qname] # n_families_total > n.fams, fam_tot -> reprseq.families.total
setkey(reprseq.families.total, qname)

# Position-wise family counts (multiplicity of hits)
# This is the base object reused for Simpson + dominance.
reprseq.family.counts <- reprseq.ecf.hits.by.pos[, .N, by = .(qname, pos, family)] # fam_counts -> reprseq.family.counts
setkey(reprseq.family.counts, qname, pos, family)

# also need n.fam.per.pos (number of unique families per position)
reprseq.fam.per.pos <- unique(reprseq.ecf.hits.by.pos)[, .(n.fam.per.pos = .N), by = .(qname, pos)] # n_fam_pos -> n.fam.per.pos, rich_pos -> reprseq.fam.per.pos
setkey(reprseq.fam.per.pos, qname, pos)

# 3) mean_eff_families (Simpson effective number) per position:
# ------------------------------------------------------------
# For each position: p_f = N_f / sum_f N_f
# Simpson effective = 1 / sum_f p_f^2
# Then mean over covered positions.
# ------------------------------------------------------------
reprseq.family.counts[, tot := sum(N), by = .(qname, pos)]
reprseq.family.counts[, p := N / tot]
reprseq.simpson.pos <- reprseq.family.counts[, .(
  eff = 1 / sum(p^2)
), by = .(qname, pos)] # simpson_pos -> reprseq.simpson.pos
setkey(reprseq.simpson.pos, qname, pos)

reprseq.eff.mean <- reprseq.simpson.pos[, .(
  mean.eff.fam = mean(eff) 
), by = qname] # mean_eff_families -> mean.eff.fam, eff_mean -> reprseq.eff.mean
setkey(reprseq.eff.mean, qname)

# 4a) segregation.overlap (segregation overlap index (0..1))
# ------------------------------------------------------------
# Intuition: families are "segregated" if their position ranges overlap little.
# Steps:
#   - For each protein-family, compute range [min(pos), max(pos)]
#   - For protein, compute average pairwise overlap fraction between family ranges
#   - segregation_index = 1 - mean_overlap
# Notes:
#   - Uses ranges (not full shapes) => robust + simple.
#   - Requires ≥2 families; else segregation_index = NA (or 1; NA is safer).
# ------------------------------------------------------------
reprseq.fam.ranges <- unique(reprseq.ecf.hits.by.pos)[, .(
  lo = min(pos),
  hi = max(pos)
), by = .(qname, family)] # fam_ranges -> reprseq.fam.ranges
setkey(reprseq.fam.ranges, qname, family)

# pairwise overlap helper for one protein
range.overlap.mean <- function(lo, hi){
  k <- length(lo)
  if (k < 2) return(NA_real_)
  ij <- combn(seq_len(k), 2)
  a <- ij[1,]; b <- ij[2,]
  inter <- pmax(0L, pmin(hi[a], hi[b]) - pmax(lo[a], lo[b]) + 1L)
  len_a <- hi[a] - lo[a] + 1L
  len_b <- hi[b] - lo[b] + 1L
  frac <- inter / pmin(len_a, len_b)   # 0..1
  mean(frac)
} # range_overlap_mean -> range.overlap.mean

# compute only where defined
reprseq.seg.overlap <- reprseq.fam.ranges[, .(
  n.fam = .N,
  mean.range.overlap = range.overlap.mean(lo, hi)
), by = qname] # n_fam -> n.fam, mean_range_overlap -> mean.range.overlap, seg -> reprseq.seg.overlap

reprseq.seg.overlap[, segregation.overlap := fifelse(
  n.fam >= 2,
  1 - mean.range.overlap,
  NA_real_
)] # segregation_index -> segregation.overlap

reprseq.seg.overlap[, c("n.fam","mean.range.overlap") := NULL]
setkey(reprseq.seg.overlap, qname)

# 4b) segregation.bary (segregation bary index (0..1))
# ------------------------------------------------------------
# Intuition: families are "segregated" if their *mass centers* (barycenters)
# are far apart along the protein, taking into account hit multiplicity.
# Steps:
#   - For each protein-position-family, count multiplicity N
#   - For each protein-family, compute barycenter mu = sum(pos*N)/sum(N)
#   - For protein, compute mean pairwise |mu_i - mu_j| / qlength
# Notes:
#   - Requires ≥2 families; else segregation_bary = NA
#   - More sensitive than range-based segregation_index to "dominant N-term vs C-term"
# ------------------------------------------------------------

# multiplicity per (qname,pos,family)
reprseq.pos.mass <- reprseq.ecf.hits.by.pos[, .N, by = .(qname, pos, family)] # fam_pos_mass -> reprseq.pos.mass

# family barycenters per protein
reprseq.bary.centers <- reprseq.pos.mass[, .(
  mu = sum(pos * N) / sum(N),
  total.mass = sum(N)
), by = .(qname, family)] # fam_centers -> reprseq.bary.centers, total_mass -> total.mass
setkey(reprseq.bary.centers, qname) 

# pairwise mean distance helper (expects vectors mu and scalar L)
mu.dist.mean <- function(mu, L){
  k <- length(mu)
  if (k < 2 || is.na(L) || L <= 0) return(NA_real_)
  ij <- combn(mu, 2)
  mean(abs(ij[1, ] - ij[2, ]) / L)
} # mu_dist_mean -> mu.dist.mean

# compute barycentric segregation per protein
reprseq.seg.bary <- reprseq.bary.centers[
  reprseq.lengths, on = "qname"
][, .(
  n.fam = .N,
  segregation.bary = mu.dist.mean(mu, length[1L])
), by = qname] # seg_bary -> reprseq.seg.bary, segregation_bary -> segregation.bary

reprseq.seg.bary[n.fam < 2, segregation.bary := NA_real_]
reprseq.seg.bary[, n.fam := NULL]
setkey(reprseq.seg.bary, qname)

# 5) n.dom.segments
# ------------------------------------------------------------
# Dominant family per position = family with max N at that position.
# Then count contiguous segments along pos among positions that are "valid".
# IMPORTANT:
# A block is a region of the protein where one evolutionary signal is stably dominant
# Hence if there is ambiguity (impurity) in the middle of a single family block,
# We define two building blocks
#
# Optional: purity filter to ignore ambiguous positions:
#   purity = max(p) at position
#   if purity < PURITY_MIN => mark as NA (doesn't contribute to segments)
# ------------------------------------------------------------
# dominant family per (qname,pos)
setorder(reprseq.family.counts, qname, pos, -N, family)  # tie-break by family for stability
dom <- reprseq.family.counts[, .SD[1L], by = .(qname, pos)]  # first row per pos = max N
dom[, purity := p]  # since p computed above from N/tot
dom[, dom.family := family] # dom -> reprseq.n.seg
dom <- dom[, .(qname, pos, dom.family, purity)]
setkey(dom, qname, pos)

PURITY_MIN <- 0.6  # tweak; set to 0 to disable
if (PURITY_MIN > 0) dom[purity < PURITY_MIN, dom.family := NA_character_]

# count segments (runs) in dom_family track, ignoring NA positions
setorder(dom, qname, pos)
dom[, prev := shift(dom.family), by = qname]
dom[, new.seg := !is.na(dom.family) & (is.na(prev) | dom.family != prev)]
reprseq.n.seg <- dom[, .(n.dom.segments = sum(new.seg)), by = qname] # segN -> reprseq.n.seg, n_dom_segments -> n.dom.segments
setkey(reprseq.n.seg, qname)

# 6) Number of ECFs per reprseq

reprseq.ecf.counts <- reprseq.ecf.hits[, .(
  n.ecf.hits = uniqueN(sname)
), by = qname] # ecf_counts -> reprseq.ecf.counts

# Assemble final table
pm <- reprseq.lengths[reprseq.covered.lengths, on="qname"]
pm[is.na(covered.length), covered.length := 0L]
pm[, coverage.fraction := fifelse(length > 0, covered.length / length, NA_real_)]

pm <- reprseq.families.total[pm, on="qname"]
pm[is.na(n.fams), n.fams := 0L]

pm <- reprseq.eff.mean[pm, on="qname"]
pm[covered.length == 0L, mean.eff.fam := NA_real_]

pm <- reprseq.seg.overlap[pm, on="qname"]
pm[n.fams < 2, segregation.overlap := NA_real_]

pm <- reprseq.seg.bary[pm, on="qname"]
pm[n.fams < 2, segregation.bary := NA_real_]

pm <- reprseq.n.seg[pm, on="qname"]
pm[is.na(n.dom.segments), n.dom.segments := 0L]

pm <- reprseq.ecf.counts[pm, on="qname"]
pm[is.na(n.ecf.hits), n.ecf.hits := 0L]

protein.metrics <- pm[, .(
  qname,
  length,
  n.ecf.hits,
  coverage.fraction,
  n.fams,
  mean.eff.fam,
  segregation.overlap,
  segregation.bary,
  n.dom.segments
)]

# # Append functional information
# pm.funct <- reprseq.functions.unique[.(protein.metrics$qname)][, .(qname = reprseq, funct, class)]
# protein.metrics <- merge(protein.metrics, pm.funct, by = "qname")
# protein.metrics[class == "", class := "unknown"]
# protein.metrics[funct == "", funct := "unknown"]

# Write to a file
fwrite(protein.metrics, file = reprseq.arch.metrics.output)

