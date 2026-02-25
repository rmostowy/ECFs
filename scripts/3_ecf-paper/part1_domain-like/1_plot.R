### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(Biostrings)))
suppressWarnings(suppressPackageStartupMessages(library(scales)))
suppressWarnings(suppressPackageStartupMessages(library(patchwork)))
# suppressWarnings(suppressPackageStartupMessages(library(tidyr)))
# suppressWarnings(suppressPackageStartupMessages(library(seqinr)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

#################################
### Set up output directory
#################################
this.part.dir <- file.path(visual.dir, "part1_domain-like")
make.dir(this.part.dir)

#################################
### Load all data:
### 1. domain tables (ECF/Pfam/ECOD)
### 2. reprseq vs domains (ECF/Pfam/ECOD)
#################################

# Load ECF tables
ecf.table <- fread(ecf.table.output)
ecf.family.table <- fread(ecf.family.table.output)

# Load Pfam table
pfam.domains <- fread(pfam.domains.output)
pfam.domains[, length := as.numeric(length)]

# Load ECOD table
ecod.domains <- fread(ecod.domains.output)

# Load reprseqs vs. domains
reprseq.vs.ecfs <- fread(reprseq.ecfs.output)
reprseq.vs.ecod <- fread(reprseq.domains.ecod.output)
reprseq.vs.pfam <- fread(reprseq.domains.pfam.output)

# Load reprseq sequences
reprseq.seqs <- readAAStringSet(reprseq.seqs.filename)


#################################
### Process loaded tables
#################################

# Filter
reprseq.vs.ecfs <- reprseq.vs.ecfs %>% 
  filter(qcov >= MIN.QCOV.FOR.DOMAIN, 
         scov >= MIN.SCOV.FOR.DOMAIN,
         eval <= MAX.EVAL.FOR.DOMAIN)

reprseq.vs.ecod <- reprseq.vs.ecod %>% 
  filter(qcov >= MIN.QCOV.FOR.DOMAIN, 
         scov >= MIN.SCOV.FOR.DOMAIN,
         eval <= MAX.EVAL.FOR.DOMAIN)

reprseq.vs.pfam <- reprseq.vs.pfam %>% 
  filter(qcov >= MIN.QCOV.FOR.DOMAIN, 
         scov >= MIN.SCOV.FOR.DOMAIN,
         eval <= MAX.EVAL.FOR.DOMAIN)

# Other
setkey(ecf.family.table, ecf)
reprseq.vs.ecfs$EF <- ecf.family.table[.(reprseq.vs.ecfs$sname)]$family

reprseq.seqs.n <- names(reprseq.seqs)
reprseq.seqs.l <- width(reprseq.seqs)
names(reprseq.seqs.l) <- reprseq.seqs.n

reprseqs.with.domains <- unique(c(reprseq.vs.ecfs$qname, reprseq.vs.ecod$qname, reprseq.vs.pfam$qname))


#################################
### Analysis 1: 
### ECF size
#################################

mean.val <- mean(ecf.table$size, na.rm = TRUE)
median.val <- median(ecf.table$size, na.rm = TRUE)

p.ecf.size <- ggplot(ecf.table, aes(x = size)) +
  geom_histogram(bins = 40, fill = "#4C78A8", color = "NA", alpha = 0.85) + 
  geom_vline(xintercept = mean.val, color = "#D62728", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = median.val, color = "#D62728", linetype = "dotted", linewidth = 1) +
  scale_x_log10() + 
  labs(
    x = "ECF size",
    y = "Count"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text  = element_text(color = "black"),
    panel.grid = element_blank()
  )
p.ecf.size

p.ecf.size.outfile <- file.path(this.part.dir, "plot1_ecf-size.png")
ggsave(p.ecf.size.outfile, p.ecf.size, width = 4, height = 3)

#################################
### Analysis 2: 
### Number of reprseqs 
### detected in ECF
#################################

reprseq.per.ecf <- reprseq.vs.ecfs[, .(no.reprseq = uniqueN(qname)), by = sname]
mean.val <- mean(reprseq.per.ecf$no.reprseq, na.rm = TRUE)
median.val <- median(reprseq.per.ecf$no.reprseq, na.rm = TRUE)

p.repr.per.ecf <- ggplot(reprseq.per.ecf, aes(x = no.reprseq)) +
  geom_histogram(bins = 40, fill = "#4C78A8", color = "NA", alpha = 0.85) + 
  geom_vline(xintercept = mean.val, color = "#D62728", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = median.val, color = "#D62728", linetype = "dotted", linewidth = 1) +
  scale_x_log10() + 
  labs(
    x = "Number of rHMMs per ECF",
    y = "Count"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text  = element_text(color = "black"),
    panel.grid = element_blank()
  )
p.repr.per.ecf

p.repr.per.ecf.outfile <- file.path(this.part.dir, "plot2_no-rhmms-per-ecf.png")
ggsave(p.repr.per.ecf.outfile, p.repr.per.ecf, width = 4, height = 3)

#################################
### Analysis 3: 
### No of ECFs per EF
#################################

# Rank the families by size
ecfs.per.family <- ecf.family.table[, .(no.ecfs = uniqueN(ecf)), by = family]
# Define bins for the number of ECFs
ecfs.per.family[, bin := cut(
  no.ecfs,
  breaks = c(-Inf, 1, 5, 10, 50, 100, 500, Inf),  # Define bin edges
  labels = c("1", "2–5", "6–10", "11–50", "51-100", "101-500", ">500"),
  right = TRUE
)]

# Count the number of families in each bin
family_counts <- ecfs.per.family[, .N, by = bin]
setorder(family_counts, bin)  # Sort by bin order

# Plot the binned data
p.no.ecf.per.ef <- ggplot(family_counts, aes(x = bin, y = N)) +
  geom_bar(stat = "identity", fill = "#4C78A8", alpha = 0.8) +
  geom_text(aes(label = N), vjust = -0.5, size = 3) +
  labs(
    x = "ECF count",
    y = "Number of EFs with this ECF count"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text  = element_text(color = "black"),
    panel.grid = element_blank()
  )

p.no.ecf.per.ef.outfile <- file.path(this.part.dir, "plot3_no-ecf-per-ef.png")
ggsave(p.no.ecf.per.ef.outfile, p.no.ecf.per.ef, width = 5, height = 4)

#################################
### Analysis 4: 
### Compare the distribution 
### of domain lengths 
### (ECF/Pfam/ECOD)
#################################

# Make reprseq lengths table and filter it
reprseq.lengths.all <- data.table(name = reprseq.seqs.n, length = reprseq.seqs.l)
reprseq.lengths.all[, db := "rHMMs"]
reprseq.lengths <- reprseq.lengths.all[name %in% reprseqs.with.domains]

# Make ecf lengths table
ecf.lengths <- ecf.table %>% select(name = ecf.name, length)
ecf.lengths[, db := "ECF"]

# Make ECOD lengths table (no filtering)
ecod.lengths <- ecod.domains %>% select(name = full.id, length)
ecod.lengths[, db := "ECOD"]

# Make PFAM lengths table (no filtering)
pfam.lengths <- pfam.domains %>% select(name = accession, length)
pfam.lengths[, db := "PFAM"]

# Combine all lengths
all.lengths <- rbind(reprseq.lengths, ecf.lengths, ecod.lengths, pfam.lengths)

# Calculate medians for each database
medians <- all.lengths[, .(median_length = median(length)), by = db]
medians$db <- factor(medians$db, levels = unique(all.lengths$db))
means <- all.lengths[, .(mean_length = mean(length)), by = db]
means$db <- factor(means$db, levels = unique(all.lengths$db))


all.lengths$db <- factor(all.lengths$db, levels = unique(all.lengths$db))
# Create the plot
p.length.distr <- ggplot(all.lengths, aes(x = length, fill = db)) +
  geom_histogram(aes(y = after_stat(density)), position = "identity", alpha = 0.6, bins = 80) +  # Histogram with transparency
  # geom_vline(data = medians, aes(xintercept = median_length, color = db),
  #            linetype = "dashed", size = 1) +  # Median line
  geom_vline(data = means, aes(xintercept = mean_length, color = db), 
             linetype = "dashed", linewidth = 1) +  # Median line
  geom_vline(data = medians, aes(xintercept = median_length, color = db), 
             linetype = "dotted", linewidth = 1) +  # Median line
  facet_grid(db ~ ., scales = "free_x") +  # Align distributions vertically
  scale_fill_manual(values = c("rHMMs" = "red", "ECF" = "#1f77b4", "PFAM" = "#ff7f0e", "ECOD" = "#2ca02c")) +  # Custom colors
  scale_color_manual(values = c("rHMMs" = "red", "ECF" = "#1f77b4", "PFAM" = "#ff7f0e", "ECOD" = "#2ca02c")) +  # Matching colors
  scale_x_log10() +
  labs(
    x = "Length",
    y = "Density",
    fill = "Database"
  ) +
  theme_bw(base_size = 14) +
  guides(color = "none", fill = "none") +
  theme(
    legend.position = NULL,
    strip.text = element_text(face = "bold")
  )
# p.length.distr

p.length.distr.outfile <- file.path(this.part.dir, "plot4_len-distr.png")
ggsave(p.length.distr.outfile, p.length.distr, width = 4, height = 7)


#################################
### Analysis 5: 
### domain rel position heatmap
#################################

# To assess whether ECFs occupy biologically plausible modular positions within proteins, 
# we mapped their normalized start and end coordinates across all representative sequences. 
# ECFs exhibit structured positional enrichment, with strong accumulation at N- and C-terminal 
# regions as well as broad internal distributions spanning substantial fractions of proteins. 
# Notably, this positional architecture closely mirrors that observed for curated Pfam domains, 
# indicating that ECFs are not randomly distributed fragments but instead behave as contiguous 
# modular units comparable to established protein domains.

# Compute relative positions in a new table
ecf.rel.pos <- copy(reprseq.vs.ecfs[prob >= 95])
ecf.rel.pos[, `:=`(
  rel_start = qstart / qlength,
  rel_end   = qend   / qlength,
  rel_len   = (qend - qstart + 1) / qlength
)]

pfam.rel.pos <- copy(reprseq.vs.pfam[prob >= 95])
pfam.rel.pos[, `:=`(
  rel_start = qstart / qlength,
  rel_end   = qend   / qlength,
  rel_len   = (qend - qstart + 1) / qlength
)]

ecod.rel.pos <- copy(reprseq.vs.ecod[prob >= 95])
ecod.rel.pos[, `:=`(
  rel_start = qstart / qlength,
  rel_end   = qend   / qlength,
  rel_len   = (qend - qstart + 1) / qlength
)]

# Define plotting function
plot.rel.pos <- function(dt, no.bins = 200, title = "2D interval map of ECF occurrences"){
  b <- ggplot_build(ggplot(dt, aes(rel_start, rel_end)) + geom_bin2d(bins = no.bins))
  upper <- quantile(b$data[[1]]$count, 0.995)  # try 0.99, 0.995, 0.999
  
  plot <- ggplot(dt, aes(rel_start, rel_end)) +
    geom_bin2d(bins = 200) +
    geom_abline(slope = 1, intercept = 0) +
    scale_fill_viridis_c(
      option = "magma",
      trans = "sqrt",
      limits = c(0, upper),
      oob = squish,
      labels = comma
    ) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    labs(
      x = "Relative start position (qstart / protein length)",
      y = "Relative end position (qend / protein length)",
      fill = "Count",
      title = "2D interval map of ECF occurrences"
    ) +
    theme_bw(base_size = 14) +
    theme(
      strip.text = element_text(face = "bold")
    )
  
  plot
}

rel.pos.ecfs <- plot.rel.pos(ecf.rel.pos)
rel.pos.pfam <- plot.rel.pos(pfam.rel.pos, title = "2D interval map of PFAM occurrences")
rel.pos.ecod <- plot.rel.pos(ecod.rel.pos, title = "2D interval map of ECOD occurrences")

rel.pos.pfam.outfile <- file.path(this.part.dir, "plot6_pfam-rel-pos.png")
ggsave(rel.pos.pfam.outfile, rel.pos.pfam, width = 6, height = 5)

rel.pos.ecod.outfile <- file.path(this.part.dir, "plot7_ecod-rel-pos.png")
ggsave(rel.pos.ecod.outfile, rel.pos.ecod, width = 6, height = 5)


#################################
### Analysis 6: 
### no of ECFs/EFs per reprseq
#################################

# Produce the ECF/EF count per reprseq table
ecf.per.reprseq.dt <- reprseq.vs.ecfs[prob>95, .(qlength = qlength[1],
                                          no.ecfs = uniqueN(sname)), 
                                      by = qname]
ef.per.reprseq.dt <- reprseq.vs.ecfs[prob>95, .(qlength = qlength[1],
                                          no.efs = uniqueN(EF)), 
                                      by = qname]
no.ecf.ef.per.reprseq.dt <- merge.data.table(ecf.per.reprseq.dt, ef.per.reprseq.dt, by = c("qname", "qlength"))

# 5 quantile bins (0.2 increments). If you want 10 bins, change by=0.1 and labels accordingly.
no.ecf.ef.per.reprseq.dt[, length_bin := cut(
  qlength,
  breaks = quantile(qlength, probs = seq(0, 1, by = 0.2), na.rm = TRUE),
  include.lowest = TRUE
)]

# --- Create efs_cat FROM no.efs (this was missing in your code)
no.ecf.ef.per.reprseq.dt[, efs_cat := cut(
  no.efs,
  breaks = c(0, 1, 2, 3, 5, 10, Inf),
  labels = c("1", "2", "3", "4–5", "6–10", "≥11"),
  right = TRUE,
  include.lowest = TRUE
)]

# Aggregate to get proportions per length bin
prop_dt <- no.ecf.ef.per.reprseq.dt[, .N, by = .(length_bin, efs_cat)]
prop_dt[, prop := N / sum(N), by = length_bin]

# Q1..Q5 labels for x axis
bin_levels <- levels(no.ecf.ef.per.reprseq.dt$length_bin)
prop_dt[, length_bin_q := factor(length_bin, levels = bin_levels,
                                 labels = paste0("Q", seq_along(bin_levels)))]

# Bin ranges for annotation on top
bin_ranges <- no.ecf.ef.per.reprseq.dt[, .(
  min_len = min(qlength, na.rm = TRUE),
  max_len = max(qlength, na.rm = TRUE)
), by = length_bin]

bin_ranges[, length_bin_q := factor(length_bin, levels = bin_levels,
                                    labels = paste0("Q", seq_along(bin_levels)))]
bin_ranges[, range_lab := paste0(comma(min_len), "–", comma(max_len), " aa")]
bin_ranges[, y := 1.02]

# Okabe–Ito palette
okabe_ito <- c(
  "1"   = "#56B4E9",
  "2"   = "#E69F00",
  "3"   = "#009E73",
  "4–5" = "#F0E442",
  "6–10"= "#0072B2",
  "≥11" = "#CC79A7"
)

# Plot
p.ef.per.reprseq <- ggplot(prop_dt, aes(length_bin_q, prop, fill = efs_cat)) +
  geom_col(width = 0.85, color = "grey20", linewidth = 0.2) +
  geom_text(
    data = bin_ranges,
    aes(x = length_bin_q, y = y, label = range_lab),
    inherit.aes = FALSE,
    size = 2.5,
    vjust = 0
  ) +
  scale_fill_manual(values = okabe_ito, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.04))) +
  labs(
    x = "Protein length (quantile bins)",
    y = "Proteins (proportion)",
    fill = "Number of EFs\nper rHMM"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

p.ef.per.reprseq.outfile <- file.path(this.part.dir, "plot8_efs-per-rhmm.png")
ggsave(p.ef.per.reprseq.outfile, p.ef.per.reprseq, width = 7, height = 5)

#################################
### Analysis 7: 
### ECF vs. domain overlap
#################################

reduce_ranges <- function(dt) {
  dt <- as.data.table(dt)
  
  # no rows or all NA → empty reduced set
  if (nrow(dt) == 0L || all(is.na(dt$qstart)) || all(is.na(dt$qend))) {
    return(data.table(qstart = integer(), qend = integer()))
  }
  
  # keep only valid ranges
  dt <- dt[!is.na(qstart) & !is.na(qend)]
  if (nrow(dt) == 0L) return(data.table(qstart = integer(), qend = integer()))
  
  # enforce start <= end
  dt[qstart > qend, c("qstart","qend") := .(qend, qstart)]
  
  setorder(dt, qstart, qend)
  
  out_start <- integer()
  out_end   <- integer()
  
  cur_start <- dt$qstart[1]
  cur_end   <- dt$qend[1]
  
  if (nrow(dt) > 1L) {
    for (i in 2:nrow(dt)) {
      if (dt$qstart[i] <= cur_end + 1L) {
        cur_end <- max(cur_end, dt$qend[i])
      } else {
        out_start <- c(out_start, cur_start)
        out_end   <- c(out_end, cur_end)
        cur_start <- dt$qstart[i]
        cur_end   <- dt$qend[i]
      }
    }
  }
  
  out_start <- c(out_start, cur_start)
  out_end   <- c(out_end, cur_end)
  
  data.table(qstart = out_start, qend = out_end)
}

covered_length <- function(dt) {
  dt <- as.data.table(dt)
  if (nrow(dt) == 0L) return(0L)
  dt <- dt[!is.na(qstart) & !is.na(qend)]
  if (nrow(dt) == 0L) return(0L)
  sum(dt$qend - dt$qstart + 1L)
}

intersection_length <- function(a, b) {
  a <- as.data.table(a)
  b <- as.data.table(b)
  
  a <- a[!is.na(qstart) & !is.na(qend)]
  b <- b[!is.na(qstart) & !is.na(qend)]
  
  if (nrow(a) == 0L || nrow(b) == 0L) return(0L)
  
  total <- 0L
  for (i in seq_len(nrow(a))) {
    for (j in seq_len(nrow(b))) {
      start <- max(a$qstart[i], b$qstart[j])
      end   <- min(a$qend[i],   b$qend[j])
      if (start <= end) total <- total + (end - start + 1L)
    }
  }
  total
}
# 
# # copy objects
# rhmm.vs.ecf <- copy(reprseq.vs.ecfs)
# rhmm.vs.pfam <- copy(reprseq.vs.pfam)
# rhmm.vs.ecod <- copy(reprseq.vs.ecod)
# 
# # index
# setkey(rhmm.vs.ecf, qname)
# setkey(rhmm.vs.pfam, qname)
# setkey(rhmm.vs.ecod, qname)
# 
# # 
# all.rhmms <- sort(unique(
#   c(rhmm.vs.pfam$qname, rhmm.vs.ecod$qname)
# ))
# 
# 
# conditional.probs.list <- lapply(1:1000, function(ri){
#   this.rhmm <- all.rhmms[ri]
#   this.rhmm.length <- as.numeric(reprseq.seqs.l[this.rhmm])
#   
#   # get unique ECF ranges
#   ecf.ranges <- rhmm.vs.ecf[.(this.rhmm)][, .(qstart, qend)]
#   ecf.red <- reduce_ranges(ecf.ranges)
#   
#   # get unique PFAM ranges
#   pfam.ranges <- rhmm.vs.pfam[.(this.rhmm)][, .(qstart, qend)]
#   pfam.red <- reduce_ranges(pfam.ranges)
#   
#   # get unique ECOD ranges
#   ecod.ranges <- rhmm.vs.ecod[.(this.rhmm)][, .(qstart, qend)]
#   ecod.red <- reduce_ranges(ecod.ranges)
#   
#   # get unique domain ranges
#   domain.red <- reduce_ranges(rbind(pfam.red, ecod.red))
#   
#   # calculate conditional probabilities
#   L <- this.rhmm.length
#   res.covered.by.both <- intersection_length(ecf.red, domain.red)
#   res.covered.by.ecf <- covered_length(ecf.red)
#   res.covered.by.domain <- covered_length(domain.red)
#   
#   res.covered.by.only.ecf <- res.covered.by.ecf - res.covered.by.both
#   res.covered.by.only.domain <- res.covered.by.domain - res.covered.by.both
#   res.covered.by.none <- this.rhmm.length - (res.covered.by.both + res.covered.by.only.ecf + res.covered.by.only.domain)
#   
#   p.ecf.given.domain <- res.covered.by.both / (res.covered.by.both + res.covered.by.only.domain)
#   if(res.covered.by.both + res.covered.by.only.domain == 0) p.ecf.given.domain <- NA
#   p.ecf.given.no.domain <- res.covered.by.only.ecf / (res.covered.by.only.ecf + res.covered.by.none)
#   if(res.covered.by.only.ecf + res.covered.by.none == 0) p.ecf.given.domain <- NA
#   
#   data.table(p.ecf.given.domain, p.ecf.given.no.domain)
#   
# })


####################################################################
# FASTER

# ---- 1) Precompute per-qname range lists once (fast) ----
mk_ranges_list <- function(dt) {
  dt <- as.data.table(dt)[, .(qname, qstart, qend)]
  dt[, .(ranges = list(.SD)), by = qname, .SDcols = c("qstart","qend")]
}

ecfL  <- mk_ranges_list(reprseq.vs.ecfs)
pfamL <- mk_ranges_list(reprseq.vs.pfam)
ecodL <- mk_ranges_list(reprseq.vs.ecod)

# Join into one table with list-columns (missing => NULL)
allL <- Reduce(function(x, y) merge(x, y, by="qname", all=TRUE),
               list(
                 setnames(ecfL,  "ranges", "ecf_ranges"),
                 setnames(pfamL, "ranges", "pfam_ranges"),
                 setnames(ecodL, "ranges", "ecod_ranges")
               ))


# If you want exactly those with PFAM or ECOD present:
allL <- allL[ lengths(pfam_ranges) > 0 | lengths(ecod_ranges) > 0 ]

# Add lengths (vectorized)
allL[, qlength := as.integer(reprseq.seqs.l[qname])]


# Add lengths (vectorized)
allL[, qlength := as.integer(reprseq.seqs.l[qname])]

# Optionally subset first 1000
# allL <- allL[1:1000]

conditional.probs <- allL[, {
  # pull range tables (NULL-safe)
  ecf.red  <- reduce_ranges(if (is.null(ecf_ranges[[1]]))  data.table(qstart=NA_integer_, qend=NA_integer_) else ecf_ranges[[1]])
  pfam.red <- reduce_ranges(if (is.null(pfam_ranges[[1]])) data.table(qstart=NA_integer_, qend=NA_integer_) else pfam_ranges[[1]])
  ecod.red <- reduce_ranges(if (is.null(ecod_ranges[[1]])) data.table(qstart=NA_integer_, qend=NA_integer_) else ecod_ranges[[1]])
  
  domain.red <- reduce_ranges(rbind(pfam.red, ecod.red))
  
  both   <- intersection_length(ecf.red, domain.red)
  ecf    <- covered_length(ecf.red)
  domain <- covered_length(domain.red)
  
  only_ecf    <- ecf - both
  only_domain <- domain - both
  none        <- qlength - (both + only_ecf + only_domain)
  
  denom1 <- both + only_domain              # = domain coverage
  denom2 <- only_ecf + none                 # = non-domain coverage
  
  p_ecf_given_domain    <- if (denom1 > 0) both / denom1 else NA_real_
  p_ecf_given_no_domain <- if (denom2 > 0) only_ecf / denom2 else NA_real_
  
  .(p.ecf.given.domain = p_ecf_given_domain,
    p.ecf.given.no.domain = p_ecf_given_no_domain,
    covered.ecf = ecf,
    covered.domain = domain,
    covered.both = both)
}, by = qname]

conditional.probs.dt <- copy(conditional.probs)
conditional.probs.dt[, delta := p.ecf.given.domain - p.ecf.given.no.domain]

# Show that ECF coverage is higher inside domain-annotated regions
p.domain.vs.ecf.cov <- ggplot(dt, aes(x = delta)) +
  geom_histogram(bins = 80) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  labs(x = expression(Delta*" = P(ECF|domain) - P(ECF|no domain)"),
       y = "Number of proteins"
       ) +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

p.domain.vs.ecf.cov.outfile <- file.path(this.part.dir, "plot9_ecf-cov-dom-reg.png")
ggsave(p.domain.vs.ecf.cov.outfile, p.domain.vs.ecf.cov, width = 6, height = 5)


####################################################

# p.length.distr # panel A
# p.ef.per.reprseq # panel B
# p.domain.vs.ecf.cov # panel C

fig2 <- 
  ( p.length.distr |
  (p.ef.per.reprseq / p.domain.vs.ecf.cov)
  )+
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(tag_levels = "A") & 
  theme(
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0.02, 0.98)
  )

fig2.outfile <- file.path(this.part.dir, "figure2.png")
ggsave(fig2.outfile, fig2, width = 10, height = 7)
