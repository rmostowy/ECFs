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
this.part.dir <- file.path(ecf.proj.dir, "part1_basic-prop")
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
  filter(prob >= MIN.PROB.FOR.DOMAIN.DETECT.IN.RHMM,
         eval <= MAX.EVAL.FOR.DOMAIN.DETECT.IN.RHMM, 
         qcov >= MIN.RHMM.COV.FOR.DOMAIN.DETECT.IN.RHMM,
         scov >= MIN.DOMAIN.COV.FOR.DOMAIN.DETECT.IN.RHMM,
         hit.length >= MIN.HIT.LENGTH.FOR.DOMAIN.DETECT.IN.RHMM
  )

reprseq.vs.ecod <- reprseq.vs.ecod %>% 
  filter(prob >= MIN.PROB.FOR.DOMAIN.DETECT.IN.RHMM,
         eval <= MAX.EVAL.FOR.DOMAIN.DETECT.IN.RHMM, 
         qcov >= MIN.RHMM.COV.FOR.DOMAIN.DETECT.IN.RHMM,
         scov >= MIN.DOMAIN.COV.FOR.DOMAIN.DETECT.IN.RHMM,
         hit.length >= MIN.HIT.LENGTH.FOR.DOMAIN.DETECT.IN.RHMM
  )

reprseq.vs.pfam <- reprseq.vs.pfam %>% 
  filter(prob >= MIN.PROB.FOR.DOMAIN.DETECT.IN.RHMM,
         eval <= MAX.EVAL.FOR.DOMAIN.DETECT.IN.RHMM, 
         qcov >= MIN.RHMM.COV.FOR.DOMAIN.DETECT.IN.RHMM,
         scov >= MIN.DOMAIN.COV.FOR.DOMAIN.DETECT.IN.RHMM,
         hit.length >= MIN.HIT.LENGTH.FOR.DOMAIN.DETECT.IN.RHMM
  )

# Other
setkey(ecf.family.table, ecf)
reprseq.vs.ecfs$EF <- ecf.family.table[.(reprseq.vs.ecfs$sname)]$family

reprseq.seqs.n <- names(reprseq.seqs)
reprseq.seqs.l <- width(reprseq.seqs)
names(reprseq.seqs.l) <- reprseq.seqs.n

reprseqs.with.domains <- unique(c(reprseq.vs.ecfs$qname, reprseq.vs.ecod$qname, reprseq.vs.pfam$qname))

#################################
### ggplot2 styling parms
#################################

ecf.fig.theme <- theme_bw(base_size = 14) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title         = element_text(face = "bold"),
    axis.text          = element_text(colour = "black"),
    legend.background  = element_blank(),
    legend.key         = element_blank()
  )


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

p.length.distr.outfile <- file.path(this.part.dir, "plot4a_len-distr.png")
ggsave(p.length.distr.outfile, p.length.distr, width = 4, height = 7)

# Plot the quantile distribution to compare the lengths at tails

q.probs <- seq(0.01, 0.99, by = 0.005)

length.quantiles <- all.lengths[
  db %in% c("ECF", "ECOD", "PFAM") & length >= 20,
  .(prob   = q.probs,
    length = quantile(length, probs = q.probs, na.rm = TRUE)),
  by = db
]

length.quantiles[, db := factor(db, levels = c("ECF", "ECOD", "PFAM"))]

p.quant.distr <- ggplot(length.quantiles, aes(x = prob, y = length, colour = db)) +
  geom_line(linewidth = 0.9) +
  scale_y_log10(
    breaks = c(10, 20, 30, 50, 100, 200, 300, 500, 1000, 2000, 5000, 3000),
    labels = comma
  ) +
  scale_x_continuous(
    breaks = c(0.01, 0.2, 0.4, 0.6, 0.8, 0.99),
    labels = c("1st", "Q1", "Q2", "Q3", "Q4", "99th")
  ) +
  # scale_colour_manual(values = db.colours, name = NULL) +
  # scale_colour_manual(values = c("ECF" = "#1f77b4", "PFAM" = "#ff7f0e", "ECOD" = "#2ca02c")) +  # Custom colors
  scale_colour_manual(values = c("ECF" = "#1f77b4", "PFAM" = "#ff7f0e", "ECOD" = "#2ca02c"),
                      name = NULL) +
  labs(x = "Percentile", y = "Length (aa)") +
  ecf.fig.theme +
  theme(legend.position = c(0.15, 0.85))

p.quant.distr.outfile <- file.path(this.part.dir, "plot4b_quant-distr.png")
ggsave(p.quant.distr.outfile, p.quant.distr, width = 5, height = 4)


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
      x = "Relative start position (qstart / rHMM length)",
      y = "Relative end position (qend / rHMM length)",
      fill = "Count",
      title = title
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

rel.pos.ecfs.outfile <- file.path(this.part.dir, "plot5_ecf-rel-pos.png")
ggsave(rel.pos.ecfs.outfile, rel.pos.ecfs, width = 6, height = 5)

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
    x = "rHMM length (quantile bins)",
    y = "rHMMs (proportion)",
    fill = "Number of EFs\nper rHMM"
  ) +
  theme_bw(base_size = 14) +
  ecf.fig.theme +
  theme(legend.position = "right")

p.ef.per.reprseq.outfile <- file.path(this.part.dir, "plot8_efs-per-rhmm.png")
ggsave(p.ef.per.reprseq.outfile, p.ef.per.reprseq, width = 7, height = 5)

#################################
### Analysis 7: 
### ECF vs. domain overlap
#################################

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
# allL <- allL[ lengths(pfam_ranges) > 0 | lengths(ecod_ranges) > 0 ]
allL <- allL[ lengths(ecf_ranges) > 0 & (lengths(pfam_ranges) > 0 | lengths(ecod_ranges) > 0) ]

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
  
  log.OR <- log(both + 0.5) + log(none + 0.5) -
    log(only_domain + 0.5) - log(only_ecf + 0.5)
  
  .(p.ecf.given.domain = p_ecf_given_domain,
    p.ecf.given.no.domain = p_ecf_given_no_domain,
    log.OR = log.OR,
    covered.ecf = ecf,
    covered.domain = domain,
    covered.both = both)
}, by = qname]

conditional.probs.dt <- copy(conditional.probs)
conditional.probs.dt[, delta := p.ecf.given.domain - p.ecf.given.no.domain]

# Show that ECF coverage is higher inside domain-annotated regions
median.log.OR <- median(conditional.probs.dt$log.OR, na.rm = TRUE)

p.domain.vs.ecf.cov <- ggplot(conditional.probs.dt, aes(x = log.OR)) +
  geom_histogram(bins = 80, fill = "#4C78A8", colour = NA, alpha = 0.85) +
  geom_vline(xintercept = 0,               linewidth = 0.5, linetype = "solid",
             colour = "black") +
  geom_vline(xintercept = median.log.OR,   linewidth = 0.5, linetype = "dashed",
             colour = "#D62728") +
  labs(
    x = "ECF enrichment within\ndomain regions (log OR)",
    y = "Number of rHMMs"
  ) +
  theme_bw(base_size = 14) +
  ecf.fig.theme

p.domain.vs.ecf.cov.outfile <- file.path(this.part.dir, "plot9_ecf-cov-dom-reg.png")
ggsave(p.domain.vs.ecf.cov.outfile, p.domain.vs.ecf.cov, width = 6, height = 5)


#################################
### Generate Figure 1
#################################

# p.length.distr # panel A
# p.ef.per.reprseq # panel B
# p.domain.vs.ecf.cov # panel C

fig1 <- 
  ( p.quant.distr | p.ef.per.reprseq | p.domain.vs.ecf.cov
  ) +
  # plot_layout(heights = c(1.2, 1)) +
  plot_annotation(tag_levels = list(c("C","D","E"))) & 
  theme(
    plot.tag = element_text(face = "bold", size = 24),
    plot.tag.position = c(0.02, 0.98)
  )

fig1.outfile <- file.path(ecf.figures.dir, "figure1.png")
ggsave(fig1.outfile, fig1, width = 14, height = 4)

#################################
### Save supplementary figures
#################################

make.dir(ecf.suppl.figures.dir)

# p.ecf.size (ECF size histogram)
# p.repr.per.ecf (No of reprseqs mapped onto each ECF)
# p.no.ecf.per.ef (No of EFs per ECF barplot)
# rel.pos.ecfs (2D interval map of ECF occurrences)
# rel.pos.pfam (2D interval map of PFAM occurrences)
# rel.pos.ecod (2D interval map of ECOD occurrences)

suppl.p.ecfs.vs.rhmms <- 
      (p.ecf.size / p.repr.per.ecf) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") & 
  theme(
    plot.tag = element_text(face = "bold", size = 14)
    # ,
    # plot.tag.position = c(0.02, 0.98)
  )
suppl.p.ecfs.vs.rhmms.outfile <- file.path(ecf.suppl.figures.dir, "suppl-1p1_ecfs-vs-rhmm.png")
ggsave(suppl.p.ecfs.vs.rhmms.outfile, suppl.p.ecfs.vs.rhmms, width = 5, height = 8)

suppl.p.no.ecf.per.ef.outfile <- file.path(ecf.suppl.figures.dir, "suppl-1p2_ecfs-vs-efs.png")
ggsave(suppl.p.no.ecf.per.ef.outfile, p.no.ecf.per.ef, width = 5, height = 4)

suppl.p.length.distr.outfile <- file.path(ecf.suppl.figures.dir, "suppl-1p3_length-distr.png")
ggsave(suppl.p.length.distr.outfile, p.length.distr, width = 4, height = 6)

 suppl.p.2d.domain.occur <- 
  (rel.pos.ecfs / rel.pos.pfam / rel.pos.ecod) +
  plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(tag_levels = "A") & 
  theme(
    plot.tag = element_text(face = "bold", size = 14)
    # ,
    # plot.tag.position = c(0.02, 0.98)
  )
suppl.p.2d.domain.occur.outfile <- file.path(ecf.suppl.figures.dir, "suppl-1p4_domain-occur.png")
ggsave(suppl.p.2d.domain.occur.outfile, suppl.p.2d.domain.occur, width = 6, height = 14)





