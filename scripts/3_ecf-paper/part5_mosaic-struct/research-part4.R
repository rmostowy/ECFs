### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(patchwork)))
suppressWarnings(suppressPackageStartupMessages(library(ggh4x)))
# suppressWarnings(suppressPackageStartupMessages(library(ComplexUpset)))
# suppressWarnings(suppressPackageStartupMessages(library(ggraph)))
# suppressWarnings(suppressPackageStartupMessages(library(broom)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

#################################
### Set up output directory
#################################
this.part.dir <- file.path(ecf.proj.dir, "part5_mosaic-struct")
make.dir(this.part.dir)

#################################
### Load all data:

#################################

# Load refseq table
reprseq.table <- fread(reprseq.table.output)
all.reprseq <- unique(reprseq.table$reprseq)
setkey(reprseq.table, reprseq)

# ECF tables
ecf.table <- fread(ecf.table.output)
ecf.family.table <- fread(ecf.family.table.output)
setkey(ecf.family.table, ecf)

# Load functional annotations
reprseq.functions.unique <- fread(reprseq.functions.unique.output)
reprseq.functions <- fread(reprseq.functions.all.output)
setkey(reprseq.functions, qname)
setkey(reprseq.functions.unique, reprseq)

# Load reprseq vs ECFs
reprseq.vs.ecfs.all <- fread(reprseq.ecfs.output)
reprseq.vs.ecfs <- reprseq.vs.ecfs.all %>%
  filter(prob >= MIN.PROB.FOR.DOMAIN.EQUIVALENCE) %>%
  # filter(prob >= 99) %>%
  filter(scov >= MIN.SCOV.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE)
setkey(reprseq.vs.ecfs, qname)
reprseq.vs.ecfs$family <- ecf.family.table[.(reprseq.vs.ecfs$sname)]$family
# reprseq.vs.ecfs[, family := ecf.family.table[.(reprseq.vs.ecfs$sname)]$family]

############################################
### 
############################################

############################################################
### Essential modularity metrics from reprseq.vs.ecfs (hits)
### Metrics:
###  1) coverage_fraction
###  2) n_families_total
###  3) mean_eff_families (Simpson effective #, per position)
###  4) segregation_index (1 - mean within-family overlap of position ranges)
###  5) n_dom_segments (dominant-family segmentation, optional purity filter)
############################################################

hits <- copy(reprseq.vs.ecfs)
setDT(hits)

stopifnot(all(c("qname","qlength","qstart","qend","family") %in% names(hits)))
hits <- hits[!is.na(qname) & !is.na(family) & !is.na(qstart) & !is.na(qend)]
hits[qstart > qend, c("qstart","qend") := .(qend, qstart)]

# One qlength per protein
QL <- hits[, .(qlength = qlength[1L]), by = qname]
setkey(QL, qname)

# ------------------------------------------------------------
# Expand to per-position table (heavy but intuitive)
# Each hit contributes 1 to each covered position (multiplicity kept).
# ------------------------------------------------------------
hit_pos <- hits[, .(pos = qstart:qend, family), by = .(qname, qstart, qend, family)]
hit_pos <- hit_pos[, .(qname, pos, family)]  # drop interval cols early

# ------------------------------------------------------------
# 1) coverage_fraction  = (# positions covered by ≥1 hit) / qlength
# ------------------------------------------------------------
covered_pos <- unique(hit_pos[, .(qname, pos)])
cov <- covered_pos[, .(covered_len = .N), by = qname]
setkey(cov, qname)

# ------------------------------------------------------------
# 2) n_families_total = distinct families in protein (over all hits)
# ------------------------------------------------------------
fam_tot <- unique(hit_pos[, .(qname, family)])[, .(n_families_total = .N), by = qname]
setkey(fam_tot, qname)

# ------------------------------------------------------------
# Position-wise family counts (multiplicity of hits)
# This is the base object reused for Simpson + dominance.
# ------------------------------------------------------------
fam_counts <- hit_pos[, .N, by = .(qname, pos, family)]
setkey(fam_counts, qname, pos, family)

# also need n_families per position (distinct, not multiplicity)
rich_pos <- unique(hit_pos)[, .(n_fam_pos = .N), by = .(qname, pos)]
setkey(rich_pos, qname, pos)

# ------------------------------------------------------------
# 3) mean_eff_families (Simpson effective number) per position:
# For each position: p_f = N_f / sum_f N_f
# Simpson effective = 1 / sum_f p_f^2
# Then mean over covered positions.
# ------------------------------------------------------------
fam_counts[, tot := sum(N), by = .(qname, pos)]
fam_counts[, p := N / tot]
simpson_pos <- fam_counts[, .(
  eff = 1 / sum(p^2)
), by = .(qname, pos)]
setkey(simpson_pos, qname, pos)

eff_mean <- simpson_pos[, .(
  mean_eff_families = mean(eff)
), by = qname]
setkey(eff_mean, qname)

# ------------------------------------------------------------
# 4) segregation_index (0..1)
# Intuition: families are "segregated" if their position ranges overlap little.
# Steps:
#   - For each protein-family, compute range [min(pos), max(pos)]
#   - For protein, compute average pairwise overlap fraction between family ranges
#   - segregation_index = 1 - mean_overlap
# Notes:
#   - Uses ranges (not full shapes) => robust + simple.
#   - Requires ≥2 families; else segregation_index = NA (or 1; NA is safer).
# ------------------------------------------------------------
fam_ranges <- unique(hit_pos)[, .(
  lo = min(pos),
  hi = max(pos)
), by = .(qname, family)]
setkey(fam_ranges, qname, family)

# pairwise overlap helper for one protein
range_overlap_mean <- function(lo, hi){
  k <- length(lo)
  if (k < 2) return(NA_real_)
  ij <- combn(seq_len(k), 2)
  a <- ij[1,]; b <- ij[2,]
  inter <- pmax(0L, pmin(hi[a], hi[b]) - pmax(lo[a], lo[b]) + 1L)
  len_a <- hi[a] - lo[a] + 1L
  len_b <- hi[b] - lo[b] + 1L
  frac <- inter / pmin(len_a, len_b)   # 0..1
  mean(frac)
}

# compute only where defined
seg <- fam_ranges[, .(
  n_fam = .N,
  mean_range_overlap = range_overlap_mean(lo, hi)
), by = qname]

seg[, segregation_index := fifelse(
  n_fam >= 2,
  1 - mean_range_overlap,
  NA_real_
)]

seg[, c("n_fam","mean_range_overlap") := NULL]
setkey(seg, qname)

# ------------------------------------------------------------
# 4b) segregation_bary (0..1)
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
fam_pos_mass <- hit_pos[, .N, by = .(qname, pos, family)]

# family barycenters per protein
fam_centers <- fam_pos_mass[, .(
  mu = sum(pos * N) / sum(N),
  total_mass = sum(N)
), by = .(qname, family)]
setkey(fam_centers, qname)

# pairwise mean distance helper (expects vectors mu and scalar L)
mu_dist_mean <- function(mu, L){
  k <- length(mu)
  if (k < 2 || is.na(L) || L <= 0) return(NA_real_)
  ij <- combn(mu, 2)
  mean(abs(ij[1, ] - ij[2, ]) / L)
}

# compute barycentric segregation per protein
seg_bary <- fam_centers[
  QL, on = "qname"
][, .(
  n_fam = .N,
  segregation_bary = mu_dist_mean(mu, qlength[1L])
), by = qname]

seg_bary[n_fam < 2, segregation_bary := NA_real_]
seg_bary[, n_fam := NULL]
setkey(seg_bary, qname)

# ------------------------------------------------------------
# 5) n_dom_segments
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
setorder(fam_counts, qname, pos, -N, family)  # tie-break by family for stability
dom <- fam_counts[, .SD[1L], by = .(qname, pos)]  # first row per pos = max N
dom[, purity := p]  # since p computed above from N/tot
dom[, dom_family := family]
dom <- dom[, .(qname, pos, dom_family, purity)]
setkey(dom, qname, pos)

PURITY_MIN <- 0.6  # tweak; set to 0 to disable
if (PURITY_MIN > 0) dom[purity < PURITY_MIN, dom_family := NA_character_]

# count segments (runs) in dom_family track, ignoring NA positions
setorder(dom, qname, pos)
dom[, prev := shift(dom_family), by = qname]
dom[, new_seg := !is.na(dom_family) & (is.na(prev) | dom_family != prev)]
segN <- dom[, .(n_dom_segments = sum(new_seg)), by = qname]
setkey(segN, qname)

# ------------------------------------------------------------
# 6) Number of ECFs per reprseq
# ------------------------------------------------------------

ecf_counts <- hits[, .(
  n_ecf_hits = uniqueN(sname)
), by = qname]

# ------------------------------------------------------------
# Assemble final table
# ------------------------------------------------------------
pm <- QL[cov, on="qname"]
pm[is.na(covered_len), covered_len := 0L]
pm[, coverage_fraction := fifelse(qlength > 0, covered_len / qlength, NA_real_)]

pm <- fam_tot[pm, on="qname"]
pm[is.na(n_families_total), n_families_total := 0L]

pm <- eff_mean[pm, on="qname"]
pm[covered_len == 0L, mean_eff_families := NA_real_]

pm <- seg[pm, on="qname"]
# if <2 families => segregation_index is NA (not defined)
pm[n_families_total < 2, segregation_index := NA_real_]

pm <- seg_bary[pm, on="qname"]
pm[n_families_total < 2, segregation_bary := NA_real_]

pm <- segN[pm, on="qname"]
pm[is.na(n_dom_segments), n_dom_segments := 0L]

pm <- ecf_counts[pm, on="qname"]
pm[is.na(n_ecf_hits), n_ecf_hits := 0L]

protein_metrics_essential <- pm[, .(
  qname,
  qlength,
  n_ecf_hits,
  coverage_fraction,
  n_families_total,
  mean_eff_families,
  segregation_index,
  segregation_bary,      # NEW
  n_dom_segments
)]

# Append functional information
pm.funct <- reprseq.functions.unique[.(protein_metrics_essential$qname)][, .(qname = reprseq, funct, class)]
protein_metrics_essential <- merge(protein_metrics_essential, pm.funct, by = "qname")
protein_metrics_essential[class == "", class := "unknown"]
protein_metrics_essential[funct == "", funct := "unknown"]

# select top N most frequent functions
N_top <- 40
top_fun <- protein_metrics_essential[, .N, by = funct][order(-N)][1:N_top]$funct
dt_top <- protein_metrics_essential[funct %chin% top_fun]
fun_order <- dt_top[, .(med = median(n_families_total)), by = funct][order(med)]$funct
dt_top[, funct := factor(funct, levels = fun_order)]


repr <- list.files('/Users/rmostowy/MGG Dropbox/Rafal Mostowy/Projects/Data/ECFs/v4/paper-ecf/part5_mosaic-struct/representative')
repr <- gsub(".png", "", repr)
protein_metrics_essential[qname %chin% repr]


# > protein_metrics_essential
# Key: <qname>
#   Index: <funct>
#   qname qlength n_ecf_hits coverage_fraction n_families_total mean_eff_families segregation_index segregation_bary n_dom_segments
# <char>   <int>      <int>             <num>            <int>             <num>             <num>            <num>          <int>
#   1: reprseq000001     536          2         0.9850746                1                 1                NA               NA              1
# 2: reprseq000005      99          4         0.7878788                1                 1                NA               NA              1
# 3: reprseq000007     109          1         0.9816514                1                 1                NA               NA              1
# 4: reprseq000008     207          6         0.3043478                1                 1                NA               NA              1
# 5: reprseq000010     140          4         0.4000000                1                 1                NA               NA              1
# ---                                                                                                                                        
#   73306: reprseq147995     170          2         0.7764706                1                 1                NA               NA              1
# 73307: reprseq147999     838          8         0.4797136                1                 1                NA               NA              1
# 73308: reprseq148000      74         13         0.9864865                1                 1                NA               NA              1
# 73309: reprseq148005      63          1         1.0000000                1                 1                NA               NA              1
# 73310: reprseq148007     250          1         0.2440000                1                 1                NA               NA              1


####### PLOTS

# H1: Most viral proteins are dominated by a single, evolutionary family

p.no.fam <- ggplot(protein_metrics_essential, aes(x = n_families_total)) +
  geom_histogram(bins=40) + 
  labs(x = "Total number of families per rHMM")

p.no.fam.high.cov <- ggplot(protein_metrics_essential[coverage_fraction >= 0.9], aes(x = n_families_total)) +
  geom_histogram(bins=40) + 
  labs(x = "Total number of families per rHMM\n(assuming min. 90% coverage")

# table(protein_metrics_essential[coverage_fraction>0.8 & n_ecf_hits > 10]$n_families_total)

p.no.fam.per.ecf.hits <- ggplot(protein_metrics_essential[coverage_fraction >= 0.], aes(x = n_ecf_hits, y = n_families_total)) +
  geom_hex(bins = 60) +
  theme_bw() +
  labs(x = "Number of ECFs per rHMM",
       y = "Number of EFs per rHMM"
       )

p.no.fam.per.length <- ggplot(protein_metrics_essential[coverage_fraction >= 0.], aes(x = qlength, y = n_families_total)) +
  geom_hex(bins = 60) +
  theme_bw() +
  labs(x = "rHMM length",
       y = "Number of EFs per rHMM"
  )

# define thresholds
thresholds <- c(1,2,3,5,10,20,30,50,75,100)
curve_dt <- rbindlist(lapply(thresholds, function(k){
  sub <- protein_metrics_essential[n_ecf_hits > k & coverage_fraction > 0.]
  sub.table <- table(sub$n_families_total)
  data.table(
    threshold = k,
    N = nrow(sub),
    mean.no.families = mean(sub$n_families_total),
    dominant.category = as.integer(names(which.max(sub.table)))
  )
}))

p.fam.div.vs.ecf.depth <- ggplot(curve_dt, aes(x = threshold)) +
  # mean number of families
  geom_line(aes(y = mean.no.families), size = 1.2, color = "black") +
  geom_point(aes(y = dominant.category), size = 5, color = "red") +
  geom_point(aes(y = mean.no.families, size = N), color = "black") +
  scale_x_continuous(trans = "log10") +
  scale_y_sqrt() +
  labs(
    x = "Minimum number of ECF hits (log scale)",
    y = "Mean number of families (black)\nDominant number of families (red)",
    title = "Family diversity increases with ECF depth"
  ) +
  theme_bw(base_size = 14)

# H2: When multiple families occur, their spatial segregation can occur 
# and depends on the underlying functional category

# 1. Explore distribution of number of families per functional category
class_order <- protein_metrics_essential[, .(med = median(n_families_total)), by = class][order(med)]$class
dt <- copy(protein_metrics_essential)
dt[, class := factor(class, levels = class_order)]

p.no.fam.per.class <- ggplot(dt, aes(class, n_families_total)) +
  geom_violin(trim = TRUE, fill = "grey80") +
  geom_boxplot(width = 0.15, outlier.size = 0.5) +
  labs(
    x = NULL,
    y = "Number of families per protein"
  ) +
  scale_y_log10() +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p.no.fam.per.funct <- ggplot(dt_top, aes(funct, n_families_total)) +
  geom_violin(trim = TRUE, fill = "grey80") +
  geom_boxplot(width = 0.15, outlier.size = 0.5) +
  labs(
    x = NULL,
    y = "Number of families per protein"
  ) +
  theme_bw(base_size = 14) +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# 2. Look at the distribution of spatial segregation in multi-family proteins per function
p.segeration.index.per.class <- ggplot(dt[n_families_total >= 2], aes(x = class, y = segregation_index)) +
  geom_violin(trim = TRUE, fill = "grey80") +
  geom_boxplot(width = 0.15, outlier.size = 0.5) +
  labs(
    x = NULL,
    y = "Segregation index"
  ) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p.segeration.index.per.funct <- ggplot(dt_top[n_families_total >= 2], aes(x = funct, y = segregation_index)) +
  geom_violin(trim = TRUE, fill = "grey80") +
  geom_boxplot(width = 0.15, outlier.size = 0.5) +
  labs(
    x = NULL,
    y = "Segregation index"
  ) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p.segerationb.index.per.class <- ggplot(dt[n_families_total >= 2], aes(x = class, y = segregation_bary)) +
  geom_violin(trim = TRUE, fill = "grey80") +
  geom_boxplot(width = 0.15, outlier.size = 0.5) +
  labs(
    x = NULL,
    y = "Segregation index"
  ) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p.segerationb.index.per.funct <- ggplot(dt_top[n_families_total >= 2], aes(x = funct, y = segregation_bary)) +
  geom_violin(trim = TRUE, fill = "grey80") +
  geom_boxplot(width = 0.15, outlier.size = 0.5) +
  labs(
    x = NULL,
    y = "Segregation index"
  ) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# H3: Spatial segregation and family diversity are partly decoupled

p.segreg.meaneff <- ggplot(
  protein_metrics_essential[n_families_total >= 2 & !is.na(segregation_index) & !is.na(mean_eff_families)],
  aes(x = mean_eff_families, y = segregation_index)
) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(trans = "sqrt", name = "Count") +
  labs(
    x = "Mean effective families (Simpson, position-wise)",
    y = "Spatial segregation index"
  ) +
  theme_bw(base_size = 14)

ggsave(file.path(this.part.dir, "h3_1_effFam_vs_segregation.png"), p.segreg.meaneff, width = 6, height = 5)

p.segreg.nofam <- ggplot(
  protein_metrics_essential[n_families_total >= 2 & !is.na(segregation_index)],
  aes(x = factor(n_families_total), y = segregation_index)
) +
  geom_violin(trim = TRUE, fill = "grey80") +
  geom_boxplot(width = 0.15, outlier.size = 0.4) +
  labs(
    x = "Total number of families in protein",
    y = "Spatial segregation index"
  ) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(this.part.dir, "h3_2_segregation_by_nFamilies.png"), p.segreg.nofam, width = 7, height = 5)


segr.by.meaneff <- ggplot(
  protein_metrics_essential[n_families_total >= 2 & !is.na(segregation_index) & !is.na(mean_eff_families)],
  aes(
    x = cut(mean_eff_families, breaks = c(1, 1.2, 1.5, 2, 3, Inf), right = FALSE),
    y = segregation_index
  )
) +
  geom_violin(trim = TRUE, fill = "grey80") +
  geom_boxplot(width = 0.15, outlier.size = 0.4) +
  labs(
    x = "Mean effective families (binned)",
    y = "Spatial segregation index"
  ) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(this.part.dir, "h3_3_segregation_by_effFamBins.png"), segr.by.meaneff, width = 7, height = 5)

# H4: Viral proteins span a continuum of architectual complexity
p1 <- ggplot(
  protein_metrics_essential[coverage_fraction > 0 & !is.na(n_dom_segments)],
  aes(x = n_dom_segments)
) +
  geom_histogram(binwidth = 1, boundary = 0, closed = "left", color = "white") +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  labs(
    x = "Number of dominant-family segments (n_dom_segments)",
    y = "Number of proteins"
  ) +
  theme_bw(base_size = 14)

ggsave(file.path(this.part.dir, "h4_1_dist_n_dom_segments.png"), p1, width = 6, height = 5)
p1

p2 <- ggplot(
  protein_metrics_essential[coverage_fraction > 0 & !is.na(n_dom_segments) & !is.na(mean_eff_families)],
  aes(x = factor(n_dom_segments), y = mean_eff_families)
) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.6) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  geom_jitter(alpha = 0.08, width = 0.15, height = 0) +
  labs(
    x = "Number of dominant-family segments (n_dom_segments)",
    y = "Mean effective number of families (Simpson, per covered position)"
  ) +
  theme_bw(base_size = 14)

ggsave(file.path(this.part.dir, "h4_2_segments_vs_mean_eff_families.png"), p2, width = 6, height = 5)
p2


p3 <- ggplot(
  protein_metrics_essential[n_families_total >= 2 & !is.na(segregation_index) & !is.na(n_dom_segments)],
  aes(x = segregation_index, y = n_dom_segments)
) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(trans = "sqrt", name = "Count") +
  labs(
    x = "Spatial segregation index",
    y = "Number of dominant-family segments (n_dom_segments)"
  ) +
  theme_bw(base_size = 14)

ggsave(file.path(this.part.dir, "h4_3_segregation_vs_segments_hex.png"), p3, width = 6, height = 5)
p3

p4 <- ggplot(
  protein_metrics_essential[coverage_fraction > 0 & !is.na(n_dom_segments) & n_ecf_hits > 0],
  aes(x = n_ecf_hits, y = n_dom_segments)
) +
  geom_hex(bins = 70) +
  scale_fill_viridis_c(trans = "sqrt", name = "Count") +
  scale_x_continuous(trans = "log10") +
  labs(
    x = "Number of ECF hits (log10 scale)",
    y = "Number of dominant-family segments (n_dom_segments)"
  ) +
  theme_bw(base_size = 14)

ggsave(file.path(this.part.dir, "h4_4_hits_vs_segments_hex.png"), p4, width = 6, height = 5)
p4

######################
### Figure ideas
######################

# Punchline 1 — Multi-family architecture is common once you have enough evidence
p.fam.div.vs.ecf.depth.outfile <- file.path(this.part.dir, "ideas-panelA.png")
ggsave(p.fam.div.vs.ecf.depth.outfile, p.fam.div.vs.ecf.depth, width = 8, height = 6.5)

# Punchline 2 — Family diversity and spatial organization are distinct axes
p.fam.div.spatial.separation <- ggplot(protein_metrics_essential, aes(x = n_families_total, y = segregation_bary)) +
  geom_bin2d(bins = 60) +
  theme_bw(base_size = 14) +
  scale_fill_viridis_c(trans = "sqrt", name = "Count") +
  labs(
    x = "Number of ECF Families in rHMM",
    y = "ECF family segregation index along protein sequence"
  )
p.fam.div.spatial.separation.outfile <- file.path(this.part.dir, "ideas-panelB.png")
ggsave(p.fam.div.spatial.separation.outfile, p.fam.div.spatial.separation, width = 8, height = 6.5)


# Punchline 3 - There are “architectural regimes” across functional classes
# dt <- copy(protein_metrics_essential)
no.functs.per.class <- 5
dt <- copy(protein_metrics_essential)
plot_dt <- rbindlist(list(
  dt[funct != "unknown" & class != "unknown",
     .N, by = .(class, funct)][
       order(class, -N)
     ][, head(.SD, no.functs.per.class), by = class][
       dt, on = .(class, funct), nomatch = 0L
     ][, panel := class],
  
  dt[funct == "unknown" | class == "unknown"][, `:=`(panel = "unknown", funct = "unknown")]
  # ,
  
  # dt[, `:=`(panel = "ALL", funct = "ALL")]
), use.names = TRUE, fill = TRUE)

plot_dt <- plot_dt[!is.na(segregation_bary)]
plot_dt[, xlab := sprintf("%s (N=%d)", funct, .N), by = .(panel, funct)]
plot_dt[, xlab := factor(xlab, levels = unique(xlab)), by = panel]

phrogs.class2 <- c(phrogs.class, "unknown")
phrog.class.labels2 <- c(phrog.class.labels, "unknown")
names(phrog.class.labels2)[phrog.class.labels2 == "unknown"] <- "unknown"

plot_dt[, panel := factor(class, levels = phrogs.class2)]

strip <- ggh4x::strip_themed(
  background_y = ggh4x::elem_list_rect(
    fill = grDevices::adjustcolor(
      phrogs.col[match(levels(plot_dt$panel), phrogs.class2)],
      alpha.f = 0.25
    ),
    colour = "black",
    linewidth = 0.4
  )
)

p.segreg.across.fam <- ggplot(
  plot_dt[!is.na(segregation_bary) & panel != "ALL"],
  aes(x = segregation_bary, y = reorder(funct, segregation_bary, FUN = median, na.rm = TRUE))
) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.18, outlier.shape = NA) +
  facet_grid2(
    rows = vars(panel),
    scales = "free_y",
    space  = "free_y",
    labeller = labeller(panel = phrog.class.labels),
    strip = strip
  ) +
  geom_vline(xintercept = mean(protein_metrics_essential$segregation_bary, na.rm = T), linetype = "dashed", color = "red", linewidth = 0.6) +
  labs(
    x = "Segregation index",
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.placement = "outside",
    strip.text.y = element_text(angle = 0, hjust = 0, face = "bold", size = 9),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 9),
    panel.spacing = unit(0.35, "lines")
  )
p.segreg.across.fam.outfile <- file.path(this.part.dir, "ideas-panelC.png")

ggsave(p.segreg.across.fam.outfile, p.segreg.across.fam, width = 6, height = 10)


# Punchline 4 (optional but powerful) — Extreme architectures exist and are interpretable exemplars

plot.ecf.distr.over.query <- function(protein_name, protein_ecf_hits_all, title = "ECF Domains and Query Alignment") {
  
  # Check if domain_data is NULL or empty
  if (is.null(protein_ecf_hits_all) || all(is.na(protein_ecf_hits_all$qstart))) {
    message("No domains detected for this protein.")
    return(NULL)
  }
  
  protein_length <- protein_ecf_hits_all$qlength[1]
  
  # sort by ECF family
  setorder(protein_ecf_hits_all, family, qstart)
  
  # Add a y-axis level for each ECF domain, reducing vertical spread for clarity
  protein_ecf_hits_all$y_level <- seq_len(nrow(protein_ecf_hits_all)) * 0.5
  
  # Define colors for families
  family_colors <- unique(protein_ecf_hits_all$family)
  no.fam <- length(family_colors)
  if(no.fam < 3) no.fam <- 3
  color_palette <- setNames(RColorBrewer::brewer.pal(no.fam, "Set3"), family_colors)
  
  # Create the plot
  p <- ggplot() +
    # Plot the query protein as a black line with the aligned segment in red
    geom_segment(aes(x = 0, xend = protein_length, y = 0, yend = 0), linewidth = 2, color = "black") +
    # geom_segment(aes(x = aln1_start, xend = aln1_end, y = 0, yend = 0), linewidth = 2, color = "red") +
    # Plot each ECF domain as a coloured segment, positioned on its y_level
    geom_segment(data = protein_ecf_hits_all, aes(x = qstart, xend = qend, y = y_level, yend = y_level, color = family),
                 linewidth = 4) +
    # Add labels for each ECF domain
    geom_text(data = protein_ecf_hits_all, aes(x = (qstart + qend) / 2, y = y_level + 0.1, label = sname),
              size = 3, hjust = 0.5, color = "black") +
    # Define custom colours for families
    scale_color_manual(values = color_palette) +
    # Customise the plot
    labs(
      x = sprintf("Position along protein: %s", protein_name),
      y = "",
      title = title,
      color = "Family"
    ) +
    # theme_classic(base_size = 12) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.margin = margin(20, 30, 20, 30)
    ) +
    coord_cartesian(xlim = c(0, protein_length))
  
  
  return(p)
}

quantile.probs <- c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
seg.quant <- quantile(protein_metrics_essential$segregation_bary, probs = quantile.probs, na.rm = T)
meaneff.quant <- quantile(protein_metrics_essential$mean_eff_families, probs = quantile.probs, na.rm = T)

# Example 1 (clean 2-block):
#   •	n_families_total == 2
# •	mean_eff_families < 1.5
# •	n_dom_segments == 2
# •	segregation_bary in top 5–10% of that subset

this.ex.dir <- file.path(this.part.dir, "vis-ecf-distr", "example1_clean2block")
make.dir(this.ex.dir)
ex1.shortlist <- protein_metrics_essential[n_families_total == 2 &
                                             mean_eff_families < 1.5 &
                                             n_dom_segments == 2 &
                                             segregation_bary >= seg.quant["90%"] &
                                             n_ecf_hits > 20]

for(index in 1:nrow(ex1.shortlist)){
  protein_name <- ex1.shortlist$qname[index]
  protein_ecf_hits_all <- reprseq.vs.ecfs[.(protein_name)]
  title <- sprintf("%s (%s)", protein_name, ex1.shortlist$funct[index])
  this.plot <- plot.ecf.distr.over.query(protein_name, protein_ecf_hits_all, title = title)
  this.plot.outfile <- file.path(this.ex.dir, sprintf("%s.png", protein_name))
  ggsave(this.plot.outfile, this.plot, width = 6, height = 5)
}

# Handpicked:
# reprseq080790 (endolysin)
# reprseq005066 (capsid assembly)
# reprseq035975 (unknown)
# reprseq040331 (type II RM system) <--

#########################################

# •	Example 2 (overlap):
#   •	n_families_total >= 2
# •	mean_eff_families high (top 10–20%)
# •	segregation_bary low (bottom 10–20%)
# •	(optional) n_dom_segments low-ish

this.ex.dir <- file.path(this.part.dir, "vis-ecf-distr", "example2_overlap")
make.dir(this.ex.dir)
ex2.shortlist <- protein_metrics_essential[n_families_total >= 2 &
                                             mean_eff_families >= meaneff.quant["90%"] &
                                             # n_dom_segments == 2 &
                                             segregation_bary <= seg.quant["10%"] &
                                             n_ecf_hits > 20]

for(index in 1:nrow(ex2.shortlist)){
  protein_name <- ex2.shortlist$qname[index]
  protein_ecf_hits_all <- reprseq.vs.ecfs[.(protein_name)]
  title <- sprintf("%s (%s)", protein_name, ex2.shortlist$funct[index])
  this.plot <- plot.ecf.distr.over.query(protein_name, protein_ecf_hits_all, title = title)
  this.plot.outfile <- file.path(this.ex.dir, sprintf("%s.png", protein_name))
  ggsave(this.plot.outfile, this.plot, width = 6, height = 5)
}

# Handpicked:
# reprseq046469 (clamp loader of DNAp)
# reprseq046736 (ATPase) <--
# reprseq055118 (tail fibre)
# reprseq059226 (DNA helicase)


###########################################

# •	Example 3 (fragmented):
#   •	n_dom_segments top 1–5%
# •	coverage decent (to avoid edge artifacts)

this.ex.dir <- file.path(this.part.dir, "vis-ecf-distr", "example3_mosaic")
make.dir(this.ex.dir)
ex3.shortlist <- protein_metrics_essential[n_families_total >= 2 &
                                             # mean_eff_families == meaneff.quant["10%"] &
                                             n_dom_segments >= 5 &
                                             segregation_bary >= seg.quant["90%"] &
                                             n_ecf_hits > 20]

for(index in 1:nrow(ex3.shortlist)){
  protein_name <- ex3.shortlist$qname[index]
  protein_ecf_hits_all <- reprseq.vs.ecfs[.(protein_name)]
  title <- sprintf("%s (%s)", protein_name, ex3.shortlist$funct[index])
  this.plot <- plot.ecf.distr.over.query(protein_name, protein_ecf_hits_all, title = title)
  this.plot.outfile <- file.path(this.ex.dir, sprintf("%s.png", protein_name))
  ggsave(this.plot.outfile, this.plot, width = 6, height = 5)
}

# Handpicked:
# reprseq004062 (tail spike)
# reprseq029825 (TLTP) <-- 
# reprseq054784 (unknown)


# •	Example 4 (asymmetric rare family):
#   •	n_families_total == 2
# •	mean_eff_families ~ 1.05–1.2
# •	but segregation_bary non-trivial
# •	(and/or mass ratio extreme if you compute it)


this.ex.dir <- file.path(this.part.dir, "vis-ecf-distr", "example4_asymm")
make.dir(this.ex.dir)
ex4.shortlist <- protein_metrics_essential[n_families_total == 2 &
                                             mean_eff_families <= 1.05 &
                                             # n_dom_segments >= 5 &
                                             segregation_bary >= seg.quant["75%"] &
                                             n_ecf_hits > 30]

for(index in 1:nrow(ex4.shortlist)){
  protein_name <- ex4.shortlist$qname[index]
  protein_ecf_hits_all <- reprseq.vs.ecfs[.(protein_name)]
  title <- sprintf("%s (%s)", protein_name, ex4.shortlist$funct[index])
  this.plot <- plot.ecf.distr.over.query(protein_name, protein_ecf_hits_all, title = title)
  this.plot.outfile <- file.path(this.ex.dir, sprintf("%s.png", protein_name))
  ggsave(this.plot.outfile, this.plot, width = 6, height = 5)
}

# Handpicked
# reprseq077019 (portal)


# Final set:
# 1. Clean 2-block: reprseq040331 (type II RM system)
# 2. Overlap: reprseq046736 (ATPase)
# 3. Fragmented: reprseq029825 (TLTP)
# 4. Asymmetric: reprseq077019 (portal)

name <- "reprseq040331"
protein_ecf_hits_all <- reprseq.vs.ecfs[.(name)]
title <- sprintf("Clean 2-block (%s)", protein_metrics_essential[qname == name]$funct)
p1 <- plot.ecf.distr.over.query(name, protein_ecf_hits_all, title = title)

name <- "reprseq046736"
protein_ecf_hits_all <- reprseq.vs.ecfs[.(name)]
title <- sprintf("Overlap (%s)", protein_metrics_essential[qname == name]$funct)
p2 <- plot.ecf.distr.over.query(name, protein_ecf_hits_all, title = title)

name <- "reprseq029825"
protein_ecf_hits_all <- reprseq.vs.ecfs[.(name)]
title <- sprintf("Fragmented (%s)", protein_metrics_essential[qname == name]$funct)
p3 <- plot.ecf.distr.over.query(name, protein_ecf_hits_all, title = title)

name <- "reprseq077019"
protein_ecf_hits_all <- reprseq.vs.ecfs[.(name)]
title <- sprintf("Asymmetric (%s)", protein_metrics_essential[qname == name]$funct)
p4 <- plot.ecf.distr.over.query(name, protein_ecf_hits_all, title = title)

p.examples <- (p1 / p2 / p3 / p4) &
  theme(
    plot.tag = element_text(face = "bold", size = 24),
    plot.tag.position = c(0.02, 0.98)
  )

p.examples.outfile <- file.path(this.part.dir, "ideas-panelD.png")
ggsave(p.examples.outfile, p.examples, width = 4, height = 11)
