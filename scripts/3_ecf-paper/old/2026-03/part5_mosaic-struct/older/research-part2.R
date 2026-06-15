### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
# suppressWarnings(suppressPackageStartupMessages(library(patchwork)))
# suppressWarnings(suppressPackageStartupMessages(library(ggh4x)))
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
# reprseq.functions <- fread(reprseq.functions.unique.output)
reprseq.functions <- fread(reprseq.functions.all.output)
setkey(reprseq.functions, qname)

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

this.reprseq <- "reprseq134177"
# # this.reprseq <- "reprseq057676"
this.reprseq.ecf.hits <- reprseq.vs.ecfs[.(this.reprseq)]
setorder(this.reprseq.ecf.hits, qname, family)

# 1. compute number of EFs per position
# 2. compute the most common EF per position
# 3. compute the 

hits <- copy(this.reprseq.ecf.hits)
L <- unique(hits$qlength)
stopifnot(length(L) == 1)

# position-level table
pos_dt <- data.table(pos = 1:L)

# expand hits into per-position table
hit_pos <- hits[, .(pos = qstart:qend, family), by = .(qname, family, qstart, qend)]

# Distinct families covering each position:
family_pos <- unique(hit_pos[, .(pos, family)])

rich_dt <- family_pos[, .(
  n_families = .N
), by = pos]

# merge with all positions
pos_dt <- merge(pos_dt, rich_dt, by = "pos", all.x = TRUE)
pos_dt[is.na(n_families), n_families := 0]

# Mean family richness
mean_family_richness <- mean(pos_dt[n_families>0]$n_families)

# coverage fraction
coverage_fraction <- mean(pos_dt$n_families > 0)

# overlap burden
union_len <- sum(pos_dt$n_families > 0)
sum_hit_len <- hits[, sum(hit.length)]
overlap_burden <- 1 - (union_len / sum_hit_len)

# positional entropy
family_counts <- hit_pos[, .N, by = .(pos, family)]
family_counts[, total := sum(N), by = pos]
family_counts[, p := N / total]
family_counts[, entropy := -p * log(p)]
entropy_dt <- family_counts[, .(
  H = sum(entropy)
), by = pos]
pos_dt <- merge(pos_dt, entropy_dt, by = "pos", all.x = TRUE)
pos_dt[is.na(H), H := 0]
mean_entropy <- mean(pos_dt$H)

# dominant family per position
dom_family <- family_counts[
  , .SD[which.max(N)], 
  by = pos
][, .(pos, family)]

pos_dt <- merge(pos_dt, dom_family, by = "pos", all.x = TRUE)

# positions with no family → NA
track <- pos_dt$family
valid_idx <- which(!is.na(track))
transitions <- sum(track[valid_idx][-1] != track[valid_idx][-length(valid_idx)])
transition_score <- transitions / length(valid_idx)

protein_metrics <- data.table(
  qname = unique(hits$qname),
  mean_family_richness = mean_family_richness,
  mean_entropy = mean_entropy,
  coverage_fraction = coverage_fraction,
  overlap_burden = overlap_burden,
  transition_score = transition_score
)
protein_metrics

#######################


hits <- copy(reprseq.vs.ecfs)
setDT(hits)

stopifnot(all(c("qname","qlength","qstart","qend","family","hit.length") %in% names(hits)))
hits <- hits[!is.na(qname) & !is.na(family) & !is.na(qstart) & !is.na(qend)]
hits[qstart > qend, c("qstart","qend") := .(qend, qstart)]

# One qlength per qname
QL <- hits[, .(qlength = qlength[1L]), by = qname]
setkey(QL, qname)

# Sum of hit.length per protein (for overlap_burden)
sum_hit_len <- hits[, .(sum_hit_len = sum(hit.length, na.rm=TRUE)), by = qname]
setkey(sum_hit_len, qname)

# hits2 <- copy(hits[1:100000])

hit_pos <- hits[, .(pos = qstart:qend, family), by = .(qname, family, qstart, qend)]
hit_pos <- hit_pos[, .(qname, pos, family)]  # drop interval cols immediately

# ------------------------------------------------------------
# Distinct families per position (per protein) => richness track
# ------------------------------------------------------------

family_pos <- unique(hit_pos)  # unique (qname,pos,family)
rich_pos <- family_pos[, .(n_families = .N), by = .(qname, pos)]
# covered length per protein and mean richness over covered positions
rich_sum <- rich_pos[, .(
  union_len = .N,
  mean_family_richness = mean(n_families)
), by = qname]
setkey(rich_sum, qname)


# coverage_fraction and overlap_burden
protein_metrics <- QL[rich_sum, on="qname"]
protein_metrics[is.na(union_len), union_len := 0L]
protein_metrics[, coverage_fraction := fifelse(qlength > 0, union_len / qlength, NA_real_)]

protein_metrics <- sum_hit_len[protein_metrics, on="qname"]
protein_metrics[is.na(sum_hit_len), sum_hit_len := 0]
protein_metrics[, overlap_burden := fifelse(sum_hit_len > 0, 1 - (union_len / sum_hit_len), NA_real_)]


# ------------------------------------------------------------
# Positional entropy using hit multiplicity per (qname,pos,family)
# (matches your original family_counts built from hit_pos)
# ------------------------------------------------------------
family_counts <- hit_pos[, .N, by = .(qname, pos, family)]
family_counts[, total := sum(N), by = .(qname, pos)]
family_counts[, p := N / total]
family_counts[, Hcomp := -p * log(p)]
entropy_pos <- family_counts[, .(H = sum(Hcomp)), by = .(qname, pos)]

# mean entropy across ALL positions; uncovered positions contribute 0
sumH <- entropy_pos[, .(sumH = sum(H)), by = qname]
setkey(sumH, qname)
protein_metrics <- sumH[protein_metrics, on="qname"]
protein_metrics[is.na(sumH), sumH := 0]
protein_metrics[, mean_entropy := fifelse(qlength > 0, sumH / qlength, NA_real_)]
protein_metrics[, sumH := NULL]


# ------------------------------------------------------------
# Dominant family per position + transition_score
# (dominant = argmax N within position; ties break by first encountered)
# ------------------------------------------------------------
# For stable tie-breaking, order before .SD[which.max()]
setorder(family_counts, qname, pos, -N, family)

dom <- family_counts[, .SD[1L], by = .(qname, pos)]
setorder(dom, qname, pos)

dom[, prev := shift(family), by = qname]
dom[, trans := !is.na(prev) & (family != prev)]
trans_sum <- dom[, .(transitions = sum(trans)), by = qname]
setkey(trans_sum, qname)

protein_metrics <- trans_sum[protein_metrics, on="qname"]
protein_metrics[is.na(transitions), transitions := 0L]
protein_metrics[, transition_score := fifelse(union_len > 0, transitions / union_len, NA_real_)]
protein_metrics[, transitions := NULL]

# Final shape
protein_metrics <- protein_metrics[, .(
  qname,
  mean_family_richness,
  mean_entropy,
  coverage_fraction,
  overlap_burden,
  transition_score
)]

protein_metrics


########################
## PLOT STATS


# pm <- as.data.table(protein_metrics)
# 
# # basic summaries
# pm[, lapply(.SD, function(x) c(
#   min=min(x), q25=quantile(x,0.25), med=median(x), q75=quantile(x,0.75), max=max(x)
# )), .SDcols = c("mean_family_richness","mean_entropy","coverage_fraction","overlap_burden","transition_score")]
# 
# # correlations (Spearman is safer for skewed stuff)
# num_cols <- c("mean_family_richness","mean_entropy","coverage_fraction","overlap_burden","transition_score")
# cor(pm[, ..num_cols], method="spearman", use="pairwise.complete.obs")
# 
# p1 <- ggplot(pm, aes(mean_family_richness)) + geom_histogram(bins=80) + theme_bw(14)
# p2 <- ggplot(pm, aes(mean_entropy)) + geom_histogram(bins=80) + theme_bw(14)
# p3 <- ggplot(pm, aes(transition_score)) + geom_histogram(bins=80) + theme_bw(14)
# p4 <- ggplot(pm, aes(overlap_burden)) + geom_histogram(bins=80) + theme_bw(14)
# p1; p2; p3; p4
# 
# 
# # define thresholds
# thr_cov_hi <- pm[, quantile(coverage_fraction, 0.75)]
# thr_olap_hi <- pm[, quantile(overlap_burden, 0.75)]
# thr_rich_hi <- pm[, quantile(mean_family_richness, 0.75)]
# thr_tr_hi   <- pm[, quantile(transition_score, 0.75)]
# 
# pm[, arch_type := fifelse(coverage_fraction < 0.5, "low_coverage",
#                           fifelse(mean_family_richness <= 1.05 & overlap_burden <= thr_olap_hi, "single_family_tiling",
#                                   fifelse(mean_family_richness >= thr_rich_hi & overlap_burden >= thr_olap_hi & transition_score < thr_tr_hi, "mixed_but_stable",
#                                           fifelse(transition_score >= thr_tr_hi & mean_family_richness >= 1.1, "segmented_or_interleaved",
#                                                   "intermediate"))))]
# 
# pm[, .N, by=arch_type][order(-N)]
# 


p1 <- ggplot(pm, aes(mean_family_richness, transition_score)) +
  geom_bin2d(bins = 100) +
  theme_bw(14)

p2 <- ggplot(pm, aes(mean_family_richness)) +
  geom_histogram(bins = 100) +
  theme_bw(14)


p3 <- ggplot(pm, aes(transition_score)) +
  geom_histogram(bins = 100) +
  theme_bw(14)

p4 <- ggplot(pm, aes(mean_family_richness, overlap_burden)) +
  geom_bin2d(bins = 100) +
  theme_bw(14)

p5 <- ggplot(pm, aes(coverage_fraction, mean_family_richness)) +
  geom_bin2d(bins = 100) +
  theme_bw(14)

p6 <- ggplot(pm, aes(mean_family_richness, transition_score)) +
  geom_density_2d(color="black") +
  theme_bw(14)


ggsave(file.path(this.part.dir, "plot1.png"), p1, width = 6, height = 5)
ggsave(file.path(this.part.dir, "plot2.png"), p2, width = 6, height = 5)
ggsave(file.path(this.part.dir, "plot3.png"), p3, width = 6, height = 5)
ggsave(file.path(this.part.dir, "plot4.png"), p4, width = 6, height = 5)
ggsave(file.path(this.part.dir, "plot5.png"), p5, width = 6, height = 5)
ggsave(file.path(this.part.dir, "plot6.png"), p6, width = 6, height = 5)


reprseq2fun <- unique(reprseq.functions[top.annotation == TRUE, .(qname, funct, class)])


pm_with_function <- merge(pm, reprseq2fun, by = 'qname', all.x = T)

ggplot(pm_with_function,
       aes(mean_family_richness, transition_score)) +
  geom_bin2d(bins = 80) +
  facet_wrap(~ class, scales = "free_y") +
  theme_bw(14) +
  labs(x = "Mean family richness",
       y = "Transition score")
