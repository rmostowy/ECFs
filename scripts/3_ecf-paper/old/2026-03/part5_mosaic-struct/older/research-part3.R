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

############################################################
### Essential modularity metrics from reprseq.vs.ecfs (hits)
### Metrics: coverage_fraction, multi_fraction,
###          mean_eff_families, fragmentation_index
############################################################

hits <- copy(reprseq.vs.ecfs)
setDT(hits)

stopifnot(all(c("qname","qlength","qstart","qend","family") %in% names(hits)))
hits <- hits[!is.na(qname) & !is.na(family) & !is.na(qstart) & !is.na(qend)]
hits[qstart > qend, c("qstart","qend") := .(qend, qstart)]

# One qlength per protein
QL <- hits[, .(qlength = qlength[1L]), by = qname]
setkey(QL, qname)

# ------------------------------------------------------------------
# Expand to per-position + keep multiplicity (each hit contributes 1)
# ------------------------------------------------------------------
# NOTE: this can be heavy. If memory becomes an issue, we can switch
# to an interval-based approach, but keep this for now as requested.
hit_pos <- hits[, .(pos = qstart:qend, family), by = .(qname, qstart, qend, family)]
hit_pos <- hit_pos[, .(qname, pos, family)]  # drop interval cols

# -------------------------
# 1) coverage_fraction
# -------------------------
covered_pos <- unique(hit_pos[, .(qname, pos)])
cov <- covered_pos[, .(union_len = .N), by = qname]
setkey(cov, qname)

# -------------------------
# 2) multi_fraction
#    positions with >=2 distinct families
# -------------------------
fam_pos <- unique(hit_pos)  # (qname,pos,family) unique
rich_pos <- fam_pos[, .(n_families = .N), by = .(qname, pos)]

multi <- rich_pos[, .(
  multi_len = sum(n_families >= 2L),
  union_len_check = .N
), by = qname]
setkey(multi, qname)

# -------------------------
# 3) mean_eff_families
#    Use family-size correction within each protein:
#    each hit from family f gets weight = 1 / nHits(protein,family)
#    Then per position compute p_f(pos) from weighted sums.
# -------------------------
# counts of hits per protein-family (before unique)
fam_hitN <- hit_pos[, .N, by = .(qname, family)]
setkey(fam_hitN, qname, family)

hit_pos <- fam_hitN[hit_pos, on = .(qname, family)]
hit_pos[, w := 1 / N]  # N = nHits for that family in that protein
hit_pos[, N := NULL]

# position-family weighted sums
w_pf <- hit_pos[, .(w_sum = sum(w)), by = .(qname, pos, family)]
w_pf[, w_tot := sum(w_sum), by = .(qname, pos)]
w_pf[, p := w_sum / w_tot]

# entropy per position => effective families = exp(H)
w_pf[, Hcomp := -p * log(p)]
Hpos <- w_pf[, .(H = sum(Hcomp)), by = .(qname, pos)]
Hpos[, eff := exp(H)]

eff_sum <- Hpos[, .(
  sum_eff = sum(eff),
  union_len = .N
), by = qname]
setkey(eff_sum, qname)

# mean effective families over covered positions
eff_sum[, mean_eff_families := sum_eff / union_len]
eff_sum[, sum_eff := NULL]

# -------------------------
# 4) fragmentation_index
#    Define dominant family per position using weighted p (above),
#    then count runs per family and summarize.
#    To reduce noise-driven switching, optionally ignore low-purity positions.
# -------------------------
# dominant family + purity
setorder(w_pf, qname, pos, -p, family)  # deterministic tie-break by family
dom <- w_pf[, .SD[1L], by = .(qname, pos)]  # first row per pos = max p
dom <- dom[, .(qname, pos, dom_family = family, purity = p)]

# OPTIONAL: drop ambiguous positions (uncomment to apply)
PURITY_MIN <- 0.6
dom[purity < PURITY_MIN, dom_family := NA_character_]

# count runs of dom_family within each protein
setorder(dom, qname, pos)
dom[, prev := shift(dom_family), by = qname]
dom[, new_run := !is.na(dom_family) & (is.na(prev) | dom_family != prev)]

runs_by_family <- dom[new_run == TRUE, .N, by = .(qname, dom_family)]
# N here = number of runs for that family in that protein
frag <- runs_by_family[, .(
  fragmentation_index = sum(N),         # total runs across families
  mean_runs_per_family = mean(N),       # optional diagnostic
  max_runs_in_family = max(N)           # optional diagnostic
), by = qname]
setkey(frag, qname)

# -------------------------
# Assemble final table
# -------------------------
pm <- QL[cov, on="qname"]
pm[is.na(union_len), union_len := 0L]
pm[, coverage_fraction := fifelse(qlength > 0, union_len / qlength, NA_real_)]

pm <- multi[pm, on="qname"]
pm[is.na(multi_len), multi_len := 0L]
pm[, multi_fraction := fifelse(union_len > 0, multi_len / union_len, NA_real_)]

pm <- eff_sum[pm, on="qname"]
# if no coverage, mean_eff_families should be NA (or 0; NA is safer)
pm[union_len == 0L, mean_eff_families := NA_real_]

pm <- frag[pm, on="qname"]
pm[is.na(fragmentation_index), fragmentation_index := 0L]

# Keep only essential outputs (plus qlength if you like)
protein_metrics_essential <- pm[, .(
  qname,
  qlength,
  coverage_fraction,
  multi_fraction,
  mean_eff_families,
  fragmentation_index
)]

protein_metrics_essential


####### PLOTS

pm <- copy(protein_metrics_essential)
setDT(pm)
setkey(pm, qname)

# merge class/function info if available
# expected columns: qname, class, funct (optional)
pmf <- merge(pm, reprseq.functions[, .(qname, class, funct)], by="qname", all.x=TRUE)

# useful transforms
pmf[, log_qlength := log10(qlength + 1)]
pmf[, frag := as.numeric(fragmentation_index)]  # for plotting

# covered positions only (optional, but helpful for some plots)
pmf_cov <- pmf[coverage_fraction > 0]

p1 <- ggplot(pmf_cov, aes(mean_eff_families, multi_fraction)) +
  geom_hex(bins = 60) +
  theme_bw() +
  labs(x = "Mean effective families per covered position",
       y = "Multi-hit fraction (stacking / overlap)",
       title = "Protein architecture space: mixing vs overlap")

ggsave(file.path(this.part.dir, "plot1_arch_space_hex.png"), p1, width = 6, height = 5)

p2 <- ggplot(pmf_cov, aes(mean_eff_families, multi_fraction)) +
  geom_hex(bins = 35) +
  facet_wrap(~ class, scales = "free_y") +
  theme_bw() +
  labs(x = "Mean effective families per covered position",
       y = "Multi-hit fraction",
       title = "Architecture space by functional class")

ggsave(file.path(this.part.dir, "plot2_arch_space_by_class.png"), p2, width = 10, height = 7)

p3 <- ggplot(pmf, aes(coverage_fraction, multi_fraction)) +
  geom_hex(bins = 60) +
  theme_bw() +
  labs(x = "Coverage fraction (positions with ≥1 ECF family)",
       y = "Multi-hit fraction",
       title = "Coverage vs overlap: are overlaps driven by sparse mapping?")

ggsave(file.path(this.part.dir, "plot3_coverage_vs_overlap.png"), p3, width = 6, height = 5)

p4 <- ggplot(pmf, aes(log_qlength, mean_eff_families)) +
  geom_hex(bins = 60) +
  theme_bw() +
  labs(x = "log10(protein length + 1)",
       y = "Mean effective families",
       title = "Length dependence of family mixing")

ggsave(file.path(this.part.dir, "plot4_length_vs_mixing.png"), p4, width = 6, height = 5)

# keep fragmentation levels not too insane (optional cap)
pmf_cov[, frag_cap := pmin(fragmentation_index, 10)]
pmf_cov[, frag_cap := factor(frag_cap, levels = sort(unique(frag_cap)))]

p5 <- ggplot(pmf_cov, aes(frag_cap, mean_eff_families)) +
  geom_boxplot(outlier.size = 0.2) +
  theme_bw() +
  labs(x = "Fragmentation index (capped at 10)",
       y = "Mean effective families",
       title = "Does fragmentation add information beyond mixing?")

ggsave(file.path(this.part.dir, "plot5_fragmentation_vs_mixing.png"), p5, width = 6, height = 5)


pm_long <- melt(
  pmf,
  id.vars = c("qname", "class"),
  measure.vars = c("coverage_fraction", "multi_fraction", "mean_eff_families"),
  variable.name = "metric",
  value.name = "value"
)

p6 <- ggplot(pm_long[!is.na(class)], aes(x = class, y = value)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.size = 0.2) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = NULL, title = "Key architecture metrics by class")

ggsave(file.path(this.part.dir, "plot6_metrics_by_class.png"), p6, width = 8, height = 9)
