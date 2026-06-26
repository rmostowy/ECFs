### load libraries -------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(igraph)
  library(ggplot2)
  library(ggraph)
  library(patchwork)
  library(broom)
  library(ggh4x)
  library(ggExtra)
}))

### load config -----------------------------------------------------------------
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

#################################
### Set up output directory
#################################
## Single output directory for the merged subsection.
this.part.dir <- file.path(ecf.proj.dir, "part2_ecfs-vs-domains")
make.dir(this.part.dir)

#################################
### Load all data:
### 1. ECF tables
### 2. ECF comparison tables (containment, former part2)
### 3. domain tables
### 4. ECF vs domain tables
### 5. reprseq vs ecfs
### 6. reprseq vs pfam
### 7. reprseq functional annotations
#################################

# ECF tables
ecf.table         <- fread(ecf.table.output)
ecf.family.table  <- fread(ecf.family.table.output)

# ECF comparison tables (containment, former part2)
ecfs.vs.ecfs      <- fread(ecfs.vs.ecfs.output)       # raw all-vs-all
ecfs.vs.ecfs.cov  <- fread(ecfs.vs.ecfs.cov.output)   # processed coverage

# domain tables
ecod.domains      <- fread(ecod.domains.output)
pfam.domains      <- fread(pfam.domains.output)

# ECF vs domain tables
ecfs.vs.ecod          <- fread(ecfs.vs.ecod.output)
ecfs.vs.pfam          <- fread(ecfs.vs.pfam.output)
ecfs.vs.ecod.hhsenser <- fread(ecfs.vs.ecod.hhsenser.output)
ecfs.vs.pfam.hhsenser <- fread(ecfs.vs.pfam.hhsenser.output)

# reprseq vs ecfs
reprseq.vs.ecfs.all   <- fread(reprseq.ecfs.output)

# reprseq vs pfam
reprseq.vs.ecod.all   <- fread(reprseq.domains.ecod.output)
reprseq.vs.pfam.all   <- fread(reprseq.domains.pfam.output)


# reprseq functional annotations
reprseq.functions.unique <- fread(reprseq.functions.unique.output)
reprseq.functions.all    <- fread(reprseq.functions.all.output)

#################################
### Process loaded tables
#################################
## Keys and lookups used throughout. Both source files keyed `ecf.family.table`
## on `ecf` and `ecf.table` on `ecf.name` -- both are set explicitly here so
## every downstream `.( )` lookup is unambiguous regardless of analysis order.
setkey(ecf.table,                  ecf.name)
setkey(ecf.family.table,           ecf)
setkey(ecod.domains,               full.id)
setkey(pfam.domains,               accession)
setkey(reprseq.functions.unique,   reprseq)

all.ecfs <- unique(ecf.table$ecf.name)
Ne <- length(all.ecfs)

## Build the directed ECF coverage network from high-confidence alignments
## (former part2 process step).
ecfs.vs.ecfs[, max.cov := pmax(qcov, scov)]
setkey(ecfs.vs.ecfs, qname, sname)

ecfs.vs.ecfs.cov$qlength <- ecf.table[.(ecfs.vs.ecfs.cov$query)]$length
ecfs.vs.ecfs.cov$slength <- ecf.table[.(ecfs.vs.ecfs.cov$subject)]$length
ecfs.vs.ecfs.cov[, qcov := pmin(qcov.min, qcov.max)]
ecfs.vs.ecfs.cov[, scov := pmin(scov.min, scov.max)]
ecfs.vs.ecfs.cov[, cov.min := pmin(scov.min, qcov.min)]
ecfs.vs.ecfs.cov[, cov.max := pmax(scov.min, qcov.min)]
eve.conf <- ecfs.vs.ecfs.cov[prob >= CONTAINMENT.PROB.THR]

ecf.network.ab <- eve.conf %>% select(query, subject, prob, qcov, scov)
ecf.network.ba <- eve.conf %>% select(query = subject, subject = query,
                                      prob, qcov = scov, scov = qcov)
ecf.network <- rbind(ecf.network.ab, ecf.network.ba)

## Attach ECOD hierarchy levels and confidence-filter the ECF-vs-domain hit
## tables (former part4 process step; used in Steps 2 and 3 below).
ecfs.vs.ecod[, x.id    := ecod.domains[.(ecfs.vs.ecod$sname)]$x_id]
ecfs.vs.ecod[, h.id    := ecod.domains[.(ecfs.vs.ecod$sname)]$h_id]
ecfs.vs.ecod[, t.id    := ecod.domains[.(ecfs.vs.ecod$sname)]$t_id]
ecfs.vs.ecod[, x.level := ecod.domains[.(ecfs.vs.ecod$sname)]$x.level]
ecfs.vs.ecod[, h.level := ecod.domains[.(ecfs.vs.ecod$sname)]$h.level]
ecfs.vs.ecod[, t.level := ecod.domains[.(ecfs.vs.ecod$sname)]$t.level]

ecfs.vs.ecod.conf <- ecfs.vs.ecod %>%
  filter(eval <= MAX.EVALUE.FOR.ECF.DOMAIN.CONFIDENT &
           prob >= MIN.PROB.FOR.ECF.DOMAIN.CONFIDENT)
ecfs.vs.pfam.conf <- ecfs.vs.pfam %>%
  filter(eval <= MAX.EVALUE.FOR.ECF.DOMAIN.CONFIDENT &
           prob >= MIN.PROB.FOR.ECF.DOMAIN.CONFIDENT)

## Confidence-filter the rHMM-vs-domain hit tables (former part4 process step;
## used in Step 3 below).
reprseq.vs.ecfs <- reprseq.vs.ecfs.all %>%
  filter(prob >= MIN.PROB.FOR.DOMAIN.DETECT.IN.RHMM,
         eval <= MAX.EVAL.FOR.DOMAIN.DETECT.IN.RHMM,
         qcov >= MIN.RHMM.COV.FOR.DOMAIN.DETECT.IN.RHMM,
         scov >= MIN.DOMAIN.COV.FOR.DOMAIN.DETECT.IN.RHMM,
         hit.length >= MIN.HIT.LENGTH.FOR.DOMAIN.DETECT.IN.RHMM)

reprseq.vs.pfam <- reprseq.vs.pfam.all %>%
  filter(prob >= MIN.PROB.FOR.DOMAIN.DETECT.IN.RHMM,
         eval <= MAX.EVAL.FOR.DOMAIN.DETECT.IN.RHMM,
         qcov >= MIN.RHMM.COV.FOR.DOMAIN.DETECT.IN.RHMM,
         scov >= MIN.DOMAIN.COV.FOR.DOMAIN.DETECT.IN.RHMM,
         hit.length >= MIN.HIT.LENGTH.FOR.DOMAIN.DETECT.IN.RHMM)

reprseq.vs.ecod <- reprseq.vs.ecod.all %>%
  filter(prob >= MIN.PROB.FOR.DOMAIN.DETECT.IN.RHMM,
         eval <= MAX.EVAL.FOR.DOMAIN.DETECT.IN.RHMM,
         qcov >= MIN.RHMM.COV.FOR.DOMAIN.DETECT.IN.RHMM,
         scov >= MIN.DOMAIN.COV.FOR.DOMAIN.DETECT.IN.RHMM,
         hit.length >= MIN.HIT.LENGTH.FOR.DOMAIN.DETECT.IN.RHMM)

## Shared plotting helpers. Not exactly "table processing", but defined once
## here and reused across every step below -- same role `base_theme` played
## in the original per-figure scripts.
base_theme <- theme_bw(base_size = 14) +
  theme(
    axis.title   = element_text(face = "bold"),
    axis.text    = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    plot.title   = element_text(face = "bold", hjust = 0)
  )

heatmap_fill <- scale_fill_viridis_c(
  option = "cividis",
  trans  = "sqrt",
  name   = "Count"
)

## Shared class ordering used by all ECF-vs-domain class plots (Steps 2 and 4).
class.plotting.order <- c("equivalent", "subdomain", "superdomain", "related", "distant", "new")
class.plotting.cols <- c(
  setNames(RColorBrewer::brewer.pal(5, "Set2"),
           rev(c("equivalent", "subdomain", "superdomain", "related", "distant"))),
  new = "grey70"
)

####################################################################
####################################################################
###  STEP 1 -- ECFs have a hierchical and nested structure
####################################################################
####################################################################

## Per-ECF connectivity / asymmetry summary.
dt.net <- copy(ecf.network)
dt.net[, subject.ef := ecf.family.table[.(subject)]$family]
ecf.no.fams <- dt.net[, .(no.fam = uniqueN(subject.ef)), by = query]

ecf.connectivity.stats <- dt.net[, .(
  mean.qcov = mean(qcov),
  mean.scov = mean(scov),
  med.qcov  = median(qcov),
  med.scov  = median(scov),
  degree    = .N
), by = query]
ecf.connectivity.stats[, family := ecf.family.table[.(ecf.connectivity.stats$query)]$family]
ecf.connectivity.stats[, length := ecf.table[.(ecf.connectivity.stats$query)]$length]
ecf.connectivity.stats[, asym_mean     := mean.qcov - mean.scov]
ecf.connectivity.stats[, asym_med      := med.qcov  - med.scov]
ecf.connectivity.stats[, strength_mean := (mean.qcov + mean.scov) / 2]
ecf.connectivity.stats[, strength_med  := (med.qcov  + med.scov)  / 2]
ecf.connectivity.stats <- merge(ecf.connectivity.stats, ecf.no.fams, by = "query")

## ---- PLOT 1: per-node median qcov/scov diagonal (p.qs.node) ----------------
## The clean per-node containment view: most ECFs sit hard against one axis,
## i.e. they are predominantly contained-in or container-of, not symmetric.
p.qs.node <- ggplot(ecf.connectivity.stats, aes(x = med.qcov, y = med.scov)) +
  geom_point(size = 0) +
  geom_bin2d(bins = 70) +
  heatmap_fill +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Median query coverage", y = "Median subject coverage") +
  base_theme
p.qs.node.marg <- ggMarginal(
  p.qs.node, type = "histogram", bins = 60,
  fill = "#16396b", colour = "white", size = 5
)
ggsave(file.path(this.part.dir, "plot1_node-qcov-scov.png"),
       p.qs.node.marg, width = 6, height = 5)

## ---- PLOT 2: raw pairwise heatmap (p.cov.heatmap) -------------------------
p.cov.heatmap <- ggplot(eve.conf, aes(x = qcov, y = scov)) +
  geom_bin2d(bins = 60) +
  heatmap_fill +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Query coverage", y = "Subject coverage") +
  base_theme
ggsave(file.path(this.part.dir, "plot2_cov-heatmap.png"),
       p.cov.heatmap, width = 6, height = 5)

## ---- PLOT 3: per-node average heatmap (p.connect.cov.heatmap) -------------
## Doubles as evidence that containment is global, not hub-driven.
p.connect.cov.heatmap <- ggplot(ecf.connectivity.stats,
                                aes(x = mean.qcov, y = mean.scov)) +
  geom_bin2d(bins = 60) +
  heatmap_fill +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Mean query coverage", y = "Mean subject coverage") +
  base_theme
ggsave(file.path(this.part.dir, "plot3_connect-cov-heatmap.png"),
       p.connect.cov.heatmap, width = 6, height = 5)

## ---- PLOT 4: asymmetry distribution (p.asym.dist) -------------------------
p.asym.dist <- ggplot(ecf.connectivity.stats, aes(x = asym_med)) +
  geom_histogram(bins = 60, fill = "#2c7fb8", color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 14) +
  labs(x = "Asymmetry (median qcov - median scov)", y = "Number of ECFs")
ggsave(file.path(this.part.dir, "plot4_asym-dist.png"),
       p.asym.dist, width = 5, height = 4)

## ---- PLOT 5: asymmetry vs length (p.asym.len) -----------------------------
## Shows that the contained-in (fragment) role concentrates among shorter ECFs.
p.asym.len <- ggplot(ecf.connectivity.stats, aes(x = asym_mean, y = length)) +
  geom_bin2d(bins = 60) +
  heatmap_fill +
  scale_y_continuous(trans = "log10") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6) +
  labs(x = "Asymmetry (mean qcov - mean scov)", y = "ECF length (aa, log10)") +
  base_theme
ggsave(file.path(this.part.dir, "plot5_asym-vs-length.png"),
       p.asym.len, width = 6, height = 5)

####################################################################
####################################################################
###  STEP 2 -- This is tail is due to ECFs contained in other ones.
###  Those fragment ECFs are shorter and connect unrelated ECFs
####################################################################
####################################################################


## ---- container-container similarity (p.cont.cont) -----------------
## Supports Step 4 (artefact control / the 11% result) but is computed here
## because it needs the containment graph. Saved with the S4 group.
containment.network <- ecf.network[qcov >= DOMAIN.CONTAINMENT.COV.THR]
containment.graph   <- graph_from_data_frame(containment.network, directed = TRUE)
containment.graph.nodes <- V(containment.graph)$name
no.nodes.cg <- length(containment.graph.nodes)


contain.neighbour.cov.max.list <- lapply(seq_len(no.nodes.cg), function(node.index) {
  out <- NULL
  this.node  <- containment.graph.nodes[node.index]
  neighbours <- neighbors(containment.graph, this.node, mode = "out")
  if (length(neighbours) >= 2) {
    neighbour.pairs    <- combn(neighbours$name, 2, simplify = FALSE)
    neighbour.pairs.dt <- as.data.table(do.call(rbind, neighbour.pairs))
    neighbour.pair.weigths <- ecfs.vs.ecfs[.(neighbour.pairs.dt)]$max.cov
    neighbour.pair.weigths <- fcoalesce(neighbour.pair.weigths, 0)
    out <- data.table(ecf = this.node, neighbour.cov.max = neighbour.pair.weigths)
  }
  out
})
contain.neighbour.cov.max <- rbindlist(contain.neighbour.cov.max.list)

contain.neighbour.no.fams.list <- lapply(seq_len(no.nodes.cg), function(node.index) {
  out <- NULL
  this.node  <- containment.graph.nodes[node.index]
  neighbours <- neighbors(containment.graph, this.node, mode = "out")
  if (length(neighbours) >= 2) {
    neighbour.fams <- ecf.family.table[.(neighbours$name)]$family
    no.fams <- length(unique(neighbour.fams))
    out <- data.table(ecf = this.node, no.fams = no.fams)
  }
  out 
})
contain.neighbour.no.fams <- rbindlist(contain.neighbour.no.fams.list)

## Decile bins on med.scov among contained ECFs (equal sample size per bin,
## not equal width) -- same logic as the `len_bin` quantile cut used for
## `ecf.ecod.domain.assignment` above.
contained.ecfs <- ecf.connectivity.stats[med.qcov >= DOMAIN.CONTAINMENT.COV.THR]

scov.bin.breaks <- quantile(contained.ecfs$med.scov, probs = seq(0, 1, by = 0.1), na.rm = TRUE)
scov.bin.breaks <- unique(scov.bin.breaks)   # guard against duplicate breaks from ties
n.scov.bins <- length(scov.bin.breaks) - 1

contained.ecfs[, scov.decile := cut(med.scov, breaks = scov.bin.breaks,
                                    include.lowest = TRUE, labels = FALSE)]
contained.ecfs <- contained.ecfs[!is.na(scov.decile)]

fragment.cov.stats <- contained.ecfs[, {
  this.range.cov <- contain.neighbour.cov.max[ecf %chin% query]$neighbour.cov.max
  mean.no.fams   <- mean(contain.neighbour.no.fams[ecf %chin% query]$no.fams)
  .(
    size           = .N,
    mean.length    = mean(length),
    prop.unrelated = lw(this.range.cov < 0.5) / length(this.range.cov),
    mean.no.fams   = mean.no.fams,
    scov.min       = min(med.scov),
    scov.max       = max(med.scov),
    scov.med       = median(med.scov)
  )
}, by = scov.decile]
setorder(fragment.cov.stats, scov.decile)

fragment.cov.stats[, scov.mid := scov.med]               # bin position for plotting
fragB <- fragment.cov.stats[scov.decile != n.scov.bins]   # drop top decile (near-equivalent matches) — same intent as the old `scov.max <= 0.9` cut

col_unrel <- "#3B6D11"

## ---- PLOT 6: Do fragments connect more diverse ECFs?

p.fragments.connect <- ggplot(fragB, aes(x = scov.mid, y = prop.unrelated)) +
  geom_line(colour = col_unrel, linewidth = 0.7) +
  geom_point(aes(size = mean.length), colour = col_unrel) +
  scale_size(range = c(2, 13), limits = range(fragB$mean.length),
             breaks = pretty(fragB$mean.length, n = 4),
             name = "Mean ECF\nsize") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8), limits = c(0.2, 0.8)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Median subject coverage  (low \u2192 more fragment like)",
       y = "Proportion of unrelated neighbours (cov < 0.5)") +
  base_theme

ggsave(file.path(this.part.dir, "plot6_fragments-connect.png"),
       p.fragments.connect, width = 7, height = 4.5)

## ---- PLOT 7: full distribution of container-container
## similarity within each scov decile -----------------------------------------
## `prop.unrelated` in the main panel collapses each bin to one number (% of
## neighbour pairs below the 0.5 "unrelated" threshold). This shows the
## underlying distribution directly -- the point being that near-zero
## (fully unrelated) pairs show up across the *whole* range of med.scov,
## not just at one extreme.
neighbour.cov.by.decile <- contained.ecfs[scov.decile != n.scov.bins, {
  data.table(neighbour.cov.max = contain.neighbour.cov.max[ecf %chin% query]$neighbour.cov.max)
}, by = scov.decile]

## reuse the same bin midpoints as fragB so the two figures read consistently
decile.labels <- fragB[, .(scov.decile, scov.mid)]
neighbour.cov.by.decile <- merge(neighbour.cov.by.decile, decile.labels, by = "scov.decile")
neighbour.cov.by.decile[, scov.decile := factor(scov.decile, levels = sort(unique(scov.decile)))]

label.order <- neighbour.cov.by.decile[, .(scov.mid = unique(scov.mid)), by = scov.decile]
setorder(label.order, scov.decile)

p.fragments.connect.dist <- ggplot(neighbour.cov.by.decile,
                                   aes(x = scov.decile, y = neighbour.cov.max)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey40") +
  geom_boxplot(fill = "#C0DD97", colour = col_unrel, alpha = 0.6,
               outlier.shape = NA, linewidth = 0.4) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.15, size = 0.6, colour = col_unrel) +
  scale_x_discrete(labels = sprintf("%.2f", label.order$scov.mid)) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.01, 0.03))) +
  labs(x = "Median subject coverage  (low \u2192 more fragment like)",
       y = "Container-container similarity\n(max coverage)") +
  base_theme +
  theme(panel.grid.minor = element_blank())

# p.fragments.connect.dist

ggsave(file.path(this.part.dir, "plot7_fragments-cov-distribution.png"),
       p.fragments.connect.dist, width = 7, height = 4.5)


# p.cont.cont <- ggplot(contain.neighbour.cov.max, aes(x = neighbour.cov.max)) +
#   geom_histogram(bins = 40, fill = "#335BAA", color = "white", linewidth = 0.2) +
#   geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.6) +
#   scale_x_continuous(breaks = seq(0, 1, 0.2), expand = c(0.03, 0.03)) +
#   coord_cartesian(xlim = c(0, 1)) +
#   labs(x = "Container-container similarity (max coverage)",
#        y = "Number of container pairs") +
#   theme_classic(base_size = 14) +
#   theme(axis.title = element_text(face = "bold"),
#         axis.text  = element_text(color = "black"),
#         plot.margin = margin(10, 10, 10, 10))
# ggsave(file.path(this.part.dir, "suppS4_container-container.png"),
#        p.cont.cont, width = 6, height = 4)

####################################################################
####################################################################
###  STEP 3 -- ECFs form a spectrum of relationships to domains
####################################################################
####################################################################

## Reusable ECF -> domain classifier (was duplicated verbatim for ECOD and
## Pfam in the part4 source; factored into one function here to avoid
## divergence).
classify.ecf.domain <- function(mt) {
  mt <- copy(mt)
  mt[, conf := (prob >= MIN.PROB.FOR.DOMAIN.RELATEDNESS)]
  mt[, is.equiv := conf & (prob >= MIN.PROB.FOR.DOMAIN.EQUIVALENCE &
                             qcov >= MIN.QCOV.FOR.DOMAIN.EQUIVALENCE &
                             scov >= MIN.SCOV.FOR.DOMAIN.EQUIVALENCE &
                             eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE)]
  mt[, is.subdom := conf & (prob >= MIN.PROB.FOR.SUBDOMAIN &
                              qcov >= MIN.QCOV.FOR.SUBDOMAIN &
                              scov <= MAX.SCOV.FOR.SUBDOMAIN &
                              eval <= MAX.EVAL.FOR.SUBDOMAIN)]
  mt[, is.superdom := conf & (prob >= MIN.PROB.FOR.SUPERDOMAIN &
                                qcov <= MAX.QCOV.FOR.SUPERDOMAIN &
                                scov >= MIN.SCOV.FOR.SUPERDOMAIN &
                                eval <= MAX.EVAL.FOR.SUPERDOMAIN)]
  mt[, {
    is.conf           <- any(conf)
    is.equiv.class    <- any(is.equiv)
    is.subdom.class   <- any(is.subdom)
    is.superdom.class <- any(is.superdom)
    class.assignment  <- "distant"
    if (is.conf) {
      class.assignment <- "related"
      if (is.equiv.class) {
        class.assignment <- "equivalent"
      } else if (is.subdom.class & !is.superdom.class) {
        class.assignment <- "subdomain"
      } else if (!is.subdom.class & is.superdom.class) {
        class.assignment <- "superdomain"
      } else if (is.subdom.class & is.superdom.class) {
        ambiguous.rows <- which(is.subdom | is.superdom)
        ambiguous.row.with.longest.hit <- ambiguous.rows[which.max(hit.length[ambiguous.rows])]
        this.is.subdomain <- is.subdom[ambiguous.row.with.longest.hit]
        class.assignment <- if (this.is.subdomain) "subdomain" else "superdomain"
      }
    }
    .(class.assignment = class.assignment)
  }, by = qname]
}

## ECOD classification by ECF length.
ecf.domain.assignment.ecod <- classify.ecf.domain(ecfs.vs.ecod)
ecf.no.domain.assignment <- all.ecfs[!all.ecfs %in% ecf.domain.assignment.ecod$qname]
ecf.domain.assignment.ecod <- rbind(ecf.domain.assignment.ecod, data.table(qname = ecf.no.domain.assignment, class.assignment = "new"))

ecf.ecod.domain.assignment <- merge(
  ecf.domain.assignment.ecod,
  ecf.table[, .(qname = ecf.name, length, size)],
  by = "qname"
)
ecf.ecod.domain.assignment[, len_bin := cut(
  length,
  breaks = round(quantile(length, probs = seq(0, 1, by = 0.1), na.rm = TRUE)),
  include.lowest = TRUE, dig.lab = 10
)]

prop_dt.ecod <- ecf.ecod.domain.assignment[, .N, by = .(len_bin, class.assignment)]
prop_dt.ecod[, prop := N / sum(N), by = len_bin]
prop_dt.ecod[, class.assignment := factor(class.assignment, levels = rev(class.plotting.order))]

## ---- PLOT 8: ECOD class by length (p.ecod.domain.cat.by.length) ------------
## Read two ways: by length -> subdomain/equivalent/superdomain spectrum;
## by row -> exact equivalence is never dominant.



p.ecod.domain.cat.by.length <- ggplot(prop_dt.ecod,
                                      aes(x = len_bin, y = prop, fill = class.assignment)) +
  geom_col(width = 0.9) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = class.plotting.cols) +
  labs(x = "ECF length (quantile bins)", y = "Proportion", fill = "ECF-ECOD\nclass") +
  theme_bw(base_size = 12) +
  guides(fill = guide_legend(reverse = TRUE, ncol = 3, byrow = TRUE)) +
  theme(legend.position = "top",
        legend.text  = element_text(size = 9),
        legend.title = element_text(size = 9, face = "bold"),
        legend.key.width  = unit(0.5, "cm"),
        legend.key.height = unit(0.5, "cm"),
        axis.title   = element_text(face = "bold"),
        axis.text.x  = element_text(angle = 45, hjust = 1))
ggsave(file.path(this.part.dir, "plot8_dom-class-ecod.png"),
       p.ecod.domain.cat.by.length, width = 6, height = 5)

## ---- PLOT 9: Pfam class by length (p.pfam.domain.cat.by.length) -----------
ecf.domain.assignment.pfam <- classify.ecf.domain(ecfs.vs.pfam)
ecf.no.domain.assignment.pfam <- all.ecfs[!all.ecfs %in% ecf.domain.assignment.pfam$qname]
ecf.domain.assignment.pfam <- rbind(ecf.domain.assignment.pfam, data.table(qname = ecf.no.domain.assignment.pfam, class.assignment = "new"))


ecf.pfam.domain.assignment <- merge(
  ecf.domain.assignment.pfam,
  ecf.table[, .(qname = ecf.name, length, size)],
  by = "qname"
)
ecf.pfam.domain.assignment[, len_bin := cut(
  length,
  breaks = quantile(length, probs = seq(0, 1, by = 0.1), na.rm = TRUE),
  include.lowest = TRUE
)]
prop_dt.pfam <- ecf.pfam.domain.assignment[, .N, by = .(len_bin, class.assignment)]
prop_dt.pfam[, prop := N / sum(N), by = len_bin]
prop_dt.pfam[, class.assignment := factor(class.assignment, levels = rev(class.plotting.order))]

p.pfam.domain.cat.by.length <- ggplot(prop_dt.pfam,
                                      aes(x = len_bin, y = prop, fill = class.assignment)) +
  geom_col(width = 0.9) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = class.plotting.cols) +
  labs(x = "ECF length (quantile bins)", y = "Proportion", fill = "Class") +
  theme_bw(base_size = 12) +
  theme(axis.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(this.part.dir, "plot9_dom-class-pfam.png"),
       p.pfam.domain.cat.by.length, width = 6, height = 5)

## ---- PLOT 10: ECOD-Pfam agreement (p.ecod.pfam.agreement) ------------------
ecod.assign <- ecf.ecod.domain.assignment[, .(qname, ecod_class = class.assignment, length, size)]
pfam.assign <- ecf.pfam.domain.assignment[, .(qname, pfam_class = class.assignment, length, size)]

all_ecf <- data.table(qname = sort(unique(c(ecod.assign$qname, pfam.assign$qname))))
all_ecf <- merge(all_ecf, ecf.table[, .(qname = ecf.name, length, size)], by = "qname")
dt.agree <- merge(all_ecf, ecod.assign, by = c("qname", "length", "size"), all.x = TRUE)
dt.agree <- merge(dt.agree, pfam.assign, by = c("qname", "length", "size"), all.x = TRUE)
dt.agree[, ecod_class := factor(ecod_class, levels = class.plotting.order)]
dt.agree[, pfam_class := factor(pfam_class, levels = class.plotting.order)]

cm <- dt.agree[, .N, by = .(ecod_class, pfam_class)]
p_cm_counts <- ggplot(cm, aes(x = pfam_class, y = ecod_class, fill = N)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(trans = "sqrt") +
  labs(x = "Pfam assignment class", y = "ECOD assignment class", fill = "Count") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
cm[, prop_row := N / sum(N), by = ecod_class]
p_cm_rowprop <- ggplot(cm, aes(x = pfam_class, y = ecod_class, fill = prop_row)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Pfam assignment class", y = "ECOD assignment class", fill = "Row %") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p.ecod.pfam.agreement <- (p_cm_counts / p_cm_rowprop) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 24),
        plot.tag.position = c(0.02, 0.98))
ggsave(file.path(this.part.dir, "plot10_ecod-pfam-agreement.png"),
       p.ecod.pfam.agreement, width = 6, height = 10)

## chi-squared agreement -> inline prose sentence (not a figure).
tab.agree <- table(dt.agree$ecod_class, dt.agree$pfam_class)
N.agree <- sum(tab.agree)
r <- nrow(tab.agree); cc <- ncol(tab.agree)
chi <- suppressWarnings(chisq.test(tab.agree))
V <- sqrt(unname(chi$statistic) / (N.agree * min(r - 1, cc - 1)))
pval_txt <- if (is.na(chi$p.value)) "NA" else if (chi$p.value < 2.2e-16) "< 2.2e-16" else
  format(chi$p.value, scientific = TRUE, digits = 2)
sent <- sprintf(
  paste0("ECOD- and Pfam-based assignment classes show a clear correspondence across ECFs ",
         "(chi-squared: X2(%d) = %.0f, p %s; Cramer's V = %.3f)."),
  unname(chi$parameter), unname(chi$statistic), pval_txt, V
)
write(sent, file = file.path(this.part.dir, "plot10-table_agreement-sentence.txt"))

## ---- plot11: ECF-domain coverage across e-value thresholds -----------------
ecf.domain.coverage <- data.table(ecf = all.ecfs)
ecf.domain.coverage[, ecod.eval0p1   := all.ecfs %in% ecfs.vs.ecod$qname]
ecf.domain.coverage[, ecod.eval0p001 := all.ecfs %in% ecfs.vs.ecod.conf$qname]
ecf.domain.coverage[, ecod.eval10    := all.ecfs %in% ecfs.vs.ecod.hhsenser$qname]
ecf.domain.coverage[, pfam.eval0p1   := all.ecfs %in% ecfs.vs.pfam$qname]
ecf.domain.coverage[, pfam.eval0p001 := all.ecfs %in% ecfs.vs.pfam.conf$qname]
ecf.domain.coverage[, pfam.eval10    := all.ecfs %in% ecfs.vs.pfam.hhsenser$qname]

assign_cat <- function(has_ecod, has_pfam) {
  fifelse(has_ecod & has_pfam, "Both (ECOD + Pfam)",
          fifelse(has_ecod & !has_pfam, "ECOD only",
                  fifelse(!has_ecod & has_pfam, "Pfam only", "Neither")))
}
long_dt <- rbindlist(list(
  data.table(ecf = ecf.domain.coverage$ecf, threshold = "1e-3",
             has_ecod = ecf.domain.coverage$ecod.eval0p001, has_pfam = ecf.domain.coverage$pfam.eval0p001),
  data.table(ecf = ecf.domain.coverage$ecf, threshold = "1e-1",
             has_ecod = ecf.domain.coverage$ecod.eval0p1,  has_pfam = ecf.domain.coverage$pfam.eval0p1),
  data.table(ecf = ecf.domain.coverage$ecf, threshold = "10",
             has_ecod = ecf.domain.coverage$ecod.eval10,   has_pfam = ecf.domain.coverage$pfam.eval10)
))
long_dt[, category := assign_cat(has_ecod, has_pfam)]
counts_dt <- long_dt[, .N, by = .(threshold, category)]
counts_dt[, threshold := factor(threshold, levels = c("1e-3", "1e-1", "10"))]
counts_dt[, category  := factor(category,
                                levels = c("Both (ECOD + Pfam)", "ECOD only", "Pfam only", "Neither"))]
cat_cols <- c("Both (ECOD + Pfam)" = "#009E73", "ECOD only" = "#0072B2",
              "Pfam only" = "#D55E00", "Neither" = "grey70")

p.ecf.domain.coverage <- ggplot(counts_dt, aes(x = threshold, y = N, fill = category)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.25) +
  scale_fill_manual(values = cat_cols, name = NULL) +
  labs(x = "E-value threshold (strict -> lenient)", y = "Number of ECFs") +
  theme_bw(base_size = 12) +
  guides(fill = guide_legend(ncol = 2)) +
  theme(legend.position = "top",
        legend.text  = element_text(size = 9),
        axis.title   = element_text(face = "bold"))
ggsave(file.path(this.part.dir, "plot11_ecf-domain-cov.png"),
       p.ecf.domain.coverage, width = 6, height = 5)

## ---- Define "new" (off-edge) ECFs and their rHMMs (shared by D and supp) ----
new.ecfs    <- ecf.domain.coverage[ecod.eval0p1 == FALSE, unique(ecf)]
reprseq.new <- unique(reprseq.vs.ecfs[sname %chin% new.ecfs, qname])

## ---- MAIN D : PHROG functional enrichment of new-ECF rHMMs (ECOD-based) -----
rs2funct <- copy(reprseq.functions.unique)
setnames(rs2funct, "reprseq", "qname")
rs2funct[, is_new  := qname %chin% reprseq.new]
rs2funct[, funct2  := fifelse(is.na(funct) | funct == "", "unknown", funct)]
rs2funct[, class2  := fifelse(is.na(class) | class == "", "unknown", class)]
rs2funct.known <- rs2funct[funct2 != "unknown"]

enrichment.phrogs <- rs2funct.known[
  , {
    a <- sum(is_new)
    b <- .N - a
    cnt <- sum(rs2funct.known$is_new) - a
    d <- (nrow(rs2funct.known) - sum(rs2funct.known$is_new)) - b
    ft <- fisher.test(matrix(c(a, b, cnt, d), nrow = 2))
    .(n_new = a, n_total = .N, odds = unname(ft$estimate), p = ft$p.value)
  },
  by = .(funct2, class2)
]
enrichment.phrogs[, `:=`(fdr = p.adjust(p, "BH"),
                         log2OR = log2(odds),
                         frac_new = n_new / n_total)]

FDR_CUT <- 0.05; NTOP <- 6; CAP <- 7
plot.ecod.dt <- enrichment.phrogs[order(fdr)][fdr <= FDR_CUT, head(.SD, NTOP), by = class2]
plot.ecod.dt[, class2 := factor(class2, levels = phrogs.class.order)]
plot.ecod.dt[, ykey := paste(funct2, class2, sep = "___")]
plot.ecod.dt[, ykey := reorder(ykey, log2OR)]
plot.ecod.dt[, log2OR_plot := log2OR]
plot.ecod.dt[is.infinite(log2OR) & log2OR > 0, log2OR_plot :=  CAP]
plot.ecod.dt[is.infinite(log2OR) & log2OR < 0, log2OR_plot := -CAP]
plot.ecod.dt[, is_inf := !is.finite(log2OR)]

strip <- ggh4x::strip_themed(
  background_y = ggh4x::elem_list_rect(
    fill = grDevices::adjustcolor(phrogs.col[match(phrogs.class.order, phrogs.class)], alpha.f = 0.25),
    colour = "black", linewidth = 0.4)
)
p.new.ecf.enrichment <- ggplot(plot.ecod.dt, aes(x = log2OR_plot, y = ykey, shape = is_inf)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
  geom_point(aes(size = n_total, colour = -log10(fdr)), alpha = 0.9) +
  ggh4x::facet_grid2(rows = vars(class2), scales = "free_y", space = "free_y",
                     labeller = labeller(class2 = phrog.class.labels), strip = strip) +
  scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
  scale_colour_gradient(low = "#D7301F", high = "darkblue", trans = "sqrt",
                        limits = c(0, max(-log10(plot.ecod.dt$fdr), na.rm = TRUE)),
                        name = expression(-log[10]("FDR"))) +
  scale_size_continuous(name = "N reprseqs", range = c(1.5, 5)) +
  scale_x_continuous(breaks = function(x) pretty(x, n = 7)) +
  labs(x = "Functional enrichment in new-ECF rHMMs\n(log2 odds ratio)", y = NULL, shape = "Infinite ORs") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        strip.placement = "outside",
        strip.text.y = element_text(angle = 0, hjust = 0, face = "bold", size = 9),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(face = "bold"),
        panel.spacing = unit(0.35, "lines"),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = .25))
ggsave(file.path(this.part.dir, "plot12_new-ecf-enrichment.png"),
       p.new.ecf.enrichment, width = 8, height = 10)

## ---- PLOT 13: same enrichment, Pfam-based (p.new.ecf.enrichment.pfams) -----
rs2pfam <- unique(reprseq.vs.pfam[, .(qname, accession)])
rs2pfam[, is_new := qname %chin% reprseq.new]

enrichment.pfams <- rs2pfam[
  , {
    a <- sum(is_new)
    b <- .N - a
    cnt <- sum(rs2pfam$is_new) - a
    d <- (nrow(rs2pfam) - sum(rs2pfam$is_new)) - b
    ft <- fisher.test(matrix(c(a, b, cnt, d), 2))
    .(n_new = a, n_total = .N, odds = unname(ft$estimate), p = ft$p.value)
  },
  by = accession
]
enrichment.pfams[, `:=`(fdr = p.adjust(p, "BH"), log2OR = log2(odds), frac_new = n_new / n_total)]
SHOW_N <- 40L
plot.pfam.dt <- enrichment.pfams[order(fdr)][fdr <= FDR_CUT, head(.SD, SHOW_N)]
CAP2 <- ceiling(max(abs(range(plot.pfam.dt$log2OR, finite = TRUE))))
plot.pfam.dt[, log2OR_plot := log2OR]
plot.pfam.dt[is.infinite(log2OR) & log2OR > 0, log2OR_plot :=  CAP2]
plot.pfam.dt[is.infinite(log2OR) & log2OR < 0, log2OR_plot := -CAP2]
plot.pfam.dt[, is_inf := !is.finite(log2OR)]
plot.pfam.dt[, funct  := pfam.domains[.(plot.pfam.dt$accession)]$description]

p.new.ecf.enrichment.pfams <- ggplot(plot.pfam.dt, aes(x = log2OR_plot, y = funct, shape = is_inf)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
  geom_point(aes(size = n_total, colour = -log10(fdr)), alpha = 0.9) +
  scale_colour_gradient(low = "#D7301F", high = "darkblue", trans = "sqrt",
                        name = expression(-log[10]("FDR"))) +
  scale_size_continuous(name = "N reprseqs", range = c(1.5, 5)) +
  scale_x_continuous(breaks = function(x) pretty(x, n = 7)) +
  labs(x = "Pfam enrichment in new-ECF rHMMs\n(log2 odds ratio)", y = NULL, shape = "Infinite ORs") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(face = "bold"),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = .25))
ggsave(file.path(this.part.dir, "plot13_new-ecf-enrichment-pfams.png"),
       p.new.ecf.enrichment.pfams, width = 8, height = 10)

####################################################################
####################################################################
###  STEP 4: ECFs link unrelated domains
####################################################################
####################################################################

## ---- SUPP: ECFs bridge multiple ECOD X-groups (p.ecod.links) ------------
## Quantitative backing for the worked example (MAIN E).
ecfs.vs.ecod.link <- ecfs.vs.ecod[prob >= MIN.PROB.FOR.ECF.DOMAIN.CONFIDENT]
ecf.ecod.linking.stats <- ecfs.vs.ecod.link[, .(
  no.unique.x = uniqueN(x.id),
  no.unique.h = uniqueN(h.id),
  no.unique.t = uniqueN(t.id)
), by = qname]
ecf.ecod.linking.stats <- merge(ecf.ecod.linking.stats,
                                ecf.ecod.domain.assignment, by = "qname", all.x = TRUE)

plot_dt.links <- ecf.ecod.linking.stats[
  class.assignment %chin% c("equivalent", "related", "subdomain", "superdomain") &
    !is.na(no.unique.x)
]
plot_dt.links[, class.assignment := factor(class.assignment, levels = class.plotting.order)]

p.ecod.links <- ggplot(plot_dt.links, aes(x = class.assignment, y = no.unique.x, fill = class.assignment)) +
  geom_violin(width = .95, trim = TRUE, scale = "width", colour = "grey25", linewidth = .35) +
  geom_point(position = position_jitter(width = .14, height = 0), alpha = .2, size = 1, colour = "grey10") +
  scale_y_continuous(breaks = sort(unique(plot_dt.links$no.unique.x))) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(5, "Set2")[c(5, 4, 3, 2)]) +
  labs(x = "ECF-ECOD class", y = "Number of ECOD\nX-groups linked") +
  theme_classic(base_size = 12) +
  guides(fill = "none") +
  theme(axis.title.x = element_text(margin = margin(t = 8), face = "bold"),
        axis.title.y = element_text(face = "bold"),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = .25))
ggsave(file.path(this.part.dir, "plot14_ecod-links.png"),
       p.ecod.links, width = 4, height = 4)

### confident hits only (you already filter upstream; keep consistent)
links <- copy(ecfs.vs.ecod.link[prob >= MIN.PROB.FOR.ECF.DOMAIN.CONFIDENT])

### --- per-ECF bridging level (nested, mutually exclusive) ---
### MIN.SECOND = minimum hits to a second group for it to "count"
MIN.SECOND <- 2L

bridge.level <- links[, {
  ## count hits per group at each level
  t.tab <- table(t.id)
  h.tab <- table(h.id)
  x.tab <- table(x.id)
  ## a level is "bridged" if >=2 groups each have >= MIN.SECOND hits
  n.t <- sum(t.tab >= MIN.SECOND)
  n.h <- sum(h.tab >= MIN.SECOND)
  n.x <- sum(x.tab >= MIN.SECOND)
  level <- if (n.x >= 2) "Bridges X-group"
  else if (n.h >= 2) "Bridges H-group"
  else if (n.t >= 2) "Bridges T-group"
  else "Within T-group"
  .(level = level, n.hits = .N)
}, by = qname]

### attach ECF-ECOD class (from your existing per-ECF assignment)
bridge.level <- merge(
  bridge.level,
  ecf.ecod.domain.assignment[, .(qname, class.assignment)],
  by = "qname", all.x = TRUE
)

### keep the meaningful classes, set orders
class.ord <- c("equivalent", "subdomain", "superdomain", "related", "distant")
level.ord <- c("Within T-group", "Bridges T-group", "Bridges H-group", "Bridges X-group")
plot.dt <- bridge.level[class.assignment %chin% class.ord]
plot.dt[, class.assignment := factor(class.assignment, levels = class.ord)]
plot.dt[, level := factor(level, levels = level.ord)]

### proportions within each class
prop.dt <- plot.dt[, .N, by = .(class.assignment, level)]
prop.dt[, prop := N / sum(N), by = class.assignment]
### bin counts for axis annotation
class.n <- plot.dt[, .(n = .N), by = class.assignment]

p.class.assignment <- ggplot(prop.dt, aes(x = class.assignment, y = prop, fill = level)) +
  geom_col(width = 0.85, colour = "white", linewidth = 0.3) +
  scale_fill_manual(
    values = c("Within T-group"  = "#D3D1C7",
               "Bridges T-group"  = "#97C459",
               "Bridges H-group"  = "#1D9E75",
               "Bridges X-group"  = "#185FA5"),
    name = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.03))) +
  labs(x = "ECF\u2013ECOD class", y = "Proportion of ECFs") +
  theme_bw(base_size = 13) +
  theme(axis.title = element_text(face = "bold"),
        legend.position = "right",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())

ggsave(file.path(this.part.dir, "plot15_class-assignment.png"),
       p.class.assignment, width = 7, height = 4)

####################################################################
####################################################################
###  STEP 5: Examples of ECFs linking diverse X
####################################################################
####################################################################

SUBDOMAIN.ECF.EXAMPLE.MIN.QCOV <- 0.9
SUBDOMAIN.ECF.EXAMPLE.MIN.PROB <- 95
SUBDOMAIN.ECF.EXAMPLE.MAX.SCOV <- 0.7
SUBDOMAIN.ECF.EXAMPLE.MIN.NO.X <- 2

##############################################################################
### STEP A -- find ECFs that bridge >=2 X-groups AND map to >=2 distinct
###           functions (the candidate pool for the worked example)
##############################################################################

ecfs.vs.ecod.link.subdomain.shortlist <- ecfs.vs.ecod.link[
  qname %chin% ecf.ecod.linking.stats[no.unique.x >= SUBDOMAIN.ECF.EXAMPLE.MIN.NO.X]$qname &
    prob >= SUBDOMAIN.ECF.EXAMPLE.MIN.PROB &
    qcov >= SUBDOMAIN.ECF.EXAMPLE.MIN.QCOV &
    scov <= SUBDOMAIN.ECF.EXAMPLE.MAX.SCOV
]

reprseq.vs.ecod.complete.domains <- reprseq.vs.ecod[scov >= SUBDOMAIN.ECF.EXAMPLE.MIN.QCOV]
reprseq.vs.ecod.best <- reprseq.vs.ecod.complete.domains[
  order(-bitscore),
  .(reprseq      = qname[1],
    rs.start     = qstart[1],
    rs.end       = qend[1],
    rs.length    = qlength[1],
    rs.bitscore  = bitscore[1],
    rs.prob      = prob[1],
    rs.qcov      = qcov[1],
    rs.scov      = scov[1]),
  by = .(sname)
]

ecf.ecod.reprseq <- merge(
  ecfs.vs.ecod.link.subdomain.shortlist,
  reprseq.vs.ecod.best,
  by = "sname",
  all.x = TRUE
)
setorder(ecf.ecod.reprseq, qname, x.id, -bitscore)

ecf.ecod.reprseq.match <- ecf.ecod.reprseq[!is.na(reprseq)]
ecf.ecod.reprseq.match[, funct := reprseq.functions.unique[.(ecf.ecod.reprseq.match$reprseq)]$funct]

all.matching.ecfs <- unique(ecf.ecod.reprseq.match$qname)
x.stats <- ecf.ecod.reprseq.match[, .(no.x = uniqueN(x.level)), by = .(qname)]
shortlist.ecfs <- unique(x.stats[no.x > 1]$qname)

ecf.ecod.reprseq.match.shortlist <- ecf.ecod.reprseq.match[qname %chin% shortlist.ecfs]
ecf.ecod.reprseq.match.shortlist[, t.name := ecod.domains[.(ecf.ecod.reprseq.match.shortlist$sname)]$t_name]

ecf.ecod.with.funct <- ecf.ecod.reprseq.match.shortlist[!is.na(funct) & funct != ""]

ecf.x.funct <- ecf.ecod.with.funct[, .(
  functions = paste(sort(unique(funct)), collapse = " / "),
  domains   = paste(sort(unique(t.name)), collapse = " / "),
  n.reprseq = uniqueN(reprseq),
  n.ecod    = uniqueN(sname)
), by = .(qname, x.id)]

ecf.pairs <- ecf.x.funct[, .(
  n.x         = uniqueN(x.id),
  n.functions = uniqueN(unlist(strsplit(functions, " / "))),
  x.groups    = paste(sort(unique(x.id)), collapse = " | "),
  funct.by.x  = paste(sprintf("X%s: %s", x.id, functions), collapse = "  ||  "),
  domain.by.x = paste(sprintf("X%s: %s", x.id, domains), collapse = "  ||  ")
), by = qname]

interesting <- ecf.pairs[n.x >= 2 & n.functions >= 2][order(-n.functions, -n.x)]

##############################################################################
### STEP B -- hand-picked examples: build candidate rows per (ECF, X-group)
##############################################################################

## Hand-picked from `interesting` above -- update as the shortlist evolves.
# ecf.examples <- c("ecf_00043", "ecf_00108", "ecf_01881", "ecf_04622", "ecf_06876")
ecf.examples <- c("ecf_00043", "ecf_00108", "ecf_01881", "ecf_04622")
`%||%` <- function(a, b) if (is.null(a)) b else a

ecf.on.reprseq <- reprseq.vs.ecfs[
  sname %chin% ecf.examples,
  .(ecf = sname, reprseq = qname,
    ecf.start = qstart, ecf.end = qend, rs.length = qlength, rs.prob = prob)
]

cand <- merge(
  ecf.ecod.reprseq.match.shortlist[
    qname %chin% ecf.examples,
    .(ecf = qname, sname, x.id, t.name, reprseq, funct)],
  ecf.on.reprseq,
  by = c("ecf", "reprseq")
)
cand <- unique(cand[!is.na(ecf.start) & !is.na(funct) & funct != ""])

dom.on.same.rs <- unique(reprseq.vs.ecod[
  , .(reprseq = qname, sname, dom.start = qstart, dom.end = qend)])

cand <- merge(cand, dom.on.same.rs, by = c("reprseq", "sname"), all.x = TRUE)

cand[, drawable := !is.na(dom.start)]
setorder(cand, ecf, x.id, -drawable, -rs.prob)


# # Option A: pick automatically
# rep.per.x <- cand[, {
#   used <- character(0)
#   x.ids <- unique(x.id)
#   ord <- x.ids[order(sapply(x.ids, function(xx) uniqueN(reprseq[x.id == xx])))]
#   
#   picked <- vector("list", length(ord))
#   for (i in seq_along(ord)) {
#     rows   <- .SD[x.id == ord[i]]
#     free   <- rows[!reprseq %chin% used]
#     chosen <- if (nrow(free) > 0) free[1] else rows[1]
#     chosen[, reused := nrow(free) == 0]
#     picked[[i]] <- chosen
#     used <- c(used, chosen$reprseq)
#   }
#   rbindlist(picked)
# }, by = ecf]
# setorder(rep.per.x, ecf, x.id)

# Option B: pick by hand
# index=3
# this.ex <- cand[ecf == ecf.examples[index]]
# View(this.ex)
example.choices <- c(
  list(c("ECOD_000010092_e1s2dA1", "ECOD_002630372_e6lhxA2")),
  list(c("ECOD_000011272_e1ttzA1", "ECOD_001108045_e3j46n1")),
  list(c("ECOD_001772495_e3g7dA1", "ECOD_000002765_e1lmb31", "ECOD_000005961_e1nr3A1")),
  list(c("ECOD_000143750_e3rgoA1", "ECOD_002130862_e6apxA3"))
)

rep.per.x.list <- lapply(1:length(example.choices), function(index){
  this.ecf <- ecf.examples[index]
  this.ecods <- example.choices[[index]]
  cand[ecf == this.ecf & sname %chin% this.ecods]  
})
rep.per.x <- rbindlist(rep.per.x.list)


# ## attach n.reprseq (how many reprseqs back this ECF/X-group/function combo,
# ## not just the single chosen representative) for the panel labels below
# rep.per.x <- merge(rep.per.x, ecf.x.funct[, .(qname, x.id, n.reprseq)],
#                    by.x = c("ecf", "x.id"), by.y = c("qname", "x.id"), all.x = TRUE)
# setorder(rep.per.x, ecf, x.id)
# 
# print(rep.per.x[, .(ecf, x.id, t.name, funct, n.reprseq, reprseq,
#                     rs.length, dom.start, dom.end, ecf.start, ecf.end,
#                     drawable, reused)])

##############################################################################
### STEP C -- publication panel: 2-4 hand-picked ECFs, 2 columns, shared
###           ECF colour legend, each ECF's X-group pair kept adjacent
##############################################################################

## Public entry point: pick 2-4 ECFs, get back one combined figure with
## `ncol` columns, one correctly-merged ECF-colour legend, light vertical
## gridlines for residue position, and no clipped domain-name titles.
##
## This is a single ggplot object (facet_wrap by column) rather than two
## separate ggplots glued together with patchwork. That matters here:
## patchwork's guide collection isn't reliable across ggplot2/patchwork
## versions when each panel only uses a subset of the fill levels (that was
## the cause of the broken/duplicated ECF legend) -- a single shared scale
## sidesteps the issue entirely instead of working around it.
draw.ecf.examples <- function(ecf.ids, ncol = 2, dom.col = "#F2C14E", ecf.cols = NULL) {
  stopifnot(length(ecf.ids) >= 1)
  if (is.null(ecf.cols)) {
    ecf.cols <- setNames(scales::hue_pal()(length(ecf.ids)), ecf.ids)
  }
  
  d <- rep.per.x[ecf %chin% ecf.ids]
  d[, ecf := factor(ecf, levels = ecf.ids)]
  
  ## assign each ECF to a column, preserving the caller's order so each
  ## ECF's X-group pair stays adjacent within its column
  split.idx <- rep(seq_len(ncol), each = ceiling(length(ecf.ids) / ncol),
                   length.out = length(ecf.ids))
  col.of.ecf <- setNames(split.idx, ecf.ids)
  d[, column := col.of.ecf[as.character(ecf)]]
  
  setorder(d, column, ecf, x.id)
  d[, y := .N - seq_len(.N) + 1, by = column]
  d[, label     := sprintf("%s\n%s", reprseq, funct)]
  d[, dom.label := sprintf("%s (X%s)", t.name, x.id)]
  
  ## per-facet y-axis breaks/labels: each column shows its own reprseq/
  ## function text rather than one global set forced onto every panel.
  y.scales <- lapply(sort(unique(d$column)), function(cc) {
    sub <- d[column == cc]
    scale_y_continuous(breaks = sub$y, labels = sub$label,
                       expand = expansion(add = 0.3))
  })
  names(y.scales) <- as.character(sort(unique(d$column)))
  
  ggplot(d) +
    geom_segment(aes(x = 1, xend = rs.length, y = y, yend = y),
                 linewidth = 3.2, colour = "grey30", lineend = "round") +
    geom_rect(data = d[!is.na(dom.start)],
              aes(xmin = dom.start, xmax = dom.end, ymin = y - 0.12, ymax = y + 0.12),
              fill = dom.col, colour = "grey35", linewidth = 0.3) +
    geom_rect(aes(xmin = ecf.start, xmax = ecf.end, ymin = y - 0.06, ymax = y + 0.06,
                  fill = ecf),
              colour = "grey15", linewidth = 0.3) +
    geom_text(data = d[!is.na(dom.start)],
              aes(x = (dom.start + dom.end) / 2, y = y + 0.40, label = dom.label),
              size = 2.5, colour = "grey20", vjust = 0, fontface = "bold") +
    facet_wrap(~ column, ncol = ncol, scales = "free") +
    ggh4x::facetted_pos_scales(y = y.scales) +
    scale_fill_manual(values = ecf.cols, name = "ECF") +
    scale_x_continuous(expand = expansion(mult = c(0.04, 0.08))) +
    coord_cartesian(clip = "off") +
    labs(x = "Residue position", y = NULL) +
    theme_classic(base_size = 12) +
    theme(panel.grid.major.x = element_line(colour = "grey85", linewidth = 0.3),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          strip.background = element_blank(),
          strip.text        = element_blank(),
          axis.text.y  = element_text(hjust = 0, colour = "grey15", lineheight = 0.9),
          axis.ticks.y = element_blank(),
          axis.line.y  = element_blank(),
          axis.title.x = element_text(face = "bold"),
          panel.spacing.x = unit(2.2, "lines"),
          plot.margin = margin(t = 18, r = 10, b = 5.5, l = 5.5))
}
# draw.ecf.examples(ecf.examples)

p.examples <- draw.ecf.examples(ecf.examples)
ggsave(file.path(this.part.dir, "plot16_examples.png"),
       p.examples, width = 10, height = 4)


# ## B4. Draw one ECF: backbone -> ECOD domain -> ECF nested inside
# draw.ecf <- function(ecf.id, dom.col = "#F2C14E", ecf.col = "#1D6F9E", title = NULL) {
#   d <- rep.per.x[ecf == ecf.id]
#   setorder(d, x.id)
#   d[, y := .N - seq_len(.N) + 1]
#   d[, label := sprintf("%s  \u2014  %s  (X%s)", reprseq, funct, x.id)]
#   
#   ggplot(d) +
#     geom_segment(aes(x = 1, xend = rs.length, y = y, yend = y),
#                  linewidth = 2.4, colour = "grey85", lineend = "round") +
#     geom_rect(data = d[!is.na(dom.start)],
#               aes(xmin = dom.start, xmax = dom.end,
#                   ymin = y - 0.26, ymax = y + 0.26),
#               fill = dom.col, colour = "grey35", linewidth = 0.3) +
#     geom_rect(aes(xmin = ecf.start, xmax = ecf.end,
#                   ymin = y - 0.16, ymax = y + 0.16),
#               fill = ecf.col, colour = "grey15", linewidth = 0.3) +
#     geom_text(data = d[!is.na(dom.start)],
#               aes(x = (dom.start + dom.end) / 2, y = y + 0.34, label = t.name),
#               size = 2.4, colour = "grey25", vjust = 0) +
#     scale_y_continuous(breaks = d$y, labels = d$label,
#                        expand = expansion(add = 0.8)) +
#     scale_x_continuous(expand = expansion(mult = c(0.01, 0.03))) +
#     labs(x = "Residue position", y = NULL, title = title %||% ecf.id) +
#     theme_minimal(base_size = 12) +
#     theme(panel.grid = element_blank(),
#           axis.text.y = element_text(hjust = 0, colour = "grey15"),
#           axis.title.x = element_text(face = "bold"),
#           plot.title = element_text(face = "bold"))
# }
# 
# index=4
# draw.ecf(ecf.examples[index], title = ecf.examples[index])
# 


################################################################################
###  ASSEMBLE MAIN FIGURE  (A -> B -> C -> D ; E added once built)
################################################################################
## Layout placeholder. Panel E (worked example) is inserted when ready.
fig.main <- (
  (p.qs.node | p.ecod.domain.cat.by.length) /
    (p.ecf.domain.coverage | p.new.ecf.enrichment)
) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16),
        plot.margin = margin(8, 8, 8, 8))

fig.main.outfile <- file.path(ecf.figures.dir, "figure2.png")
ggsave(fig.main.outfile, fig.main, width = 12, height = 11)

################################################################################
###  NOTES ON WHAT WAS CUT / RELOCATED (intentionally NOT generated here)
################################################################################
## CUT (former part2):
##   p.asym.fam, p.asym.deg, p.glm.model + GLM table  -> IRR ~ 1.1 neighbour-
##     family-diversity result. Reduce to a single supplementary sentence.
##   p.asym.str (asymmetry vs strength), p.str.len (strength vs length) -> sanity
##     checks, dropped.
## RELOCATE (former part4):
##   p.domain.class.per.fam (pA/pB/pC/pD: per-family domain-class diversity) ->
##     concerns ECF *families*, belongs with the architecture figure's supplement
##     (or cut). Not generated in this subsection.
################################################################################
