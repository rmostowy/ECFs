### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(patchwork)))
suppressWarnings(suppressPackageStartupMessages(library(ggh4x)))
suppressWarnings(suppressPackageStartupMessages(library(ggrepel)))
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
this.part.dir <- file.path(ecf.proj.dir, "part4_arch-regimes")
make.dir(this.part.dir)

#################################
### Load all data:
### 1. ECF tables
### 2. reprseq lengths
### 3. reprseq annotations
### 4. reprseq vs ECF table
### 5. reprseq arch metrics
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

# Load reprseq architecture metrics
reprseq.arch.metrics <- fread(reprseq.arch.metrics.output)

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

# append functions onto protein metrics
pm.funct <- reprseq.functions.unique[.(reprseq.arch.metrics$qname)][, .(qname = reprseq, funct, class)]
reprseq.arch.metrics <- merge(reprseq.arch.metrics, pm.funct, by = "qname")
reprseq.arch.metrics[class == "", class := "unknown"]
reprseq.arch.metrics[funct == "", funct := "unknown"]

base_theme <- theme_bw(base_size = 14) +
  theme(
    axis.title   = element_text(face = "bold"),
    axis.text    = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    plot.title   = element_text(face = "bold", hjust = 0)
  )

###########################################
### 1. Investigate Number of EFs  
### as a function of coverage
###########################################

p.div.cov <- ggplot(reprseq.arch.metrics[coverage.fraction > 0], aes(x = n.ecf.hits, y = n.fams)) +
  geom_bin2d(bins = 60) +
  scale_fill_viridis_c(trans = "log10", name = "Count") +
  scale_x_continuous(trans = "log10") +
#  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  base_theme +
  labs(
    x = "Number of ECF hits (log10 scale)",
    y = "EF diversity"
  )

ggsave(file.path(this.part.dir, "plot1_div-cov.png"),
       p.div.cov, width = 6, height = 5)
  
p.div.cov.short <- ggplot(reprseq.arch.metrics[coverage.fraction > 0 & length < 1000], aes(x = n.ecf.hits, y = n.fams)) +
  geom_bin2d(bins = 60) +
  scale_fill_viridis_c(trans = "log10", name = "Count") +
  # geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  scale_x_continuous(trans = "log10") +
  base_theme +
  labs(
    x = "Number of ECF hits (log10 scale)",
    y = "EF diversity"
  )

ggsave(file.path(this.part.dir, "plot2_div-cov-short.png"),
       p.div.cov.short, width = 6, height = 5)

###########################################
### 2. Investigate Number of EFs  
### vs. segregation
###########################################

segregation.range <- range(reprseq.arch.metrics$n.fams)
segregation.range.min <- floor(min(segregation.range))
segregation.range.max <- floor(max(segregation.range))

# Punchline 2 — Family diversity and spatial organization are distinct axes
p.fam.div.spatial.separation.bary <- ggplot(reprseq.arch.metrics, aes(x = n.fams, y = segregation.bary)) +
  geom_bin2d(bins = 60) +
  theme_bw(base_size = 14) +
  scale_fill_viridis_c(trans = "log10", name = "Count") +
  scale_x_sqrt(breaks = c(1,2,3,4,5,10,15,20,30,40)) +
  labs(
    x = "EF diversity",
    y = "EF segregation"
  ) + 
  base_theme
ggsave(file.path(this.part.dir, "plot3_div-segr.png"),
       p.fam.div.spatial.separation.bary, width = 6, height = 5)


###########################################
### 3a. Eight-example grid for Panel C,
### organised as 4 quadrants x 2 examples,
### matching Panel B's diversity/segregation axes
###########################################

MIN.HITS.FOR.EXAMPLE <- 20

low.div.range  <- c(2, 2)
high.div.range <- c(4, 15)   # deliberately stops short of the >20 tail

arch <- reprseq.arch.metrics[n.ecf.hits > MIN.HITS.FOR.EXAMPLE & !is.na(segregation.bary)]
arch[, fam.evenness := mean.eff.fam / n.fams]

# segregation cutoffs computed WITHIN each diversity band, not globally —
# this is what keeps "high diversity, high segregation" from being an empty set
seg.q.low.div  <- quantile(arch[n.fams %between% low.div.range]$segregation.bary,
                           probs = c(0.25, 0.75), na.rm = TRUE)
seg.q.high.div <- quantile(arch[n.fams %between% high.div.range]$segregation.bary,
                           probs = c(0.25, 0.75), na.rm = TRUE)

quadrant.defs <- list(
  "low_div-high_seg"  = arch[n.fams %between% low.div.range  & segregation.bary >= seg.q.low.div["75%"]],
  "low_div-low_seg"   = arch[n.fams %between% low.div.range  & segregation.bary <= seg.q.low.div["25%"]],
  "high_div-high_seg" = arch[n.fams %between% high.div.range & segregation.bary >= seg.q.high.div["75%"]],
  "high_div-low_seg"  = arch[n.fams %between% high.div.range & segregation.bary <= seg.q.high.div["25%"]]
)

# quadrants where "low segregation" should mean genuine co-occurring families,
# not one dominant family plus a stray hit from another — rank by evenness
# and keep the best-supported cases instead of sampling broadly
EVENNESS.RANKED.QUADRANTS <- c("low_div-low_seg", "high_div-low_seg")

MAX.CANDIDATES.PER.QUADRANT <- 100
set.seed(1234)

for (qname.label in names(quadrant.defs)) {
  shortlist <- quadrant.defs[[qname.label]]
  
  if (qname.label %in% EVENNESS.RANKED.QUADRANTS) {
    shortlist <- shortlist[order(-fam.evenness)]
    shortlist <- head(shortlist, MAX.CANDIDATES.PER.QUADRANT)
  } else {
    shortlist <- shortlist[order(n.fams)]
    if (nrow(shortlist) > MAX.CANDIDATES.PER.QUADRANT) {
      shortlist <- shortlist[order(qname)][sample(.N, MAX.CANDIDATES.PER.QUADRANT)]
      shortlist <- shortlist[order(n.fams)]
    }
  }
  
  this.ex.dir <- file.path(this.part.dir, "vis-ecf-distr", paste0("panelC_", qname.label))
  make.dir(this.ex.dir)
  
  for (i in seq_len(nrow(shortlist))) {
    protein_name <- shortlist$qname[i]
    protein_filename <- file.path(this.ex.dir, sprintf("%s.png", protein_name))
    if(!file.exists(protein_filename)){
      protein.ecf.hits_all <- reprseq.vs.ecfs[.(protein_name)]
      title <- sprintf(
        "%s (%s)\nn.fams = %d, seg = %.2f, evenness = %.2f",
        protein_name, shortlist$funct[i], shortlist$n.fams[i], shortlist$segregation.bary[i], shortlist$fam.evenness[i]
      )
      this.plot <- plot.ecf.distr.over.query(protein_name, protein.ecf.hits_all, title = title)
      ggsave(protein_filename, this.plot, width = 6, height = 5)  
    }
  }
}

###########################################
### 3b. Assemble Panel C from hand-picked
### examples — edit the qnames below once
### you've inspected the candidate folders
###########################################

# EDIT THIS — replace each placeholder with a qname from the matching
# panelC_<quadrant> folder. Layout mirrors Panel B: top row = high segregation,
# bottom row = low segregation; left two columns = low diversity, right two = high diversity.
panel.c.picks <- c(
  # top row — high segregation
  "reprseq005982",   # low diversity, high segregation
  "reprseq042119",   # low diversity, high segregation
  "reprseq007270",  # high diversity, high segregation
  "reprseq038339",  # high diversity, high segregation
  # bottom row — low segregation
  "reprseq049478",    # low diversity, low segregation
  "reprseq074693",    # low diversity, low segregation
  "reprseq034239",   # high diversity, low segregation
  "reprseq043569"    # high diversity, low segregation
)

quadrant.labels <- c(
  "Low diversity, high segregation", "Low diversity, high segregation",
  "High diversity, high segregation", "High diversity, high segregation",
  "Low diversity, low segregation", "Low diversity, low segregation",
  "High diversity, low segregation", "High diversity, low segregation"
)

final.picks <- data.table(qname = panel.c.picks, quadrant = quadrant.labels)
final.picks <- merge(final.picks, reprseq.arch.metrics[, .(qname, funct, n.fams, segregation.bary)],
                     by = "qname", sort = FALSE)


###############

make.strip <- function(label, fill = "grey88") {
  ggplot() +
    annotate("text", x = 0, y = 0, label = label, fontface = "bold", size = 3.8, lineheight = 0.9) +
    xlim(-1, 1) + ylim(-1, 1) +
    theme_void() +
    theme(plot.background = element_rect(fill = fill, colour = "black", linewidth = 0.4))
}

strips <- list(
  "low_div-high_seg"  = make.strip("Low diversity, high segregation"),
  "high_div-high_seg" = make.strip("High diversity, high segregation"),
  "low_div-low_seg"   = make.strip("Low diversity, low segregation"),
  "high_div-low_seg"  = make.strip("High diversity, low segregation")
)

# shortened per-protein titles — the quadrant phrase now lives in the strip, not repeated 8 times
panel.c.plots <- lapply(seq_len(nrow(final.picks)), function(i) {
  row <- final.picks[i]
  protein.ecf.hits_all <- reprseq.vs.ecfs[.(row$qname)]
  title <- sprintf("%s\nn.fams = %d, seg = %.2f", row$funct, row$n.fams, row$segregation.bary)
  plot.ecf.distr.over.query(row$qname, protein.ecf.hits_all, title = title)
})

# final.picks order assumed: 1-2 low_div-high_seg, 3-4 high_div-high_seg, 5-6 low_div-low_seg, 7-8 high_div-low_seg
row.high.seg <- (
  ((strips[["low_div-high_seg"]] | strips[["high_div-high_seg"]]) + plot_layout(widths = c(2, 2))) /
    (panel.c.plots[[1]] | panel.c.plots[[2]] | panel.c.plots[[3]] | panel.c.plots[[4]])
) + plot_layout(heights = c(0.15, 1))

row.low.seg <- (
  ((strips[["low_div-low_seg"]] | strips[["high_div-low_seg"]]) + plot_layout(widths = c(2, 2))) /
    (panel.c.plots[[5]] | panel.c.plots[[6]] | panel.c.plots[[7]] | panel.c.plots[[8]])
) + plot_layout(heights = c(0.15, 1))

design <- "
AAAABBBB
CCDDEEFF
GGGGHHHH
IIJJKKLL
"

p.examples <- (
  strips[["low_div-high_seg"]] + strips[["high_div-high_seg"]] +
    panel.c.plots[[1]] + panel.c.plots[[2]] + panel.c.plots[[3]] + panel.c.plots[[4]] +
    strips[["low_div-low_seg"]] + strips[["high_div-low_seg"]] +
    panel.c.plots[[5]] + panel.c.plots[[6]] + panel.c.plots[[7]] + panel.c.plots[[8]]
) + plot_layout(design = design, heights = c(0.2, 3, 0.2, 3))

# p.examples

ggsave(file.path(this.part.dir, "plot4_examples.pdf"), p.examples, width = 14, height = 8)

###########################################
### 4. Investigate AHO of per-family ECFs
### depending on the function
###########################################

# Compute per-family AHO for each reprseq
reprseq.vs.ecfs.append <- merge(reprseq.vs.ecfs, reprseq.vs.ecfs[, .N, by = .(qname, family)], by = c("qname", "family"))
reprseq.vs.ecfs.filt <- reprseq.vs.ecfs.append[N >= MIN.NO.ECF.PER.FAMILY.FOR.AHO]

reprseq.ecf.overlap.from <- reprseq.vs.ecfs.filt %>% select(qname, sname.ecf1 = sname, qstart.ecf1 = qstart, qend.ecf1 = qend, family)
reprseq.ecf.overlap.to <- reprseq.vs.ecfs.filt %>% select(qname, sname.ecf2 = sname, qstart.ecf2 = qstart, qend.ecf2 = qend, family)
reprseq.ecf.overlaps <- reprseq.ecf.overlap.from[reprseq.ecf.overlap.to, on = .(qname, family), allow.cartesian = TRUE, nomatch = 0L]
reprseq.ecf.overlaps <- reprseq.ecf.overlaps[sname.ecf1 < sname.ecf2]
reprseq.ecf.overlaps[, aho := calculate_overlap_jaccard(qstart.ecf1, qend.ecf1, qstart.ecf2, qend.ecf2)]
reprseq.ef.aho.stats <- reprseq.ecf.overlaps[, .(mean.aho = mean(aho)), by = .(qname, family)]

# produce AHO stats for the NF most frequent functions per class
NF <- 6

fun2prot <- unique(reprseq.functions[, .(class, funct, qname)])

# count proteins per function and pick top NF per class
top_fun <- fun2prot[, .(n_proteins = uniqueN(qname)), by = .(class, funct)]
setorder(top_fun, class, -n_proteins)
top_fun <- top_fun[, head(.SD, NF), by = class]
setkey(top_fun, class, funct)

# map (class, funct) -> proteins (qname), then attach AHO rows (qname,family)
top_fun2prot <- fun2prot[top_fun, on = .(class, funct), nomatch = 0L]
setkey(top_fun2prot, qname)

# join to AHO (this matches your example: take all AHO rows for proteins in the function)
top_aho <- reprseq.ef.aho.stats[top_fun2prot, on = "qname", allow.cartesian = TRUE]
top_aho <- top_aho[!is.na(family)]

# summarize AHO distribution per (class, funct)
top.aho.stats <- top_aho[, .(
  n_proteins = uniqueN(qname),
  n_pf = .N,  # protein-family observations
  mean_aho = mean(mean.aho, na.rm = TRUE),
  median_aho = median(mean.aho, na.rm = TRUE),
  q25 = quantile(mean.aho, 0.25, na.rm = TRUE, type = 2),
  q75 = quantile(mean.aho, 0.75, na.rm = TRUE, type = 2)
), by = .(class, funct)]
setkey(top.aho.stats, class, funct)

# keep protein counts from top_fun (optional sanity)
top.aho.stats <- top_fun[top.aho.stats, on = .(class, funct)]

# Order functions within each class by mean_aho
top.aho.stats[, funct_ord := reorder(funct, mean_aho), by = class]

strip <- ggh4x::strip_themed(
  background_y = ggh4x::elem_list_rect(
    fill = grDevices::adjustcolor(
      phrogs.col[match(phrogs.class.order, phrogs.class)],
      alpha.f = 0.25
    ),
    colour = "black",
    linewidth = 0.4
  )
)

p.aho.top.funct <- ggplot(top.aho.stats, aes(x = mean_aho, y = funct_ord)) +
  geom_errorbarh(aes(xmin = q25, xmax = q75), height = 0.2, linewidth = 0.35) +
  geom_point(size = 2) +
  ggh4x::facet_grid2(
    rows = vars(class),
    scales = "free_y",
    space  = "free_y",
    labeller = labeller(class = phrog.class.labels),
    strip = strip
  ) +
  # facet_grid(rows = vars(class), scales = "free_y", space = "free_y") +
  theme_bw(base_size = 12) +
  scale_x_continuous(limits = c(0,1)) +
  labs(x = "Mean within-ECF-family AHO",
       y = NULL
  ) +
  theme(panel.grid.minor = element_blank(),
        strip.placement = "outside",
        strip.text.y = element_text(angle = 0, hjust = 0, face = "bold", size = 9)
  )

ggsave(file.path(this.part.dir, "plot5_aho-funct.pdf"), p.aho.top.funct, width = 7, height = 8)

p.aho.dist <- ggplot(reprseq.ef.aho.stats, aes(x = mean.aho)) +
  geom_histogram(bins = 60, fill = "#2c7fb8", color = "white") +
  geom_vline(xintercept = median(reprseq.ef.aho.stats$mean.aho), linetype = "dashed", col = "red", linewidth = 2) +
  theme_bw(base_size = 14) +
  labs(
    x = "Within-ECF-family AHO",
    y = "Frequency"
  )
ggsave(file.path(this.part.dir, "plot6_aho-distr.pdf"), p.aho.dist, width = 6, height = 5)


###########################################
### 5. Correlation structure of 
### architecture metrics
###########################################

corr.vars <- c(
  "length",
  "n.ecf.hits",
  "coverage.fraction",
  "n.fams",
  "mean.eff.fam",
  "n.dom.segments",
  "segregation.bary",
  "segregation.overlap"
)

corr.dt <- copy(reprseq.arch.metrics[, ..corr.vars])
corr.mat <- suppressWarnings(cor(corr.dt, use = "pairwise.complete.obs", method = "spearman"))

corr.long <- as.data.table(as.table(corr.mat))
setnames(corr.long, c("V1", "V2", "N"), c("x", "y", "rho"))
corr.long$x <- factor(corr.long$x, levels = rev(corr.vars))
corr.long$y <- factor(corr.long$y, levels = rev(corr.vars))

p.sup.arch.corr <- ggplot(corr.long, aes(x = x, y = y, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 3) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, limits = c(-1, 1)) +
  labs(x = NULL, y = NULL, fill = "Spearman\nrho") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_blank()
  )

ggsave(file.path(this.part.dir, "plot7_corr-matrix.pdf"), p.sup.arch.corr, width = 6, height = 5)


###########################################
### 6. Pairwise density grid —
### lower triangle = joint density (bin2d),
### upper triangle = Spearman rho,
### diagonal = variable name
###########################################

make.pairs.grid <- function(dt, vars, corr.method = "spearman") {
  n <- length(vars)
  corr.mat <- suppressWarnings(cor(dt[, ..vars], use = "pairwise.complete.obs", method = corr.method))
  rho.pal <- scales::col_numeric(c("#2166ac", "white", "#b2182b"), domain = c(-1, 1))
  
  panels <- list()
  for (row.i in seq_len(n)) {
    for (col.i in seq_len(n)) {
      v.row <- vars[row.i]; v.col <- vars[col.i]
      
      if (row.i == col.i) {
        # diagonal: just the variable name — standard pairs-plot convention
        p <- ggplot() +
          annotate("text", x = 0, y = 0, label = v.row, size = 3.2, fontface = "bold") +
          xlim(-1, 1) + ylim(-1, 1) + theme_void() +
          theme(plot.background = element_rect(fill = "grey92", colour = "white"))
        
      } else if (row.i > col.i) {
        # lower triangle: joint density, not raw scatter — avoids overplotting at this N
        p <- ggplot(dt, aes(x = .data[[v.col]], y = .data[[v.row]])) +
          geom_bin2d(bins = 40) +
          scale_fill_viridis_c(trans = "log10", guide = "none") +
          theme_void() +
          theme(plot.background = element_rect(fill = "white", colour = "grey80", linewidth = 0.2))
        
      } else {
        # upper triangle: the correlation coefficient, background-tinted to match the heatmap palette
        rho.val <- corr.mat[v.row, v.col]
        p <- ggplot() +
          annotate("text", x = 0, y = 0, label = sprintf("%.2f", rho.val), size = 3.5,
                   color = if (abs(rho.val) > 0.6) "white" else "black") +
          xlim(-1, 1) + ylim(-1, 1) + theme_void() +
          theme(plot.background = element_rect(fill = rho.pal(rho.val), colour = "white"))
      }
      panels[[length(panels) + 1]] <- p
    }
  }
  wrap_plots(panels, nrow = n, ncol = n)
}

p.sup.arch.pairs <- make.pairs.grid(reprseq.arch.metrics, corr.vars)
ggsave(file.path(this.part.dir, "plot8_corr-matrix-grid.pdf"), p.sup.arch.pairs, width = 10, height = 10)

###########################################
### 7. PCA of the five architecture metrics
###########################################

### 7a. Prepare data ------------------------------------------------------
arch.pca.vars <- c("n.fams", "mean.eff.fam", "n.dom.segments", "segregation.bary", "segregation.overlap")
pca.dt <- na.omit(copy(reprseq.arch.metrics[, ..arch.pca.vars]))

# rank-transform to stay Spearman-consistent with the rest of the paper —
# prcomp() has no built-in rank-based mode, so this is done manually before fitting
pca.dt.ranked <- pca.dt[, lapply(.SD, rank)]

### 7b. Fit PCA -------------------------------------------------------------
arch.pca <- prcomp(pca.dt.ranked, scale. = TRUE)
var.explained <- arch.pca$sdev^2 / sum(arch.pca$sdev^2)
cum.var.explained <- cumsum(var.explained)
loadings.pc1.pc2 <- arch.pca$rotation[, 1:2]

### 7c. Check: do n.fams / segregation.bary track PC1 / PC2? ----------------
pc.scores <- as.data.table(arch.pca$x[, 1:2])
check.dt <- cbind(pca.dt, pc.scores)

pca.tracking.cor <- list(
  n.fams_vs_PC1            = cor(check.dt$n.fams, check.dt$PC1, method = "spearman"),
  segregation.bary_vs_PC2  = cor(check.dt$segregation.bary, check.dt$PC2, method = "spearman")
)
# print(pca.tracking.cor)

### 7d. Scoping check (not plotted) ------------------------------------------
# Confirms why length / n.ecf.hits / coverage.fraction are excluded from the
# main PCA: including them merges depth/size with complexity into PC1 rather
# than leaving two clean architecture axes. Reported as a sentence in text,
# not as a figure — see partial-correlation analysis for the confound story.
arch.pca.vars.full <- c(arch.pca.vars, "n.ecf.hits", "length", "coverage.fraction")
pca.dt.full <- na.omit(reprseq.arch.metrics[, ..arch.pca.vars.full])
pca.dt.full.ranked <- pca.dt.full[, lapply(.SD, rank)]
arch.pca.full <- prcomp(pca.dt.full.ranked, scale. = TRUE)

var.explained.full <- (arch.pca.full$sdev^2 / sum(arch.pca.full$sdev^2))[1:4]
# print(round(var.explained.full, 3))


### 7e. Plots

### PLOT 1: Scree plot ------------------------------------------------------------
scree.dt <- data.table(
  PC = factor(paste0("PC", seq_along(var.explained)), levels = paste0("PC", seq_along(var.explained))),
  var.explained = var.explained,
  cum.var = cum.var.explained
)

pc1.pc2.cum.label <- scales::percent(cum.var.explained[2], accuracy = 1)

p.pca.scree <- ggplot(scree.dt, aes(x = PC)) +
  geom_col(aes(y = var.explained), fill = "#2a78d6", width = 0.6) +
  geom_line(aes(y = cum.var, group = 1), color = "#1baf7a", linewidth = 0.8) +
  geom_point(aes(y = cum.var), color = "#1baf7a", size = 2.5) +
  geom_text(aes(y = var.explained, label = scales::percent(var.explained, accuracy = 1)),
            vjust = -0.6, size = 3, color = "#2a78d6") +
  geom_text(aes(y = cum.var, label = scales::percent(cum.var, accuracy = 1)),
            vjust = -1.2, size = 3, color = "#1baf7a") +
  scale_y_continuous(labels = scales::percent, name = "Variance explained", limits = c(0, 1.08)) +
  labs(
    x = NULL,
    title = "A. Two components capture most architectural variation",
    subtitle = sprintf(
      "PCA on five architecture metrics (rank-transformed). PC1 + PC2 together explain %s of total variance\n(blue bars: per-component; green line: cumulative).",
      pc1.pc2.cum.label
    )
  ) +
  theme_classic(base_size = 12) +
  theme(plot.subtitle = element_text(size = 9, color = "grey30"))

### PLOT 2: Loadings ----------------------------------------------------------------
loadings.dt <- as.data.table(loadings.pc1.pc2, keep.rownames = "metric")
loadings.long <- melt(loadings.dt, id.vars = "metric", variable.name = "PC", value.name = "loading")
loadings.long[, metric := factor(metric, levels = rev(arch.pca.vars))]

# dynamic facet labels carrying each component's variance share, so the reader
# doesn't have to cross-reference panel A to know how much each axis matters
pc.facet.labels <- c(
  PC1 = sprintf("PC1 (%s of variance)", scales::percent(var.explained[1], accuracy = 1)),
  PC2 = sprintf("PC2 (%s of variance)", scales::percent(var.explained[2], accuracy = 1))
)

p.pca.loadings <- ggplot(loadings.long, aes(x = metric, y = loading, fill = loading)) +
  geom_col() +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  facet_wrap(~PC, ncol = 1, labeller = labeller(PC = pc.facet.labels)) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, guide = "none") +
  coord_flip() +
  labs(
    x = NULL, y = "Loading (contribution to component, signed)",
    title = "B. PC1 reflects overall diversity; PC2 contrasts segregation against overlap-driven diversity",
    subtitle = "Red = positive contribution, blue = negative. On PC2, segregation.bary and mean.eff.fam load in\nopposite directions, consistent with overlapping families inflating mean.eff.fam while reducing segregation."
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.subtitle = element_text(size = 9, color = "grey30"),
    strip.text = element_text(face = "bold")
  )

### PLOT 3: Tracking ------------------------------------------------------------------
rho.fam <- pca.tracking.cor$n.fams_vs_PC1
rho.seg <- pca.tracking.cor$segregation.bary_vs_PC2

p.pca.tracking.fam <- ggplot(check.dt[n.fams <= 15], aes(x = factor(n.fams), y = PC1)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2, fill = "#2a78d6", alpha = 0.3) +
  labs(
    x = "n.fams", y = "PC1 score",
    subtitle = sprintf("n.fams vs. PC1: Spearman \u03c1 = %.2f", rho.fam)
  ) +
  theme_classic(base_size = 12) +
  theme(plot.subtitle = element_text(size = 9.5, face = "italic"))

p.pca.tracking.seg <- ggplot(check.dt, aes(x = segregation.bary, y = PC2)) +
  geom_bin2d(bins = 60) +
  scale_fill_viridis_c(trans = "log10", name = "Count") +
  labs(
    x = "segregation.bary", y = "PC2 score",
    subtitle = sprintf("segregation.bary vs. PC2: Spearman \u03c1 = %.2f", rho.seg)
  ) +
  theme_classic(base_size = 12) +
  theme(plot.subtitle = element_text(size = 9.5, face = "italic"))

p.pca.tracking <- (p.pca.tracking.fam | p.pca.tracking.seg) +
  plot_annotation(
    title = "C. n.fams and segregation.bary closely track PC1 and PC2",
    subtitle = "The two metrics used throughout the main text are strong proxies for the formal principal components,\njustifying their use as representative measures of diversity and segregation.",
    theme = theme(
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 9, color = "grey30")
    )
  )

### 7d. Save ------------------------------------------------------------------
ggsave(file.path(this.part.dir, "plot9_pca-scree.pdf"), p.pca.scree, width = 6.5, height = 5.5)
ggsave(file.path(this.part.dir, "plot10_pca-loadings.pdf"), p.pca.loadings, width = 6.5, height = 6)
ggsave(file.path(this.part.dir, "plot11_pca-tracking.pdf"), p.pca.tracking, width = 10, height = 5.5)


###########################################
### X. Assemble Panel C from hand-picked
### examples — edit the qnames below once
### you've inspected the candidate folders
###########################################

TOP.N.FUNCTS <- 20
MIN.N.FOR.FUNCT.SUMMARY <- 30

# add IQR alongside mean/median — same NA-correct split as before:
# diversity computed over every protein with this function, segregation only over n.fams >= 2
diversity.by.funct <- reprseq.arch.metrics[
  !is.na(n.fams) & class != "unknown" & funct != "unknown",
  .(n_proteins_total = .N,
    mean.fam = mean(n.fams, na.rm = TRUE),
    med.fam  = as.double(median(n.fams, na.rm = TRUE)),
    q25.fam  = quantile(n.fams, 0.25, na.rm = TRUE, type = 2),
    q75.fam  = quantile(n.fams, 0.75, na.rm = TRUE, type = 2)),
  by = .(class, funct)
]

segregation.by.funct <- reprseq.arch.metrics[
  !is.na(segregation.bary) & class != "unknown" & funct != "unknown",
  .(n_proteins_multi = .N,
    mean.seg = mean(segregation.bary, na.rm = TRUE),
    med.seg  = median(segregation.bary, na.rm = TRUE),
    q25.seg  = quantile(segregation.bary, 0.25, na.rm = TRUE, type = 2),
    q75.seg  = quantile(segregation.bary, 0.75, na.rm = TRUE, type = 2)),
  by = .(class, funct)
]

function.2d.summary <- merge(diversity.by.funct, segregation.by.funct, by = c("class", "funct"))
function.2d.summary <- function.2d.summary[
  n_proteins_total >= MIN.N.FOR.FUNCT.SUMMARY & n_proteins_multi >= MIN.N.FOR.FUNCT.SUMMARY
]

# reference lines computed from the full qualifying set, not just the subset you end up showing —
# keeps the reference stable if you later change TOP.N.FUNCTS
ref.fam <- median(function.2d.summary$mean.fam)
ref.seg <- median(function.2d.summary$mean.seg)

# "most numerous" = most proteins behind the function, full stop
func.subset <- function.2d.summary[order(-n_proteins_total)][1:min(TOP.N.FUNCTS, .N)]
func.subset[, class := factor(class, levels = phrogs.class.order)]

class.colors <- setNames(phrogs.col, phrogs.class)

p.func.labeled <- ggplot(func.subset, aes(x = mean.fam, y = mean.seg)) +
  geom_vline(xintercept = ref.fam, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_hline(yintercept = ref.seg, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_point(aes(size = n_proteins_total, color = class), alpha = 0.85) +
  geom_text_repel(aes(label = funct), size = 3.2, color = "black",
                  max.overlaps = 30, box.padding = 0.4, segment.size = 0.3) +
  scale_color_manual(values = class.colors, labels = phrog.class.labels, name = "PHROG class") +
  scale_size_continuous(range = c(2, 8), name = "Proteins (n)") +
  labs(x = "Mean EF diversity", y = "Mean EF segregation") +
  base_theme

# p.func.labeled

p.func.labeled.errbars <- p.func.labeled +
  geom_errorbarh(aes(xmin = q25.fam, xmax = q75.fam, color = class), height = 0, alpha = 0.7, linewidth = 0.4) +
  geom_errorbar(aes(ymin = q25.seg, ymax = q75.seg, color = class), width = 0, alpha = 0.7, linewidth = 0.4)

ggsave(file.path(this.part.dir, "plot12_functions-on-2d.pdf"), p.func.labeled.errbars, width = 9, height = 6)


#################################
### Generate Figure 4
#################################

# p.div.cov                          # panel A
# p.fam.div.spatial.separation.bary  # panel B
# p.examples                         # panel C
# p.func.labeled.errbars             # panel D

design <- "
AB
CC
DD
"

fig4 <-
  ( p.div.cov + p.fam.div.spatial.separation.bary +
      wrap_elements(full = p.examples) +
      p.func.labeled.errbars
  ) +
  plot_layout(design = design, heights = c(1, 1.8, 1.3)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 24),
    plot.tag.position = c(0.02, 0.98)
  )

fig4.outfile <- file.path(ecf.figures.dir, "figure4.png")
ggsave(fig4.outfile, fig4, width = 12, height = 17)

#################################
### Save supplementary figures
#################################

make.dir(ecf.suppl.figures.dir)

# p.div.cov.short (diversity vs. coverage, restricted to proteins <1000 aa —
# confirms the long-tail pattern in panel A isn't driven by length alone)
suppl.p.div.cov.short.outfile <- file.path(ecf.suppl.figures.dir, "suppl-4p1_div-cov-short.png")
ggsave(suppl.p.div.cov.short.outfile, p.div.cov.short, width = 6, height = 5)

# AHO: global distribution (A), then per-function breakdown (B) —
# reordered from script sequence to match the narrative (overview before detail)
suppl.p.aho <-
  (p.aho.dist / p.aho.top.funct) +
  plot_layout(heights = c(1, 2.5)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))
suppl.p.aho.outfile <- file.path(ecf.suppl.figures.dir, "suppl-4p2_aho.png")
ggsave(suppl.p.aho.outfile, suppl.p.aho, width = 8, height = 11)

# Correlation structure: pairwise density grid only (matrix dropped per review —
# the grid's diagonal cells already carry variable names, no outer tagging needed)
suppl.p.corr.outfile <- file.path(ecf.suppl.figures.dir, "suppl-4p3_corr-structure.png")
ggsave(suppl.p.corr.outfile, p.sup.arch.pairs, width = 9, height = 9)

# PCA trio: scree, loadings, tracking.
# Each plot's title is already baked into its own text ("A. ...", "B. ...", "C. ...").
# p.pca.tracking is wrapped so its internal title survives being stacked here —
# without wrap_elements(), a title set via plot_annotation() on a nested patch
# is silently dropped when merged into a larger composite.
suppl.p.pca <- p.pca.scree / p.pca.loadings / wrap_elements(full = p.pca.tracking)
suppl.p.pca.outfile <- file.path(ecf.suppl.figures.dir, "suppl-4p4_pca.png")
ggsave(suppl.p.pca.outfile, suppl.p.pca, width = 8, height = 16)

