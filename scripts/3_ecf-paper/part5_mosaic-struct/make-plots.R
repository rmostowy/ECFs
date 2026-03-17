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


###########################################
### 1. Investigate AHO of per-family ECFs
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

p.aho.dist <- ggplot(reprseq.ef.aho.stats, aes(x = mean.aho)) +
  geom_histogram(bins = 60, fill = "#2c7fb8", color = "white") +
  geom_vline(xintercept = median(reprseq.ef.aho.stats$mean.aho), linetype = "dashed", col = "red", size = 2) +
  theme_bw(base_size = 14) +
  labs(
    x = "Within-ECF-family AHO",
    y = "Frequency"
  )

low.prob.thr <- 0.05
high.prob.thr <- 0.95
low.thr <- quantile(reprseq.ef.aho.stats$mean.aho, probs = c(low.prob.thr, high.prob.thr))[1]
high.thr <- quantile(reprseq.ef.aho.stats$mean.aho, probs = c(low.prob.thr, high.prob.thr))[2]

low.thr.reprseqs <- unique(reprseq.ef.aho.stats[mean.aho < low.thr]$qname)
high.thr.reprseqs <- unique(reprseq.ef.aho.stats[mean.aho > high.thr]$qname)

my.table(fun2prot[qname %chin% low.thr.reprseqs]$funct)[1:20]
my.table(fun2prot[qname %chin% high.thr.reprseqs]$funct)[1:20]


###########################################
### 2. Analyse ECF Family diversity per 
### ECF coverage
###########################################

# define thresholds
thresholds <- c(1,2,3,5,10,20,30,50,75,100)
curve_dt <- rbindlist(lapply(thresholds, function(k){
  sub <- reprseq.arch.metrics[n.ecf.hits > k & coverage.fraction > 0.]
  sub.table <- table(sub$n.fams)
  data.table(
    threshold = k,
    N = nrow(sub),
    mean.no.families = mean(sub$n.fams),
    dominant.category = as.integer(names(which.max(sub.table)))
  )
}))

p.fam.div.vs.ecf.depth <- ggplot(curve_dt, aes(x = threshold)) +
  # mean number of families
  geom_line(aes(y = mean.no.families), size = 1.2, color = "black") +
  geom_point(aes(y = dominant.category), size = 3, color = "red", shape = "triangle") +
  geom_point(aes(y = mean.no.families, size = N), color = "black") +
  scale_x_continuous(trans = "log10") +
  scale_y_sqrt(breaks = function(x) pretty(x, n = 15)) +
  labs(
    x = "Minimum number of ECF hits (log scale)",
    y = "Mean number of families (black)\nDominant number of families (red)"
  ) +
  theme_bw(base_size = 14) + 
  theme(
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.margin = margin(2,2,2,2)
  )

###########################################
### 3. Analyse ECF Family diversity vs
### family segregation
###########################################


segregation.range <- range(reprseq.arch.metrics$n.fams)
segregation.range.min <- floor(min(segregation.range))
segregation.range.max <- floor(max(segregation.range))

# Punchline 2 — Family diversity and spatial organization are distinct axes
p.fam.div.spatial.separation.bary <- ggplot(reprseq.arch.metrics, aes(x = n.fams, y = segregation.bary)) +
  geom_bin2d(bins = 60) +
  theme_bw(base_size = 14) +
  scale_fill_viridis_c(trans = "sqrt", name = "Count") +
  scale_x_sqrt(breaks = c(1,2,3,4,5,10,15,20,30,40)) +
  labs(
    x = "Number of ECF Families in rHMM",
    y = "Segregation index (Bary)\nalong protein sequence"
  ) + 
  theme(
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

p.fam.div.spatial.separation.overlap <- ggplot(reprseq.arch.metrics, aes(x = n.fams, y = segregation.overlap)) +
  geom_bin2d(bins = 60) +
  theme_bw(base_size = 14) +
  scale_fill_viridis_c(trans = "sqrt", name = "Count") +
  scale_x_sqrt(breaks = c(1,2,3,4,5,10,15,20,30,40)) +
  labs(
    x = "Number of ECF Families in rHMM",
    y = "ECF family segregation index\nalong protein sequence"
  ) + 
  theme(
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

###########################################
### 4. Visualise archtectural regimes
### across functional classes
###########################################



# Punchline 3 - There are “architectural regimes” across functional classes
no.functs.per.class <- 5
# dt <- copy(reprseq.arch.metrics)
function.vs.arch <- rbindlist(list(
  reprseq.arch.metrics[funct != "unknown" & class != "unknown",
     .N, by = .(class, funct)][
       order(class, -N)
     ][, head(.SD, no.functs.per.class), by = class][
       reprseq.arch.metrics, on = .(class, funct), nomatch = 0L
     ][, panel := class],
  
  reprseq.arch.metrics[funct == "unknown" | class == "unknown"][, `:=`(panel = "unknown", funct = "unknown")]
  # ,
  
  # dt[, `:=`(panel = "ALL", funct = "ALL")]
), use.names = TRUE, fill = TRUE)

function.vs.arch <- function.vs.arch[!is.na(segregation.bary)]
function.vs.arch[, xlab := sprintf("%s (N=%d)", funct, .N), by = .(panel, funct)]
function.vs.arch[, xlab := factor(xlab, levels = unique(xlab)), by = panel]

phrogs.class2 <- c(phrogs.class, "unknown")
phrog.class.labels2 <- c(phrog.class.labels, "unknown")
names(phrog.class.labels2)[phrog.class.labels2 == "unknown"] <- "unknown"

function.vs.arch[, panel := factor(class, levels = phrogs.class2)]

# strip <- ggh4x::strip_themed(
#   background_y = ggh4x::elem_list_rect(
#     fill = grDevices::adjustcolor(
#       phrogs.col[match(levels(function.vs.arch$panel), phrogs.class2)],
#       alpha.f = 0.25
#     ),
#     colour = "black",
#     linewidth = 0.4
#   )
# )
# 
# p.segreg.across.fam <- ggplot(
#   function.vs.arch[!is.na(segregation.bary) & panel != "ALL"],
#   aes(x = segregation.bary, y = reorder(funct, segregation.bary, FUN = median, na.rm = TRUE))
# ) +
#   geom_violin(scale = "width", trim = TRUE) +
#   geom_boxplot(width = 0.18, outlier.shape = NA) +
#   facet_grid2(
#     rows = vars(panel),
#     scales = "free_y",
#     space  = "free_y",
#     labeller = labeller(panel = phrog.class.labels),
#     strip = strip
#   ) +
#   geom_vline(xintercept = mean(reprseq.arch.metrics$segregation.bary, na.rm = T), linetype = "dashed", color = "red", linewidth = 0.6) +
#   labs(
#     x = "Segregation index",
#     y = NULL
#   ) +
#   theme_classic(base_size = 12) +
#   theme(
#     strip.placement = "outside",
#     strip.text.y = element_text(angle = 0, hjust = 0, face = "bold", size = 9),
#     axis.text.y = element_text(size = 8),
#     axis.text.x = element_text(size = 9),
#     panel.spacing = unit(0.35, "lines")
#   )

strip <- ggh4x::strip_themed(
  background_x = ggh4x::elem_list_rect(
    fill = grDevices::adjustcolor(
      phrogs.col[match(levels(function.vs.arch$panel), phrogs.class2)],
      alpha.f = 0.25
    ),
    colour = "black",
    linewidth = 0.4
  )
)

phrogs.class2 <- c(phrogs.class, "unknown")
phrog.class.labels2 <- c(phrog.class.labels, "unknown")
names(phrog.class.labels2)[phrog.class.labels2 == "unknown"] <- "unknown"

p.segreg.across.fam <- ggplot(
  function.vs.arch[!is.na(segregation.bary) & panel != "ALL"],
  aes(y = segregation.bary, x = reorder(funct, segregation.bary, FUN = median, na.rm = TRUE))
) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.18, outlier.shape = NA) +
  facet_grid2(
    cols = vars(panel),
    scales = "free_x",
    space  = "free_x",
    # labeller = labeller(panel = phrog.class.labels),
    labeller = labeller(panel = as_labeller(phrog.class.labels2, label_wrap_gen(width = 20))),
    strip = strip
  ) +
  geom_hline(yintercept = mean(reprseq.arch.metrics$segregation.bary, na.rm = T), linetype = "dashed", color = "red", linewidth = 0.6) +
  labs(
    y = "Segregation index",
    x = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
  strip.placement = "outside",
  strip.text.x = element_text(angle = 0, hjust = 0, face = "bold", size = 9),
  axis.text.y = element_text(size = 8),
  axis.text.x = element_text(angle = 45, size = 9, hjust = 1),
  panel.spacing = unit(0.35, "lines")
)

# labeller = labeller(panel = as_labeller(phrog.class.labels, label_wrap_gen(width = 20)))

###########################################
### 5. Show representative examples
### of architectural regimes
###########################################

quantile.probs <- c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
seg.quant <- quantile(reprseq.arch.metrics$segregation.bary, probs = quantile.probs, na.rm = T)
meaneff.quant <- quantile(reprseq.arch.metrics$mean.eff.fam, probs = quantile.probs, na.rm = T)

# Example 1 (clean 2-block):
#   •	n.fams == 2
# •	mean.eff.fam < 1.5
# •	n.dom.segments == 2
# •	segregation.bary in top 5–10% of that subset

this.ex.dir <- file.path(this.part.dir, "vis-ecf-distr", "example1_clean2block")
make.dir(this.ex.dir)
ex1.shortlist <- reprseq.arch.metrics[n.fams == 2 &
                                        mean.eff.fam < 1.5 &
                                        n.dom.segments == 2 &
                                        segregation.bary >= seg.quant["90%"] &
                                        n.ecf.hits > 20]

for(index in 1:nrow(ex1.shortlist)){
  protein_name <- ex1.shortlist$qname[index]
  protein.ecf.hits.all <- reprseq.vs.ecfs[.(protein_name)]
  title <- sprintf("%s (%s)", protein_name, ex1.shortlist$funct[index])
  this.plot <- plot.ecf.distr.over.query(protein_name, protein.ecf.hits.all, title = title)
  this.plot.outfile <- file.path(this.ex.dir, sprintf("%s.png", protein_name))
  ggsave(this.plot.outfile, this.plot, width = 6, height = 5)
}

# Handpicked:
# reprseq012506 (portal)
# reprseq035975 (unknown)
# reprseq097228 (DNA methyltransferase) <--
# reprseq080790 (endonuclease)

#########################################

# •	Example 2 (overlap):
#   •	n.fams >= 2
# •	mean.eff.fam high (top 10–20%)
# •	segregation.bary low (bottom 10–20%)
# •	(optional) n.dom.segments low-ish

this.ex.dir <- file.path(this.part.dir, "vis-ecf-distr", "example2_overlap")
make.dir(this.ex.dir)
ex2.shortlist <- reprseq.arch.metrics[n.fams >= 2 &
                                        mean.eff.fam >= meaneff.quant["90%"] &
                                        # n.dom.segments == 2 &
                                        segregation.bary <= seg.quant["10%"] &
                                        n.ecf.hits > 20]

for(index in 1:nrow(ex2.shortlist)){
  # cat(index,"\n")
  protein_name <- ex2.shortlist$qname[index]
  protein.ecf.hits.all <- reprseq.vs.ecfs[.(protein_name)]
  title <- sprintf("%s (%s)", protein_name, ex2.shortlist$funct[index])
  this.plot <- plot.ecf.distr.over.query(protein_name, protein.ecf.hits.all, title = title)
  this.plot.outfile <- file.path(this.ex.dir, sprintf("%s.png", protein_name))
  ggsave(this.plot.outfile, this.plot, width = 6, height = 5)
}

# Handpicked:
# reprseq010149 (unknown)
# reprseq046736 (ATPase) <--
# reprseq055118 (tail fibre)
# reprseq065324 (DNA methyltransferase)


###########################################

# •	Example 3 (fragmented):
#   •	n.dom.segments top 1–5%
# •	coverage decent (to avoid edge artifacts)

this.ex.dir <- file.path(this.part.dir, "vis-ecf-distr", "example3_mosaic")
make.dir(this.ex.dir)
ex3.shortlist <- reprseq.arch.metrics[n.fams >= 2 &
                                        # mean.eff.fam == meaneff.quant["10%"] &
                                        n.dom.segments >= 5 &
                                        segregation.bary >= seg.quant["90%"] &
                                        n.ecf.hits > 20]

for(index in 1:nrow(ex3.shortlist)){
  protein_name <- ex3.shortlist$qname[index]
  protein.ecf.hits_all <- reprseq.vs.ecfs[.(protein_name)]
  title <- sprintf("%s (%s)", protein_name, ex3.shortlist$funct[index])
  this.plot <- plot.ecf.distr.over.query(protein_name, protein.ecf.hits_all, title = title)
  this.plot.outfile <- file.path(this.ex.dir, sprintf("%s.png", protein_name))
  ggsave(this.plot.outfile, this.plot, width = 6, height = 5)
}

# Handpicked:
# reprseq016545 (unknown)
# reprseq029825 (TLTP) <-- 
# reprseq054784 (unknown)
# reprseq141653 (phosphoesterase)

################################################

# Final set:
# 1. Clean 2-block: 
# reprseq012506 (portal)
# reprseq035975 (unknown)
# reprseq097228 (DNA methyltransferase) <--
# reprseq080790 (endonuclease)

name <- "reprseq012506"
protein.ecf.hits_all <- reprseq.vs.ecfs[.(name)]
title <- sprintf("Clean 2-block\n(%s)", reprseq.arch.metrics[qname == name]$funct)
p1a <- plot.ecf.distr.over.query(name, protein.ecf.hits_all, title = title)

name <- "reprseq035975"
protein.ecf.hits_all <- reprseq.vs.ecfs[.(name)]
title <- sprintf("Clean 2-block\n(%s)", reprseq.arch.metrics[qname == name]$funct)
p1b <- plot.ecf.distr.over.query(name, protein.ecf.hits_all, title = title)

name <- "reprseq097228"
protein.ecf.hits_all <- reprseq.vs.ecfs[.(name)]
title <- sprintf("Clean 2-block\n(%s)", reprseq.arch.metrics[qname == name]$funct)
p1c <- plot.ecf.distr.over.query(name, protein.ecf.hits_all, title = title)

name <- "reprseq080790"
protein.ecf.hits_all <- reprseq.vs.ecfs[.(name)]
title <- sprintf("Clean 2-block\n(%s)", reprseq.arch.metrics[qname == name]$funct)
p1d <- plot.ecf.distr.over.query(name, protein.ecf.hits_all, title = title)

p1 <- (p1a | p1b | p1c | p1d)

# Final set:
# 2. Overlap
# reprseq010149 (unknown)
# reprseq046736 (ATPase) <--
# reprseq055118 (tail fibre)
# reprseq065324 (DNA methyltransferase)

name <- "reprseq010149"
protein.ecf.hits_all <- reprseq.vs.ecfs[.(name)]
title <- sprintf("Overlap\n(%s)", reprseq.arch.metrics[qname == name]$funct)
p2a <- plot.ecf.distr.over.query(name, protein.ecf.hits_all, title = title)

name <- "reprseq046736"
protein.ecf.hits_all <- reprseq.vs.ecfs[.(name)]
title <- sprintf("Overlap\n(%s)", reprseq.arch.metrics[qname == name]$funct)
p2b <- plot.ecf.distr.over.query(name, protein.ecf.hits_all, title = title)

name <- "reprseq055118"
protein.ecf.hits_all <- reprseq.vs.ecfs[.(name)]
title <- sprintf("Overlap\n(%s)", reprseq.arch.metrics[qname == name]$funct)
p2c <- plot.ecf.distr.over.query(name, protein.ecf.hits_all, title = title)

name <- "reprseq065324"
protein.ecf.hits_all <- reprseq.vs.ecfs[.(name)]
title <- sprintf("Overlap\n(%s)", reprseq.arch.metrics[qname == name]$funct)
p2d <- plot.ecf.distr.over.query(name, protein.ecf.hits_all, title = title)

p2 <- (p2a | p2b | p2c | p2d)

# Final set:
# 3. Fragmented
# reprseq016545 (unknown)
# reprseq029825 (TLTP) <-- 
# reprseq054784 (unknown)
# reprseq141653 (phosphoesterase)

name <- "reprseq016545"
protein.ecf.hits_all <- reprseq.vs.ecfs[.(name)]
title <- sprintf("Fragmented\n(%s)", reprseq.arch.metrics[qname == name]$funct)
p3a <- plot.ecf.distr.over.query(name, protein.ecf.hits_all, title = title)

name <- "reprseq029825"
protein.ecf.hits_all <- reprseq.vs.ecfs[.(name)]
title <- sprintf("Fragmented\n(%s)", reprseq.arch.metrics[qname == name]$funct)
p3b <- plot.ecf.distr.over.query(name, protein.ecf.hits_all, title = title)

name <- "reprseq054784"
protein.ecf.hits_all <- reprseq.vs.ecfs[.(name)]
title <- sprintf("Fragmented\n(%s)", reprseq.arch.metrics[qname == name]$funct)
p3c <- plot.ecf.distr.over.query(name, protein.ecf.hits_all, title = title)

name <- "reprseq141653"
protein.ecf.hits_all <- reprseq.vs.ecfs[.(name)]
title <- sprintf("Fragmented\n(%s)", reprseq.arch.metrics[qname == name]$funct)
p3d <- plot.ecf.distr.over.query(name, protein.ecf.hits_all, title = title)

p3 <- (p3a | p3b | p3c | p3d)


p.examples <- (p1 / p2 / p3) &
  theme(
    plot.tag = element_text(face = "bold", size = 24),
    plot.tag.position = c(0.02, 0.98)
  )

###########################################
### 6. Correlation structure of 
### architecture metrics
###########################################

corr.vars <- c(
  "length",
  "n.ecf.hits",
  "coverage.fraction",
  "n.fams",
  "mean.eff.fam",
  "segregation.bary",
  "segregation.overlap",
  "n.dom.segments"
)

corr.dt <- copy(reprseq.arch.metrics[, ..corr.vars])
corr.mat <- suppressWarnings(cor(corr.dt, use = "pairwise.complete.obs", method = "spearman"))


corr.long <- as.data.table(as.table(corr.mat))
setnames(corr.long, c("V1", "V2", "N"), c("x", "y", "rho"))

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

###########################################
### 7. ECF depth versus architecture metrics
###########################################

p.depth.2a <- ggplot(reprseq.arch.metrics[coverage.fraction > 0], aes(x = n.ecf.hits, y = n.fams)) +
  geom_bin2d(bins = 60) +
  scale_fill_viridis_c(trans = "sqrt", name = "Count") +
  scale_x_continuous(trans = "log10") +
  theme_bw(base_size = 12) +
  labs(
    x = "Number of ECF hits (log10 scale)",
    y = "Number of ECF families"
  )

p.depth.2b <- ggplot(reprseq.arch.metrics[coverage.fraction > 0 & !is.na(mean.eff.fam)], aes(x = n.ecf.hits, y = mean.eff.fam)) +
  geom_bin2d(bins = 60) +
  scale_fill_viridis_c(trans = "sqrt", name = "Count") +
  scale_x_continuous(trans = "log10") +
  theme_bw(base_size = 12) +
  labs(
    x = "Number of ECF hits (log10 scale)",
    y = "Mean effective number of families"
  )

p.depth.2c <- ggplot(reprseq.arch.metrics[coverage.fraction > 0 & !is.na(n.dom.segments)], aes(x = n.ecf.hits, y = n.dom.segments)) +
  geom_bin2d(bins = 60) +
  scale_fill_viridis_c(trans = "sqrt", name = "Count") +
  scale_x_continuous(trans = "log10") +
  theme_bw(base_size = 12) +
  labs(
    x = "Number of ECF hits (log10 scale)",
    y = "Number of dominant-family segments"
  )

p.depth <- (p.depth.2a / p.depth.2b / p.depth.2c) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    plot.tag.position = c(0.02, 0.98)
  )

###########################################
### 8. Correlation structure of 
### architecture metrics
###########################################

class.plot.dt <- copy(reprseq.arch.metrics[class != "unknown"])
class.plot.dt[, class := factor(class, levels = phrogs.class.order)]

vars <- c("n.fams", "mean.eff.fam", "n.dom.segments")
quantile.probs <- c(0, 0.5, 0.75, 0.9, 0.95, 0.99)
res <- class.plot.dt[
  , lapply(.SD, function(x)
    quantile(x, probs = quantile.probs, na.rm = TRUE)
  ),
  by = class,
  .SDcols = vars
]

strip.y <- ggh4x::strip_themed(
  background_x = ggh4x::elem_list_rect(
    fill = grDevices::adjustcolor(
      phrogs.col[match(levels(class.plot.dt$class), phrogs.class)],
      alpha.f = 0.25
    ),
    colour = "black",
    linewidth = 0.4
  )
)

make.class.metric.plot <- function(metric, ylab){
  ggplot(class.plot.dt[!is.na(get(metric))], aes(x = class, y = get(metric))) +
    geom_violin(trim = TRUE, fill = "grey85", colour = "grey20", linewidth = 0.3) +
    geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.3) +
    ggh4x::facet_grid2(
      cols = vars(class),
      scales = "free_x",
      space = "free_x",
      labeller = labeller(class = as_labeller(phrog.class.labels, label_wrap_gen(width = 20))),
      strip = strip.y
    ) +
    labs(x = NULL, y = ylab) +
    theme_classic(base_size = 12) +
    theme(
      strip.placement = "outside",
      strip.text.x = element_text(angle = 0, hjust = 0, face = "bold", size = 9),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing = unit(0.35, "lines")
    )
}

p.arch.class.3a <- make.class.metric.plot("n.fams", "Number of ECF families")
p.arch.class.3b <- make.class.metric.plot("mean.eff.fam", "Mean effective number of families")
p.arch.class.3c <- make.class.metric.plot("n.dom.segments", "Number of dominant-family segments")

p.arch.class <- (p.arch.class.3a / p.arch.class.3b / p.arch.class.3c) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    plot.tag.position = c(0.02, 0.98)
  )

##################

# Figure 6:
p.fam.div.vs.ecf.depth # PanelA
p.fam.div.spatial.separation.bary # PanelB
p.segreg.across.fam # PanelC
p.examples # panelD

fig <- (
  (p.fam.div.vs.ecf.depth | p.fam.div.spatial.separation.bary) /
  p.segreg.across.fam /
    wrap_elements(p.examples)
) + 
  plot_layout(heights =  c(1.2, 1, 3.5),
              widths = c(1, 1)) +   # middle number controls spacing
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 24),
    plot.tag.position = c(0.1, 1.05),
    plot.margin =  margin(10, 10, 10, 10)
  )
fig.outfile <- file.path(ecf.figures.dir, "figure5.png")
ggsave(fig.outfile, fig, width = 14, height = 18)

### Supplementary Figures
p.aho.top.funct # Within-ECF-family mean AHO in top function
p.aho.dist # Distribution of all Within-ECF-family AHO values
p.fam.div.spatial.separation.overlap # Support for Panel B using another segregation metric
p.sup.arch.corr # Correlation matrix of architectures
p.depth # ECF depth vs arch stats
p.arch.class # Distribution of archtecture metrics across classes

ggsave(file.path(ecf.suppl.figures.dir, "suppl-5p1_within-fam-aho.png"),
       p.aho.top.funct, width = 6, height = 8)

ggsave(file.path(ecf.suppl.figures.dir, "suppl-5p2_aho-distr.png"),
       p.aho.dist, width = 5, height = 4)

ggsave(file.path(ecf.suppl.figures.dir, "suppl-5p3_fam-div-seg-overlap.png"),
       p.fam.div.spatial.separation.overlap, width = 6, height = 5)

ggsave(file.path(ecf.suppl.figures.dir, "suppl-5p4_corr-arch.png"),
       p.sup.arch.corr, width = 7, height = 6)

ggsave(file.path(ecf.suppl.figures.dir, "suppl-5p5_depth-vs-arch.png"),
       p.depth, width = 6, height = 13)

ggsave(file.path(ecf.suppl.figures.dir, "suppl-5p6_corr-arch.png"),
       p.arch.class, width = 12, height = 16)

