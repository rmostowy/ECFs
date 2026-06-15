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
reprseq.functions <- fread(reprseq.functions.all.output)
reprseq.functions.unique <- fread(reprseq.functions.unique.output)
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
### Plot ECF distributions in 
### most common reprseq functions
############################################

# Define plotting function
# Function to visualise ECF domains and alignment on the query
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

# define output dir
ecfs.per.funct.dir <- file.path(this.part.dir, "vis-ecf-distr")
make.dir(ecfs.per.funct.dir)

set.seed(1234)

# 1. Get 20 most common functions, and then for each function, 
# randomly draw 'sample.size.per.funct' reprseqs from the reprseq.vs.ecfs table
top.functs <- names(my.table(reprseq.functions[top.annotation == TRUE]$funct)[1:20])
sample.size.per.funct <- 10
reprseqs.to.plot <- lapply(1:length(top.functs), function(fi){
  this.funct <- top.functs[fi]
  this.funct.reprseqs.all <- reprseq.functions[funct == this.funct]$qname
  this.funct.reprseqs.with.ecf <- this.funct.reprseqs.all[this.funct.reprseqs.all %in% reprseq.vs.ecfs$qname]
  this.funct.reprseqs.to.plot <- sample(this.funct.reprseqs.with.ecf, sample.size.per.funct)
})
for(fi in 1:length(top.functs)){
  this.funct <- top.functs[fi]
  this.funct.reprseqs <- reprseqs.to.plot[[fi]]
  this.funct.dir <- file.path(ecfs.per.funct.dir, sprintf("funct_%*.2d", 2, fi))
  make.dir(this.funct.dir)
  for(ri in 1:sample.size.per.funct){
    this.reprseq <- this.funct.reprseqs[ri]
    this.reprseq.ecf.hits <- reprseq.vs.ecfs[.(this.reprseq)]
    title <- sprintf("%s (%s)", this.reprseq, this.funct)
    this.reprseq.plot <- plot.ecf.distr.over.query(this.reprseq, this.reprseq.ecf.hits, title = title)
    this.reprseq.plot.outfile <- file.path(this.funct.dir, sprintf("%s.png", this.reprseq))
    ggsave(this.reprseq.plot.outfile, this.reprseq.plot, width = 6, height = 5)
  }
}

# 2. Draw 'sample.from.all.reprseqs' random reprseqs (those with annot and ECFs)
sample.from.all.reprseqs <- 100

reprseqs.with.funct.and.ecf <- sort(unique(intersect(reprseq.functions$qname, reprseq.vs.ecfs$qname)))
random.sample <- sort(sample(reprseqs.with.funct.and.ecf, sample.from.all.reprseqs))
this.funct.all.dir <- file.path(ecfs.per.funct.dir, "all")
make.dir(this.funct.all.dir)
for(ri in 1:sample.from.all.reprseqs){
  this.reprseq <- random.sample[ri]
  this.reprseq.funct <- reprseq.functions.unique[reprseq == this.reprseq]$funct
  this.reprseq.ecf.hits <- reprseq.vs.ecfs[.(this.reprseq)]
  title <- sprintf("%s (%s)", this.reprseq, this.reprseq.funct)
  this.reprseq.plot <- plot.ecf.distr.over.query(this.reprseq, this.reprseq.ecf.hits, title = title)
  this.reprseq.plot.outfile <- file.path(this.funct.all.dir, sprintf("%s.png", this.reprseq))
  ggsave(this.reprseq.plot.outfile, this.reprseq.plot, width = 6, height = 5)
}




############################################
### Compute within-EF AHO per reprseq 
############################################

# #######
# ##|||||
# ## TRIAL RUN
# 
# this.reprseq <- "reprseq134177"
# # this.reprseq <- "reprseq057676"
# this.reprseq.ecf.hits <- reprseq.vs.ecfs[.(this.reprseq)]
# 
# reprseq.ecf1.overlap <- this.reprseq.ecf.hits %>% select(qname, sname.ecf1 = sname, qstart.ecf1 = qstart, qend.ecf1 = qend, family)
# reprseq.ecf2.overlap <- this.reprseq.ecf.hits %>% select(qname, sname.ecf2 = sname, qstart.ecf2 = qstart, qend.ecf2 = qend, family)
# this.reprseq.ecf.overlaps <- reprseq.ecf1.overlap[reprseq.ecf2.overlap, on = .(qname, family), allow.cartesian = TRUE, nomatch = 0L]
# this.reprseq.ecf.overlaps <- this.reprseq.ecf.overlaps[sname.ecf1 < sname.ecf2]
# # this.reprseq.ecf.overlaps[, `:=`(
# #   overlap.start = pmax(qstart.ecf1, qstart.ecf2),
# #   overlap.end = pmin(qend.ecf1, qend.ecf2),
# #   overlap.length = pmax(0, overlap.end - overlap.start + 1)
# # )]
# # overlap.length <- pmax(0, overlap.end - overlap.start + 1)
# this.reprseq.ecf.overlaps[, aho := calculate_overlap_jaccard(qstart.ecf1, qend.ecf1, qstart.ecf2, qend.ecf2)]
# 
# ##|||||
# #######

MIN.NO.ECF.PER.FAMILY=3

# Append stats about number of ecfs per family in each reprseq
reprseq.vs.ecfs.append <- merge(reprseq.vs.ecfs, reprseq.vs.ecfs[, .N, by = .(qname, family)], by = c("qname", "family"))
reprseq.vs.ecfs.filt <- reprseq.vs.ecfs.append[N >= MIN.NO.ECF.PER.FAMILY]

reprseq.ecf.overlap.from <- reprseq.vs.ecfs.filt %>% select(qname, sname.ecf1 = sname, qstart.ecf1 = qstart, qend.ecf1 = qend, family)
reprseq.ecf.overlap.to <- reprseq.vs.ecfs.filt %>% select(qname, sname.ecf2 = sname, qstart.ecf2 = qstart, qend.ecf2 = qend, family)
reprseq.ecf.overlaps <- reprseq.ecf.overlap.from[reprseq.ecf.overlap.to, on = .(qname, family), allow.cartesian = TRUE, nomatch = 0L]
reprseq.ecf.overlaps <- reprseq.ecf.overlaps[sname.ecf1 < sname.ecf2]
reprseq.ecf.overlaps[, aho := calculate_overlap_jaccard(qstart.ecf1, qend.ecf1, qstart.ecf2, qend.ecf2)]
reprseq.ef.aho.stats <- reprseq.ecf.overlaps[, .(mean.aho = mean(aho)), by = .(qname, family)]


# example.function <- "endonuclease"
# 
# # 1. find all reprseqs assigned to this function
# all.reprseqs.with.that.function <- reprseq.functions[funct == "endonuclease"]$qname
# # 2. aho stats
# aho.stats <- reprseq.ef.aho.stats[qname %chin% all.reprseqs.with.that.function]$mean.aho

NF <- 10

# 0) proteins per function (within class) — unique qname like your example
fun2prot <- unique(reprseq.functions[, .(class, funct, qname)])

# 1) count proteins per function and pick top NF per class
top_fun <- fun2prot[, .(n_proteins = uniqueN(qname)), by = .(class, funct)]
setorder(top_fun, class, -n_proteins)
top_fun <- top_fun[, head(.SD, NF), by = class]
setkey(top_fun, class, funct)

# 2) map (class, funct) -> proteins (qname), then attach AHO rows (qname,family)
top_fun2prot <- fun2prot[top_fun, on = .(class, funct), nomatch = 0L]
setkey(top_fun2prot, qname)

# join to AHO (this matches your example: take all AHO rows for proteins in the function)
top_aho <- reprseq.ef.aho.stats[top_fun2prot, on = "qname", allow.cartesian = TRUE]

# top_aho now has: qname, family, mean.aho, class, funct

# 3) summarize AHO distribution per (class, funct)
top_fun_aho_stats <- top_aho[, .(
  n_proteins = uniqueN(qname),
  n_pf = .N,  # protein-family observations
  mean_aho = mean(mean.aho, na.rm = TRUE),
  median_aho = median(mean.aho, na.rm = TRUE),
  q25 = quantile(mean.aho, 0.25, na.rm = TRUE, type = 2),
  q75 = quantile(mean.aho, 0.75, na.rm = TRUE, type = 2)
), by = .(class, funct)]

# keep protein counts from top_fun (optional sanity)
setkey(top_fun_aho_stats, class, funct)
top_fun_aho_stats <- top_fun[top_fun_aho_stats, on = .(class, funct)]
top_fun_aho_stats[]

#########################
### PLOT
#########################

# Order functions within each class by mean_aho
top_fun_aho_stats[, funct_ord := reorder(funct, mean_aho), by = class]

p1 <- ggplot(top_fun_aho_stats, aes(x = mean_aho, y = funct_ord)) +
  geom_errorbarh(aes(xmin = q25, xmax = q75), height = 0.2, linewidth = 0.35) +
  geom_point(size = 2) +
  facet_grid(rows = vars(class), scales = "free_y", space = "free_y") +
  theme_bw(base_size = 12) +
  scale_x_continuous(limits = c(0,1)) +
  labs(x = "Mean AHO (across protein–family rows)",
       y = NULL,
       title = "AHO by predicted function (top functions per class)") +
  theme(panel.grid.minor = element_blank(),
        strip.text.y = element_text(angle = 0, face = "bold"))
p1


##############

# median per (class, funct)
ord <- top_aho[, .(med = median(mean.aho, na.rm=TRUE)), by=.(class, funct)]
setorder(ord, class, -med)

# add rank within class and join back
ord[, rank := seq_len(.N), by = class]
setkey(ord, class, funct)

setkey(top_aho, class, funct)
top_aho <- ord[top_aho]  # adds med + rank

# now define funct_ord per class using rank
top_aho[, funct_ord := factor(funct, levels = ord[class == .BY$class][order(rank), funct]), by = class]

# # build per-class levels table
# lev <- ord[order(class, rank), .(funct_levels = list(funct)), by = class]
# 
# # join levels then set factor (still avoids list-cols in by= for heavy ops; lev has only #classes rows)
# setkey(lev, class)
# top_aho <- lev[top_aho, on="class"]
# 
# top_aho[, funct_ord := factor(funct, levels = funct_levels[[1L]]), by = class]
# top_aho[, funct_levels := NULL]

p2 <- ggplot(top_aho, aes(x = mean.aho, y = funct_ord)) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.3) +
  facet_grid(rows = vars(class), scales = "free_y", space = "free_y") +
  theme_bw(base_size = 12) +
  labs(x = "AHO (protein–family)",
       y = NULL,
       title = "Distribution of AHO values per function (top functions per class)") +
  theme(panel.grid.minor = element_blank(),
        strip.text.y = element_text(angle = 0, face = "bold"))


# reprseq.ecf1.overlap <- this.reprseq.ecf.hits %>% select(qname, sname.ecf1 = sname, qstart.ecf1 = qstart, qend.ecf1 = qend, family.ecf1 = family)
# reprseq.ecf2.overlap <- this.reprseq.ecf.hits %>% select(qname, sname.ecf2 = sname, qstart.ecf2 = qstart, qend.ecf2 = qend, family.ecf2 = family)
# reprseq.ecf1.overlap[reprseq.ecf2.overlap, on = .(qname, family), allow.cartesian = TRUE, nomatch = 0L]
# 
# 
# ecfs.ecod2.overlap <- ecfs.full.vs.ecod %>% select(qname, sname.ecod2 = sname, t.id.ecod2 = t.id,
#                                                    sstart.ecod2 = sstart, send.ecod2 = send, slength.ecod2 = slength,
#                                                    scov.ecod2 = scov, pident.ecod2 = pident)
# shared.ecf.table <- ecfs.ecod1.overlap[ecfs.ecod2.overlap, on = .(qname), allow.cartesian = TRUE]
# 
# this.reprseq.ecf.hits[this.reprseq.ecf.hits, ]
# 






this.reprseq <- "reprseq000014"
reprseq.vs.ecfs[.(this.reprseq)]

protein_name <- this.reprseq
protein_ecf_hits_all <- reprseq.vs.ecfs[.(protein_name)]
title <- this.reprseq
plot.ecf.distr.over.query(protein_name, protein_ecf_hits_all, title = title)

