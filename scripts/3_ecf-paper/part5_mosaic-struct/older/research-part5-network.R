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

protein_metrics_essential <- fread("protein-metrics-tmp.txt")

# -----------------------------
# 0) Basic hygiene
# -----------------------------
stopifnot(all(c("qname","n_ecf_hits","coverage_fraction","n_families_total",
                "mean_eff_families","segregation_bary","n_dom_segments",
                "funct","class") %in% names(protein_metrics_essential)))

stopifnot(all(c("qname","family") %in% names(reprseq.vs.ecfs)))

# Families per protein (distinct)
fam_by_prot <- unique(reprseq.vs.ecfs[, .(qname, family)])
setkey(fam_by_prot, qname)

# -----------------------------
# 1) Define "mosaic modular" proteins
#    (simple + tunable thresholds)
# -----------------------------
# Recommend: require decent coverage + >=2 families + >=2 segments + some separation
# Use quantiles so it adapts to your dataset.

cov_min   <- 0.6
hits_min  <- 5
fam_min   <- 2
seg_min_q <- 0.75  # "high separation" = top quartile of segregation_bary
seg_min   <- quantile(protein_metrics_essential[!is.na(segregation_bary)]$segregation_bary,
                      probs = seg_min_q, na.rm = TRUE)
me_max   <- 1.1
segN_min  <- 2

mosaic_proteins <- protein_metrics_essential[
  coverage_fraction >= cov_min &
    n_ecf_hits >= hits_min &
    n_families_total >= fam_min &
    n_dom_segments >= segN_min &
    mean_eff_families <= me_max &
    !is.na(segregation_bary) & segregation_bary >= seg_min,
  .(qname, funct, class, n_ecf_hits, n_families_total, mean_eff_families,
    segregation_bary, n_dom_segments)
]
setkey(mosaic_proteins, qname)


# -----------------------------
# 2) Family co-occurrence in mosaic proteins
# -----------------------------
fam_mosaic <- fam_by_prot[mosaic_proteins, on="qname", nomatch=0]
# count in how many mosaic proteins each family appears
fam_occ <- fam_mosaic[, .(n_proteins = uniqueN(qname)), by = family]
setorder(fam_occ, -n_proteins)

# co-occurrence pairs within each protein
fam_list <- fam_mosaic[, .(fams = list(sort(unique(family)))), by = qname]
pair_dt <- fam_list[, {
  f <- fams[[1]]
  if (length(f) < 2) return(NULL)
  ij <- t(combn(f, 2))
  data.table(f1 = ij[,1], f2 = ij[,2])
}, by = qname][, qname := NULL]

# edge weights: number of mosaic proteins where pair co-occurs
edges <- pair_dt[, .(w = .N), by = .(f1, f2)]
setorder(edges, -w)

# Family–family co-occurrence edges (node = family)
edges <- fam_list[
  lengths(fams) >= 2,
  {
    v <- sort(unique(fams[[1]]))
    pairs <- t(combn(v, 2))
    data.table(from = pairs[, 1], to = pairs[, 2])
  },
  by = qname
][
  , .(weight = .N), by = .(from, to)
][
  order(-weight, from, to)
]

# Save for Cytoscape
fwrite(edges, file.path(this.part.dir, "ecf_family_cooccurrence_edges.tsv"), sep = "\t")

# Optional: node table
nodes <- data.table(id = sort(unique(c(edges$from, edges$to))))
ecf.classification.table <- fread(ecf.classification.table.output)
ecf.classification.table.ef <- merge(ecf.classification.table, ecf.family.table, by = "ecf")
ef_1 <- ecf.classification.table.ef[!duplicated(family)]
nodes.out <- merge(nodes, ef_1, by.x = "id", by.y = "family", all.x = T)
fwrite(nodes.out, file.path(this.part.dir, "ecf_family_cooccurrence_nodes.tsv"), sep = "\t")





