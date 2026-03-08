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
reprseq.vs.ecod.all <- fread(reprseq.domains.ecod.output)
reprseq.vs.pfam.all <- fread(reprseq.domains.pfam.output)


reprseq.vs.ecfs <- reprseq.vs.ecfs.all %>%
  filter(prob >= MIN.PROB.FOR.DOMAIN.EQUIVALENCE) %>%
  # filter(prob >= 99) %>%
  filter(scov >= MIN.SCOV.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE)
setkey(reprseq.vs.ecfs, qname)
reprseq.vs.ecfs$family <- ecf.family.table[.(reprseq.vs.ecfs$sname)]$family
# reprseq.vs.ecfs[, family := ecf.family.table[.(reprseq.vs.ecfs$sname)]$family]

reprseq.vs.ecod <- reprseq.vs.ecod.all %>% 
  filter(
         eval <= MAX.EVAL.FOR.DOMAIN.DETECT.IN.RHMM, 
        
         hit.length >= MIN.HIT.LENGTH.FOR.DOMAIN.DETECT.IN.RHMM
  )
reprseq.vs.pfam <- reprseq.vs.pfam.all %>% 
  filter(
         
         hit.length >= MIN.HIT.LENGTH.FOR.DOMAIN.DETECT.IN.RHMM
  )

reprseq.with.ecf.hits <- unique(reprseq.vs.ecfs$qname)
reprseq.with.ecod.hits <- unique(reprseq.vs.ecod$qname)
reprseq.with.pfam.hits <- unique(reprseq.vs.pfam$qname)

############################################
### 
############################################

# All representative sequences
all.reprseq <- unique(reprseq.table$reprseq)
N <- length(all.reprseq)

# Build annotation table
annot.dt <- data.table(
  qname = all.reprseq
)
annot.dt[, ECF  := qname %in% reprseq.with.ecf.hits]
annot.dt[, Pfam := qname %in% reprseq.with.pfam.hits]
annot.dt[, ECOD := qname %in% reprseq.with.ecod.hits]

# list all unique family combinations per protein
reprseq.fam <- unique(reprseq.vs.ecfs[, .(qname, family)])
setkey(reprseq.fam, qname)
reprseq.fam.comb <- reprseq.fam[, .(fams = list(sort(unique(family)))), by = qname]
reprseq.fam.comb[, fam_key := vapply(fams, paste, collapse = "|", FUN.VALUE = character(1))]

# Find domain-dark proteins
reprseq.dd.dt <- annot.dt[ECF == T & Pfam == F & ECOD == F]
reprseq.dd <- reprseq.dd.dt$qname

# List all unique family combinations in domain-dark proteins (domain-dark FCs)
reprseq.fam.comb.dd <- reprseq.fam.comb[qname %chin% reprseq.dd]

# Link all proteins with domain-dark FCs
reprseq.fam.comb.dd.all <- reprseq.fam.comb[fam_key %chin% reprseq.fam.comb.dd$fam_key]

# ---- parameters
TOP_N <- 20L

# attach top functional annotation
dd.funct <- reprseq.functions.unique[
  .(reprseq.fam.comb.dd.all$qname),
  .(qname = reprseq, funct, class)
]

dd.dt <- merge(reprseq.fam.comb.dd.all, dd.funct, by = "qname", all.x = TRUE)
dd.dt[is.na(funct) | funct == "", funct := "unknown"]
dd.dt[is.na(class) | class == "", class := "unknown"]

# ---- summarize family combinations
famcomp.dd.summary <- dd.dt[, .(
  n_proteins = .N,
  qnames = paste(qname, collapse = ";"),
  funct_counts = paste(
    .SD[, .N, by = funct][order(-N, funct),
                          paste0(funct, " (", N, ")")],
    collapse = "; "
  ),
  class_counts = paste(
    .SD[, .N, by = class][order(-N, class),
                          paste0(class, " (", N, ")")],
    collapse = "; "
  )
), by = fam_key]


# recover fams as list-column if you want it
famcomp.dd.summary[, fams := strsplit(fam_key, "|", fixed = TRUE)]

# sort by frequency
setorder(famcomp.dd.summary, -n_proteins)

# top N
famcomp.dd.top <- famcomp.dd.summary[1:min(.N, TOP_N)]

TOP_FUN <- 40

top_keys <- famcomp.dd.summary[order(-n_proteins)][1:min(.N, TOP_N)]$fam_key

hm.funct <- dd.dt[
  fam_key %chin% top_keys & funct != "unknown",
  .N,
  by = .(fam_key, funct)
]

# keep top functions overall for readability
top_functs <- hm.funct[, .(tot = sum(N)), by = funct][order(-tot)][1:min(.N, TOP_FUN)]$funct
hm.funct <- hm.funct[funct %chin% top_functs]

# order rows by frequency in famcomp.dd.summary
hm.funct[, fam_key := factor(fam_key, levels = famcomp.dd.summary[order(-n_proteins)][fam_key %chin% top_keys]$fam_key)]

# order columns by abundance
hm.funct[, funct := factor(funct, levels = hm.funct[, .(tot = sum(N)), by = funct][order(tot)]$funct)]

p.hm.funct <- ggplot(hm.funct, aes(x = funct, y = fam_key, fill = N)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(trans = "sqrt", name = "Proteins") +
  labs(x = "Predicted function", y = "ECF family combination") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_blank()
  )

p.hm.funct


# 
# 
# 
# reprseq.dd.ef <- fam_by_prot[reprseq.dd, on="qname", nomatch=0]
# # count in how many mosaic proteins each family appears
# reprseq.dd.fam.list <- reprseq.dd.ef[, .(fams = list(sort(unique(family)))), by = qname]
# 
# 
# 
# 
# reprseq.dd.dt <- annot.dt[ECF == T & Pfam == F & ECOD == F]
# reprseq.dd <- reprseq.dd.dt$qname
# 
# # Families per protein (distinct)
# fam_by_prot <- unique(reprseq.vs.ecfs[, .(qname, family)])
# setkey(fam_by_prot, qname)
# 
# reprseq.dd.ef <- fam_by_prot[reprseq.dd, on="qname", nomatch=0]
# # count in how many mosaic proteins each family appears
# reprseq.dd.fam.list <- reprseq.dd.ef[, .(fams = list(sort(unique(family)))), by = qname]
# 
# reprseq.dd.fam.funct <- cbind(reprseq.dd.fam.list, reprseq.functions.unique[.(reprseq.dd.fam.list$qname), .(funct, class)])
# reprseq.dd.fam.funct[funct == "", funct := "unknown"]
# reprseq.dd.fam.funct[class == "", class := "unknown"]
# 
# unique(reprseq.dd.fam.funct$fams)





########################################################################################
########################################################################################
########################################################################################

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





