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
this.part.dir <- file.path(ecf.proj.dir, "part7_funct-rescue")
make.dir(this.part.dir)

#################################
### Load all data:
### 1. reprseq data
### 2. ECF tables
### 3. reprseq functions
### 4. reprseq vs. ECFs
### 5. reprseq vs. domains
#################################

# Load refseq table
reprseq.table <- fread(reprseq.table.output)

# ECF tables
ecf.table <- fread(ecf.table.output)
ecf.family.table <- fread(ecf.family.table.output)

# Load functional annotations
reprseq.functions.unique <- fread(reprseq.functions.unique.output)
reprseq.functions <- fread(reprseq.functions.all.output)

# Load reprseq vs ECFs
reprseq.vs.ecfs.all <- fread(reprseq.ecfs.output)
reprseq.vs.ecod.all <- fread(reprseq.domains.ecod.output)
reprseq.vs.pfam.all <- fread(reprseq.domains.pfam.output)

#################################
### Process loaded tables
#################################

all.reprseq <- unique(reprseq.table$reprseq)
setkey(reprseq.table, reprseq)

setkey(ecf.family.table, ecf)

setkey(reprseq.functions, qname)
setkey(reprseq.functions.unique, reprseq)

reprseq.vs.ecfs <- reprseq.vs.ecfs.all %>% 
  filter(prob >= MIN.PROB.FOR.DOMAIN.DETECT.IN.RHMM,
         eval <= MAX.EVAL.FOR.DOMAIN.DETECT.IN.RHMM, 
         qcov >= MIN.RHMM.COV.FOR.DOMAIN.DETECT.IN.RHMM,
         scov >= MIN.DOMAIN.COV.FOR.DOMAIN.DETECT.IN.RHMM,
         hit.length >= MIN.HIT.LENGTH.FOR.DOMAIN.DETECT.IN.RHMM
  )
reprseq.vs.ecfs$family <- ecf.family.table[.(reprseq.vs.ecfs$sname)]$family
setkey(reprseq.vs.ecfs, qname)
# reprseq.vs.ecfs[, family := ecf.family.table[.(reprseq.vs.ecfs$sname)]$family]

reprseq.vs.ecod <- reprseq.vs.ecod.all %>% 
  filter(
    eval <= MAX.EVAL.FOR.DOMAIN.DETECT.IN.RHMM, 
    hit.length >= MIN.HIT.LENGTH.FOR.DOMAIN.DETECT.IN.RHMM
  )
reprseq.vs.pfam <- reprseq.vs.pfam.all %>% 
  filter(
    eval <= MAX.EVAL.FOR.DOMAIN.DETECT.IN.RHMM, 
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

############################
### 1. functional enrichment in DD proteins
############################

# =============================================================================
# FUNCTIONAL ENRICHMENT OF ECF FAMILIES IN DOMAIN-DARK PROTEINS
# =============================================================================

MIN.SCORE        <- 0.5
MIN.DD.WITH.FAM  <- 5L    # min domain-dark proteins carrying the family with that function
MIN.DD.TOTAL.FUNCT <- 20L # min domain-dark proteins with that function overall
MAX.QVAL         <- 0.01  # FDR threshold
TOP_FUNCT        <- 25L   # top functions by domain-dark protein count to consider

# -----------------------------------------------------------------------------
# Step 1: Build the domain-dark protein × family × function table
# -----------------------------------------------------------------------------

# Best function per domain-dark protein (highest score, top.annotation only)
dd.funct <- reprseq.functions[
  qname %chin% reprseq.dd &
    top.annotation == TRUE  &
    score >= MIN.SCORE
][
  order(qname, -score)
][
  , .SD[1], by = qname   # best hit per protein
][
  , .(qname, funct, class, score)
]

# Families present in each domain-dark protein (one row per protein-family pair)
dd.fam <- reprseq.fam[qname %chin% reprseq.dd]

# Number of domain-dark proteins with each function (for filtering)
dd.funct.counts <- dd.funct[, .(n_dd_funct = .N), by = funct]
top.functs <- dd.funct.counts[
  order(-n_dd_funct)
][1:min(.N, TOP_FUNCT)]$funct

# Restrict to top functions
dd.funct.top <- dd.funct#[funct %chin% top.functs]

# Total domain-dark proteins with any top function annotation
N.dd.annotated <- uniqueN(dd.funct.top$qname)
N.dd.total     <- length(reprseq.dd)

# Restrict to only annotated domain-dark proteins
reprseq.dd.annotated <- dd.funct.top$qname  # only those with a top function

# Recompute N.dd.total as annotated only
N.dd.annotated.only <- length(unique(reprseq.dd.annotated))

# -----------------------------------------------------------------------------
# Step 2: For each ECF family × function pair, build 2x2 contingency table
# and run Fisher's exact test
# -----------------------------------------------------------------------------

# All domain-dark proteins: annotate with function (NA if unannotated)
dd.base <- data.table(qname = reprseq.dd)
dd.base <- merge(dd.base, dd.funct.top[, .(qname, funct)], 
                 by = "qname", all.x = TRUE)
# dd.base <- dd.base[!is.na(funct)]

# Get all families present in domain-dark proteins
all.dd.families <- dd.fam[, unique(family)]

# For each family, which domain-dark proteins carry it?
fam.members <- dd.fam[, .(qname = list(qname)), by = family]
setkey(fam.members, family)

# # Run enrichment for each family × function pair
# enrich.results <- rbindlist(lapply(all.dd.families, function(fam) {
#   
#   # Proteins carrying this family
#   has.fam    <- fam.members[fam]$qname[[1]]
#   n.has.fam  <- length(has.fam)
#   if (n.has.fam < MIN.DD.WITH.FAM) return(NULL)
#   
#   rbindlist(lapply(top.functs, function(fn) {
#     
#     # 2x2 table:
#     #               has function   no function
#     # has family:       a               b
#     # no family:        c               d
#     
#     a <- sum(has.fam %chin% dd.funct.top[funct == fn]$qname)
#     if (a < MIN.DD.WITH.FAM) return(NULL)
#     
#     b <- n.has.fam - a
#     c <- dd.funct.counts[funct == fn]$n_dd_funct - a
#     d <- N.dd.total - a - b - c
#     # d = N.dd.annotated.only - a - b - c
#     if (c < 0 || d < 0) return(NULL)
#     
#     mat <- matrix(c(a, b, c, d), nrow = 2)
#     ft  <- fisher.test(mat, alternative = "greater")
#     
#     data.table(
#       family    = fam,
#       funct     = fn,
#       n_dd_fam_funct  = a,   # domain-dark with both family AND function
#       n_dd_fam        = n.has.fam,
#       n_dd_funct      = a + c,
#       odds.ratio      = ft$estimate,
#       pval            = ft$p.value,
#       frac_funct_in_fam = a / n.has.fam  # of proteins with this family, fraction with this function
#     )
#   }))
# }))

# -----------------------------------------------------------------------------
# Step 2: For each ECF family x function pair, build 2x2 contingency table
# and run Fisher's exact test
#
# UPDATE: restrict to domain-dark proteins with a confident functional call
# (present in dd.funct.top) before building the table. Pooling "different
# function" and "no confident annotation" into b conflates functional
# specificity with annotation completeness -- family membership predicted
# annotation completeness directly for EF0137 (OR ~ 47, p < 2.2e-16), which
# would otherwise inflate odds ratios for families that are easier to
# annotate for reasons unrelated to the function being tested.
# -----------------------------------------------------------------------------

enrich.results <- rbindlist(lapply(all.dd.families, function(fam) {
  
  # All domain-dark proteins carrying this family
  has.fam.all   <- fam.members[fam]$qname[[1]]
  n.has.fam.all <- length(has.fam.all)
  
  # Restrict to those with a confident top-25 functional call --
  # this is the population the 2x2 table is now built from
  has.fam   <- has.fam.all[has.fam.all %chin% dd.funct.top$qname]
  n.has.fam <- length(has.fam)
  
  if (n.has.fam < MIN.DD.WITH.FAM) return(NULL)
  
  rbindlist(lapply(top.functs, function(fn) {
    
    # 2x2 table (annotated domain-dark proteins only):
    #               has function   no function
    # has family:       a               b
    # no family:        c               d
    
    a <- sum(has.fam %chin% dd.funct.top[funct == fn]$qname)
    if (a < MIN.DD.WITH.FAM) return(NULL)
    
    b <- n.has.fam - a
    c <- dd.funct.counts[funct == fn]$n_dd_funct - a
    d <- N.dd.annotated.only - a - b - c
    if (c < 0 || d < 0) return(NULL)
    
    mat <- matrix(c(a, b, c, d), nrow = 2)
    ft  <- fisher.test(mat, alternative = "greater")
    
    data.table(
      family             = fam,
      funct              = fn,
      n_dd_fam_funct     = a,                          # annotated DD proteins: family AND function
      n_dd_fam           = n.has.fam,                  # annotated DD proteins with this family
      n_dd_fam_total     = n.has.fam.all,               # ALL DD proteins with this family (annotated or not)
      frac_fam_annotated = n.has.fam / n.has.fam.all,    # diagnostic: annotation completeness for this family
      n_dd_funct         = a + c,
      odds.ratio         = ft$estimate,
      pval               = ft$p.value,
      frac_funct_in_fam  = a / n.has.fam                # of ANNOTATED proteins with this family, fraction with this function
    )
  }))
}))

# -----------------------------------------------------------------------------
# Step 3: Multiple testing correction and filtering
# -----------------------------------------------------------------------------

enrich.results[, qval := p.adjust(pval, method = "BH")]
# enrich.results[, qval := p.adjust(pval, method = "bonferroni")]
enrich.sig <- enrich.results[
  qval        <= MAX.QVAL     &
    odds.ratio  >  1            &
    n_dd_fam_funct >= MIN.DD.WITH.FAM
]
setorder(enrich.sig, funct, -odds.ratio)

# -----------------------------------------------------------------------------
# Step 4: Summarise per function — how many significant ECF families?
# and what is the functional diversity (unique families per function)
# -----------------------------------------------------------------------------

# Total unique ECF families per function across ALL domain-dark proteins
fam.diversity.per.funct <- merge(
  dd.fam, 
  dd.funct.top[, .(qname, funct)], 
  by = "qname"
)[, .(
  n_unique_fams = uniqueN(family),
  n_dd_proteins = uniqueN(qname)
), by = funct]

# Significant enriched families per function
sig.per.funct <- enrich.sig[, .(
  n_sig_fams       = .N,
  top_family       = family[which.max(odds.ratio)],
  top_odds         = max(odds.ratio),
  top_frac         = frac_funct_in_fam[which.max(odds.ratio)]
), by = funct]

funct.summary <- merge(fam.diversity.per.funct, sig.per.funct, 
                       by = "funct", all.x = TRUE)
funct.summary[is.na(n_sig_fams), n_sig_fams := 0]
setorder(funct.summary, -n_unique_fams)

# -----------------------------------------------------------------------------
# Step 5: Visualisation — two-panel figure
# Panel A: function diversity vs rescue (scatter)
# Panel B: top enriched families per function (dot plot)
# -----------------------------------------------------------------------------

# --- Panel A: scatter showing ECF family diversity vs fraction domain-dark
# with bubble size = n_dd_proteins

p.scatter <- ggplot(
  funct.summary,
  aes(x = n_unique_fams, 
      y = n_sig_fams,
      size = n_dd_proteins,
      label = funct)
) +
  geom_point(aes(color = n_sig_fams > 0), alpha = 0.8) +
  ggrepel::geom_text_repel(
    size = 2.8, max.overlaps = 15,
    segment.color = "grey70"
  ) +
  scale_color_manual(
    values = c("FALSE" = "grey70", "TRUE" = "#4E79A7"),
    labels = c("No enriched ECF families", "Has enriched ECF families"),
    name   = NULL
  ) +
  scale_size_continuous(
    name   = "Domain-dark\nproteins",
    range  = c(2, 10)
  ) +
  labs(
    x = "Unique ECF families in domain-dark proteins",
    y = "Significantly enriched ECF families"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        axis.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# --- Panel B: dot plot of top enriched families per function
# Show top 3 families per function, sized by n_dd_fam_funct, 
# colored by odds ratio

top3.per.funct <- enrich.sig[
  , .SD[order(-odds.ratio)][1:min(.N, 3L)],
  by = funct
]

# Order functions by ECF family diversity (links to Panel A)
funct.order.b <- funct.summary[order(-n_unique_fams)]$funct
top3.per.funct[, funct := factor(funct, levels = funct.order.b)]

# Order families within each function by odds ratio
top3.per.funct[, family.rank := frank(-odds.ratio, ties.method = "first"), 
               by = funct]
top3.per.funct[, family.label := paste0(family, "\n(n=", n_dd_fam_funct, ")")]

p.dotplot <- ggplot(
  top3.per.funct,
  aes(x     = factor(family.rank),
      y     = funct,
      size  = n_dd_fam_funct,
      color = log2(odds.ratio))
) +
  geom_point(alpha = 0.85) +
  geom_text(
    aes(label = family),
    size = 2.0, vjust = -1.2, color = "black"
  ) +
  scale_color_gradient2(
    low      = "steelblue",
    mid      = "white",
    high     = "#E15759",
    midpoint = 1,
    name     = "log2(Odds ratio)"
  ) +
  scale_size_continuous(
    name  = "Domain-dark\nproteins",
    range = c(2, 10)
  ) +
  scale_x_discrete(labels = c("1" = "Top", "2" = "2nd", "3" = "3rd")) +
  labs(
    x     = "Enrichment rank",
    y     = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position   = "right",
    axis.title = element_text(face = "bold"),
    axis.text.y       = element_text(size = 9)
  )

# --- Combine
p.combined <- p.scatter / p.dotplot +
  plot_layout(heights = c(0.8, 1.2)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 24)
  )

# print(p.combined)

fig.outfile <- file.path(ecf.figures.dir, "figure5.png")
ggsave(
  fig.outfile,
  p.combined,
  width = 12, height = 20, units = "cm", scale = 1.5
)

# # --- Save enrichment table
# fwrite(enrich.sig,     file.path(this.part.dir, "ecf_funct_enrichment.csv"))
# fwrite(funct.summary,  file.path(this.part.dir, "ecf_funct_summary.csv"))


