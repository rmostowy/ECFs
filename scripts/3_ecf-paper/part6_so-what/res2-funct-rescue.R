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
this.part.dir <- file.path(ecf.proj.dir, "part6_so-what")
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
dd.funct.top <- dd.funct[funct %chin% top.functs]

# Total domain-dark proteins with any top function annotation
N.dd.annotated <- uniqueN(dd.funct.top$qname)
N.dd.total     <- length(reprseq.dd)

# -----------------------------------------------------------------------------
# Step 2: For each ECF family × function pair, build 2x2 contingency table
# and run Fisher's exact test
# -----------------------------------------------------------------------------

# All domain-dark proteins: annotate with function (NA if unannotated)
dd.base <- data.table(qname = reprseq.dd)
dd.base <- merge(dd.base, dd.funct.top[, .(qname, funct)], 
                 by = "qname", all.x = TRUE)

# Get all families present in domain-dark proteins
all.dd.families <- dd.fam[, unique(family)]

# For each family, which domain-dark proteins carry it?
fam.members <- dd.fam[, .(qname = list(qname)), by = family]
setkey(fam.members, family)

# Run enrichment for each family × function pair
enrich.results <- rbindlist(lapply(all.dd.families, function(fam) {
  
  # Proteins carrying this family
  has.fam    <- fam.members[fam]$qname[[1]]
  n.has.fam  <- length(has.fam)
  if (n.has.fam < MIN.DD.WITH.FAM) return(NULL)
  
  rbindlist(lapply(top.functs, function(fn) {
    
    # 2x2 table:
    #               has function   no function
    # has family:       a               b
    # no family:        c               d
    
    a <- sum(has.fam %chin% dd.funct.top[funct == fn]$qname)
    if (a < MIN.DD.WITH.FAM) return(NULL)
    
    b <- n.has.fam - a
    c <- dd.funct.counts[funct == fn]$n_dd_funct - a
    d <- N.dd.total - a - b - c
    if (c < 0 || d < 0) return(NULL)
    
    mat <- matrix(c(a, b, c, d), nrow = 2)
    ft  <- fisher.test(mat, alternative = "greater")
    
    data.table(
      family    = fam,
      funct     = fn,
      n_dd_fam_funct  = a,   # domain-dark with both family AND function
      n_dd_fam        = n.has.fam,
      n_dd_funct      = a + c,
      odds.ratio      = ft$estimate,
      pval            = ft$p.value,
      frac_funct_in_fam = a / n.has.fam  # of proteins with this family, fraction with this function
    )
  }))
}))

# -----------------------------------------------------------------------------
# Step 3: Multiple testing correction and filtering
# -----------------------------------------------------------------------------

enrich.results[, qval := p.adjust(pval, method = "BH")]
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
    y = "Significantly enriched ECF families",
    title = "A"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
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
    y     = NULL,
    title = "B"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position   = "right",
    axis.text.y       = element_text(size = 9)
  )

# --- Combine
p.combined <- p.scatter + p.dotplot +
  plot_layout(widths = c(1.2, 1))

print(p.combined)

ggsave(
  file.path(this.part.dir, "functional_rescue_enrichment.pdf"),
  p.combined,
  width = 32, height = 14, units = "cm", scale = 1.5
)

# --- Save enrichment table
fwrite(enrich.sig,     file.path(this.part.dir, "ecf_funct_enrichment.csv"))
fwrite(funct.summary,  file.path(this.part.dir, "ecf_funct_summary.csv"))




# ############################
# ### 1. Filter by top EF combinations
# ############################
# 
# # ---- parameters
# TOP_N       <- 25L    # top ECF family combinations by frequency
# MIN_SCORE   <- 0.5    # minimum PHROG annotation score
# MIN_PROTEINS_IN_CELL <- 3  # hide cells with very few proteins
# 
# # ---- Use reprseq.functions (all hits) filtered to top annotations
# # This recovers more functional signal than unique-only
# funct.dt <- reprseq.functions[
#   top.annotation == TRUE & score >= MIN_SCORE,
#   .(qname, funct, class, score)
# ]
# 
# # For each protein keep the highest-scoring function
# # (in case of multiple hits to same function, take best)
# setorder(funct.dt, qname, funct, -score)
# funct.dt <- funct.dt[, .SD[1], by = .(qname, funct)]
# 
# # Best single function per protein for row ordering
# funct.best <- funct.dt[, .SD[which.max(score)], by = qname]
# setnames(funct.best, c("funct", "class", "score"), 
#          c("best.funct", "best.class", "best.score"))
# setkey(funct.best, qname)
# 
# # ---- Select top N family combinations (by protein count)
# top_keys <- famcomp.dd.summary[
#   order(-n_proteins)
# ][1:min(.N, TOP_N)]$fam_key
# 
# # ---- Build heatmap data
# # Join all proteins with those family combinations to all their functions
# hm.funct <- dd.dt[fam_key %chin% top_keys][
#   funct.dt, 
#   on = .(qname),
#   nomatch = 0
# ][
#   , .(N = .N), by = .(fam_key, funct, class)
# ]
# 
# # ---- Filter to top.annotation functions only
# # (already enforced above via funct.dt filter)
# 
# # ---- Keep only functions appearing in >= 2 family combinations
# # avoids singleton columns that add noise
# funct.in.multiple.keys <- hm.funct[
#   , .(n_keys = uniqueN(fam_key)), by = funct
# ][n_keys >= 2]$funct
# 
# hm.funct <- hm.funct[funct %chin% funct.in.multiple.keys]
# 
# # ---- Suppress very small cells
# hm.funct[N < MIN_PROTEINS_IN_CELL, N := NA]
# 
# # ---- Compute functional purity per fam_key
# # (proportion of proteins in dominant function)
# # used for row ordering - most coherent combinations first
# purity.dt <- hm.funct[!is.na(N), .(
#   dom_N     = max(N),
#   total_N   = sum(N),
#   dom_funct = funct[which.max(N)]
# ), by = fam_key]
# purity.dt[, purity := dom_N / total_N]
# setorder(purity.dt, -purity, -total_N)
# 
# # ---- Clean fam_key labels for y axis
# # Replace pipe with " + " and wrap at 3 families per line
# clean.fam.label <- function(key) {
#   fams <- strsplit(key, "|", fixed = TRUE)[[1]]
#   if (length(fams) <= 3) {
#     paste(fams, collapse = " + ")
#   } else {
#     # wrap: first 3 on line 1, rest on line 2
#     paste0(
#       paste(fams[1:3], collapse = " + "),
#       "\n+ ",
#       paste(fams[4:length(fams)], collapse = " + ")
#     )
#   }
# }
# 
# hm.funct[, fam_label := sapply(as.character(fam_key), clean.fam.label)]
# purity.dt[, fam_label := sapply(as.character(fam_key), clean.fam.label)]
# 
# # ---- Factor levels
# # Rows: ordered by purity (most coherent on top)
# row.order <- purity.dt$fam_label
# hm.funct[, fam_label := factor(fam_label, levels = rev(row.order))]
# 
# # Columns: ordered by PHROG class then total abundance
# # so related functions cluster together
# col.order <- hm.funct[!is.na(N), .(
#   tot   = sum(N, na.rm = TRUE),
#   class = class[which.max(N)]
# ), by = funct][order(class, -tot)]$funct
# 
# hm.funct[, funct := factor(funct, levels = col.order)]
# 
# # ---- Add purity annotation for right-side bar
# purity.plot.dt <- purity.dt[, .(
#   fam_label = factor(fam_label, levels = rev(row.order)),
#   purity,
#   total_N
# )]
# 
# # ---- Main heatmap
# p.hm <- ggplot(hm.funct, aes(x = funct, y = fam_label, fill = N)) +
#   geom_tile(color = "white", linewidth = 0.3) +
#   geom_text(
#     data = hm.funct[!is.na(N) & N >= 10],
#     aes(label = N),
#     size = 2.2, color = "white", fontface = "bold"
#   ) +
#   scale_fill_viridis_c(
#     trans   = "sqrt",
#     name    = "Proteins",
#     na.value = "grey93",
#     option  = "D"
#   ) +
#   labs(
#     x = "Predicted function",
#     y = "ECF family combination"
#   ) +
#   theme_bw(base_size = 11) +
#   theme(
#     axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
#     axis.text.y     = element_text(size = 8),
#     axis.line       = element_blank(),
#     legend.position = "right",
#     panel.grid      = element_blank()
#   )
# 
# # ---- Purity sidebar
# p.purity <- ggplot(purity.plot.dt, 
#                    aes(x = purity, y = fam_label)) +
#   geom_col(fill = "steelblue", alpha = 0.7) +
#   geom_text(
#     aes(label = sprintf("%.0f%%", purity * 100)),
#     hjust = -0.1, size = 2.5
#   ) +
#   scale_x_continuous(
#     limits = c(0, 1.15),
#     breaks = c(0, 0.5, 1),
#     labels = c("0", "0.5", "1")
#   ) +
#   labs(x = "Purity", y = NULL) +
#   theme_bw(base_size = 11) +
#   theme(
#     axis.text.y  = element_blank(),
#     axis.ticks.y = element_blank(),
#     panel.grid.minor = element_blank()
#   )
# 
# # ---- Combine
# p.combined <- p.hm + p.purity +
#   plot_layout(widths = c(4, 1))
# 
# print(p.combined)
# 
# ggsave(
#   file.path(this.part.dir, "functional_rescue_heatmap.pdf"),
#   p.combined,
#   width = 28, height = 14, units = "cm", scale = 1.5
# )
# 
# 
# ############################
# ### 2. Filter by top functs
# ############################
# 
# 
# # ---- Revised strategy: start from top functions, not top EF combinations
# 
# TOP_FUNCT    <- 20L   # top functions to show
# MIN_PURITY   <- 0.4   # minimum purity to keep a row
# MIN_PROTEINS <- 10L   # minimum proteins in dominant function per row
# 
# # ---- Step 1: Select top functions by total protein count
# # across ALL domain-dark proteins (not just top combinations)
# funct.totals <- dd.dt[
#   funct != "unknown" & !is.na(funct)
# ][
#   funct.dt, on = .(qname), nomatch = 0
# ][
#   , .(total = .N), by = funct
# ][order(-total)]
# 
# top_functs <- funct.totals[1:min(.N, TOP_FUNCT)]$funct
# 
# # ---- Step 2: For each top function, find the ECF family combinations
# # that are most associated with it
# # Use ALL proteins with that function (not just domain-dark)
# # but flag whether they are domain-dark
# hm.data <- dd.dt[funct != "unknown"][
#   funct.dt, on = .(qname), nomatch = 0
# ][
#   funct %chin% top_functs
# ][
#   , .(
#     n_proteins  = .N,
#     n_dd        = sum(qname %chin% reprseq.dd)  # how many are domain-dark
#   ), by = .(fam_key, funct)
# ]
# 
# # ---- Step 3: For each function, keep only combinations
# # where >= MIN_PROTEINS proteins have that function
# hm.data <- hm.data[n_proteins >= MIN_PROTEINS]
# 
# # ---- Step 4: Filter out single-family rows that appear in many functions
# # (these are the endonuclease-dominated EF0001-type rows)
# # Keep a fam_key-funct pair only if the fam_key is somewhat specific to that function
# fam_key.specificity <- hm.data[, .(
#   n_functs    = uniqueN(funct),
#   top_funct   = funct[which.max(n_proteins)],
#   purity      = max(n_proteins) / sum(n_proteins)
# ), by = fam_key]
# 
# # Remove combinations that appear across too many functions
# # (low specificity = not informative)
# specific.keys <- fam_key.specificity[
#   n_functs <= 3 & purity >= MIN_PURITY
# ]$fam_key
# 
# hm.data <- hm.data[fam_key %chin% specific.keys]
# 
# # ---- Step 5: Limit to top N combinations per function
# # to avoid one function dominating with hundreds of rows
# TOP_KEYS_PER_FUNCT <- 4L
# top.per.funct <- hm.data[
#   order(funct, -n_proteins),
#   .SD[1:min(.N, TOP_KEYS_PER_FUNCT)],
#   by = funct
# ]
# hm.data <- hm.data[
#   paste(fam_key, funct) %chin% paste(top.per.funct$fam_key, top.per.funct$funct)
# ]
# 
# # ---- Step 6: Clean labels and reorder
# hm.data[, fam_label := sapply(as.character(fam_key), clean.fam.label)]
# 
# # Order columns: functions by total protein count (most frequent right)
# funct.order <- funct.totals[funct %chin% hm.data$funct][order(total)]$funct
# hm.data[, funct := factor(funct, levels = funct.order)]
# 
# # Order rows: by number of distinct functions they appear in (most specific first)
# # then by max protein count
# row.order <- hm.data[, .(
#   n_functs  = uniqueN(funct),
#   max_n     = max(n_proteins),
#   n_fams    = lengths(strsplit(fam_key[1], "|", fixed = TRUE))
# ), by = fam_label][order(-n_fams, n_functs, -max_n)]$fam_label
# 
# hm.data[, fam_label := factor(fam_label, levels = row.order)]
# 
# # ---- Flag multi-family combinations for visual emphasis
# hm.data[, n_fams := lengths(strsplit(as.character(fam_key), "|", fixed = TRUE))]
# hm.data[, combo_type := ifelse(n_fams > 1, 
#                                "Multi-family", 
#                                "Single family")]
# 
# # ---- Purity per row
# purity.sidebar <- hm.data[, .(
#   purity  = max(n_proteins) / sum(n_proteins),
#   total_n = sum(n_proteins),
#   n_fams  = n_fams[1]
# ), by = fam_label]
# purity.sidebar[, fam_label := factor(fam_label, levels = row.order)]
# 
# # ---- Plot
# p.hm <- ggplot(hm.data, aes(x = funct, y = fam_label)) +
#   geom_tile(aes(fill = n_proteins), 
#             color = "white", linewidth = 0.3) +
#   geom_text(
#     data = hm.data[n_proteins >= 15],
#     aes(label = n_proteins),
#     size = 2.2, color = "white", fontface = "bold"
#   ) +
#   # Highlight multi-family rows with a left border marker
#   geom_tile(
#     data = hm.data[n_fams > 1],
#     aes(x = funct, y = fam_label),
#     fill = NA, color = "#E15759", linewidth = 0.6
#   ) +
#   scale_fill_viridis_c(
#     trans    = "sqrt",
#     name     = "Proteins",
#     na.value = "grey93",
#     limits   = c(1, 500),
#     oob      = scales::squish
#   ) +
#   facet_grid(combo_type ~ ., scales = "free_y", space = "free_y") +
#   labs(
#     x = "Predicted function",
#     y = "ECF family combination"
#   ) +
#   theme_bw(base_size = 11) +
#   theme(
#     axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
#     axis.text.y      = element_text(size = 8),
#     axis.line        = element_blank(),
#     # panel.grid       = element_blank(),
#     strip.text       = element_text(face = "bold", size = 9),
#     strip.background = element_rect(fill = "grey95"),
#     legend.position  = "right"
#   )
# 
# p.purity <- ggplot(purity.sidebar, aes(x = purity, y = fam_label)) +
#   geom_col(aes(fill = factor(n_fams > 1)), alpha = 0.8) +
#   scale_fill_manual(
#     values = c("FALSE" = "steelblue", "TRUE" = "#E15759"),
#     guide  = "none"
#   ) +
#   geom_text(
#     aes(label = sprintf("%.0f%%", purity * 100)),
#     hjust = -0.1, size = 2.3
#   ) +
#   scale_x_continuous(limits = c(0, 1.2), breaks = c(0, 0.5, 1)) +
#   facet_grid(
#     factor(n_fams > 1, labels = c("Single family", "Multi-family")) ~ .,
#     scales = "free_y", space = "free_y"
#   ) +
#   labs(x = "Purity", y = NULL) +
#   theme_bw(base_size = 11) +
#   theme(
#     axis.text.y      = element_blank(),
#     axis.ticks.y     = element_blank(),
#     panel.grid.minor = element_blank(),
#     strip.text       = element_blank(),
#     strip.background = element_blank()
#   )
# 
# p.combined <- p.hm + p.purity +
#   plot_layout(widths = c(5, 1))
# 
# print(p.combined)
# 
# ggsave(
#   file.path(this.part.dir, "functional_rescue_heatmap_v2.pdf"),
#   p.combined,
#   width = 16, height = 16, units = "cm", scale = 1.5
# )
# 
# 
# 
# # # ---- parameters
# # TOP_N <- 20L
# # 
# # # attach top functional annotation
# # dd.funct <- reprseq.functions.unique[
# #   .(reprseq.fam.comb.dd.all$qname),
# #   .(qname = reprseq, funct, class)
# # ]
# # 
# # dd.dt <- merge(reprseq.fam.comb.dd.all, dd.funct, by = "qname", all.x = TRUE)
# # dd.dt[is.na(funct) | funct == "", funct := "unknown"]
# # dd.dt[is.na(class) | class == "", class := "unknown"]
# # 
# # # ---- summarize family combinations
# # famcomp.dd.summary <- dd.dt[, .(
# #   n_proteins = .N,
# #   qnames = paste(qname, collapse = ";"),
# #   funct_counts = paste(
# #     .SD[, .N, by = funct][order(-N, funct),
# #                           paste0(funct, " (", N, ")")],
# #     collapse = "; "
# #   ),
# #   class_counts = paste(
# #     .SD[, .N, by = class][order(-N, class),
# #                           paste0(class, " (", N, ")")],
# #     collapse = "; "
# #   )
# # ), by = fam_key]
# # 
# # 
# # # recover fams as list-column if you want it
# # famcomp.dd.summary[, fams := strsplit(fam_key, "|", fixed = TRUE)]
# # 
# # # sort by frequency
# # setorder(famcomp.dd.summary, -n_proteins)
# # 
# # # top N
# # famcomp.dd.top <- famcomp.dd.summary[1:min(.N, TOP_N)]
# # 
# # TOP_FUN <- 40
# # 
# # top_keys <- famcomp.dd.summary[order(-n_proteins)][1:min(.N, TOP_N)]$fam_key
# # 
# # hm.funct <- dd.dt[
# #   fam_key %chin% top_keys & funct != "unknown",
# #   .N,
# #   by = .(fam_key, funct)
# # ]
# # 
# # # keep top functions overall for readability
# # top_functs <- hm.funct[, .(tot = sum(N)), by = funct][order(-tot)][1:min(.N, TOP_FUN)]$funct
# # hm.funct <- hm.funct[funct %chin% top_functs]
# # 
# # # order rows by frequency in famcomp.dd.summary
# # hm.funct[, fam_key := factor(fam_key, levels = famcomp.dd.summary[order(-n_proteins)][fam_key %chin% top_keys]$fam_key)]
# # 
# # # order columns by abundance
# # hm.funct[, funct := factor(funct, levels = hm.funct[, .(tot = sum(N)), by = funct][order(tot)]$funct)]
# # 
# # p.hm.funct <- ggplot(hm.funct, aes(x = funct, y = fam_key, fill = N)) +
# #   geom_tile(color = "white", linewidth = 0.2) +
# #   scale_fill_viridis_c(trans = "sqrt", name = "Proteins") +
# #   labs(x = "Predicted function", y = "ECF family combination") +
# #   theme_bw(base_size = 12) +
# #   theme(
# #     axis.text.x = element_text(angle = 45, hjust = 1),
# #     axis.line = element_blank()
# #   )
# # 
# # p.hm.funct
# # 
# # 
# # # 
# # # 
# # # 
# # # reprseq.dd.ef <- fam_by_prot[reprseq.dd, on="qname", nomatch=0]
# # # # count in how many mosaic proteins each family appears
# # # reprseq.dd.fam.list <- reprseq.dd.ef[, .(fams = list(sort(unique(family)))), by = qname]
# # # 
# # # 
# # # 
# # # 
# # # reprseq.dd.dt <- annot.dt[ECF == T & Pfam == F & ECOD == F]
# # # reprseq.dd <- reprseq.dd.dt$qname
# # # 
# # # # Families per protein (distinct)
# # # fam_by_prot <- unique(reprseq.vs.ecfs[, .(qname, family)])
# # # setkey(fam_by_prot, qname)
# # # 
# # # reprseq.dd.ef <- fam_by_prot[reprseq.dd, on="qname", nomatch=0]
# # # # count in how many mosaic proteins each family appears
# # # reprseq.dd.fam.list <- reprseq.dd.ef[, .(fams = list(sort(unique(family)))), by = qname]
# # # 
# # # reprseq.dd.fam.funct <- cbind(reprseq.dd.fam.list, reprseq.functions.unique[.(reprseq.dd.fam.list$qname), .(funct, class)])
# # # reprseq.dd.fam.funct[funct == "", funct := "unknown"]
# # # reprseq.dd.fam.funct[class == "", class := "unknown"]
# # # 
# # # unique(reprseq.dd.fam.funct$fams)
# # 
# # 
# # 
# # 
# # 
# # ########################################################################################
# # ########################################################################################
# # ########################################################################################
# # 
# # # -----------------------------
# # # 1) Define "mosaic modular" proteins
# # #    (simple + tunable thresholds)
# # # -----------------------------
# # # Recommend: require decent coverage + >=2 families + >=2 segments + some separation
# # # Use quantiles so it adapts to your dataset.
# # 
# # cov_min   <- 0.6
# # hits_min  <- 5
# # fam_min   <- 2
# # seg_min_q <- 0.75  # "high separation" = top quartile of segregation_bary
# # seg_min   <- quantile(protein_metrics_essential[!is.na(segregation_bary)]$segregation_bary,
# #                       probs = seg_min_q, na.rm = TRUE)
# # me_max   <- 1.1
# # segN_min  <- 2
# # 
# # mosaic_proteins <- protein_metrics_essential[
# #   coverage_fraction >= cov_min &
# #     n_ecf_hits >= hits_min &
# #     n_families_total >= fam_min &
# #     n_dom_segments >= segN_min &
# #     mean_eff_families <= me_max &
# #     !is.na(segregation_bary) & segregation_bary >= seg_min,
# #   .(qname, funct, class, n_ecf_hits, n_families_total, mean_eff_families,
# #     segregation_bary, n_dom_segments)
# # ]
# # setkey(mosaic_proteins, qname)
# # 
# # 
# # # -----------------------------
# # # 2) Family co-occurrence in mosaic proteins
# # # -----------------------------
# # fam_mosaic <- fam_by_prot[mosaic_proteins, on="qname", nomatch=0]
# # # count in how many mosaic proteins each family appears
# # fam_occ <- fam_mosaic[, .(n_proteins = uniqueN(qname)), by = family]
# # setorder(fam_occ, -n_proteins)
# # 
# # # co-occurrence pairs within each protein
# # fam_list <- fam_mosaic[, .(fams = list(sort(unique(family)))), by = qname]
# # pair_dt <- fam_list[, {
# #   f <- fams[[1]]
# #   if (length(f) < 2) return(NULL)
# #   ij <- t(combn(f, 2))
# #   data.table(f1 = ij[,1], f2 = ij[,2])
# # }, by = qname][, qname := NULL]
# # 
# # # edge weights: number of mosaic proteins where pair co-occurs
# # edges <- pair_dt[, .(w = .N), by = .(f1, f2)]
# # setorder(edges, -w)
# # 
# # # Family–family co-occurrence edges (node = family)
# # edges <- fam_list[
# #   lengths(fams) >= 2,
# #   {
# #     v <- sort(unique(fams[[1]]))
# #     pairs <- t(combn(v, 2))
# #     data.table(from = pairs[, 1], to = pairs[, 2])
# #   },
# #   by = qname
# # ][
# #   , .(weight = .N), by = .(from, to)
# # ][
# #   order(-weight, from, to)
# # ]
# # 
# # # Save for Cytoscape
# # fwrite(edges, file.path(this.part.dir, "ecf_family_cooccurrence_edges.tsv"), sep = "\t")
# # 
# # # Optional: node table
# # nodes <- data.table(id = sort(unique(c(edges$from, edges$to))))
# # ecf.classification.table <- fread(ecf.classification.table.output)
# # ecf.classification.table.ef <- merge(ecf.classification.table, ecf.family.table, by = "ecf")
# # ef_1 <- ecf.classification.table.ef[!duplicated(family)]
# # nodes.out <- merge(nodes, ef_1, by.x = "id", by.y = "family", all.x = T)
# # fwrite(nodes.out, file.path(this.part.dir, "ecf_family_cooccurrence_nodes.tsv"), sep = "\t")
# # 
# # 
# # 
# # 
# # 
# 
# 
# 
# 
# 
# 
# 
