### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(patchwork)))
suppressWarnings(suppressPackageStartupMessages(library(ggnewscale)))
# suppressWarnings(suppressPackageStartupMessages(library(ggh4x)))
# suppressWarnings(suppressPackageStartupMessages(library(Biostrings)))
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
this.part.dir <- file.path(ecf.proj.dir, "part6_rbp-struct")
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

ecod.domains <- fread(ecod.domains.output)
pfam.domains <- fread(pfam.domains.output)

# reprseq lengths
reprseq.lengths <- fread(reprseq.lengths.filename)

# Load functional annotations
reprseq.functions.unique <- fread(reprseq.functions.unique.output)
reprseq.functions <- fread(reprseq.functions.all.output)

# Load reprseq vs ECFs
reprseq.vs.ecfs.all <- fread(reprseq.ecfs.output)

# Load reprseq vs ECOD
reprseq.vs.ecod.all <- fread(reprseq.domains.ecod.output)

# Load reprseq vs ECOD
reprseq.vs.pfam.all <- fread(reprseq.domains.pfam.output)

# Load ECFs vs ECOF
ecfs.vs.ecod <- fread(ecfs.vs.ecod.output)

# Protein arch metrics
reprseq.arch.metrics <- fread(reprseq.arch.metrics.output)


#################################
### Process loaded tables
#################################

setkey(ecf.family.table, ecf)
setkey(ecod.domains, full.id)
setkey(pfam.domains, accession)

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

reprseq.vs.ecod <- reprseq.vs.ecod.all %>% 
  filter(prob >= MIN.PROB.FOR.DOMAIN.DETECT.IN.RHMM,
         eval <= MAX.EVAL.FOR.DOMAIN.DETECT.IN.RHMM, 
         qcov >= MIN.RHMM.COV.FOR.DOMAIN.DETECT.IN.RHMM,
         scov >= MIN.DOMAIN.COV.FOR.DOMAIN.DETECT.IN.RHMM,
         hit.length >= MIN.HIT.LENGTH.FOR.DOMAIN.DETECT.IN.RHMM
  )
reprseq.vs.ecod <- cbind(reprseq.vs.ecod, t.id = ecod.domains[.(reprseq.vs.ecod$sname)]$t_id)
# sort ECOD hits
setorder(reprseq.vs.ecod, qname, t.id, -bitscore)
reprseq.vs.ecod.top <- reprseq.vs.ecod[, .SD[1], by = .(qname, t.id)]

reprseq.vs.pfam <- reprseq.vs.pfam.all %>% 
  filter(prob >= MIN.PROB.FOR.DOMAIN.DETECT.IN.RHMM,
         eval <= MAX.EVAL.FOR.DOMAIN.DETECT.IN.RHMM, 
         qcov >= MIN.RHMM.COV.FOR.DOMAIN.DETECT.IN.RHMM,
         scov >= MIN.DOMAIN.COV.FOR.DOMAIN.DETECT.IN.RHMM,
         hit.length >= MIN.HIT.LENGTH.FOR.DOMAIN.DETECT.IN.RHMM
  )
reprseq.vs.pfam <- cbind(reprseq.vs.pfam, descr = pfam.domains[.(reprseq.vs.pfam$accession)]$description)
setorder(reprseq.vs.pfam, qname, accession, -bitscore)
reprseq.vs.pfam.top <- reprseq.vs.pfam[, .SD[1], by = .(qname, accession)]

setkey(reprseq.arch.metrics, qname)

#############################################################
### 1. Plot EF distributions in three representative 
### RBPs, find consensus and make a figure panel
#############################################################

# =============================================================================
# PARAMETERS
# =============================================================================

focal.proteins <- c("reprseq019028", "reprseq031110", "reprseq057378", 
                    "reprseq012341", "reprseq125186")
prob.threshold <- MIN.PROB.FOR.DOMAIN.DETECT.IN.RHMM

# ECF family colours - fixed palette for known families
family.colors <- c(
  EF0007 = "#4E79A7", EF0013 = "#F28E2B", EF0033 = "#E15759",
  EF0065 = "#76B7B2", EF0079 = "#59A14F", EF0080 = "#EDC948",
  EF0212 = "#B07AA1", EF0261 = "#FF9DA7", EF0635 = "#9C755F",
  EF2655 = "#BAB0AC", EF2175 = "#D4A6C8", EF0072 = "#82C0CC"
)

# ECOD T-level colours - distinct muted palette to avoid confusion with ECFs
# Collect all T-levels observed in focal proteins first
focal.ecod.tids <- reprseq.vs.ecod[
  qname %chin% focal.proteins, unique(t.id)
]
focal.ecod.tids <- focal.ecod.tids[!is.na(focal.ecod.tids)]

ecod.base.colors <- c(
  "#1B7837", "#762A83", "#C51B7D", "#D6604D", "#4D4D4D",
  "#878787", "#D9D9D9", "#A6DBA0", "#C2A5CF", "#F4A582"
)
ecod.colors <- setNames(
  ecod.base.colors[seq_along(focal.ecod.tids)],
  focal.ecod.tids
)

# Y positions: ECF blocks above backbone, ECOD blocks below
# backbone at y = 0, ECF blocks at y > 0, ECOD blocks at y < -2
BACKBONE.Y  <-  0
ECOD.YMIN   <- -2.8
ECOD.YMAX   <- -2.0
ECF.YOFFSET <-  0.2   # gap between backbone and ECF blocks

# =============================================================================
# PROCESS PROTEINS
# =============================================================================

assign.gantt.levels <- function(fc) {
  fc <- copy(fc)  # avoid modifying .SD in place
  setorder(fc, map.start)
  fc[, y.level := 0L]
  level.ends <- c()
  
  for (i in seq_len(nrow(fc))) {
    placed <- FALSE
    for (lev in seq_along(level.ends)) {
      if (fc[i, map.start] > level.ends[lev] + 10) {
        set(fc, i, "y.level", lev)
        level.ends[lev] <- fc[i, map.end]
        placed <- TRUE
        break
      }
    }
    if (!placed) {
      new.lev <- length(level.ends) + 1L
      set(fc, i, "y.level", new.lev)
      level.ends[new.lev] <- fc[i, map.end]
    }
  }
  fc
}

process.protein <- function(protein.id, ecf.data) {
  hits        <- ecf.data[qname == protein.id]
  prot.length <- hits[, max(qlength)]
  
  # Denoised family boundaries
  family.coords <- hits[, .(
    map.start = as.integer(quantile(qstart, 0.25)),
    map.end   = as.integer(quantile(qend,   0.75)),
    n.ecfs    = .N
  ), by = family]
  family.coords[, protein := protein.id]
  
  # Per-position family assignment using denoised boundaries
  position.counts <- family.coords[, .(
    position = seq(map.start, map.end)
  ), by = .(family, n.ecfs)][,
                             .(weighted.count = sum(n.ecfs)), by = .(family, position)
  ]
  
  position.assignment <- position.counts[
    position.counts[, .I[which.max(weighted.count)], by = position]$V1,
    .(position, family)
  ]
  position.assignment[, protein := protein.id]
  
  # Structure blocks for PyMOL
  setorder(position.assignment, position)
  position.assignment[, block := rleid(family)]
  structure.blocks <- position.assignment[, .(
    protein   = protein.id,
    family    = family[1],
    res.start = min(position),
    res.end   = max(position)
  ), by = block][, block := NULL]
  
  list(
    family.coords       = family.coords,
    position.assignment = position.assignment,
    structure.blocks    = structure.blocks,
    prot.length         = prot.length
  )
}

# =============================================================================
# PROCESS ECF DATA (unchanged from existing code)
# =============================================================================

results <- lapply(focal.proteins, process.protein,
                  ecf.data = reprseq.vs.ecfs)
names(results) <- focal.proteins

# =============================================================================
# PROCESS ECOD DATA
# =============================================================================

process.ecod <- function(protein.id, ecod.data) {
  hits <- ecod.data[qname == protein.id & !is.na(t.id)]
  if (nrow(hits) == 0) return(data.table())
  
  # One block per T-level: use denoised boundaries same as ECFs
  ecod.coords <- hits[, .(
    map.start = as.integer(quantile(qstart, 0.25)),
    map.end   = as.integer(quantile(qend,   0.75)),
    n.hits    = .N
  ), by = t.id]
  ecod.coords[, protein := protein.id]
  ecod.coords
}

all.ecod.coords <- rbindlist(
  lapply(focal.proteins, process.ecod, ecod.data = reprseq.vs.ecod.top)
)
ecod.t.names <- sapply(1:nrow(all.ecod.coords), function(k) ecod.domains[t_id == all.ecod.coords$t.id[k]]$t_name[1])
all.ecod.coords$t.name <- ecod.t.names

ecod.legend.labels <- all.ecod.coords[, .(
  label = paste0(t.id[1], " (", t.name[1], ")")
), by = t.id][, setNames(label, t.id)]

# =============================================================================
# PROCESS PFAM DATA
# =============================================================================

process.pfam <- function(protein.id, pfam.data) {
  hits <- pfam.data[qname == protein.id & !is.na(accession)]
  if (nrow(hits) == 0) return(data.table())
  
  # One block per T-level: use denoised boundaries same as ECFs
  pfam.coords <- hits[, .(
    map.start = as.integer(quantile(qstart, 0.25)),
    map.end   = as.integer(quantile(qend,   0.75)),
    n.hits    = .N
  ), by = accession]
  pfam.coords[, protein := protein.id]
  pfam.coords
}

all.pfam.coords <- rbindlist(
  lapply(focal.proteins, process.pfam, pfam.data = reprseq.vs.pfam.top)
)
all.pfam.coords$description <- pfam.domains[.(all.pfam.coords$accession)]$description

pfam.legend.labels <- all.pfam.coords[, .(
  label = paste0(accession[1], " (", description[1], ")")
), by = accession][, setNames(label, accession)]


# =============================================================================
# PREPARE PLOTTING DATA
# =============================================================================

plotting.proteins <- focal.proteins[focal.proteins != "reprseq057378"]   # all three

all.fc <- rbindlist(lapply(results, `[[`, "family.coords"))
all.pa <- rbindlist(lapply(results, `[[`, "position.assignment"))
all.fc <- all.fc[protein %chin% plotting.proteins]
all.pa <- all.pa[protein %chin% plotting.proteins]

# Gantt stacking for ECF blocks - start above backbone
all.fc <- all.fc[, assign.gantt.levels(.SD), by = protein]
all.fc[, y.level := y.level + ECF.YOFFSET]   # shift up from backbone
all.fc[, y.half  := 0.15 + 0.35 * sqrt(n.ecfs / max(n.ecfs)), by = protein]

# Backbone segments from position assignment
setorder(all.pa, protein, position)
all.pa[, block := rleid(family), by = protein]
backbone.segments <- all.pa[, .(
  x      = min(position),
  xend   = max(position),
  family = family[1]
), by = .(protein, block)]

# Protein lengths
prot.lengths <- data.table(
  protein     = focal.proteins,
  prot.length = sapply(results, `[[`, "prot.length")
)[protein %chin% plotting.proteins]

all.ecod.coords <- all.ecod.coords[protein %chin% plotting.proteins]

# =============================================================================
# COLOUR VALIDATION
# Check which ECF families appear in plot but are missing from palette
# =============================================================================

missing.ecf <- setdiff(all.fc$family, names(family.colors))
if (length(missing.ecf) > 0) {
  warning("ECF families missing from palette: ", 
          paste(missing.ecf, collapse = ", "))
}

missing.ecod <- setdiff(all.ecod.coords$t.id, names(ecod.colors))
if (length(missing.ecod) > 0) {
  warning("ECOD T-levels missing from palette: ", 
          paste(missing.ecod, collapse = ", "))
}

# =============================================================================
# PLOT
# =============================================================================

# Add track labels directly into prot.lengths so they facet correctly
track.labels <- prot.lengths[, .(
  protein,
  label = "consensus",
  y     = BACKBONE.Y
)]
track.labels <- rbind(
  track.labels,
  prot.lengths[, .(
    protein,
    label = "ECOD",
    y     = (ECOD.YMIN + ECOD.YMAX) / 2
  )]
)

# Compute the actual y range across all proteins dynamically
y.max <- all.fc[, max(y.level + y.half)] + 0.3
y.min <- ECOD.YMIN - 0.4

p <- ggplot() +
  
  # --- Grey full-length backbone underlay ---
  geom_segment(
    data = prot.lengths,
    aes(x = 1, xend = prot.length, 
        y = BACKBONE.Y, yend = BACKBONE.Y),
    linewidth = 3, color = "grey80"
  ) +
  
  # --- Coloured consensus backbone (ECF position assignment) ---
  geom_segment(
    data = backbone.segments,
    aes(x = x, xend = xend,
        y = BACKBONE.Y, yend = BACKBONE.Y,
        color = family),
    linewidth = 3
  ) +
  scale_color_manual(
    values = family.colors, 
    guide  = "none"
  ) +
  new_scale_color() +
  
  # --- ECF family blocks (above backbone) ---
  geom_rect(
    data = all.fc,
    aes(xmin = map.start, xmax = map.end,
        ymin  = y.level - y.half,
        ymax  = y.level + y.half,
        fill  = family),
    alpha = 0.85, color = "white", linewidth = 0.3
  ) +
  
  # --- ECF labels on wide boxes ---
  geom_text(
    data = all.fc[map.end - map.start > 40],
    aes(x     = (map.start + map.end) / 2,
        y     = y.level,
        label = sprintf("%s (n=%d)", family, n.ecfs)),
    size = 2.5, fontface = "bold", color = "black"
  ) +
  
  scale_fill_manual(
    values = family.colors,
    name   = "ECF family",
    drop   = FALSE
  ) +
  new_scale_fill() +
  
  # --- ECOD blocks (below backbone) ---
  geom_rect(
    data = all.ecod.coords,
    aes(xmin = map.start, xmax = map.end,
        ymin  = ECOD.YMIN,
        ymax  = ECOD.YMAX,
        fill  = t.id),
    alpha = 0.85, color = "white", linewidth = 0.3
  ) +
  
  # --- ECOD labels on wide boxes ---
  geom_text(
    data = all.ecod.coords[map.end - map.start > 80],
    aes(x     = (map.start + map.end) / 2,
        y     = (ECOD.YMIN + ECOD.YMAX) / 2,
        label = t.id),
    size = 2.2, color = "white", fontface = "bold"
  ) +
  
  scale_fill_manual(
    values = ecod.colors,
    name   = "ECOD T-level",
    labels = ecod.legend.labels,
    drop   = TRUE   # only show T-levels present in data
  ) +
  
  # --- Track labels on y axis ---
  # Add "ECF" and "ECOD" text annotations per facet
  # Then in the ggplot chain, replace the two geom_text calls with:
  geom_text(
    data = track.labels,
    aes(x = 0, y = y, label = label),
    hjust = 1.05, size = 2.5, color = "grey50"
  ) +
  # --- Scales and facets ---
  scale_x_continuous(
    breaks = seq(0, 1200, by = 200),
    expand = c(0.08, 0)   # extra left margin for track labels
  ) +
  scale_y_continuous(
    limits = c(y.min, y.max),
    expand = c(0, 0)
  ) +
  facet_wrap(~ protein, nrow = 1, scales = "free_x") +
  labs(x = "Position (aa)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.text         = element_text(face = "bold", size = 10),
    legend.position    = "bottom",
    legend.box         = "vertical"   # stack ECF and ECOD legends
  )

print(p)

fig.outfile <- file.path(ecf.figures.dir, "figure6A.pdf")
ggsave(
  fig.outfile, p,
  width = 22, height = 9, units = "cm", scale = 1.5
)


