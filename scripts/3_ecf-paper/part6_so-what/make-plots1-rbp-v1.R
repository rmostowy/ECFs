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
this.part.dir <- file.path(ecf.proj.dir, "part6_so-what")
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

# reprseq lengths
reprseq.lengths <- fread(reprseq.lengths.filename)

# Load functional annotations
reprseq.functions.unique <- fread(reprseq.functions.unique.output)
reprseq.functions <- fread(reprseq.functions.all.output)

# Load reprseq vs ECFs
reprseq.vs.ecfs.all <- fread(reprseq.ecfs.output)

# Load reprseq vs ECOD
reprseq.vs.ecod.all <- fread(reprseq.domains.ecod.output)

# Load ECFs vs ECOF
ecfs.vs.ecod <- fread(ecfs.vs.ecod.output)

# Protein arch metrics
reprseq.arch.metrics <- fread(reprseq.arch.metrics.output)


#################################
### Process loaded tables
#################################

setkey(ecf.family.table, ecf)

setkey(ecod.domains, full.id)

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

family.colors <- c(
  EF0007 = "#4E79A7", EF0013 = "#F28E2B", EF0033 = "#E15759",
  EF0065 = "#76B7B2", EF0079 = "#59A14F", EF0080 = "#EDC948",
  EF0212 = "#B07AA1", EF0261 = "#FF9DA7", EF0635 = "#9C755F",
  EF2655 = "#BAB0AC", EF2175 = "#D4A6C8", EF0072 = "#82C0CC"
)

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

results <- lapply(focal.proteins, process.protein, 
                  ecf.data = reprseq.vs.ecfs)
names(results) <- focal.proteins

# =============================================================================
# PREPARE PLOTTING DATA
# =============================================================================

plotting.proteins <- focal.proteins[c(1,2,3)]

all.fc <- rbindlist(lapply(results, `[[`, "family.coords"))
all.pa <- rbindlist(lapply(results, `[[`, "position.assignment"))
all.fc <- all.fc[protein %chin% plotting.proteins]
all.pa <- all.pa[protein %chin% plotting.proteins]

# Gantt y-level assignment and height scaling per protein
all.fc <- all.fc[, assign.gantt.levels(.SD), by = protein]
all.fc[, y.half := 0.15 + 0.35 * sqrt(n.ecfs / max(n.ecfs)), by = protein]

# Backbone segments from position assignment
setorder(all.pa, protein, position)
all.pa[, block := rleid(family), by = protein]
backbone.segments <- all.pa[, .(
  x      = min(position),
  xend   = max(position),
  family = family[1]
), by = .(protein, block)]

# Grey underlay per protein
prot.lengths <- data.table(
  protein     = focal.proteins,
  prot.length = sapply(results, `[[`, "prot.length")
)
prot.lengths <- prot.lengths[protein %chin% plotting.proteins]


# =============================================================================
# PLOT
# =============================================================================

p <- ggplot() +
  # Grey full-length backbone underlay
  geom_segment(
    data = prot.lengths,
    aes(x = 1, xend = prot.length, y = -0.8, yend = -0.8),
    linewidth = 3, color = "grey80"
  ) +
  # Coloured consensus backbone
  geom_segment(
    data = backbone.segments,
    aes(x = x, xend = xend, y = -0.8, yend = -0.8, color = family),
    linewidth = 3
  ) +
  scale_color_manual(values = family.colors, guide = "none") +
  new_scale_color() +
  # Family blocks
  geom_rect(
    data = all.fc,
    aes(xmin = map.start, xmax = map.end,
        ymin  = y.level - y.half,
        ymax  = y.level + y.half,
        fill  = family),
    alpha = 0.85, color = "white", linewidth = 0.3
  ) +
  # Labels on wide enough boxes only
  geom_text(
    data = all.fc[map.end - map.start > 80],
    aes(x     = (map.start + map.end) / 2,
        y     = y.level,
        label = sprintf("%s (n=%d)", family, n.ecfs)),
    size = 2.5, fontface = "bold", color = "black"
  ) +
  scale_fill_manual(values = family.colors, name = "ECF family", drop = FALSE) +
  scale_x_continuous(breaks = seq(0, 1200, by = 200), expand = c(0.01, 0)) +
  scale_y_continuous(expand = c(0.3, 0.3)) +
  facet_wrap(~ protein, nrow = 1, scales = "free_x") +
  labs(x = "Position (aa)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.text         = element_text(face = "bold", size = 10),
    legend.position    = "bottom"
  )

ggsave(file.path(this.part.dir, "ecf_rbp_panels.pdf"), p,
       width = 22, height = 7, units = "cm", scale = 1.5)

# =============================================================================
# PYMOL TABLES
# =============================================================================

# Helper to convert hex to PyMOL RGB string
hex.to.pymol.rgb <- function(hex) {
  rgb <- col2rgb(hex) / 255
  sprintf("[%.4f, %.4f, %.4f]", rgb[1], rgb[2], rgb[3])
}

all.structure.blocks <- rbindlist(lapply(results, `[[`, "structure.blocks"))
all.structure.blocks[, pymol.selection := paste0("resi ", res.start, "-", res.end)]

pymol.table <- all.structure.blocks[, .(
  selection  = paste(pymol.selection, collapse = " or "),
  hex.color  = family.colors[family[1]],
  pymol.rgb  = hex.to.pymol.rgb(family.colors[family[1]])
), by = .(protein, family)]

fwrite(pymol.table, file.path(this.part.dir, "ecf_pymol_selections.csv"))
fwrite(pymol.table,          file.path(this.part.dir, "ecf_pymol_selections.csv"))


# # --- Parameters ---
# focal.proteins <- c("reprseq019028", "reprseq031110", "reprseq057378")
# 
# # --- Consistent family colour palette across all proteins ---
# # Built from all families observed across all focal proteins
# all.families <- reprseq.vs.ecfs[
#   qname %in% focal.proteins, 
#   unique(family)
# ]
# 
# # Define colours - extend as needed
# base.colors <- c(
#   "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
#   "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
#   "#D37295", "#FABFD2", "#8CD17D", "#B6992D", "#499894"
# )
# family.colors <- setNames(
#   base.colors[seq_along(all.families)], 
#   sort(all.families)
# )
# 
# # --- Function to process one protein ---
# process.protein <- function(protein.id, ecf.data) {
#   
#   hits <- ecf.data[qname == protein.id]
#   prot.length <- hits[, max(qlength)]
#   
#   # Denoised family boundaries
#   family.coords <- hits[, .(
#     map.start = as.integer(quantile(qstart, 0.25)),
#     map.end   = as.integer(quantile(qend,   0.75)),
#     n.ecfs    = .N
#   ), by = family]
#   
#   # Per-position family assignment using denoised boundaries
#   position.counts <- family.coords[, {
#     data.table(position = seq(map.start, map.end))
#   }, by = .(family, n.ecfs)
#   ][, .(weighted.count = sum(n.ecfs)), by = .(family, position)]
#   
#   position.assignment <- position.counts[
#     position.counts[, .I[which.max(weighted.count)], by = position]$V1
#   ][, .(position, family, weighted.count)]
#   
#   # Collapse to contiguous blocks for structure coloring
#   setorder(position.assignment, position)
#   position.assignment[, block := rleid(family)]
#   
#   structure.blocks <- position.assignment[, .(
#     protein   = protein.id,
#     family    = family[1],
#     res.start = min(position),
#     res.end   = max(position)
#   ), by = block][, block := NULL]
#   
#   list(
#     family.coords     = family.coords[, protein := protein.id],
#     structure.blocks  = structure.blocks,
#     prot.length       = prot.length,
#     position.assignment = position.assignment[, protein := protein.id]
#   )
# }
# 
# # --- Process all three proteins ---
# results <- lapply(focal.proteins, process.protein, 
#                   ecf.data = reprseq.vs.ecfs
#                   )
# names(results) <- focal.proteins
# 
# # --- Function to build one panel ---
# make.panel <- function(res, protein.id, fam.colors) {
#   
#   fc          <- res$family.coords
#   pa          <- res$position.assignment
#   prot.length <- res$prot.length
#   
#   # Y stacking: order families by start position
#   # Width of rect scales with n.ecfs (normalised to max)
#   fc[, n.ecfs.norm := n.ecfs / max(n.ecfs)]
#   fc[, y.level := frankv(map.start, ties.method = "first")]
#   
#   # Half-height scales with normalised ECF count
#   # so more-supported families get taller boxes
#   fc[, y.half := 0.45 * n.ecfs.norm]
#   
#   # Consensus coloring on the query sequence
#   # Colour the backbone by position assignment
#   # Build segments of contiguous same-family positions
#   setorder(pa, position)
#   pa[, block := rleid(family)]
#   
#   backbone.segments <- pa[, .(
#     x    = min(position),
#     xend = max(position),
#     family = family[1]
#   ), by = block]
#   
#   p <- ggplot(fc) +
#     # Coloured consensus backbone
#     geom_segment(
#       data = backbone.segments,
#       aes(x = x, xend = xend, 
#           y = 0, yend = 0, 
#           color = family),
#       linewidth = 4,
#       show.legend = FALSE
#     ) +
#     scale_color_manual(values = fam.colors) +
#     # Family blocks - width scaled by n.ecfs
#     ggnewscale::new_scale_color() +
#     geom_rect(aes(
#       xmin = map.start,
#       xmax = map.end,
#       ymin = y.level - y.half,
#       ymax = y.level + y.half,
#       fill = family
#     ), alpha = 0.85, color = "white", linewidth = 0.3) +
#     # Family labels inside rects
#     geom_text(aes(
#       x     = (map.start + map.end) / 2,
#       y     = y.level,
#       label = paste0(family, "\n(n=", n.ecfs, ")")
#     ), size = 2.5, fontface = "bold", color = "white") +
#     scale_fill_manual(values = fam.colors, drop = FALSE) +
#     scale_x_continuous(
#       limits = c(0, prot.length + 10),
#       breaks = seq(0, prot.length, by = 200),
#       expand = c(0.01, 0)
#     ) +
#     labs(
#       title = protein.id,
#       x     = "Position (aa)",
#       y     = NULL,
#       fill  = "ECF family"
#     ) +
#     theme_minimal(base_size = 11) +
#     theme(
#       axis.text.y        = element_blank(),
#       axis.ticks.y       = element_blank(),
#       panel.grid.major.y = element_blank(),
#       panel.grid.minor.y = element_blank(),
#       plot.title         = element_text(face = "bold", size = 10),
#       legend.position    = "none"  # shared legend added below
#     )
#   
#   p
# }
# 
# # --- Build panels ---
# panels <- mapply(
#   make.panel,
#   res        = results,
#   protein.id = focal.proteins,
#   MoreArgs   = list(fam.colors = family.colors),
#   SIMPLIFY   = FALSE
# )
# 
# # --- Shared legend ---
# # Build a dummy plot just for the legend
# legend.data <- data.table(
#   family = names(family.colors),
#   x = 1, y = seq_along(family.colors)
# )
# 
# legend.plot <- ggplot(legend.data, aes(x, y, fill = family)) +
#   geom_tile() +
#   scale_fill_manual(values = family.colors) +
#   theme_void() +
#   theme(legend.position = "right",
#         legend.title = element_text(face = "bold"))
# 
# shared.legend <- cowplot::get_legend(legend.plot)
# 
# # --- Combine panels ---
# combined <- (panels[[1]] | panels[[2]] | panels[[3]]) +
#   plot_layout(guides = "collect") &
#   theme(legend.position = "bottom")
# 
# print(combined)
# 
# # Save figure
# ggsave("ecf_rbp_panels.pdf", combined, 
#        width = 18, height = 6, units = "cm", scale = 1.5)
# 
# # --- Save PyMOL coordinate table ---
# all.structure.blocks <- rbindlist(
#   lapply(results, `[[`, "structure.blocks")
# )
# 
# # Add PyMOL selection string per row
# all.structure.blocks[, pymol.selection := paste0("resi ", res.start, "-", res.end)]
# 
# # Collapsed per protein + family for direct PyMOL input
# pymol.table <- all.structure.blocks[, .(
#   selection = paste(pymol.selection, collapse = " or ")
# ), by = .(protein, family)]
# 
# # Add hex colour for convenience
# pymol.table[, hex.color := family.colors[family]]
# 
# fwrite(all.structure.blocks, "ecf_structure_blocks.csv")
# fwrite(pymol.table,          "ecf_pymol_selections.csv")
# 
# cat("\n=== PyMOL selections ===\n")
# print(pymol.table)
# 
# 
# ################################
# 
# 
# make.panel <- function(res, protein.id, fam.colors) {
#   
#   fc          <- res$family.coords
#   pa          <- res$position.assignment
#   prot.length <- res$prot.length
#   
#   # --- Layout logic ---
#   # Assign y levels using a greedy interval packing algorithm
#   # so non-overlapping families share y levels (Gantt-style)
#   setorder(fc, map.start)
#   fc[, y.level := 0L]
#   
#   # Track end position of last family at each y level
#   level.ends <- c()
#   
#   for (i in seq_len(nrow(fc))) {
#     placed <- FALSE
#     for (lev in seq_along(level.ends)) {
#       if (fc[i, map.start] > level.ends[lev] + 10) {
#         fc[i, y.level := lev]
#         level.ends[lev] <- fc[i, map.end]
#         placed <- TRUE
#         break
#       }
#     }
#     if (!placed) {
#       new.lev <- length(level.ends) + 1
#       fc[i, y.level := new.lev]
#       level.ends[new.lev] <- fc[i, map.end]
#     }
#   }
#   
#   # Scale half-height more aggressively
#   fc[, n.ecfs.norm := sqrt(n.ecfs / max(n.ecfs))]  # sqrt for better visual range
#   fc[, y.half := 0.15 + 0.35 * n.ecfs.norm]        # min height 0.15, max 0.50
#   
#   # --- Backbone segments from position assignment ---
#   setorder(pa, position)
#   pa[, block := rleid(family)]
#   backbone.segments <- pa[, .(
#     x      = min(position),
#     xend   = max(position),
#     family = family[1]
#   ), by = block]
#   
#   # Y position for backbone: just below all boxes
#   backbone.y <- -0.8
#   
#   # --- Build plot ---
#   p <- ggplot() +
#     # Grey full-length backbone
#     annotate("segment",
#              x = 1, xend = prot.length,
#              y = backbone.y, yend = backbone.y,
#              linewidth = 3, color = "grey80") +
#     # Coloured consensus backbone below boxes
#     geom_segment(
#       data = backbone.segments,
#       aes(x = x, xend = xend,
#           y = backbone.y, yend = backbone.y,
#           color = family),
#       linewidth = 3
#     ) +
#     scale_color_manual(values = fam.colors, guide = "none") +
#     # Family blocks
#     ggnewscale::new_scale_color() +
#     geom_rect(
#       data = fc,
#       aes(xmin = map.start, xmax = map.end,
#           ymin = y.level - y.half,
#           ymax = y.level + y.half,
#           fill = family),
#       alpha = 0.85, color = "white", linewidth = 0.3
#     ) +
#     # Labels only if box is wide enough
#     geom_text(
#       data = fc[(map.end - map.start) > 80],
#       aes(x     = (map.start + map.end) / 2,
#           y     = y.level,
#           label = paste0(family, "\nn=", n.ecfs)),
#       size = 2.5, fontface = "bold", color = "white"
#     ) +
#     scale_fill_manual(values = fam.colors, drop = FALSE) +
#     scale_x_continuous(
#       limits = c(0, prot.length + 10),
#       breaks = seq(0, prot.length, by = 200),
#       expand = c(0.01, 0)
#     ) +
#     scale_y_continuous(expand = c(0.3, 0.3)) +
#     labs(
#       title = protein.id,
#       x     = "Position (aa)",
#       y     = NULL,
#       fill  = "ECF family"
#     ) +
#     theme_minimal(base_size = 11) +
#     theme(
#       axis.text.y        = element_blank(),
#       axis.ticks.y       = element_blank(),
#       panel.grid.major.y = element_blank(),
#       panel.grid.minor.y = element_blank(),
#       plot.title         = element_text(face = "bold", size = 10),
#       legend.position    = "none"
#     )
#   
#   p
# }
# 
# # --- Rebuild panels ---
# panels <- mapply(
#   make.panel,
#   res        = results,
#   protein.id = focal.proteins,
#   MoreArgs   = list(fam.colors = family.colors),
#   SIMPLIFY   = FALSE
# )
# 
# # --- Clean shared legend from first panel ---
# legend.panel <- panels[[1]] + 
#   theme(legend.position = "right")
# 
# shared.legend <- cowplot::get_legend(
#   ggplot(results[[1]]$family.coords, 
#          aes(x = map.start, y = 1, fill = family)) +
#     geom_tile() +
#     scale_fill_manual(values = family.colors, name = "ECF family") +
#     theme_void() +
#     theme(legend.position = "right",
#           legend.title = element_text(face = "bold", size = 10),
#           legend.text  = element_text(size = 9))
# )
# 
# # --- Combine ---
# combined <- cowplot::plot_grid(
#   cowplot::plot_grid(
#     panels[[1]], panels[[2]], panels[[3]],
#     nrow    = 1,
#     labels  = c("A", "B", "C"),
#     align   = "h",
#     axis    = "tb"
#   ),
#   shared.legend,
#   nrow         = 1,
#   rel_widths   = c(1, 0.12)
# )
# 
# print(combined)
# 
# ggsave("ecf_rbp_panels_v2.pdf", combined,
#        width = 22, height = 7, units = "cm", scale = 1.5)
# 
# 
# 
# ####################################
# 
# 
# # Combine all family.coords into one data.table
# all.fc <- rbindlist(lapply(results, `[[`, "family.coords"))
# all.pa <- rbindlist(lapply(results, `[[`, "position.assignment"))
# 
# # Gantt y-level assignment per protein
# all.fc[, y.level := 0L]
# all.fc[, y.half  := 0.15 + 0.35 * sqrt(n.ecfs / max(n.ecfs)), by = protein]
# 
# for (prot in focal.proteins) {
#   idx        <- which(all.fc$protein == prot)
#   sub        <- all.fc[idx][order(map.start)]
#   level.ends <- c()
#   
#   for (i in seq_len(nrow(sub))) {
#     placed <- FALSE
#     for (lev in seq_along(level.ends)) {
#       if (sub[i, map.start] > level.ends[lev] + 10) {
#         sub[i, y.level := lev]
#         level.ends[lev] <- sub[i, map.end]
#         placed <- TRUE
#         break
#       }
#     }
#     if (!placed) {
#       new.lev <- length(level.ends) + 1
#       sub[i, y.level := new.lev]
#       level.ends[new.lev] <- sub[i, map.end]
#     }
#   }
#   all.fc[idx] <- sub[order(map.start)]
# }
# 
# # Backbone segments
# setorder(all.pa, protein, position)
# all.pa[, block := rleid(family), by = protein]
# backbone.segments <- all.pa[, .(
#   x      = min(position),
#   xend   = max(position),
#   family = family[1]
# ), by = .(protein, block)]
# 
# # Protein lengths for grey underlay
# prot.lengths <- data.table(
#   protein    = focal.proteins,
#   prot.length = sapply(results, `[[`, "prot.length")
# )
# 
# p <- ggplot() +
#   # Grey full-length underlay
#   geom_segment(
#     data = prot.lengths,
#     aes(x = 1, xend = prot.length, y = -0.8, yend = -0.8),
#     linewidth = 3, color = "grey80"
#   ) +
#   # Coloured consensus backbone
#   geom_segment(
#     data = backbone.segments,
#     aes(x = x, xend = xend,
#         y = -0.8, yend = -0.8,
#         color = family),
#     linewidth = 3
#   ) +
#   scale_color_manual(values = family.colors, guide = "none") +
#   ggnewscale::new_scale_color() +
#   # Family blocks
#   geom_rect(
#     data = all.fc,
#     aes(xmin = map.start, xmax = map.end,
#         ymin  = y.level - y.half,
#         ymax  = y.level + y.half,
#         fill  = family),
#     alpha = 0.85, color = "white", linewidth = 0.3
#   ) +
#   # Labels on wide enough boxes only
#   geom_text(
#     data = all.fc[(map.end - map.start) > 80],
#     aes(x     = (map.start + map.end) / 2,
#         y     = y.level,
#         label = sprintf("%s (n=%d)", family, n.ecfs)),
#     size = 2.5, fontface = "bold", color = "black"
#   ) +
#   scale_fill_manual(values = family.colors, 
#                     name   = "ECF family",
#                     drop   = FALSE) +
#   scale_x_continuous(
#     breaks = seq(0, 1200, by = 200),
#     expand = c(0.01, 0)
#   ) +
#   scale_y_continuous(expand = c(0.3, 0.3)) +
#   facet_wrap(~ protein, nrow = 1, scales = "free_x") +
#   labs(x = "Position (aa)", y = NULL) +
#   theme_minimal(base_size = 11) +
#   theme(
#     axis.text.y        = element_blank(),
#     axis.ticks.y       = element_blank(),
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor.y = element_blank(),
#     strip.text         = element_text(face = "bold", size = 10),
#     legend.position    = "bottom"
#   )
# 
# print(p)
# 
# ggsave("ecf_rbp_panels_v3.pdf", p,
#        width = 22, height = 7, units = "cm", scale = 1.5)
