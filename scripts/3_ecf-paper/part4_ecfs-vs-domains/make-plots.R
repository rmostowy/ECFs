### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(patchwork)))
# suppressWarnings(suppressPackageStartupMessages(library(ggh4x)))
suppressWarnings(suppressPackageStartupMessages(library(ComplexUpset)))
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
this.part.dir <- file.path(ecf.proj.dir, "part4_ecfs-vs-domains")
make.dir(this.part.dir)

#################################
### Load all data:
### 1. ECF tables
### 2. domain tables
### 3. ecf vs domain tables
### 4. reprseq vs ecfs
### 5. 
### 6. 
### 7. 
#################################

# ECF tables
ecf.table <- fread(ecf.table.output)
ecf.family.table <- fread(ecf.family.table.output)

# domain tables
ecod.domains <- fread(ecod.domains.output)
pfam.domains <- fread(pfam.domains.output)

# ecf vs domains tables
ecfs.vs.ecod <- fread(ecfs.vs.ecod.output)
ecfs.vs.pfam <- fread(ecfs.vs.pfam.output)
ecfs.vs.ecod.hhsenser <- fread(ecfs.vs.ecod.hhsenser.output)
ecfs.vs.pfam.hhsenser <- fread(ecfs.vs.pfam.hhsenser.output)

# reprseq vs ecfs
reprseq.vs.ecfs.all <- fread(reprseq.ecfs.output)

# functional annotations
reprseq.functions.unique <- fread(reprseq.functions.unique.output)
reprseq.functions.all <- fread(reprseq.functions.all.output)

#################################
### Process loaded tables
#################################

all.ecfs <- unique(ecf.table$ecf.name)
Ne <- length(all.ecfs)

ecfs.vs.ecod.conf <- ecfs.vs.ecod %>%
  filter(eval <= MAX.EVALUE.FOR.ECF.DOMAIN.CONFIDENT & prob >= MIN.PROB.FOR.ECF.DOMAIN.CONFIDENT)
ecfs.vs.pfam.conf <- ecfs.vs.pfam %>%
  filter(eval <= MAX.EVALUE.FOR.ECF.DOMAIN.CONFIDENT & prob >= MIN.PROB.FOR.ECF.DOMAIN.CONFIDENT)

reprseq.vs.ecfs <- reprseq.vs.ecfs.all %>% 
  filter(prob >= MIN.PROB.FOR.DOMAIN.DETECT.IN.RHMM,
         eval <= MAX.EVAL.FOR.DOMAIN.DETECT.IN.RHMM, 
         qcov >= MIN.RHMM.COV.FOR.DOMAIN.DETECT.IN.RHMM,
         scov >= MIN.DOMAIN.COV.FOR.DOMAIN.DETECT.IN.RHMM,
         hit.length >= MIN.HIT.LENGTH.FOR.DOMAIN.DETECT.IN.RHMM
  )

#################################
### 1. Compare how many ECFs hit
### domain databases
#################################

ecf.domain.coverage <- data.table(ecf = all.ecfs)

# ECOD
ecf.domain.coverage[, ecod.eval0p1 := all.ecfs %in% ecfs.vs.ecod$qname]
ecf.domain.coverage[, ecod.eval0p001 := all.ecfs %in% ecfs.vs.ecod.conf$qname]
ecf.domain.coverage[, ecod.eval10 := all.ecfs %in% ecfs.vs.ecod.hhsenser$qname]
# PFAM
ecf.domain.coverage[, pfam.eval0p1 := all.ecfs %in% ecfs.vs.pfam$qname]
ecf.domain.coverage[, pfam.eval0p001 := all.ecfs %in% ecfs.vs.pfam.conf$qname]
ecf.domain.coverage[, pfam.eval10 := all.ecfs %in% ecfs.vs.pfam.hhsenser$qname]

# dt <- copy(ecf.domain.coverage)

# Helper: assign the 4-category label for a given pair of booleans
assign_cat <- function(has_ecod, has_pfam) {
  fifelse(
    has_ecod & has_pfam, "Both (ECOD + Pfam)",
    fifelse(
      has_ecod & !has_pfam, "ECOD only",
      fifelse(!has_ecod & has_pfam, "Pfam only", "Neither")
    )
  )
}

# Build a long table: one row per ECF × threshold
# (thresholds ordered: strict -> lenient)
long_list <- list(
  data.table(ecf = ecf.domain.coverage$ecf,
             threshold = "1e-3",
             has_ecod = ecf.domain.coverage$ecod.eval0p001,
             has_pfam = ecf.domain.coverage$pfam.eval0p001),
  data.table(ecf = ecf.domain.coverage$ecf,
             threshold = "1e-1",
             has_ecod = ecf.domain.coverage$ecod.eval0p1,
             has_pfam = ecf.domain.coverage$pfam.eval0p1),
  data.table(ecf = ecf.domain.coverage$ecf,
             threshold = "10",
             has_ecod = ecf.domain.coverage$ecod.eval10,
             has_pfam = ecf.domain.coverage$pfam.eval10)
)
long_dt <- rbindlist(long_list)

# Category + counts
long_dt[, category := assign_cat(has_ecod, has_pfam)]

counts_dt <- long_dt[, .N, by = .(threshold, category)]

# Factor ordering for nicer plot
counts_dt[, threshold := factor(threshold, levels = c("1e-3", "1e-1", "10"))]
counts_dt[, category := factor(
  category,
  levels = c("Both (ECOD + Pfam)", "ECOD only", "Pfam only", "Neither")
)]

# (Optional) a clean, colorblind-friendly palette
cat_cols <- c(
  "Both (ECOD + Pfam)" = "#009E73",
  "ECOD only"          = "#0072B2",
  "Pfam only"          = "#D55E00",
  "Neither"            = "grey70"
)

# Plot
p.ecf.domain.coverage <- ggplot(counts_dt, aes(x = threshold, y = N, fill = category)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.25) +
  scale_fill_manual(values = cat_cols, name = NULL) +
  labs(
    x = "E-value threshold (strict → lenient)",
    y = "Number of ECFs"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )
p.ecf.domain.coverage.outfile <- file.path(this.part.dir, "plot1_ecf-domain-cov.png")
ggsave(p.ecf.domain.coverage.outfile, p.ecf.domain.coverage, width = 6, height = 5)


#################################
### 2. 
### 
#################################



#################################
### 3. 
### 
#################################


##########################################################
##########################################################
##########################################################
##########################################################

reprseq.vs.ecfs <- reprseq.vs.ecfs.all %>% filter(qcov >= MIN.QCOV.FOR.DOMAIN, scov >= MIN.SCOV.FOR.DOMAIN)
reprseq.vs.ecod <- reprseq.vs.ecod.all %>% filter(qcov >= MIN.QCOV.FOR.DOMAIN, scov >= MIN.SCOV.FOR.DOMAIN)
reprseq.vs.pfam <- reprseq.vs.pfam.all %>% filter(qcov >= MIN.QCOV.FOR.DOMAIN, scov >= MIN.SCOV.FOR.DOMAIN)

reprseq.with.ecf.hits <- unique(reprseq.vs.ecfs$qname)
reprseq.with.ecod.hits <- unique(reprseq.vs.ecod$qname)
reprseq.with.pfam.hits <- unique(reprseq.vs.pfam$qname)
reprseq.with.domain.hits <- unique(c(reprseq.with.ecf.hits, reprseq.with.ecod.hits, reprseq.with.pfam.hits))

top.annotations <- reprseq.functions.all[reprseq.functions.all$top.annotation] %>% select(funct, class, annotation) %>% distinct()
setkey(reprseq.functions.all, annotation)

setkey(metadata.table, accession)

ps.vs.ecf <- ps.vs.ecf.all %>% filter(qcov >= MIN.QCOV.FOR.DOMAIN, scov >= MIN.SCOV.FOR.DOMAIN)
ps.vs.ecod <- ps.vs.ecod.all %>% filter(qcov >= MIN.QCOV.FOR.DOMAIN, scov >= MIN.SCOV.FOR.DOMAIN)
ps.vs.pfam <- ps.vs.pfam.all %>% filter(qcov >= MIN.QCOV.FOR.DOMAIN, scov >= MIN.SCOV.FOR.DOMAIN)

#################################
### 1. Compare annotation
### coverage between three 
### domain db in RefSeq (UpSet)
#################################

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
annot.dt[, combo := paste0(as.integer(ECF), as.integer(Pfam), as.integer(ECOD))]
combo_counts <- annot.dt[, .N, by = combo][order(-N)]
upset.df <- as.data.frame(annot.dt[, .(ECF, Pfam, ECOD)])

# A minimal, clean theme set
upset_themes <- list(
  "overall_sizes" = theme_bw(base_size = 8) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank()
          ),
  "intersections_matrix" = theme_bw(base_size = 8) +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    ),
  "Intersection size" = theme_bw(base_size = 8) +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    )
)

p.upset <- upset(
  upset.df,
  intersect = c("ECF", "Pfam", "ECOD"),
  mode = "exclusive_intersection",          # exclusive intersections (standard UpSet)
  sort_intersections_by = "cardinality",    # largest intersections first
  min_size = 1,                             # keep all
  width_ratio = 0.40,                       # space for set sizes vs intersections
  stripes = "white",                        # cleaner background
  set_sizes = (
    upset_set_size() +
      labs(x = "Proteins annotated")   # <- your new label
  ),
  base_annotations = list(
    "Intersection size" = intersection_size(
      text = list(size = 3),
      counts = TRUE
    )
  ),
  themes = upset_themes
)
# p.upset

# Save
p.upset.outfile <- file.path(this.part.dir, "plot1_upset-reprseq-domain-cov.png")
ggsave(p.upset.outfile, p.upset, width = 5.5, height = 4.2)

#################################
### 2. Compare the functional 
### coverage between the three 
### domain db in RefSeq
#################################

# Prepare functional coverage per domain table for plotting
funct.cov.list <- lapply(1:nrow(top.annotations), function(annotation.index){
  this.annotation <- top.annotations$annotation[annotation.index]
  this.funct <- top.annotations$funct[annotation.index]
  this.class <- top.annotations$class[annotation.index]
  this.reprseq.functions.all <- reprseq.functions.all[.(this.annotation)]
  all.reprseqs.with.this.annotation <- unique(this.reprseq.functions.all$qname)
  reprseqs.with.any.hit <- all.reprseqs.with.this.annotation %in% reprseq.with.domain.hits
  reprseqs.with.ecf.hit <- all.reprseqs.with.this.annotation %in% reprseq.with.ecf.hits
  reprseqs.with.ecod.hit <- all.reprseqs.with.this.annotation %in% reprseq.with.ecod.hits
  reprseqs.with.pfam.hit <- all.reprseqs.with.this.annotation %in% reprseq.with.pfam.hits
  
  data.table(funct = this.funct, class = this.class, annotation = this.annotation,
             no.reprseq = length(all.reprseqs.with.this.annotation),
             no.reprseq.with.any = length(reprseqs.with.any.hit),
             no.reprseq.with.ecf = lw(reprseqs.with.ecf.hit),
             no.reprseq.with.ecod = lw(reprseqs.with.ecod.hit),
             no.reprseq.with.pfam = lw(reprseqs.with.pfam.hit)
  )  
})
funct.cov <- rbindlist(funct.cov.list)
funct.cov[, `:=`(ecf_coverage = no.reprseq.with.ecf / no.reprseq,
                 ecod_coverage = no.reprseq.with.ecod / no.reprseq,
                 pfam_coverage = no.reprseq.with.pfam / no.reprseq)]
funct.cov <- funct.cov %>% filter(no.reprseq >= MIN.NUM.PROT.PER.ANNOT)

# Prepare the data for plotting
melted_data <- melt(funct.cov, id.vars = c("funct", "class"), 
                    measure.vars = c("ecf_coverage", "ecod_coverage", "pfam_coverage"),
                    variable.name = "database", value.name = "coverage")


# # Define the colours for the domain coverages
# coverage_colours <- c("ecf_coverage" = "#08306B", "ecod_coverage" = "#006D2C", "pfam_coverage" = "#800026")
# # Set the order of the facets
# melted_data$class <- factor(melted_data$class, levels = phrogs.class.order)
# # Count the number of functions per class
# n_funct <- melted_data[, .N, by = class]
# # Define strip backgrounds using ggh4x
# strip <- strip_themed(
#   background_x = elem_list_rect(fill = phrogs.col[match(levels(melted_data$class), phrogs.class)])
# )
# 
# # Create the ggplot object with customized strips and varying plot heights
# p.domain.cov <- ggplot(melted_data) +
#   geom_col(aes(x = funct, y = coverage, fill = database), 
#            position = "dodge", width = 0.7, linewidth = 0.2) +
#   scale_fill_manual(values = coverage_colours) +
#   theme(strip.text.y = element_text(face = "bold", size = 7, angle = 360)) +
#   labs(x = "Function", y = "Domain Coverage") +  
#   facet_grid2(cols = vars(class), scales = "free_x",  
#               labeller = labeller(class = phrog.class.labels),
#               strip = strip, axes = "x", space = "free") +  # Apply the `strip` object
#   guides(color = "none") + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),  
#         strip.placement = "outside",
#         panel.spacing = unit(0.1, "lines"))
# 
# p.domain.cov.outfile <- file.path(this.part.dir, "plot2_phrogs-domain-cov.png")
# ggsave(p.domain.cov.outfile, p.domain.cov, width = 14, height = 6)



########################
### GPT's version
########################

### after you create `melted_data` and set factor levels for class

# Pretty DB labels + consistent factor order
melted_data[, database := factor(
  database,
  levels = c("ecf_coverage", "ecod_coverage", "pfam_coverage"),
  labels = c("ECF", "ECOD", "Pfam")
)]

# Order functions within each class by ECF coverage (you can switch to max/mean if you want)
ord <- melted_data[database == "ECF",
                   .(ord_key = coverage),
                   by = .(class, funct)]
setorder(ord, class, -ord_key)
ord[, funct_ord := factor(funct, levels = unique(funct)), by = class]

melted_data <- merge(melted_data, ord[, .(class, funct, funct_ord)], by = c("class", "funct"))
melted_data[, funct := funct_ord]
melted_data[, funct_ord := NULL]

# Colorblind-safe palette (Okabe–Ito-ish; strong but not neon)
coverage_colours2 <- c(
  "ECF"  = "#0072B2",
  "ECOD" = "#009E73",
  "Pfam" = "#D55E00"
)

# (optional) show only top N functions per class (keeps plot readable)
TOP_N <- 10
top_funct <- melted_data[database == "ECF",
                         .SD[order(-coverage)][1:min(.N, TOP_N)],
                         by = class]$funct
melted_plot <- melted_data[funct %in% top_funct]
melted_plot[, class := factor(class, levels = phrogs.class)]


strip <- ggh4x::strip_themed(
  background_y = ggh4x::elem_list_rect(
    fill = grDevices::adjustcolor(
      phrogs.col[match(levels(melted_plot$class), phrogs.class)],
      alpha.f = 0.25
    ),
    colour = "black",
    linewidth = 0.4
  )
)

p.domain.cov.pretty <- ggplot(melted_plot, aes(x = coverage, y = funct, fill = database)) +
  geom_col(position = position_dodge(width = 0.75),
           width = .9, color = "white", linewidth = 0.2) +
  facet_grid2(
    rows = vars(class),
    scales = "free_y",
    space = "free_y",
    labeller = labeller(class = phrog.class.labels),
    strip = strip
  ) +
  scale_fill_manual(values = coverage_colours2, name = "Database") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Fraction of proteins with ≥1 hit",
       y = NULL) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    strip.placement = "outside",
    strip.text.y = element_text(angle = 0, hjust = 0, face = "bold", size = 9),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 9),
    panel.spacing = unit(0.35, "lines")
  )
p.domain.cov.pretty

p.domain.cov.outfile <- file.path(this.part.dir, "plot2_phrogs-domain-cov.png")
ggsave(p.domain.cov.outfile, p.domain.cov.pretty, width = 7.5, height = 9, dpi = 300)

#################################
### 3. Plot only per functional
### class, with total and unknown
#################################

all.phrogs.classes <- names(phrog.class.labels)
# annot.dt

# Funct cov per class
per.class.cov.list <- lapply(1:length(all.phrogs.classes), function(ci){
  this.class <- all.phrogs.classes[ci]
  reprseqs.with.this.class <- unique(reprseq.functions.all[class == this.class]$qname)
  this.class.annot.dt <- annot.dt[qname %chin% reprseqs.with.this.class]
  this.class.no.hits <- apply(this.class.annot.dt[,.(ECF, Pfam, ECOD)], 2, lw)
  this.class.funct.cov <- this.class.no.hits/length(reprseqs.with.this.class)
  data.table(category = this.class,
             ECF = this.class.funct.cov["ECF"],
             Pfam = this.class.funct.cov["Pfam"],
             ECOD = this.class.funct.cov["ECOD"])
})
per.class.cov <- rbindlist(per.class.cov.list)

# Funct cov in proteins without funct
reprseqs.no.funct <- all.reprseq[!all.reprseq %in% reprseq.functions.all$qname]

no.funct.annot.dt <- annot.dt[qname %chin% reprseqs.no.funct]
no.funct.no.hits <- apply(no.funct.annot.dt[,.(ECF, Pfam, ECOD)], 2, lw)
no.funct.funct.cov <- no.funct.no.hits/length(reprseqs.no.funct)

domain.cov.dt <- rbind(per.class.cov,
                       data.table(category = "unknown",
                                  ECF = no.funct.funct.cov["ECF"],
                                  Pfam = no.funct.funct.cov["Pfam"],
                                  ECOD = no.funct.funct.cov["ECOD"])
)

# Total domain coverage
total.no.hits <- apply(annot.dt[,.(ECF, Pfam, ECOD)], 2, lw)
total.funct.cov <- total.no.hits/nrow(annot.dt)

# Combine
domain.cov.dt <- rbind(domain.cov.dt,
                       data.table(category = "total",
                                  ECF = total.funct.cov["ECF"],
                                  Pfam = total.funct.cov["Pfam"],
                                  ECOD = total.funct.cov["ECOD"])
)


# ---- domain.cov.dt -> long format (dt_long) ----
domain.cov.dt.long <- melt(
  domain.cov.dt,
  id.vars = "category",
  variable.name = "database",
  value.name = "coverage"
)


# Make database a factor to control order + match your palette names
domain.cov.dt.long[, database := factor(database, levels = c("ECF", "ECOD", "Pfam"))]

# Order categories (optional: keep your preferred order)
cat_order <- c(
  "head and packaging",
  "connector",
  "tail",
  "DNA, RNA and nucleotide metabolism",
  "integration and excision",
  "moron, auxiliary metabolic gene and host takeover",
  "transcription regulation",
  "lysis",
  "antidefence",
  "other",
  "unknown",
  "total"
)
domain.cov.dt.long[, category := factor(category, levels = rev(cat_order))]

# phrog.class.labels
p.domain.cov.cat <- ggplot(domain.cov.dt.long, aes(x = coverage, y = category, fill = database)) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.7,
    color = "white",
    linewidth = 0.2
  ) +
  scale_fill_manual(values = coverage_colours2, name = "Database") +
  scale_x_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_discrete(labels = phrog.class.labels) +   # <- THIS
  labs(
    x = "Fraction of proteins with ≥1 hit",
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 9),
    panel.spacing = unit(0.35, "lines")
  )
# p.domain.cov.cat

p.domain.cov.cat.outfile <- file.path(this.part.dir, "plot3_phrogs-domain-cov-class.png")
ggsave(p.domain.cov.cat.outfile, p.domain.cov.cat, width = 4.5, height = 5.5)

#################################
### 4. Compare annotation
### coverage between three 
### domain db in PhageScope (UpSet)
#################################

phagescope.metadata[, ecf.hit := MGG_Protein_ID %in% ps.vs.ecf$qname]
phagescope.metadata[, ecod.hit := MGG_Protein_ID %in% ps.vs.ecod$qname]
phagescope.metadata[, pfam.hit := MGG_Protein_ID %in% ps.vs.pfam$qname]

protein.vs.domains <- phagescope.metadata[, c("MGG_Protein_ID", "ecf.hit", "ecod.hit", "pfam.hit")]

protein.vs.category <- phagescope.metadata[
  , .(category = unlist(strsplit(Protein_classification, ";", fixed = TRUE))),
  by = .(MGG_Protein_ID)
]

protein.vs.category[category %in% c("unsorted", "hypothetical")]$category <- "unknown"
unique.categories <- unique(protein.vs.category$category)


df.upset <- protein.vs.domains[, .(
  ECF  = as.integer(as.logical(ecf.hit)),
  ECOD = as.integer(as.logical(ecod.hit)),
  Pfam = as.integer(as.logical(pfam.hit))
)] |> as.data.frame()


p.upset.phagescope <- upset(
  df.upset,
  intersect = colnames(df.upset),
  mode = "exclusive_intersection",          # exclusive intersections (standard UpSet)
  sort_intersections_by = "cardinality",    # largest intersections first
  min_size = 1,                             # keep all
  width_ratio = 0.40,                       # space for set sizes vs intersections
  stripes = "white",                        # cleaner background
  set_sizes = (
    upset_set_size() +
      labs(x = "Proteins annotated")   # <- your new label
  ),
  themes = upset_themes
)
p.upset.phagescope

p.upset.phagescope.outfile <- file.path(this.part.dir, "plot4_upset-ps-domain-cov.png")
ggsave(p.upset.phagescope.outfile, p.upset.phagescope, width = 5.5, height = 4.2)

#################################
### 4. Compare annotation
### coverage between functional
### categories in PhageScope
#################################

functional.cov.domains.table.list <- lapply(1:length(unique.categories), function(index){
  this.category <- unique.categories[index]
  this.category.proteins <- protein.vs.category[category %chin% this.category]$MGG_Protein_ID
  this.category.proteins.no <- length(this.category.proteins)
  this.category.proteins.vs.domains <- protein.vs.domains[MGG_Protein_ID %chin% this.category.proteins]
  this.category.ecf.coverage <- length(which(this.category.proteins.vs.domains$ecf.hit))
  this.category.ecod.coverage <- length(which(this.category.proteins.vs.domains$ecod.hit))
  this.category.pfam.coverage <- length(which(this.category.proteins.vs.domains$pfam.hit))
  this.category.out <- data.table(category = this.category,
                                  ECF = this.category.ecf.coverage/this.category.proteins.no,
                                  ECOD = this.category.ecod.coverage/this.category.proteins.no,
                                  Pfam = this.category.pfam.coverage/this.category.proteins.no
  )
  return(this.category.out)
})
functional.cov.domains.table <- rbindlist(functional.cov.domains.table.list)

ps.total.number.of.hits <- apply(protein.vs.domains[, .(ecf.hit, ecod.hit, pfam.hit)], 2, lw)
ps.total.cov <- ps.total.number.of.hits/nrow(protein.vs.domains)
functional.cov.domains.table <- rbind(
  functional.cov.domains.table,
  data.table(category = "total", 
             ECF = ps.total.cov["ecf.hit"],
             ECOD = ps.total.cov["ecod.hit"],
             Pfam = ps.total.cov["pfam.hit"])
)

dt_long <- melt(
  functional.cov.domains.table,
  id.vars = "category",
  variable.name = "database",
  value.name = "coverage"
)

# functional.cov.domains.table$category
# [1] "unknown"      "integration"  "replication"  "assembly"     "regulation"   "infection"    "lysis"        "tRNA_related" "packaging"   
# [10] "immune"       "total"    

cat_order <- c(
  "lysis",
  "integration",
  "replication",
  "tRNA_related",
  "regulation",
  "packaging",
  "assembly",
  "infection",
  "immune",
  "unknown",
  "total"
)
dt_long[, category := factor(category, levels = rev(cat_order))]

p.phagescope.cov <- ggplot(dt_long, aes(x = coverage, y = category, fill = database)) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.7,
    color = "white",
    linewidth = 0.2
  ) +
  scale_fill_manual(values = coverage_colours2, name = "Database") +
  scale_x_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Fraction of proteins with ≥1 hit",
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9),
    panel.spacing = unit(0.35, "lines")
  )
p.phagescope.cov

p.phagescope.cov.outfile <- file.path(this.part.dir, "plot5_ps-domain-cov.png")
ggsave(p.phagescope.cov.outfile, p.phagescope.cov, width = 4.5, height = 5.5)

#################################
### Generate Figure 4
#################################

# p.upset # panel A
# p.domain.cov.cat # panel B
# p.upset.phagescope # panel C
# p.phagescope.cov # panel D


# wrap the ComplexUpset objects so they count as ONE panel each
p_upset_onepanel <- wrap_elements(full = p.upset)
p_upset_ps_onepanel <- wrap_elements(full = p.upset.phagescope)

fig4 <- (
  (p_upset_onepanel | p.domain.cov.cat) /
    (p_upset_ps_onepanel | p.phagescope.cov)
) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 24),
    plot.tag.position = c(0.02, 0.98)
  )

fig4.outfile <- file.path(ecf.figures.dir, "figure4.png")
ggsave(fig4.outfile, fig4, width = 12, height = 10)

#################################
### Generate Supplementary Figs
#################################

p.domain.cov.pretty.outfile <- file.path(ecf.suppl.figures.dir, "suppl-4p1_refseq-dom-cov-per-funct.png")
ggsave(p.domain.cov.pretty.outfile, p.domain.cov.pretty,  width = 7.5, height = 9)
