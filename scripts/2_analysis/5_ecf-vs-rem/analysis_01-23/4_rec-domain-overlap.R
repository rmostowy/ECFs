### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(stringr)))
# suppressWarnings(suppressPackageStartupMessages(library(seqinr)))
# suppressWarnings(suppressPackageStartupMessages(library(Biostrings)))
# suppressWarnings(suppressPackageStartupMessages(library(parallel)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

source("functions.R")

# Define output folder
table.dir <- "tables"
output.dir <- file.path(visual.dir, "part5_rem", "recs-vs-domains")
make.dir(output.dir)

# Load refseq table
reprseq.table <- fread(reprseq.table.output)
all.reprseq <- unique(reprseq.table$reprseq)

# Load recent recombination table
rem.table <- fread(rem.table.output)
rem.table[, cov.max.full := pmax(qcov.full, scov.full, na.rm = T)]

# Load recombination domain overlap table
rec.domain.overlap.table.outfile <- file.path(table.dir, "rem-domain-overlaps.txt")
rec.domain.overlaps <- fread(rec.domain.overlap.table.outfile)
rec.domain.overlaps[, row.id := .I]

# Load relative classes
relative.classes.filename <- file.path(table.dir, "relative-classes.txt")
relative.classes <- fread(relative.classes.filename)

# Load ECOD domain table
ecod.domains <- fread(ecod.domains.output)
setkey(ecod.domains, full.id)

# Load PFAM domain table
pfam.domains <- fread(pfam.domains.output)
setkey(pfam.domains, accession)

# Load ECF families
ecf.family.table <- fread(ecf.family.table.output)
setkey(ecf.family.table, ecf)

# Load ECF table
ecf.class.table <- fread(ecf.classification.table.output)
setkey(ecf.class.table, ecf)

# Translate ECFs to ECODs
ecf.class.table$family <- ecf.family.table[.(ecf.class.table$ecf)]$family
ecf2ecod <- ecf.class.table %>% select(ecf, family, t.id)
ecod.tid2tname <- unique(ecod.domains %>% select(t.id = t_id, t.name = t_name))
setkey(ecod.tid2tname, t.id)
ecf2ecod$t.name <- ecod.tid2tname[.(ecf2ecod$t.id)]$t.name
ecf2ecod[is.na(t.name), t.name := "unknown"]
setorder(ecf2ecod, family)

# Load functional annotations
reprseq.functions <- fread(reprseq.functions.unique.output)
setkey(reprseq.functions, reprseq)

# Load reprseq vs ECFs
reprseq.ecfs.all <- fread(reprseq.ecfs.output)
reprseq.ecfs <- reprseq.ecfs.all %>%
  filter(prob >= MIN.PROB.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(scov >= MIN.SCOV.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE)
reprseq.ecfs$EF <- ecf.family.table[.(reprseq.ecfs$sname)]$family
setkey(reprseq.ecfs, qname)

# Load reprseq vs ECFs relaxed
reprseq.ecfs.relaxed <- fread(reprseq.ecfs.all.output)

# Load reprseq vs ECOD
reprseq.ecod.all <- fread(reprseq.domains.ecod.output)
reprseq.ecod <- reprseq.ecod.all %>%
  filter(prob >= MIN.PROB.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(scov >= MIN.SCOV.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE)
reprseq.ecod$t.id <- ecod.domains[.(reprseq.ecod$sname)]$t_id
setkey(reprseq.ecod, qname)

# Load reprseq vs Pfam
reprseq.pfam.all <- fread(reprseq.domains.pfam.output)
reprseq.pfam <- reprseq.pfam.all %>%
  filter(prob >= MIN.PROB.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(scov >= MIN.SCOV.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE)
setkey(reprseq.pfam, qname)

#######################################################
## Compile a table of best overlapping domains with reprseqs
#######################################################

# Make a table to best rr-overlapping ECODs per reprseq
rows.from <- rec.domain.overlaps %>% select(reprseq = protein1.name, domain = protein1.ecod.best, score = protein1.ecod.score)
rows.to <- rec.domain.overlaps %>% select(reprseq = protein2.name, domain = protein2.ecod.best, score = protein2.ecod.score)
protein.rec.ecod.overlaps <- rbind(rows.from, rows.to)
protein.rec.ecod.overlaps[is.na(score), score := 0]
setorder(protein.rec.ecod.overlaps, reprseq, -score)
protein.rec.ecod.overlaps.best <- protein.rec.ecod.overlaps[, .SD[1], by = .(reprseq)]
protein.rec.ecod.overlaps.best$id <- ecod.domains[.(protein.rec.ecod.overlaps.best$domain)]$t_id
protein.rec.ecod.overlaps.best$name <- ecod.domains[.(protein.rec.ecod.overlaps.best$domain)]$t_name

# Make a table to best rr-overlapping PFAMs per reprseq
rows.from <- rec.domain.overlaps %>% select(reprseq = protein1.name, domain = protein1.pfam.best, score = protein1.pfam.score)
rows.to <- rec.domain.overlaps %>% select(reprseq = protein2.name, domain = protein2.pfam.best, score = protein2.pfam.score)
protein.rec.pfam.overlaps <- rbind(rows.from, rows.to)
protein.rec.pfam.overlaps$domain <- sapply(protein.rec.pfam.overlaps$domain, get.pfam.accession, USE.NAMES = F)
protein.rec.pfam.overlaps[is.na(score), score := 0]
setorder(protein.rec.pfam.overlaps, reprseq, -score)
protein.rec.pfam.overlaps.best <- protein.rec.pfam.overlaps[, .SD[1], by = .(reprseq)]
protein.rec.pfam.overlaps.best$id <- pfam.domains[.(protein.rec.pfam.overlaps.best$domain)]$id
protein.rec.pfam.overlaps.best$name <- pfam.domains[.(protein.rec.pfam.overlaps.best$domain)]$description

# Make a table to best rr-overlapping ECFs per reprseq
rows.from <- rec.domain.overlaps %>% select(reprseq = protein1.name, domain = protein1.ecf.best, score = protein1.ecf.score)
rows.to <- rec.domain.overlaps %>% select(reprseq = protein2.name, domain = protein2.ecf.best, score = protein2.ecf.score)
protein.rec.ecf.overlaps <- rbind(rows.from, rows.to)
protein.rec.ecf.overlaps[is.na(score), score := 0]
setorder(protein.rec.ecf.overlaps, reprseq, -score)
protein.rec.ecf.overlaps.best <- protein.rec.ecf.overlaps[, .SD[1], by = .(reprseq)]
protein.rec.ecf.overlaps.best$id <- ecf.class.table[.(protein.rec.ecf.overlaps.best$domain)]$family
protein.rec.ecf.overlaps.best$name <- protein.rec.ecf.overlaps.best$id

# Merge the table
protein.rec.domain.overlaps.best <- rbind(data.table("ECOD", protein.rec.ecod.overlaps.best),
      data.table("PFAM", protein.rec.pfam.overlaps.best),
      data.table("ECF", protein.rec.ecf.overlaps.best)
)
setnames(protein.rec.domain.overlaps.best, "V1", "database")

#######################################################
## Visualise overlaps of recombinations with domains
#######################################################

# Visualise the distribution of overlaps
protein.rec.domain.overlaps.best.hist <- protein.rec.domain.overlaps.best[score>0]
protein.rec.domain.overlaps.best.hist$database <- factor(protein.rec.domain.overlaps.best.hist$database, levels = unique(protein.rec.domain.overlaps.best$database))
p.score.overlap.histogram <- protein.rec.domain.overlaps.best.hist %>%
  ggplot(aes(x = score, fill = database)) +
  geom_histogram(binwidth = 0.01, color = "black", alpha = 0.7) +
  facet_grid(database ~ .) +  # One histogram per database
  theme_minimal() +
  scale_fill_manual(values = c("ECOD" = "blue", "PFAM" = "red", "ECF" = "green")) +
  labs(title = "Score Distribution Across Databases (Excluding 0 Overlap)",
       x = "Overlap Score",
       y = "Frequency",
       fill = "Database") +
  theme(legend.position = "none")
p.score.overlap.histogram.outfile <- file.path(output.dir, "plot1_rec-vs-domains.png")
ggsave(p.score.overlap.histogram.outfile, p.score.overlap.histogram, width = 9, height = 12)

# Visualise three categories of overlaps
overlap_categories <- protein.rec.domain.overlaps.best[, .(
  zero_overlap = sum(score == 0),
  partial_overlap = sum(score > 0 & score <= 0.9),
  complete_overlap = sum(score > 0.9)
), by = database]
overlap_categories_melted <- melt(overlap_categories, id.vars = "database", 
                                  variable.name = "OverlapType", value.name = "Count")

# Define colors for categories
overlap_colors <- c("zero_overlap" = "#F7FBFF", 
                    "partial_overlap" = "#9ECAE1", 
                    "complete_overlap" = "#08306B")

# Plot stacked proportions of overlap categories
p.overlap <- ggplot(overlap_categories_melted, aes(x = database, y = Count, fill = OverlapType)) +
  geom_bar(stat = "identity", position = "fill") +  # Stacked proportion bar plot
  scale_fill_manual(values = overlap_colors) +
  theme_minimal() +
  labs(title = "Proportion of Zero, Partial, and Complete Overlaps",
       x = "Database",
       y = "Proportion",
       fill = "Overlap Type") +
  theme(legend.position = "right") +
  scale_y_continuous(labels = scales::percent)

# Save the plot
p.overlap.outfile <- file.path(output.dir, "plot2_overlap-categories.png")
ggsave(p.overlap.outfile, p.overlap, width = 5, height = 4)


#######################################################
## Visualise overlap score distributions per most frequent domains
#######################################################

plot.overlaps.per.domain <- function(db.name, min.no.overlaps = 20){
  this.db.overlaps <- protein.rec.domain.overlaps.best[database == db.name & score > 0]
  overlap_counts <- this.db.overlaps[!is.na(id), .N, by = .(id, name)]
  top_id_counts <- overlap_counts[N >= min.no.overlaps]
  this.db.overlaps[, label := str_wrap(name, width = 30)]
  this.db.overlaps$id <- factor(this.db.overlaps$id, levels = top_id_counts$id)
  
  ggplot(this.db.overlaps[id %in% top_id_counts$id], 
         aes(x = factor(id, levels = top_id_counts$id), y = score, color = id)) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2) +  
    geom_boxplot(outlier.shape = NA, alpha = 0.3, width = 0.5, color = "black") +  
    scale_x_discrete(labels = setNames(this.db.overlaps$label, 
                                       this.db.overlaps$id)) +  # Use formatted labels
    labs(title = "Score Distribution for Top domain IDs",
         x = "IF and LABEL",
         y = "Score") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10),  
          legend.position = "none")  
  
}

p.overlaps.per.domain.ecod <- plot.overlaps.per.domain('ECOD')
p.overlaps.per.domain.ecod.outfile <- file.path(output.dir, "plot3_ecod-score-distr.png")
ggsave(p.overlaps.per.domain.ecod.outfile, p.overlaps.per.domain.ecod, width = 15, height = 7)

p.overlaps.per.domain.pfam <- plot.overlaps.per.domain('PFAM')
p.overlaps.per.domain.pfam.outfile <- file.path(output.dir, "plot4_pfam-score-distr.png")
ggsave(p.overlaps.per.domain.pfam.outfile, p.overlaps.per.domain.pfam, width = 15, height = 7)

