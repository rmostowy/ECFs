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

# Define output folder
output.dir <- file.path(visual.dir, "analysis_02-10", "figures", "part1_ecod-comp")
make.dir(output.dir)

# Load refseq table
reprseq.table <- fread(reprseq.table.filename)
all.reprseq <- unique(reprseq.table$reprseq)

# Calculate NCBI overlaps per reprseq
max.overlap.per.reprseq <- reprseq.table[, .(reprseq, ncbi.overlap)]
setorder(max.overlap.per.reprseq, reprseq, -ncbi.overlap)
max.overlap.per.reprseq <- max.overlap.per.reprseq[, .SD[1], by = reprseq]
setkey(max.overlap.per.reprseq, reprseq)

# Load recent recombination table
recent.rec.table <- fread("rececent-rec-table-sub.txt")
recent.rec.table[, cov.max.full := pmax(qcov.full, scov.full, na.rm = T)]

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

ecf.class.table$family <- ecf.family.table[.(ecf.class.table$ecf)]$family
ecf2ecod <- ecf.class.table %>% select(ecf, family, t.id)
ecf2ecod$t.name <- sapply(1:nrow(ecf2ecod), function(k){
  out <- "unknown"
  this.t.id <- ecf2ecod$t.id[k]
  if(this.t.id != "unknown") out <- ecod.domains[t_id == this.t.id]$t_name[1]  
  return(out)
})
setorder(ecf2ecod, family)

# Load functional annotations
reprseq.functions <- fread(reprseq.functions.unique.output)
setkey(reprseq.functions, reprseq)

# Load recombination domain overlap table
rec.domain.overlaps <- fread("rec-domain-overlaps.txt")

rows.from <- rec.domain.overlaps %>% select(reprseq = protein1.name, domain = protein1.ecod.best, score = protein1.ecod.score)
rows.to <- rec.domain.overlaps %>% select(reprseq = protein2.name, domain = protein2.ecod.best, score = protein2.ecod.score)

protein.rec.ecod.overlaps <- rbind(rows.from, rows.to)
protein.rec.ecod.overlaps[is.na(score), score := 0]
setorder(protein.rec.ecod.overlaps, reprseq, -score)
protein.rec.ecod.overlaps.best <- protein.rec.ecod.overlaps[, .SD[1], by = .(reprseq)]

###############################################
#### Distribution of overlap scores
###############################################

p.ecod.rec.overlaps <- ggplot(protein.rec.ecod.overlaps.best, aes(x = score)) + geom_histogram(fill = "darkgreen", bins = 100)
p.ecod.rec.overlaps <- p.ecod.rec.overlaps + ggtitle("Histogram of overlap scores with ECOD per reprseq")
p.ecod.rec.overlaps.outfile <- file.path(output.dir, "overlap-score-distr.png")
ggsave(p.ecod.rec.overlaps.outfile, p.ecod.rec.overlaps, width = 6, height = 4)


protein.rec.ecod.overlaps.best$t.id <- ecod.domains[.(protein.rec.ecod.overlaps.best$domain)]$t_id
protein.rec.ecod.overlaps.best$t.name <- ecod.domains[.(protein.rec.ecod.overlaps.best$domain)]$t_name

###############################################
#### Distribution of functions
###############################################

# number of functions to plot
no.top.functions <- 20

# Count occurrences of each TID, filtering out NA values
tid_counts <- protein.rec.ecod.overlaps.best[!is.na(t.id), .N, by = .(t.id, t.name)]
setorder(tid_counts, -N)
top_tid_counts <- tid_counts[1:no.top.functions]

# Plot
p.funct.distribution <- ggplot(top_tid_counts, aes(x = reorder(t.id, -N), y = N, fill = t.id)) +
  geom_bar(stat = "identity", show.legend = FALSE) + 
  geom_text(aes(label = t.name), angle = 90, hjust = -0.1, vjust = 0.5, size = 3.5, color = "black") + 
  scale_y_continuous(expand = expansion(mult = c(0, 1))) + 
  labs(title = "Top ECOD TIDs by Frequency",
       x = "TID",
       y = "Frequency") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Keep TIDs readable
p.funct.distribution.outfile <- file.path(output.dir, "funct-distr.png")
ggsave(p.funct.distribution.outfile, p.funct.distribution, width = 6, height = 8)

###############################################
#### Distribution of overlap scores per functions
###############################################

# Create a new label column with TID and TNAME combined, inserting a line break after ~30 characters
# protein.rec.ecod.overlaps.best[, label := paste0(t.id, "\n", str_wrap(t.name, width = 30))]
min.no.overlaps <- 10

protein.rec.ecod.overlaps.best.nonzero <- protein.rec.ecod.overlaps.best[score > 0]
overlap_counts <- protein.rec.ecod.overlaps.best.nonzero[!is.na(t.id), .N, by = .(t.id, t.name)]
top_tid_counts <- overlap_counts[N >= min.no.overlaps]

protein.rec.ecod.overlaps.best.nonzero[, label := str_wrap(t.name, width = 30)]

# Ensure the order remains the same as before
protein.rec.ecod.overlaps.best.nonzero$t.id <- factor(protein.rec.ecod.overlaps.best.nonzero$t.id, 
                                              levels = top_tid_counts$t.id)

# Jitter plot with updated labels
p.overlap.per.funct.distribution <- ggplot(protein.rec.ecod.overlaps.best.nonzero[t.id %in% top_tid_counts$t.id], 
                                           aes(x = factor(t.id, levels = top_tid_counts$t.id), y = score, color = t.id)) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +  
  geom_boxplot(outlier.shape = NA, alpha = 0.3, width = 0.5, color = "black") +  
  scale_x_discrete(labels = setNames(protein.rec.ecod.overlaps.best.nonzero$label, 
                                     protein.rec.ecod.overlaps.best.nonzero$t.id)) +  # Use formatted labels
  labs(title = "Score Distribution for Top ECOD TIDs",
       x = "TID and Function",
       y = "Score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10),  
        legend.position = "none")  

p.overlap.per.funct.distribution.outfile <- file.path(output.dir, "overlap-per-funct-distr.png")
ggsave(p.overlap.per.funct.distribution.outfile, p.overlap.per.funct.distribution, width = 8, height = 9)

