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
output.dir <- file.path(visual.dir, "analysis_02-10", "figures", "part3_ecf-comp")
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
  filter(eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE)
reprseq.ecod$t.id <- ecod.domains[.(reprseq.ecod$sname)]$t_id
setkey(reprseq.ecod, qname)

# Load reprseq vs Pfam
reprseq.pfam.all <- fread(reprseq.domains.pfam.output)
reprseq.pfam <- reprseq.pfam.all %>%
  filter(prob >= MIN.PROB.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE)
setkey(reprseq.pfam, qname)

# Load ECFs vs ECOD
ecfs.vs.ecod <- fread(ecfs.vs.ecod.output)
setkey(ecfs.vs.ecod, qname)

# Load recombination domain overlap table
rec.domain.overlaps <- fread("rec-domain-overlaps.txt")
rec.domain.overlaps[, pair.index := 1:nrow(rec.domain.overlaps)]

rows.from <- rec.domain.overlaps %>% select(reprseq = protein1.name, domain = protein1.ecf.best, score = protein1.ecf.score)
rows.to <- rec.domain.overlaps %>% select(reprseq = protein2.name, domain = protein2.ecf.best, score = protein2.ecf.score)

protein.rec.ecf.overlaps <- rbind(rows.from, rows.to)
protein.rec.ecf.overlaps[is.na(score), score := 0]
setorder(protein.rec.ecf.overlaps, reprseq, -score)
protein.rec.ecf.overlaps.best <- protein.rec.ecf.overlaps[, .SD[1], by = .(reprseq)]

###############################################
#### Distribution of overlap scores
###############################################

p.rec.overlaps <- ggplot(protein.rec.ecf.overlaps.best, aes(x = score)) + 
  geom_histogram(fill = "darkgreen", bins = 100) + ggtitle("Histogram of overlap scores with ECFs per reprseq")
p.rec.overlaps.outfile <- file.path(output.dir, "overlap-score-distr.png")
ggsave(p.rec.overlaps.outfile, p.rec.overlaps, width = 6, height = 4)

p.rec.overlaps.nonzero <- ggplot(protein.rec.ecf.overlaps.best[score>0], aes(x = score)) + 
  geom_histogram(fill = "darkgreen", bins = 100) + ggtitle("Histogram of non-zero overlap scores with ECFs per reprseq")
p.rec.overlaps.nonzero.outfile <- file.path(output.dir, "overlap-score-distr-nz.png")
ggsave(p.rec.overlaps.nonzero.outfile, p.rec.overlaps.nonzero, width = 6, height = 4)

protein.rec.ecf.overlaps.best$family <- ecf.family.table[.(protein.rec.ecf.overlaps.best$domain)]$family
protein.rec.ecf.overlaps.best.nonzero <- protein.rec.ecf.overlaps.best[score>0]

###############################################
#### ECFs associated with recent-rec
###############################################
ecf.frequency <- protein.rec.ecf.overlaps.best.nonzero[, .N, by = domain]
setorder(ecf.frequency, -N)
common.ecfs <- ecf.frequency[N>=5]$domain

avg_ecf_scores <- protein.rec.ecf.overlaps.best.nonzero[domain %in% common.ecfs, 
                                                          .(avg_score = mean(score, na.rm = TRUE)), 
                                                          by = .(domain, family)]
rr.associated.ecfs <- sort(avg_ecf_scores[avg_score>0.9]$domain)
rr.associated.ecfs.reprseq.scov <- sapply(1:length(rr.associated.ecfs), function(k) mean(reprseq.ecfs[sname == rr.associated.ecfs[k]]$scov))

rr.associated.ecfs.metadata <- ecf.class.table[.(rr.associated.ecfs)]
# View()

index=14
this.ecf <- rr.associated.ecfs[index]
this.ecf.ecod.hits <- ecfs.vs.ecod[.(this.ecf)]
this.ecf.ecod.hits$t.id <- ecod.domains[.(this.ecf.ecod.hits$sname)]$t_id
this.ecf.ecod.hits$t.name <- ecod.domains[.(this.ecf.ecod.hits$sname)]$t_name
no.hits <- nrow(this.ecf.ecod.hits[!is.na(qstart)])
if(no.hits > 0){
  hist(this.ecf.ecod.hits$scov, xlim = c(0,1))
  my.table(this.ecf.ecod.hits$t.id)
  my.table(this.ecf.ecod.hits$t.name)
} else{
  cat("No hits!\n")
}






scov.distribution.associated <- reprseq.ecfs.relaxed[sname %in% rr.associated.ecfs] %>% select(sname, scov)
scov.distribution.nonassociated <- reprseq.ecfs.relaxed[!sname %in% rr.associated.ecfs] %>% select(sname, scov)

scov.associated.vs.other <- rbind(data.table(class = "rr.associated", scov.distribution.associated),
                                  data.table(class = "other", scov.distribution.nonassociated)
)

plot.scov.rr.associated.ecfs <-  ggplot(scov.associated.vs.other, aes(x = scov, fill = class)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(title = "SCOV Density Plot",
       x = "SCOV",
       y = "Density",
       fill = "Class") +
  scale_fill_manual(values = c("rr.associated" = "blue", "other" = "red"))
plot.scov.rr.associated.pfams.outfile <- file.path(output.dir, "scov-rr-associated-pfams.png")
ggsave(plot.scov.rr.associated.pfams.outfile, plot.scov.rr.associated.pfams, width = 8, height = 9)


####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################


###############################################
#### Distribution of functions
###############################################

# number of functions to plot
no.top.functions <- 20

# Count occurrences of each TID, filtering out NA values
tid_counts <- protein.rec.pfam.overlaps.best[!is.na(id), .N, by = .(id, name)]
setorder(tid_counts, -N)
top_tid_counts <- tid_counts[1:no.top.functions]

# Plot
p.funct.distribution <- ggplot(top_tid_counts, aes(x = reorder(id, -N), y = N, fill = id)) +
  geom_bar(stat = "identity", show.legend = FALSE) + 
  geom_text(aes(label = name), angle = 90, hjust = -0.1, vjust = 0.5, size = 3.5, color = "black") + 
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

protein.rec.pfam.overlaps.best.nonzero <- protein.rec.pfam.overlaps.best[score > 0]
overlap_counts <- protein.rec.pfam.overlaps.best.nonzero[!is.na(id), .N, by = .(id, name)]
top_tid_counts <- overlap_counts[N >= min.no.overlaps]

protein.rec.pfam.overlaps.best.nonzero[, label := str_wrap(name, width = 30)]

# Ensure the order remains the same as before
protein.rec.pfam.overlaps.best.nonzero$id <- factor(protein.rec.pfam.overlaps.best.nonzero$id, 
                                              levels = top_tid_counts$id)

# Jitter plot with updated labels
p.overlap.per.funct.distribution <- ggplot(protein.rec.pfam.overlaps.best.nonzero[id %in% top_tid_counts$id],
                                           aes(x = factor(id, levels = top_tid_counts$id), y = score, color = id)) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +  
  geom_boxplot(outlier.shape = NA, alpha = 0.3, width = 0.5, color = "black") +  
  scale_x_discrete(labels = setNames(protein.rec.pfam.overlaps.best.nonzero$label, 
                                     protein.rec.pfam.overlaps.best.nonzero$id)) +  # Use formatted labels
  labs(title = "Score Distribution for Top pfam TIDs",
       x = "TID and Function",
       y = "Score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10),  
        legend.position = "none")  

p.overlap.per.funct.distribution.outfile <- file.path(output.dir, "overlap-per-funct-distr.png")
ggsave(p.overlap.per.funct.distribution.outfile, p.overlap.per.funct.distribution, width = 8, height = 9)

###############################################
#### PFAMs associated with recent-rec
###############################################
avg_pfam_scores <- protein.rec.pfam.overlaps.best.nonzero[, 
                                                          .(avg_score = mean(score, na.rm = TRUE)), 
                                                          by = .(domain, name)]
rr.associated.pfams <- sort(avg_pfam_scores[avg_score>0.9]$domain)

scov.distribution.associated <- reprseq.pfam[accession %in% rr.associated.pfams] %>% select(accession, scov)
scov.distribution.nonassociated <- reprseq.pfam[!accession %in% rr.associated.pfams] %>% select(accession, scov)

scov.associated.vs.other <- rbind(data.table(class = "rr.associated", scov.distribution.associated),
                                  data.table(class = "other", scov.distribution.nonassociated)
)

plot.scov.rr.associated.pfams <-  ggplot(scov.associated.vs.other, aes(x = scov, fill = class)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(title = "SCOV Density Plot",
       x = "SCOV",
       y = "Density",
       fill = "Class") +
  scale_fill_manual(values = c("rr.associated" = "blue", "other" = "red"))
plot.scov.rr.associated.pfams.outfile <- file.path(output.dir, "scov-rr-associated-pfams.png")
ggsave(plot.scov.rr.associated.pfams.outfile, plot.scov.rr.associated.pfams, width = 8, height = 9)



# scov.distribution.nonassociated <- 
 


# # rr.associated.pfams <- sort(unique(protein.rec.pfam.overlaps.best.nonzero[score>0.9]$domain))
# 
# rr.associated.pfams.pairs <- lapply(1:length(rr.associated.pfams), function(index){
#   this.pfam <- rr.associated.pfams[index]
#   this.pfam.descr <- pfam.domains[.(this.pfam)]$description
#   this.pfam.rr <- rec.domain.overlaps[grepl(this.pfam, rec.domain.overlaps$protein1.pfam.best) | grepl(this.pfam, rec.domain.overlaps$protein2.pfam.best)]
#   this.pfam.rr <- this.pfam.rr[protein1.pfam.score >0 | protein2.pfam.score>0]
#   this.pfam.rr$pair.index
# })
# rr.associated.pfams.pairs <- sort(unlist(rr.associated.pfams.pairs))
# 
# 
# rr.swapped.ecods.list <- lapply(1:length(rr.associated.pfams.pairs), function(index){
#   this.pair.index <- rr.associated.pfams.pairs[index]
#   this.pair.rr <- recent.rec.table[this.pair.index]
#   pair.qname.nonoverlap.ecods <- pair.sname.nonoverlap.ecods <- NA
#   
#   pair.qname <- this.pair.rr$qname
#   pair.qstart <- this.pair.rr$qstart
#   pair.qend <- this.pair.rr$qend
#   pair.qname.ecod.hits <- reprseq.ecod[.(pair.qname)]
#   pair.qname.ecod.hits[, rr.overlap.jaccard := calculate_overlap_jaccard(pair.qstart, pair.qend, qstart, qend)]
#   pair.qname.ecod.hits.nonoverlap <- pair.qname.ecod.hits[rr.overlap.jaccard==0]
#   no.non.overlapping.ecods.query <- nrow(pair.qname.ecod.hits.nonoverlap)
#   if(no.non.overlapping.ecods.query>0) pair.qname.nonoverlap.ecods <- paste(unique(pair.qname.ecod.hits.nonoverlap$t.id), collapse = "|")
#   
#   pair.sname <- this.pair.rr$sname
#   pair.sstart <- this.pair.rr$sstart
#   pair.send <- this.pair.rr$send
#   pair.sname.ecod.hits <- reprseq.ecod[.(pair.sname)]
#   pair.sname.ecod.hits[, rr.overlap.jaccard := calculate_overlap_jaccard(pair.qstart, pair.qend, qstart, qend)]
#   pair.sname.ecod.hits.nonoverlap <- pair.sname.ecod.hits[rr.overlap.jaccard==0]
#   no.non.overlapping.ecods.subject <- nrow(pair.sname.ecod.hits.nonoverlap)
#   if(no.non.overlapping.ecods.subject>0) pair.sname.nonoverlap.ecods <- paste(unique(pair.sname.ecod.hits.nonoverlap$t.id), collapse = "|")
#   
#   data.table(pair.index = this.pair.index, qname = pair.qname, sname = pair.sname,
#              query.nonoverlap.ecods = pair.qname.nonoverlap.ecods,
#              subject.nonoverlap.ecods = pair.sname.nonoverlap.ecods)  
# })
# rr.swapped.ecods <- rbindlist(rr.swapped.ecods.list)
# 
# 
# 
# no.non.overlapping.ecods.query
# no.non.overlapping.ecods.subject
# 
# 
