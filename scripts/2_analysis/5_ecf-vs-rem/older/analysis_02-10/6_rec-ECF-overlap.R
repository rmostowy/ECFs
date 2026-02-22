### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(stringr)))
# suppressWarnings(suppressPackageStartupMessages(library(seqinr)))
suppressWarnings(suppressPackageStartupMessages(library(Biostrings)))
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
setkey(reprseq.table, reprseq)
all.reprseq <- unique(reprseq.table$reprseq)

# Calculate NCBI overlaps per reprseq
max.overlap.per.reprseq <- reprseq.table[, .(reprseq, ncbi.overlap)]
setorder(max.overlap.per.reprseq, reprseq, -ncbi.overlap)
max.overlap.per.reprseq <- max.overlap.per.reprseq[, .SD[1], by = reprseq]
setkey(max.overlap.per.reprseq, reprseq)

# Load recent recombination table
recent.rec.table <- fread("rececent-rec-table-sub.txt")
recent.rec.table[, cov.max.full := pmax(qcov.full, scov.full, na.rm = T)]

# Load ECF table
ecf.table <- fread(ecf.table.output)
setkey(ecf.table, ecf.name)

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
common.ecfs <- ecf.frequency[N>=3]$domain

avg_ecf_scores <- protein.rec.ecf.overlaps.best.nonzero[domain %in% common.ecfs, 
                                                          .(avg_score = mean(score, na.rm = TRUE)), 
                                                          by = .(domain, family)]
rr.associated.ecfs <- sort(avg_ecf_scores[avg_score>0.9]$domain)
# rr.associated.ecfs.reprseq.scov <- sapply(1:length(rr.associated.ecfs), function(k) mean(reprseq.ecfs[sname == rr.associated.ecfs[k]]$scov))

ecf.associated.pairs <- lapply(1:length(rr.associated.ecfs), function(index){
  this.ecf <- rr.associated.ecfs[index]
  this.ecf.pair.indices <- sort(union(rec.domain.overlaps[protein1.ecf.best == this.ecf]$pair.index, rec.domain.overlaps[protein2.ecf.best == this.ecf]$pair.index))  
  this.ecf.pair.indices
})
ecf.associated.pairs.n <- sapply(ecf.associated.pairs, length)

rr.associated.ecf.metadata <- cbind(ecf.class.table[.(rr.associated.ecfs)], no.pairs = ecf.associated.pairs.n)
rr.associated.ecf.metadata$ecf.length <- ecf.table[.(rr.associated.ecf.metadata$ecf)]$length

setorder(rr.associated.ecf.metadata, family)

###############################################
#### Extract exchanged sequences per ECF
###############################################

# Load reprseq sequences
reprseq.seqs <- readAAStringSet(reprseq.seqs.filename)
sapply(reprseq.seqs, length)
reprseq.seqs.n <- names(reprseq.seqs)
reprseq.seqs.l <- nchar(reprseq.seqs)
reprseq.lengths <- data.table(name = reprseq.seqs.n, length = reprseq.seqs.l)
setkey(reprseq.lengths, name)

subseq.outpath <- file.path(output.dir, "ecf-subseqs")
make.dir(subseq.outpath)

ecf2find <- "ecf_00084"
index <- which(rr.associated.ecfs == ecf2find)

for(index in 5:length(rr.associated.ecfs)){
  cat(index,"\n")
  this.ecf <- rr.associated.ecfs[index]
  this.ecf.associated.reprseqs <- rbind(rec.domain.overlaps[protein1.ecf.best == this.ecf, .(qname = protein1.name, pair.index)], 
                                        rec.domain.overlaps[protein2.ecf.best == this.ecf, .(qname = protein2.name, pair.index)])
  setorder(this.ecf.associated.reprseqs, qname, pair.index)
  this.ecf.associated.reprseqs.unique <- this.ecf.associated.reprseqs[, .SD[1], by = qname]
  
  ######
  ##||||

  coords.list <- lapply(1:nrow(this.ecf.associated.reprseqs), function(k){
    pi <- this.ecf.associated.reprseqs$pair.index[k]
    pi.qname <- this.ecf.associated.reprseqs$qname[k]
    rt <- recent.rec.table[pi]
    is.query <- pi.qname == rt$qname
    is.subject <- pi.qname == rt$sname
    if(is.query){
      rstart <- rt$qstart
      rend <- rt$qend
    } else if(is.subject){
      rstart <- rt$sstart
      rend <- rt$send
    }
    pi.ecf <- reprseq.ecfs[qname == pi.qname & sname == this.ecf]
    data.table(pair.index = pi, reprseq = pi.qname, rstart, rend, 
               ecf.start = pi.ecf$qstart, ecf.end = pi.ecf$qend)  
  })
  coords <- rbindlist(coords.list)
  
  ##||||
  ######
  
  
  
  this.ecf.reprseq.hits.list <- lapply(1:nrow(this.ecf.associated.reprseqs.unique), function(i){
    this.pair.index <- this.ecf.associated.reprseqs.unique$pair.index[i]
    this.pair.reprseq <- this.ecf.associated.reprseqs.unique$qname[i]
    this.pair.rec.table <- recent.rec.table[this.pair.index]
    is.query <- this.pair.reprseq == this.pair.rec.table$qname
    is.subject <- this.pair.reprseq == this.pair.rec.table$sname
    if(is.query){
      start <- this.pair.rec.table$qstart
      end <- this.pair.rec.table$qend
    } else if(is.subject){
      start <- this.pair.rec.table$sstart
      end <- this.pair.rec.table$send
    } else{
      cat("reprseq missing!\n")
    }
    data.table(name = this.pair.reprseq, start, end)  
  })
  this.ecf.reprseq.hits <- rbindlist(this.ecf.reprseq.hits.list)
  setkey(this.ecf.reprseq.hits, name)
  
  this.ecf.reprseq.hits$length <- reprseq.lengths[.(this.ecf.reprseq.hits$name)]$length
  this.ecf.reprseq.seqs <- reprseq.seqs[reprseq.seqs.n %in% this.ecf.reprseq.hits$name]
  
  # Extract subsequences using Biostrings' efficient subseq function
  extracted.subseqs <- this.ecf.reprseq.hits[, {
    # Extract full sequence from AAStringSet
    full_seq <- as.character(this.ecf.reprseq.seqs[[name]])
    
    # Extract the sub-sequence using `subseq()`
    sub_seq <- subseq(full_seq, start = start, end = end)
    
    # Return as a list to keep structure
    .(sub_seq)
  }, by = name]
  
  this.ecf.reprseq.subseqs <- AAStringSet(extracted.subseqs$sub_seq)
  names(this.ecf.reprseq.subseqs) <- this.ecf.reprseq.hits$name
  this.ecf.outfile <- file.path(subseq.outpath, sprintf("%s.fasta", this.ecf))
  writeXStringSet(this.ecf.reprseq.subseqs, this.ecf.outfile)
  
  # Find all ECOD domains associated with those reprseqs
  this.ecf.associated.reprseq.ecod.hits <- reprseq.ecod[qname %in% this.ecf.associated.reprseqs.unique$qname]
  this.ecf.associated.reprseq.ecod.hits.best <- this.ecf.associated.reprseq.ecod.hits[
    order(-scov),  # Sort by scov in descending order
    .SD[1],  # Select the first (best) hit per group
    by = .(qname, t.id)
  ]
  best.ecod.hits <- this.ecf.associated.reprseq.ecod.hits.best[, .(qname, sname, qstart, qend, t.id)]
  best.ecod.hits$t.name <- ecod.domains[.(best.ecod.hits$sname)]$t_name
  setnames(best.ecod.hits, c("name", "domain", "start", "end", "id", "description"))
  setorder(best.ecod.hits, name, start)
  
  # Find all PFAM domains associated with those reprseqs
  this.ecf.associated.reprseq.pfam.hits <- reprseq.pfam[qname %in% this.ecf.associated.reprseqs.unique$qname & scov >= MIN.SCOV.FOR.DOMAIN]
  this.ecf.associated.reprseq.pfam.hits$id <- pfam.domains[.(this.ecf.associated.reprseq.pfam.hits$accession)]$id
  this.ecf.associated.reprseq.pfam.hits$descr <- pfam.domains[.(this.ecf.associated.reprseq.pfam.hits$accession)]$description
  best.pfam.hits <- this.ecf.associated.reprseq.pfam.hits[, .(qname, sname, qstart, qend, id, descr)]
  setnames(best.pfam.hits, c("name", "domain", "start", "end", "id", "description"))
  setorder(best.pfam.hits, name, start)
  
  ecod.data.outfile <- file.path(subseq.outpath, sprintf("%s-ecod.txt", this.ecf))
  pfam.data.outfile <- file.path(subseq.outpath, sprintf("%s-pfam.txt", this.ecf))
  fragment.data.outfile <- file.path(subseq.outpath, sprintf("%s-frag.txt", this.ecf))
  
  fwrite(best.ecod.hits, file = ecod.data.outfile)
  fwrite(best.pfam.hits, file = pfam.data.outfile)
  fwrite(this.ecf.reprseq.hits, file = fragment.data.outfile)
  
}

###############################################
#### Extract exchanged sequences per EFs
###############################################

subseq.outpath <- file.path(output.dir, "EF-subseqs")
make.dir(subseq.outpath)

all.efs <- unique(rr.associated.ecf.metadata$family)


for(index in 1:length(all.efs)){
  cat(index,"\n")
  this.ef <- all.efs[index]
  this.ef.ecfs <- rr.associated.ecf.metadata[family == this.ef]$ecf
  
  this.ef.associated.reprseqs <- rbindlist(lapply(this.ef.ecfs, function(ecf){
    rbind(rec.domain.overlaps[protein1.ecf.best == ecf, .(qname = protein1.name, pair.index)], 
          rec.domain.overlaps[protein2.ecf.best == ecf, .(qname = protein2.name, pair.index)])
  }))
  setorder(this.ef.associated.reprseqs, qname, pair.index)
  this.ef.associated.reprseqs.unique <- this.ef.associated.reprseqs[, .SD[1], by = qname]
  
  
  this.ef.reprseq.hits.list <- lapply(1:nrow(this.ef.associated.reprseqs.unique), function(i){
    this.pair.index <- this.ef.associated.reprseqs.unique$pair.index[i]
    this.pair.reprseq <- this.ef.associated.reprseqs.unique$qname[i]
    this.pair.rec.table <- recent.rec.table[this.pair.index]
    is.query <- this.pair.reprseq == this.pair.rec.table$qname
    is.subject <- this.pair.reprseq == this.pair.rec.table$sname
    if(is.query){
      start <- this.pair.rec.table$qstart
      end <- this.pair.rec.table$qend
    } else if(is.subject){
      start <- this.pair.rec.table$sstart
      end <- this.pair.rec.table$send
    } else{
      cat("reprseq missing!\n")
    }
    data.table(name = this.pair.reprseq, start, end)  
  })
  this.ef.reprseq.hits <- rbindlist(this.ef.reprseq.hits.list)
  setkey(this.ef.reprseq.hits, name)
  
  this.ef.reprseq.hits$length <- reprseq.lengths[.(this.ef.reprseq.hits$name)]$length
  this.ef.reprseq.seqs <- reprseq.seqs[reprseq.seqs.n %in% this.ef.reprseq.hits$name]
  
  # Extract subsequences using Biostrings' efficient subseq function
  extracted.subseqs <- this.ef.reprseq.hits[, {
    # Extract full sequence from AAStringSet
    full_seq <- as.character(this.ef.reprseq.seqs[[name]])
    
    # Extract the sub-sequence using `subseq()`
    sub_seq <- subseq(full_seq, start = start, end = end)
    
    # Return as a list to keep structure
    .(sub_seq)
  }, by = name]
  
  this.ef.reprseq.subseqs <- AAStringSet(extracted.subseqs$sub_seq)
  names(this.ef.reprseq.subseqs) <- this.ef.reprseq.hits$name
  this.ef.outfile <- file.path(subseq.outpath, sprintf("%s.fasta", this.ef))
  writeXStringSet(this.ef.reprseq.subseqs, this.ef.outfile)
  
  # Find all ECOD domains associated with those reprseqs
  this.ef.associated.reprseq.ecod.hits <- reprseq.ecod[qname %in% this.ef.associated.reprseqs.unique$qname]
  this.ef.associated.reprseq.ecod.hits.best <- this.ef.associated.reprseq.ecod.hits[
    order(-scov),  # Sort by scov in descending order
    .SD[1],  # Select the first (best) hit per group
    by = .(qname, t.id)
  ]
  best.ecod.hits <- this.ef.associated.reprseq.ecod.hits.best[, .(qname, sname, qstart, qend, t.id)]
  best.ecod.hits$t.name <- ecod.domains[.(best.ecod.hits$sname)]$t_name
  setnames(best.ecod.hits, c("name", "domain", "start", "end", "id", "description"))
  setorder(best.ecod.hits, name, start)
  
  # Find all PFAM domains associated with those reprseqs
  this.ef.associated.reprseq.pfam.hits <- reprseq.pfam[qname %in% this.ef.associated.reprseqs.unique$qname & scov >= MIN.SCOV.FOR.DOMAIN]
  this.ef.associated.reprseq.pfam.hits$id <- pfam.domains[.(this.ef.associated.reprseq.pfam.hits$accession)]$id
  this.ef.associated.reprseq.pfam.hits$descr <- pfam.domains[.(this.ef.associated.reprseq.pfam.hits$accession)]$description
  best.pfam.hits <- this.ef.associated.reprseq.pfam.hits[, .(qname, sname, qstart, qend, id, descr)]
  setnames(best.pfam.hits, c("name", "domain", "start", "end", "id", "description"))
  setorder(best.pfam.hits, name, start)
  
  ecod.data.outfile <- file.path(subseq.outpath, sprintf("%s-ecod.txt", this.ef))
  pfam.data.outfile <- file.path(subseq.outpath, sprintf("%s-pfam.txt", this.ef))
  fragment.data.outfile <- file.path(subseq.outpath, sprintf("%s-frag.txt", this.ef))
  
  fwrite(best.ecod.hits, file = ecod.data.outfile)
  fwrite(best.pfam.hits, file = pfam.data.outfile)
  fwrite(this.ecf.reprseq.hits, file = fragment.data.outfile)
}



# this.ecf.reprseq.hits <- reprseq.ecfs[sname == this.ecf & scov >= 0.9]
# reprseq.seqs.n %in% this.ecf.reprseq.hits$qname









# rec.domain.overlaps[protein1.ecf.best == this.ecf]



###############################################################################################
###############################################################################################
###############################################################################################

rr.associated.ecfs.metadata <- ecf.class.table[.(rr.associated.ecfs)]
setorder(rr.associated.ecfs.metadata, family)
View(rr.associated.ecfs.metadata)

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
