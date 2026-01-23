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
table.dir <- "tables"
output.dir <- file.path(visual.dir, "part5_rem", "rem-equiv-ecfs")
make.dir(output.dir)

# Load refseq table
reprseq.table <- fread(reprseq.table.output)
setkey(reprseq.table, reprseq)
all.reprseq <- unique(reprseq.table$reprseq)

# Load recent recombination table
rem.table <- fread(rem.table.output)
rem.table[, cov.max.full := pmax(qcov.full, scov.full, na.rm = T)]

# Load recombination domain overlap table
rem.domain.overlap.table.outfile <- file.path(table.dir, "rem-domain-overlaps.txt")
rem.domain.overlaps <- fread(rem.domain.overlap.table.outfile)

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

# Load ECFs vs ECOD
ecfs.vs.ecod <- fread(ecfs.vs.ecod.output)
setkey(ecfs.vs.ecod, qname)

##################################################
### Select REM-Equivalent ECFs (RE-ECFs)
##################################################

# Get ECF recombination overlap table
ECF.OVERLAP.SCORE.THR <- 0.9
ECF.MIN.OVERLAP.SCORE.THR <- 0.9
MIN.NO.PROTEINS.WITH.HIGH.OVERLAPPING.ECF <- 3
RE.ECF.NO.PAIRS.THR <- 20

rows.from <- rem.domain.overlaps %>% select(reprseq = protein1.name, domain = protein1.ecf.best, score = protein1.ecf.score)
rows.to <- rem.domain.overlaps %>% select(reprseq = protein2.name, domain = protein2.ecf.best, score = protein2.ecf.score)
protein.rem.ecf.overlaps <- rbind(rows.from, rows.to)
protein.rem.ecf.overlaps <- protein.rem.ecf.overlaps[is.na(score), score := 0]
protein.rem.ecf.overlaps.nonzero <- protein.rem.ecf.overlaps[score > 0]

# To select recombination-equivalent ECFs, we need to do two things:
# a) consider only ECFs which qualify as RE (REM-equivalent) in at least MIN.NO.PROTEINS.WITH.HIGH.OVERLAPPING.ECF proteins
# b) consider only ECFs which have an average overlap score of at least 0.9

# Hence we first need to count the number of reprseq per ECF with high score
protein.rem.ecf.overlaps.high <- protein.rem.ecf.overlaps[score >= ECF.OVERLAP.SCORE.THR]
ecf.vs.reprseq.high <- protein.rem.ecf.overlaps.high %>% select(domain, reprseq)
ecf.high.overlap.freq <- ecf.vs.reprseq.high[, .N, by = domain]
re.ecf.shortlist <- ecf.high.overlap.freq[N >= MIN.NO.PROTEINS.WITH.HIGH.OVERLAPPING.ECF]$domain

# Then we need to calculate mean overlap for those ECFs for pre-filtering overlap scores
filtered_data <- protein.rem.ecf.overlaps.nonzero[domain %in% re.ecf.shortlist]
mean.overlap.per.ecf <- filtered_data[, .(mean.score = mean(score, na.rm = TRUE)), by = domain]

# Then we obtain RE ECFs
re.ecfs <- sort(mean.overlap.per.ecf[mean.score >= ECF.MIN.OVERLAP.SCORE.THR]$domain)
re.ecf.associated.pairs <- lapply(1:length(re.ecfs), function(index){
  this.ecf <- re.ecfs[index]
  this.ecf.rem.indices <- union(rem.domain.overlaps[protein1.ecf.best == this.ecf]$this.rr.pair.index,
                                rem.domain.overlaps[protein2.ecf.best == this.ecf]$this.rr.pair.index
  )
  this.ecf.rem.indices <- sort(this.ecf.rem.indices)
  this.ecf.rem.indices
})
re.ecf.associated.pairs.n <- sapply(re.ecf.associated.pairs, length)

re.ecfs.metadata <- cbind(ecf.class.table[.(re.ecfs)], no.pairs = re.ecf.associated.pairs.n)
re.ecfs.metadata$ecf.length <- ecf.table[.(re.ecfs)]$length
# setorder(re.ecfs.metadata, family)

# Load reprseq sequences
reprseq.seqs <- readAAStringSet(reprseq.seqs.filename)
reprseq.seqs.n <- names(reprseq.seqs)
reprseq.seqs.l <- nchar(reprseq.seqs)
reprseq.lengths <- data.table(name = reprseq.seqs.n, length = reprseq.seqs.l)
setkey(reprseq.lengths, name)

subseq.outpath <- file.path(output.dir, "ecf-subseqs")
make.dir(subseq.outpath)

re.ecf.common <- re.ecfs.metadata[no.pairs >= RE.ECF.NO.PAIRS.THR]
re.ecf.common.outfile <- file.path(output.dir, "re-ecf-common.txt")
fwrite(re.ecf.common, file = re.ecf.common.outfile)

# ecf2find <- "ecf_00084"
# index <- which(rr.associated.ecfs == ecf2find)
iq.extensions2remove <- c("bionj", "ckp.gz", "iqtree", "log", "mldist")

for(index in 1:nrow(re.ecf.common)){
  cat(index,"\n")
  this.ecf <- re.ecf.common[index]$ecf
  this.ecf.associated.reprseqs <- rbind(rem.domain.overlaps[protein1.ecf.best == this.ecf, .(qname = protein1.name, rem.index)], 
                                        rem.domain.overlaps[protein2.ecf.best == this.ecf, .(qname = protein2.name, rem.index)])
  setorder(this.ecf.associated.reprseqs, qname, rem.index)
  this.ecf.associated.reprseqs.unique <- this.ecf.associated.reprseqs[, .SD[1], by = qname]
  
  this.ecf.reprseq.hits.list <- lapply(1:nrow(this.ecf.associated.reprseqs.unique), function(i){
    this.rem.index <- this.ecf.associated.reprseqs.unique$rem.index[i]
    this.rem.reprseq <- this.ecf.associated.reprseqs.unique$qname[i]
    this.pair.rem.table <- rem.table[this.rem.index]
    is.query <- this.rem.reprseq == this.pair.rem.table$qname
    is.subject <- this.rem.reprseq == this.pair.rem.table$sname
    if(is.query){
      start <- this.pair.rem.table$qstart
      end <- this.pair.rem.table$qend
    } else if(is.subject){
      start <- this.pair.rem.table$sstart
      end <- this.pair.rem.table$send
    } else{
      cat("reprseq missing!\n")
    }
    data.table(name = this.rem.reprseq, start, end)  
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
  
  this.ecf.aln.outfile <- file.path(subseq.outpath, sprintf("%s.aln", this.ecf))
  clustalo.cmd <- paste("clustalo -i", this.ecf.outfile, "-o", this.ecf.aln.outfile, "--force --auto")
  cat("Aligning ECF sequences...\n")
  system(clustalo.cmd)
  cat("Done!\n\n")

  cat("Running IQ-tree...\n")
  iqtree.cmd <- sprintf("iqtree3 -s %s -m LG -T 8 --redo --quiet --fast", this.ecf.aln.outfile)
  system(iqtree.cmd)
  cat("Done!\n\n")
  
  iqfiles2rm <- sprintf("%s.%s", this.ecf.aln.outfile, iq.extensions2remove)
  delete.dir(iqfiles2rm)
  
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
  this.ecf.associated.reprseq.pfam.hits <- reprseq.pfam[qname %in% this.ecf.associated.reprseqs.unique$qname]
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

# 
# 
# 
# ################
# 
# ecf.vs.ecf.cov <- fread(ecfs.vs.ecfs.cov.output)
# 
# 
# # Load reprseq sequences
# reprseq.seqs <- readAAStringSet(reprseq.seqs.filename)
# reprseq.seqs.n <- names(reprseq.seqs)
# reprseq.seqs.l <- nchar(reprseq.seqs)
# reprseq.lengths <- data.table(name = reprseq.seqs.n, length = reprseq.seqs.l)
# setkey(reprseq.lengths, name)
# 
# subseq.outpath <- file.path(output.dir, "EF-subseqs")
# make.dir(subseq.outpath)
# 
# all.efs <- unique(re.ecfs.metadata$family)
# 
# for(index in 1:length(all.efs)){
#   cat(index,"\n")
#   this.ef <- all.efs[index]
#   this.ef.ecfs <- re.ecfs.metadata[family == this.ef]$ecf
#   
#   ecf.vs.ecf.cov[query %chin% this.ef.ecfs & subject %chin% this.ef.ecfs]
# 
#   this.ef.associated.reprseqs <- rbindlist(lapply(this.ef.ecfs, function(ecf){
#     rbind(rec.domain.overlaps[protein1.ecf.best == ecf, .(qname = protein1.name, row.id)],
#           rec.domain.overlaps[protein2.ecf.best == ecf, .(qname = protein2.name, row.id)])
#   }))
#   setorder(this.ef.associated.reprseqs, qname, row.id)
#   this.ef.associated.reprseqs.unique <- this.ef.associated.reprseqs[, .SD[1], by = qname]
# 
# 
#   this.ef.reprseq.hits.list <- lapply(1:nrow(this.ef.associated.reprseqs.unique), function(i){
#     this.row.id <- this.ef.associated.reprseqs.unique$row.id[i]
#     this.pair.reprseq <- this.ef.associated.reprseqs.unique$qname[i]
#     this.pair.rec.table <- recent.rec.table[this.row.id]
#     is.query <- this.pair.reprseq == this.pair.rec.table$qname
#     is.subject <- this.pair.reprseq == this.pair.rec.table$sname
#     if(is.query){
#       start <- this.pair.rec.table$qstart
#       end <- this.pair.rec.table$qend
#     } else if(is.subject){
#       start <- this.pair.rec.table$sstart
#       end <- this.pair.rec.table$send
#     } else{
#       cat("reprseq missing!\n")
#     }
#     data.table(name = this.pair.reprseq, start, end)
#   })
#   this.ef.reprseq.hits <- rbindlist(this.ef.reprseq.hits.list)
#   setkey(this.ef.reprseq.hits, name)
# 
#   this.ef.reprseq.hits$length <- reprseq.lengths[.(this.ef.reprseq.hits$name)]$length
#   this.ef.reprseq.seqs <- reprseq.seqs[reprseq.seqs.n %in% this.ef.reprseq.hits$name]
# 
#   # Extract subsequences using Biostrings' efficient subseq function
#   extracted.subseqs <- this.ef.reprseq.hits[, {
#     # Extract full sequence from AAStringSet
#     full_seq <- as.character(this.ef.reprseq.seqs[[name]])
# 
#     # Extract the sub-sequence using `subseq()`
#     sub_seq <- subseq(full_seq, start = start, end = end)
# 
#     # Return as a list to keep structure
#     .(sub_seq)
#   }, by = name]
# 
#   this.ef.reprseq.subseqs <- AAStringSet(extracted.subseqs$sub_seq)
#   names(this.ef.reprseq.subseqs) <- this.ef.reprseq.hits$name
#   this.ef.outfile <- file.path(subseq.outpath, sprintf("%s.fasta", this.ef))
#   writeXStringSet(this.ef.reprseq.subseqs, this.ef.outfile)
# 
#   # Find all ECOD domains associated with those reprseqs
#   this.ef.associated.reprseq.ecod.hits <- reprseq.ecod[qname %in% this.ef.associated.reprseqs.unique$qname]
#   this.ef.associated.reprseq.ecod.hits.best <- this.ef.associated.reprseq.ecod.hits[
#     order(-scov),  # Sort by scov in descending order
#     .SD[1],  # Select the first (best) hit per group
#     by = .(qname, t.id)
#   ]
#   best.ecod.hits <- this.ef.associated.reprseq.ecod.hits.best[, .(qname, sname, qstart, qend, t.id)]
#   best.ecod.hits$t.name <- ecod.domains[.(best.ecod.hits$sname)]$t_name
#   setnames(best.ecod.hits, c("name", "domain", "start", "end", "id", "description"))
#   setorder(best.ecod.hits, name, start)
# 
#   # Find all PFAM domains associated with those reprseqs
#   this.ef.associated.reprseq.pfam.hits <- reprseq.pfam[qname %in% this.ef.associated.reprseqs.unique$qname & scov >= MIN.SCOV.FOR.DOMAIN]
#   this.ef.associated.reprseq.pfam.hits$id <- pfam.domains[.(this.ef.associated.reprseq.pfam.hits$accession)]$id
#   this.ef.associated.reprseq.pfam.hits$descr <- pfam.domains[.(this.ef.associated.reprseq.pfam.hits$accession)]$description
#   best.pfam.hits <- this.ef.associated.reprseq.pfam.hits[, .(qname, sname, qstart, qend, id, descr)]
#   setnames(best.pfam.hits, c("name", "domain", "start", "end", "id", "description"))
#   setorder(best.pfam.hits, name, start)
# 
#   ecod.data.outfile <- file.path(subseq.outpath, sprintf("%s-ecod.txt", this.ef))
#   pfam.data.outfile <- file.path(subseq.outpath, sprintf("%s-pfam.txt", this.ef))
#   fragment.data.outfile <- file.path(subseq.outpath, sprintf("%s-frag.txt", this.ef))
# 
#   fwrite(best.ecod.hits, file = ecod.data.outfile)
#   fwrite(best.pfam.hits, file = pfam.data.outfile)
#   fwrite(this.ef.reprseq.hits, file = fragment.data.outfile)
# }
# 
# 
# 
# 
# 