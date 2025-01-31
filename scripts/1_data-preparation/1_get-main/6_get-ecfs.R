### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(Biostrings)))
suppressWarnings(suppressPackageStartupMessages(library(stringr)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################


#########################
# 1. Load reprseq and ECF data
#########################

# load reprseq data
reprseq.table <- fread(reprseq.table.filename)
all.reprseq <- unique(reprseq.table$reprseq)

# make an ECF table
ecf.a3m.list <- list.files(ecfs.a3m.dir, full.names = T, recursive = F)
ecf.full.table.list <- lapply(1:length(ecf.a3m.list), function(ecf.index){
  this.ecf.file.path <- ecf.a3m.list[ecf.index]
  this.ecf.file.name <- basename(this.ecf.file.path)
  this.ecf.name = sub("_final\\.a3m$", "", this.ecf.file.name)
  this.ecf.a3m <- readAAStringSet(this.ecf.file.path)
  this.ecf.repr <- names(this.ecf.a3m)[1]
  this.ecf.seqs <- as.character(this.ecf.a3m)
  this.ecf.length <- nchar(this.ecf.seqs[1])
  this.ecf.size <- length(this.ecf.seqs)
  this.ecf.total.seq.length <- sum(nchar(this.ecf.seqs))
  this.ecf.total.no.letters <- sum(nchar(gsub("-", "", this.ecf.seqs)))
  
  # vector which counts the number of X characters per sequence in the alignment
  no.x.per.seq <- str_count(this.ecf.seqs, "X")
  # total number of sequencing with X
  no.seqs.with.x <- lw(no.x.per.seq>0)
  # maximum X count per any sequence
  x.count.max <- max(no.x.per.seq)
  # X count in the reference sequence
  x.count.ref <- no.x.per.seq[1]
  # total X count
  x.count.total <- sum(no.x.per.seq)
  
  data.table(ecf.name = this.ecf.name, er.name = this.ecf.repr,
             length = this.ecf.length, size = this.ecf.size, 
             total.seq.length = this.ecf.total.seq.length, total.no.aa = this.ecf.total.no.letters - x.count.total,
             no.seqs.with.x, x.count.max, x.count.ref, x.count.total)
})
ecf.full.table <- rbindlist(ecf.full.table.list)

# extract and edit name table
ecf.name.table <- ecf.full.table[, .(ecf.name, er.name)]
ecf.name.table[, c("reprseq", "suffix") := tstrsplit(er.name, "_")]
ecf.name.table[, suffix := NULL]
setkey(ecf.name.table, er.name)
# ecf.name.table <- ecf.name.table[reprseq %chin% all.reprseq]
# fwrite(ecf.name.table, file = ecf.name.table.output)

# load reprseq to ecf mapping
reprseq.vs.ecfs.all <- fread(ecf.hits.filename)
reprseq.vs.ecfs.all[, hit.length := qend - qstart + 1]
reprseq.vs.ecfs.all[, qcov := round((qend - qstart + 1)/qlength, 4)]
reprseq.vs.ecfs.all[, scov := round((send - sstart + 1)/slength, 4)]
reprseq.vs.ecfs.all$sname <- ecf.name.table[.(reprseq.vs.ecfs.all$sname)]$ecf.name
reprseq.vs.ecfs.all <- reprseq.vs.ecfs.all[qname %in% all.reprseq]

# filter reprseq vs ecf data
reprseq.vs.ecfs.all.filtered <- reprseq.vs.ecfs.all %>%
  filter(eval <= MAX.EVAL.FOR.ECF.FILTERING) %>% 
  filter(prob >= MIN.PROB.FOR.ECF.FILTERING) %>%
  filter(scov >= MIN.SCOV.FOR.ECF.FILTERING) %>%
  filter(hit.length >= MIN.HIT.LENGTH.FOR.ECF.FILTERING)


#########################
# 2. Filter ECFs
#########################

# A) Identify ECFs with more than 10 unknown characters 
# in the reference sequence (this requirement should
# become redundant once Krzysztof fixes the pipeline)
has.excessive.x.in.ref <- ecf.full.table$x.count.ref > ACCEPTED.NUM.UNKNOWN.AA

# B) Has an informative alignment
# MIN.PERCENTAGE.OF.INFORMATIVE.BASES <- 0.5
# MIN.NO.REPRSEQ.PER.ECF <- 1
ecf.full.table[, prop.informative.bases := total.no.aa/total.seq.length]
has.informative.alignment <- ecf.full.table$prop.informative.bases >= MIN.PERCENTAGE.OF.INFORMATIVE.BASES

# C) Identify ECFs which do not hit any reprseq
ecf.reprseq.map <- reprseq.vs.ecfs.all.filtered[, .(ecf = sname, reprseq = qname)]
reprseq.per.ecf.table <- ecf.reprseq.map[, .N, by = ecf]
reprseq.per.ecf.table.filtered <- reprseq.per.ecf.table[N >= MIN.NO.REPRSEQ.PER.ECF]
is.found.in.enough.reprseqs <- ecf.full.table$ecf.name %in% reprseq.per.ecf.table.filtered$ecf

# Make a list of all ECFs which pass these filters
ecf.passes.filter <- !has.excessive.x.in.ref & has.informative.alignment & is.found.in.enough.reprseqs
ecf.full.table.filtered <- ecf.full.table[ecf.passes.filter]
all.ecfs <- ecf.full.table.filtered$ecf.name
N <- length(all.ecfs)

ecf.full.table.filtered[ , prop.informative.bases := NULL]
fwrite(ecf.full.table.filtered, file = ecf.table.output)

#####################################
# 3. Process reprseq to ECF mapping
#####################################

# remove fishy reprseqs and ecfs
reprseq.vs.ecfs.all <- reprseq.vs.ecfs.all[qname %chin% all.reprseq]
reprseq.vs.ecfs.all <- reprseq.vs.ecfs.all[sname %chin% all.ecfs]

# filter ECF hits
reprseq.vs.ecfs <- reprseq.vs.ecfs.all %>%
  filter(eval <= MAX.EVAL.FOR.DOMAIN) %>% 
  filter(prob >= MIN.PROB.FOR.DOMAIN) %>%
  filter(qcov >= MIN.QCOV.FOR.DOMAIN) %>%
  filter(scov >= MIN.SCOV.FOR.DOMAIN) %>%
  filter(hit.length >= MIN.DOMAIN.LENGTH)

# write to file
fwrite(reprseq.vs.ecfs, file = reprseq.ecfs.output)

#####################################
# 4. Process ECF to ECF mapping
#####################################

ecfs.vs.ecfs.all <- fread(ecf.hits.vs.ecf.filename)
ecfs.vs.ecfs.all[, hit.length := qend - qstart + 1]
ecfs.vs.ecfs.all[, qcov := round((qend - qstart + 1)/qlength, 4)]
ecfs.vs.ecfs.all[, scov := round((send - sstart + 1)/slength, 4)]
ecfs.vs.ecfs.all[, qname := gsub("_final", "", qname)]
ecfs.vs.ecfs.all$sname <- ecf.name.table[.(ecfs.vs.ecfs.all$sname)]$ecf.name
ecfs.vs.ecfs.all <- ecfs.vs.ecfs.all[!is.na(sname)]
ecfs.vs.ecfs.all <- ecfs.vs.ecfs.all[qname %chin% all.ecfs]
ecfs.vs.ecfs.all <- ecfs.vs.ecfs.all[sname %chin% all.ecfs]

# filter ECF-to-ECF hits
ecfs.vs.ecfs <- ecfs.vs.ecfs.all %>%
  filter(eval <= MAX.EVAL.FOR.ECF) %>% 
  filter(prob >= MIN.PROB.FOR.ECF) %>%
  filter(hit.length >= MIN.ECF.HIT.LENGTH | qname == sname) # don't filter by hit length if it's a self-hit

# write to file
fwrite(ecfs.vs.ecfs, file = ecfs.vs.ecfs.output)


#####################################
# 5. Process ECF to ECOD mapping
#####################################

ecfs.vs.ecod.all <- fread(ecf.hits.vs.ecod.filename)
ecfs.vs.ecod.all[, hit.length := qend - qstart + 1]
ecfs.vs.ecod.all[, qcov := round((qend - qstart + 1)/qlength, 4)]
ecfs.vs.ecod.all[, scov := round((send - sstart + 1)/slength, 4)]
ecfs.vs.ecod.all[, qname := gsub("_final", "", qname)]
ecfs.vs.ecod.all <- ecfs.vs.ecod.all[qname %chin% all.ecfs]

# filter ECF-to-ECF hits
ecfs.vs.ecod <- ecfs.vs.ecod.all %>%
  filter(eval <= MAX.EVAL.FOR.ECF) %>% 
  filter(prob >= MIN.PROB.FOR.ECF) %>%
  filter(hit.length >= MIN.ECF.HIT.LENGTH)


# write to file
fwrite(ecfs.vs.ecod, file = ecfs.vs.ecod.output)

#####################################
# 6. Process ECF to PFAM mapping
#####################################

ecfs.vs.pfam.all <- fread(ecf.hits.vs.pfam.filename)
ecfs.vs.pfam.all[, hit.length := qend - qstart + 1]
ecfs.vs.pfam.all[, qcov := round((qend - qstart + 1)/qlength, 4)]
ecfs.vs.pfam.all[, scov := round((send - sstart + 1)/slength, 4)]
ecfs.vs.pfam.all[, qname := gsub("_final", "", qname)]
ecfs.vs.pfam.all <- ecfs.vs.pfam.all[qname %chin% all.ecfs]

# filter ECF-to-ECF hits
ecfs.vs.pfam <- ecfs.vs.pfam.all %>%
  filter(eval <= MAX.EVAL.FOR.ECF) %>% 
  filter(prob >= MIN.PROB.FOR.ECF) %>%
  filter(hit.length >= MIN.ECF.HIT.LENGTH)

# write to file
fwrite(ecfs.vs.pfam, file = ecfs.vs.pfam.output)


#####################################
# 7. Calculate ECF qcov/scov
#####################################

use.cores <- 8
no.iterations <- 20
ecfs.vs.ecfs.cov <- ecfs.vs.ecfs %>% filter(prob >= MIN.PROB.FOR.ECF.QCOV.SCOV)

pair.table <- get.all.pairs(ecfs.vs.ecfs.cov)
no.pairs <- nrow(pair.table)

pairs.per.iteration <- floor(no.pairs/no.iterations)
start.end.per.iteration <- lapply(1:no.iterations, function(index){
  end.pair <- pairs.per.iteration*index
  start.pair <- end.pair - pairs.per.iteration + 1
  if(index == no.iterations) end.pair <- no.pairs
  return(c(start.pair, end.pair))
})

pb <- txtProgressBar(min = 0, max = no.iterations, style = 3, initial = 0)
pair.table.hits.list <- lapply(1:no.iterations, function(iteration.index){
  this.iteration.pairs.index <- start.end.per.iteration[[iteration.index]][1]:start.end.per.iteration[[iteration.index]][2]
  this.iteration.pairs <- pair.table[this.iteration.pairs.index]
  this.iteration.hits <- get.pair.table.hits(this.iteration.pairs, ecfs.vs.ecfs.cov, no.cores = use.cores)
  this.iteration.output <- cbind(this.iteration.pairs, this.iteration.hits)
  setTxtProgressBar(pb, iteration.index)
  return(this.iteration.output)
})
close(pb)

pair.table.hits <- rbindlist(pair.table.hits.list)
new.col.names <- c("query", "subject", "prob", "scov.min", "scov.max", "qcov.min", "qcov.max", "pident", "hit.length")
setnames(pair.table.hits, new = new.col.names)
pair.table.hits[, prob := prob * 100]
pair.table.hits[, pident := pident * 100]

setorder(pair.table.hits, query, subject)

fwrite(pair.table.hits, file = ecfs.vs.ecfs.cov.output)


