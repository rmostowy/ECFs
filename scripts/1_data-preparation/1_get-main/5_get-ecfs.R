### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

#########################
# 1. ECF name conversion
#########################

# load input data
reprseq.table <- fread(reprseq.table.filename)
all.reprseq <- unique(reprseq.table$reprseq)

# extract representative sequences for each ecf (Krzysztof's code is now wrong)
ecf.a3m.list <- list.files(ecfs.a3m.dir, full.names = T, recursive = F)
ecf.a3m.repr <- sapply(ecf.a3m.list, function(file) readLines(file, n = 1))
ecf.a3m.repr <- gsub(">", "", ecf.a3m.repr)
ecf.a3m.list.basename <- basename(ecf.a3m.list)
ecf.name.table <- data.table(ecf.name = sub("_final\\.a3m$", "", ecf.a3m.list.basename),
                             er.name = ecf.a3m.repr)
# Up until now, the code we used was the following (commented out); 
# we can reintroduce it when Krzysztof fixes his
# # load ECF name table
# ecf.name.table <- fread(ecf.name.table.filename, header = F)
# setnames(ecf.name.table, c("er.name", "ecf.name"))
# ecf.name.table[, c("reprseq", "suffix") := tstrsplit(er.name, "_")]
# ecf.name.table[, suffix := NULL]
ecf.name.table[, c("reprseq", "suffix") := tstrsplit(er.name, "_")]
ecf.name.table[, suffix := NULL]
ecf.name.table <- ecf.name.table[reprseq %chin% all.reprseq]
fwrite(ecf.name.table, file = ecf.name.table.output)
setkey(ecf.name.table, er.name)

# henceforth, ECFs created only from qECFs from corrupted reprseq are removed
# (even if they were found subsequently in other reprseqs)
all.ecfs <- unique(ecf.name.table$ecf.name)
Ne <- length(all.ecfs)

# Read the file efficiently
ecf.vs.qecf.raw <- fread(ecf.vs.qecf.filename, header = FALSE, sep = "\n")
# Separate the 'ecf' and 'qcfs' parts
ecf.vs.qecf.raw[, c("ecf", "qcfs") := tstrsplit(V1, ",")]
# Split the 'qcfs' part into individual elements and process
ecf.vs.qecf <- ecf.vs.qecf.raw[, .(qecf = unlist(strsplit(qcfs, " "))), by = ecf]
ecf.vs.qecf[, qecf := gsub("ecf_", "", qecf)]
# Remove all ECFs which are not in the 'all.ecfs' vector
ecf.vs.qecf <- ecf.vs.qecf[ecf %in% all.ecfs]
ecf.vs.qecf[, c("reprseq", "suffix") := tstrsplit(qecf, "_")]
ecf.vs.qecf[, suffix := NULL]
fwrite(ecf.vs.qecf, file = ecf.vs.qecf.output)

all.ecfs <- unique(ecf.vs.qecf$ecf)

#####################################
# 2. Process reprseq to ECF mapping
#####################################

# load reprseq to ecf mapping
reprseq.vs.ecfs.all <- fread(ecf.hits.filename)
reprseq.vs.ecfs.all[, hit.length := qend - qstart + 1]
reprseq.vs.ecfs.all[, qcov := round((qend - qstart + 1)/qlength, 4)]
reprseq.vs.ecfs.all[, scov := round((send - sstart + 1)/slength, 4)]
# remove fishy reprseqs and ecfs
reprseq.vs.ecfs.all <- reprseq.vs.ecfs.all[qname %chin% all.reprseq]
reprseq.vs.ecfs.all <- reprseq.vs.ecfs.all[sname %chin% ecf.name.table$er.name]
reprseq.vs.ecfs.all$sname <- ecf.name.table[.(reprseq.vs.ecfs.all$sname)]$ecf.name


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
# 2. Process ECF to ECF mapping
#####################################

ecfs.vs.ecfs.all <- fread(ecf.hits.vs.ecf.filename)
ecfs.vs.ecfs.all[, hit.length := qend - qstart + 1]
ecfs.vs.ecfs.all[, qcov := round((qend - qstart + 1)/qlength, 4)]
ecfs.vs.ecfs.all[, scov := round((send - sstart + 1)/slength, 4)]
ecfs.vs.ecfs.all[, qname := gsub("_final", "", qname)]
ecfs.vs.ecfs.all$sname <- ecf.name.table[.(ecfs.vs.ecfs.all$sname)]$ecf.name
ecfs.vs.ecfs.all <- ecfs.vs.ecfs.all[!is.na(sname)]
ecfs.vs.ecfs.all <- ecfs.vs.ecfs.all[qname %chin% all.ecfs]

# filter ECF-to-ECF hits
ecfs.vs.ecfs <- ecfs.vs.ecfs.all %>%
  filter(eval <= MAX.EVAL.FOR.ECF) %>% 
  filter(prob >= MIN.PROB.FOR.ECF) %>%
  filter(hit.length >= MIN.ECF.HIT.LENGTH)

# write to file
fwrite(ecfs.vs.ecfs, file = ecfs.vs.ecfs.output)

#####################################
# 3. Process ECF to ECOD mapping
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
# 4. Process ECF to PFAM mapping
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
# 5. Calculate ECF qcov/scov
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


