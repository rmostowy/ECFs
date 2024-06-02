### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(readxl)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

# load reprseq tables
reprseq.table <- fread(reprseq.table.filename)
all.reprseq <- unique(reprseq.table$reprseq)

# load phrogs tables
phrogs.functions <- fread(phrogs.functions.output)
phrogs.count <- fread(phrogs.count.output)
phrog.hhpred.hits <- fread(phrog.hhpred.hits.filename)
phrog.hhpred.hits[, hit.length := qend - qstart + 1]
phrog.hhpred.hits[, qcov := round(hit.length/qlength*100)]
phrog.hhpred.hits[, scov := round((send - sstart + 1)/slength*100)]

# load antidefence tables
antidefence.hhpred.hits <- fread(antidefence.hhpred.hits.filename)
antidefence.hhpred.hits[, hit.length := send - sstart + 1]
antidefence.hhpred.hits[, qcov := round((qend - qstart + 1)/qlength*100)]
antidefence.hhpred.hits[, scov := round(hit.length/slength*100)]
antidefence.table <- fread(antidefence.table.output)

# filter mapping tables
phrog.hhpred.hits.filtered <- phrog.hhpred.hits %>%
  filter(hit.length >= MINIMUM.ALIGNMENT.LENGTH.FOR.ANNOTATION) %>%
  filter(prob >= DEFAULT.MINIMUM.PROB.FOR.ANNOTATION) %>%
  filter(qcov >= DEFAULT.MINIMUM.COV.FOR.ANNOTATION) %>%
  filter(scov >= DEFAULT.MINIMUM.COV.FOR.ANNOTATION) %>%
  filter((qname %in% all.reprseq)) %>%
  select(qname, sname, qcov, scov, pident, eval, prob)

antidefence.hhpred.hits.filtered <- antidefence.hhpred.hits %>%
  filter(hit.length >= MINIMUM.ALIGNMENT.LENGTH.FOR.ANNOTATION) %>%
  filter(prob >= DEFAULT.MINIMUM.PROB.FOR.ANNOTATION) %>%
  filter(qcov >= DEFAULT.MINIMUM.COV.FOR.ANNOTATION) %>%
  filter(scov >= DEFAULT.MINIMUM.COV.FOR.ANNOTATION) %>%
  filter((sname %in% all.reprseq)) %>%
  select(qname, sname, qcov, scov, pident, eval, prob)

# merge functional information
phrog.hhpred.hits.filtered$phrog <- as.numeric(gsub("phrog_", "", phrog.hhpred.hits.filtered$sname))
phrog.hhpred.hits.filtered.annot <- merge(phrog.hhpred.hits.filtered, phrogs.functions, by = "phrog")
antidefence.hhpred.hits.filtered.annot <- merge(antidefence.hhpred.hits.filtered, antidefence.table, by.x = "qname", by.y = "accession", all.x = TRUE)

functional.annot.phrogs <- phrog.hhpred.hits.filtered.annot %>% 
  select(qname, sname, qcov, scov, pident, eval, prob, funct, class)
functional.annot.phrogs$annotation <- sprintf("%s|%s", functional.annot.phrogs$funct, functional.annot.phrogs$class)

functional.annot.antidef <- antidefence.hhpred.hits.filtered.annot %>% 
  select(qname = sname, sname = qname, qcov = scov, scov = qcov, pident, eval, prob, funct = annotation)
functional.annot.antidef$class <- "antidefence"
functional.annot.antidef$annotation <- sprintf("%s|%s", functional.annot.antidef$funct, functional.annot.antidef$class)

reprseq.functional.annot.all <- rbind(functional.annot.phrogs, functional.annot.antidef)
reprseq.functional.annot.all[, score := qcov*scov/100]
# setkey(reprseq.functional.annot, qname)

# identify top PHROGs functional annotations:
# 1) linked to PHROGs with a minimum number of sequences MIN.NUM.PHROG.SEQ.PER.ANNOT
# 2) linked to a minimum number of MIN.NUM.PROT.PER.ANNOT reprseq

functions.largest.phrogs <- phrogs.count[phrogs.count$no.seqs >= MIN.NUM.PHROG.SEQ.PER.ANNOT]
functions.freq.reprseq <- functional.annot.phrogs[, .N, by = annotation]
functions.freq.reprseq <- functions.freq.reprseq[functions.freq.reprseq$N >= MIN.NUM.PROT.PER.ANNOT]
top.annotations <- intersect(functions.largest.phrogs$annotation, functions.freq.reprseq$annotation)
non.top.annotations <- phrogs.count$annotation[!phrogs.count$annotation %in% top.annotations]
ad.annotations <- unique(functional.annot.antidef$annotation)

reprseq.functional.annot.all$annotation.cat <- "NT"
reprseq.functional.annot.all$annotation.cat[reprseq.functional.annot.all$annotation %in% top.annotations] <- "T"
reprseq.functional.annot.all$annotation.cat[reprseq.functional.annot.all$annotation %in% ad.annotations] <- "AD"
setkey(reprseq.functional.annot.all, qname)

unique.reprseq <- unique(reprseq.functional.annot.all$qname)

index=227
this.reprseq <- unique.reprseq[index]
reprseq.functional.annot.all[.(this.reprseq)]




# top.annotations.reprseq <- unique(reprseq.functional.annot.all$qname[reprseq.functional.annot.all$annotation %in% top.annotations])
# non.top.annotations.reprseq <- unique(reprseq.functional.annot.all$qname[reprseq.functional.annot.all$annotation %in% non.top.annotations])
# ad.annotations.reprseq <- unique(reprseq.functional.annot.all$qname[reprseq.functional.annot.all$annotation %in% ad.annotations])
# 
# x <- sort(intersect(intersect(top.annotations.reprseq, non.top.annotations.reprseq), ad.annotations.reprseq))
# 
# this.x <- reprseq.functional.annot.all[reprseq.functional.annot.all$qname == x[6]]
# setorder(this.x, -score)
# this.x




# functional.annot.phrogs.top <- functional.annot.phrogs[functional.annot.phrogs$annotation %in% top.annotations]









# 
# 
# functional.annot.phrogs
# 
# 
# no.reprseq.per.function <- reprseq.functional.annot.all[ , .N, by = funct]
# top.functions.by.no.reprseq <- no.reprseq.per.function[no.reprseq.per.function$N >= MIN.NUM.PROT.PER.ANNOT]
# 
# # top.functional.hits <- 


# choose the strongest annotation per reprseq
reprseq.functional.annot.best <- copy(reprseq.functional.annot.all)
reprseq.functional.annot.best[, cov := qcov*scov/100]

setorder(reprseq.functional.annot.best, qname, -cov, -prob, -pident, eval)
reprseq.functional.annot.best <- reprseq.functional.annot.best[ , .SD[1], by = qname]

# Add a column to indicate if the annotation was unique or not
reprseq.functional.annot.all[, num.annotations := .N, by = qname]
reprseq.functional.annot.best <- reprseq.functional.annot.best %>%
  merge(reprseq.functional.annot.all[, .(qname, num.annotations)], by = "qname")

reprseq.functional.annot.best.sub <- reprseq.functional.annot.best %>%
  select(reprseq = qname, funct, class, num.annotations)

reprseq.functional.annot.unique <- merge(data.table(reprseq = all.reprseq), reprseq.functional.annot.best.sub, by = "reprseq", all.x = T)

# How to assign function?
# BY DEFAULT
# 1. Rank all functions per reprseq by score (eg, by scov*qcov)
# 2. Pick top function
# ----
# FOR TOP PHROG FUNCTIONS
# 1. Rank all functions per reprseq by score (eg, by scov*qcov)
# 2. If functions are only non-top PHROGs, assign to "rare-phrogs"
# 3. If functions are only top PHROGs or only top AD, no problem.
# 4. If functions are (top and nontop) OR (AD and nontop), also no problem
# 5. 

# T, NT, AD -> pick best (T, AD)
# T, NT -> T
# T, AD -> pick best
# NT, AD -> AD
# T -> T
# NT -> rare
# AD -> AD






# TODO
# 1. make a list of most common PHROG functions
# 2. think how to quantify common antidefence functions (maybe there's no point)
# 3. think how to assign unique functions per reprseq -> if by score, what about common functions?






