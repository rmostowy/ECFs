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
phrog.hhpred.hits <- phrog.hhpred.hits[qname %in% all.reprseq]

# load antidefence tables
antidefence.hhpred.hits <- fread(antidefence.hhpred.hits.filename)
antidefence.hhpred.hits[, hit.length := send - sstart + 1]
antidefence.hhpred.hits[, qcov := round((qend - qstart + 1)/qlength*100)]
antidefence.hhpred.hits[, scov := round(hit.length/slength*100)]
antidefence.hhpred.hits <- antidefence.hhpred.hits[sname %in% all.reprseq]
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
functional.annot.phrogs <- functional.annot.phrogs[!functional.annot.phrogs$funct %in% UNSPECIFIC.ANNOTATIONS]

functional.annot.antidef <- antidefence.hhpred.hits.filtered.annot %>% 
  select(qname = sname, sname = qname, qcov = scov, scov = qcov, pident, eval, prob, funct = annotation)
functional.annot.antidef$class <- "antidefence"
functional.annot.antidef$annotation <- sprintf("%s|%s", functional.annot.antidef$funct, functional.annot.antidef$class)

reprseq.functional.annot.all <- rbind(functional.annot.phrogs, functional.annot.antidef)
reprseq.functional.annot.all[, score := qcov*scov/10000]
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

reprseq.functional.annot.all$top.annotation <- FALSE
reprseq.functional.annot.all$top.annotation[reprseq.functional.annot.all$annotation %in% top.annotations] <- TRUE
reprseq.functional.annot.all$top.annotation[reprseq.functional.annot.all$annotation %in% ad.annotations] <- TRUE
setkey(reprseq.functional.annot.all, qname)

unique.annotations.count <- reprseq.functional.annot.all[, .(no.unique.annotations = uniqueN(funct)), by = qname]

reprseq.functional.annot.best <- copy(reprseq.functional.annot.all)
setorder(reprseq.functional.annot.best, qname, -score, -pident, eval)
reprseq.functional.annot.best <- reprseq.functional.annot.best[ , .SD[1], by = qname]

reprseq.functional.annot.best.full <- merge(reprseq.functional.annot.best, unique.annotations.count, by = "qname")

reprseq.functional.annot.best.full.sub <- reprseq.functional.annot.best.full %>%
  select(reprseq = qname, score, funct, class, annotation, top.annotation, no.unique.annotations)

reprseq.functional.annot.unique <- merge(data.table(reprseq = all.reprseq), reprseq.functional.annot.best.full.sub, by = "reprseq", all.x = T)
reprseq.functional.annot.unique$no.unique.annotations[is.na(reprseq.functional.annot.unique$no.unique.annotations)] <- 0

fwrite(reprseq.functional.annot.all, file = reprseq.functions.all.output)
fwrite(reprseq.functional.annot.unique, file = reprseq.functions.unique.output)

