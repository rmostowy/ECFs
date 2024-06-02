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
  filter((qname %in% all.reprseq)) %>%
  select(qname, sname, qcov, scov, pident, eval, prob)

antidefence.hhpred.hits.filtered <- antidefence.hhpred.hits %>%
  filter(hit.length >= MINIMUM.ALIGNMENT.LENGTH.FOR.ANNOTATION) %>%
  filter(prob >= DEFAULT.MINIMUM.PROB.FOR.ANNOTATION) %>%
  filter(scov >= DEFAULT.MINIMUM.COV.FOR.ANNOTATION) %>%
  filter((sname %in% all.reprseq)) %>%
  select(qname, sname, qcov, scov, pident, eval, prob)

# merge functional information
phrog.hhpred.hits.filtered$phrog <- as.numeric(gsub("phrog_", "", phrog.hhpred.hits.filtered$sname))
phrog.hhpred.hits.filtered.annot <- merge(phrog.hhpred.hits.filtered, phrogs.functions, by = "phrog")
antidefence.hhpred.hits.filtered.annot <- merge(antidefence.hhpred.hits.filtered, antidefence.table, by.x = "qname", by.y = "accession", all.x = TRUE)

functional.annot.phrogs <- phrog.hhpred.hits.filtered.annot %>% 
  select(qname, sname, qcov, scov, pident, eval, prob, funct, class)

functional.annot.antidef <- antidefence.hhpred.hits.filtered.annot %>% 
  select(qname = sname, sname = qname, qcov = scov, scov = qcov, pident, eval, prob, funct = annotation)
functional.annot.antidef$class <- "antidefence"

reprseq.functional.annot.all <- rbind(functional.annot.phrogs, functional.annot.antidef)
# setkey(reprseq.functional.annot, qname)

# identify top annotations
top.phrogs.functions <- phrogs.count[phrogs.count$no.seqs >= MIN.NUM.PHROG.SEQ.PER.ANNOT]

no.reprseq.per.function <- reprseq.functional.annot.all[ , .N, by = funct]
top.functions.by.no.reprseq <- no.reprseq.per.function[no.reprseq.per.function$N >= MIN.NUM.PROT.PER.ANNOT]

# top.functional.hits <- 


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



# TODO
# 1. make a list of most common PHROG functions
# 2. think how to quantify common antidefence functions (maybe there's no point)
# 3. think how to assign unique functions per reprseq -> if by score, what about common functions?






