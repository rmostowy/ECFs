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

reprseq.functional.annot <- rbind(functional.annot.phrogs, functional.annot.antidef)
setkey(reprseq.functional.annot, qname)

reprseq.functional.annot.unique <- copy(reprseq.functional.annot)
reprseq.functional.annot.unique[, cov := qcov*scov/100]

setorder(reprseq.functional.annot.unique, qname, -cov, -pident, eval)
reprseq.functional.annot.unique <- reprseq.functional.annot.unique[ , .SD[1], by=qname]

reprseq.functional.annot.unique <- reprseq.functional.annot.unique %>%
  select(reprseq = qname, funct, class)

# TODO
# 1. consider adding annotation indices to reprseq.functional.annot (but what about antidefence?)
# 2. expand reprseq.functional.annot.unique with unknowns
# 3. write them to files






