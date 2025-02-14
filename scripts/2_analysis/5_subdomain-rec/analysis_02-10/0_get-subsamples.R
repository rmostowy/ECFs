### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
# suppressWarnings(suppressPackageStartupMessages(library(seqinr)))
# suppressWarnings(suppressPackageStartupMessages(library(Biostrings)))
# suppressWarnings(suppressPackageStartupMessages(library(parallel)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################
set.seed(1234)

recent.rec.table <- fread(reprseq.recent.rec.table.output)
n <- nrow(recent.rec.table)

prop.subsample <- 0.1
no.rows2sample <- floor(n*prop.subsample)
rows2sample <- sort(sample(n, no.rows2sample))
recent.rec.table.sub <- recent.rec.table[rows2sample]
fwrite(recent.rec.table.sub, file = "rececent-rec-table-sub.txt")

prop.subsample <- 0.2
no.rows2sample <- floor(n*prop.subsample)
rows2sample <- sort(sample(n, no.rows2sample))
recent.rec.table.sub <- recent.rec.table[rows2sample]
fwrite(recent.rec.table.sub, file = "rececent-rec-table-sub-20p.txt")

prop.subsample <- 0.3
no.rows2sample <- floor(n*prop.subsample)
rows2sample <- sort(sample(n, no.rows2sample))
recent.rec.table.sub <- recent.rec.table[rows2sample]
fwrite(recent.rec.table.sub, file = "rececent-rec-table-sub-30p.txt")

prop.subsample <- 0.4
no.rows2sample <- floor(n*prop.subsample)
rows2sample <- sort(sample(n, no.rows2sample))
recent.rec.table.sub <- recent.rec.table[rows2sample]
fwrite(recent.rec.table.sub, file = "rececent-rec-table-sub-40p.txt")

prop.subsample <- 0.5
no.rows2sample <- floor(n*prop.subsample)
rows2sample <- sort(sample(n, no.rows2sample))
recent.rec.table.sub <- recent.rec.table[rows2sample]
fwrite(recent.rec.table.sub, file = "rececent-rec-table-sub-50p.txt")

prop.subsample <- 0.7
no.rows2sample <- floor(n*prop.subsample)
rows2sample <- sort(sample(n, no.rows2sample))
recent.rec.table.sub <- recent.rec.table[rows2sample]
fwrite(recent.rec.table.sub, file = "rececent-rec-table-sub-70p.txt")
