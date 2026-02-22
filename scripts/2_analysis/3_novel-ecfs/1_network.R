### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(igraph)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

# Load ECF tables
ecf.table <- fread(ecf.table.output)
setkey(ecf.table, ecf.name)
all.ecfs <- unique(ecf.table$ecf.name)
# ecf.lengths <- ecf.table$length

MAX.ECF.LENGTH <- 800
ECF.MIN.COV.THR <- 0.1
ECF.MAX.COV.THR <- 0.9
ecf4network <- ecf.table[length <= MAX.ECF.LENGTH]$ecf.name

# Calculate coverage
ecfs.vs.ecfs.cov <- fread(ecfs.vs.ecfs.cov.output)
ecfs.vs.ecfs.cov$qlength <- ecf.table[.(ecfs.vs.ecfs.cov$query)]$length
ecfs.vs.ecfs.cov$slength <- ecf.table[.(ecfs.vs.ecfs.cov$subject)]$length
ecfs.vs.ecfs.cov[, qcov := pmin(qcov.min, qcov.max)]
ecfs.vs.ecfs.cov[, scov := pmin(scov.min, scov.max)]
ecfs.vs.ecfs.cov[, cov.min := pmin(scov.min, qcov.min)]
ecfs.vs.ecfs.cov[, cov.max := pmax(scov.min, qcov.min)]

ecfs.vs.ecfs.cov.filtered <- ecfs.vs.ecfs.cov[query %chin% ecf4network & subject %chin% ecf4network]
ecfs.vs.ecfs.cov.filtered <- ecfs.vs.ecfs.cov.filtered[cov.min >= ECF.MIN.COV.THR & cov.max >= ECF.MAX.COV.THR]
# ecfs.vs.ecfs.cov <- ecfs.vs.ecfs.cov[!(ecfs.vs.ecfs.cov$qlength*ECF.LENGTH.DIFF.THR < ecfs.vs.ecfs.cov$slength | ecfs.vs.ecfs.cov$qlength/ECF.LENGTH.DIFF.THR > ecfs.vs.ecfs.cov$slength)]
# ecf.network.full.self <- data.table(query = all.ecfs, subject = all.ecfs, prob = 100, pident = 100, hit.length = ecf.lengths$length, qcov = 1, scov = 1, cov.min = 1, cov.max = 1)
# ecf.network.full <- rbind(ecf.network.full.nonself, ecf.network.full.self)

# # filter ECF pairs with an extreme length difference
# filter by min/max coverage

# # create an ECF network
# ecf.network.full.nonself <- ecfs.vs.ecfs.cov %>%
#   select(query, subject, prob, pident, hit.length, qcov, scov, cov.min, cov.max)

ecf.network <- ecfs.vs.ecfs.cov.filtered %>%
  select(query, subject, weight = cov.max)
setkey(ecf.network, query, subject)

ecf.network.outfile <- file.path(ecf.network.dir, "ecf-vis-network.txt")
fwrite(ecf.network, file = ecf.network.outfile, sep = "\t")
