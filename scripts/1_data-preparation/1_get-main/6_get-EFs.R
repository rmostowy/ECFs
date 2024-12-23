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

# Load ECF table
ecfs.vs.qecfs <- fread(ecf.vs.qecf.output)
ecfs.vs.qecfs[, c("reprseq", "suffix") := tstrsplit(qecf, "_")]
ecfs.vs.qecfs[, suffix := NULL]
setkey(ecfs.vs.qecfs, ecf)
all.ecfs <- unique(ecfs.vs.qecfs$ecf)

# Calculate ECF lengths from a3m files
ecf.a3m.list <- list.files(ecfs.a3m.dir, full.names = T, recursive = F)
ecf.a3m.list.basename <- basename(ecf.a3m.list)
ecf.a3m.seqs <- sapply(ecf.a3m.list, function(file) readLines(file, n = 2)[2])
ecf.lengths <- data.table(name = sub("_final\\.a3m$", "", ecf.a3m.list.basename),
                             length = sapply(ecf.a3m.seqs, nchar))
ecf.lengths <- ecf.lengths[name %in% all.ecfs]
setkey(ecf.lengths, name)
fwrite(ecf.lengths, file = ecf.lengths.output)

# Calculate coverage
ecfs.vs.ecfs.cov <- fread(ecfs.vs.ecfs.cov.output)
ecfs.vs.ecfs.cov$qlength <- ecf.lengths[.(ecfs.vs.ecfs.cov$query)]$length
ecfs.vs.ecfs.cov$slength <- ecf.lengths[.(ecfs.vs.ecfs.cov$subject)]$length
ecfs.vs.ecfs.cov[, scov := pmin(scov.min, scov.max)]
ecfs.vs.ecfs.cov[, qcov := pmin(qcov.min, qcov.max)]
ecfs.vs.ecfs.cov[, cov.min := pmin(scov.min, qcov.min)]
ecfs.vs.ecfs.cov[, cov.max := pmax(scov.min, qcov.min)]

# filter ECF pairs with an extreme length difference
ecfs.vs.ecfs.cov <- ecfs.vs.ecfs.cov[!(ecfs.vs.ecfs.cov$qlength*ECF.LENGTH.DIFF.THR < ecfs.vs.ecfs.cov$slength | ecfs.vs.ecfs.cov$qlength/ECF.LENGTH.DIFF.THR > ecfs.vs.ecfs.cov$slength)]
# filter by min/max coverage
ecfs.vs.ecfs.cov <- ecfs.vs.ecfs.cov[ecfs.vs.ecfs.cov$cov.min >= ECF.MIN.COV.THR & ecfs.vs.ecfs.cov$cov.max >= ECF.MAX.COV.THR]

# create an ECF network
ecf.network.full.nonself <- ecfs.vs.ecfs.cov %>% 
  select(query, subject, prob, pident, qcov, scov, cov.min, cov.max)
ecf.network.full.self <- data.table(query = all.ecfs, subject = all.ecfs, prob = 100, pident = 100, qcov = 1, scov = 1, cov.min = 1, cov.max = 1)
ecf.network.full <- rbind(ecf.network.full.nonself, ecf.network.full.self)

# create a MCL input from the network
ecf.network.mcl <- ecf.network.full %>% select(from = query, to = subject)

# create a temp mcl directory
mcl.run.dir <- file.path(data.tmp.dir, "mcl")
mcl.tmp.dir <- file.path(mcl.run.dir, "tmp")
make.dir(mcl.run.dir)
make.dir(mcl.tmp.dir)

# run MCL
input.filename <- file.path(mcl.run.dir, "mcl-input.txt")
fwrite(ecf.network.mcl, file = input.filename, sep = " ")
output.filename <- file.path(mcl.run.dir, "mcl-output.txt")
run.mcl(input.filename, output.filename, inflation.factor = 2, tmp.dir = mcl.tmp.dir)

# process MCL output
mcl.output <- readLines(output.filename)
Nf <- length(mcl.output)

mcl.output.dt <- data.table(index = 1:Nf, mcl_output = mcl.output)
mcl.vs.ecf <- mcl.output.dt[, .(ecf = unlist(tstrsplit(mcl_output, "\t", fixed = TRUE))), by = index]
# Create family names
mcl.vs.ecf[, ecf.family := sprintf("EF%.*d", nchar(Nf), index)]
ecf.family.table <- mcl.vs.ecf[, .(ecf.family, ecf)]

# Write to file
fwrite(ecf.family.table, file = ecf.family.table.output)
fwrite(ecf.network.full, file = ecf.family.network.output)

