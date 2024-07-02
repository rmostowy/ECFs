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
ecfs.vs.qecfs$reprseq <- substr(ecfs.vs.qecfs$qecf, 1, regexpr("_", ecfs.vs.qecfs$qecf)-1)
setkey(ecfs.vs.qecfs, ecf)
all.ecfs <- unique(ecfs.vs.qecfs$ecf)

# Extract ECF lengths from ECF-ECF comparison
ecfs.vs.ecfs.all <- fread(ecfs.vs.ecfs.output)
ecf.lengths <- rbind(
  ecfs.vs.ecfs.all %>% select(name = qname, length = qlength),
  ecfs.vs.ecfs.all %>% select(name = sname, length = slength)
)
ecf.lengths <- unique(ecf.lengths)
setorder(ecf.lengths, name)
setkey(ecf.lengths, name)

# Calculate coverage
ecfs.vs.ecfs.cov <- fread(ecfs.vs.ecfs.cov.output)
ecfs.vs.ecfs.cov$qlength <- ecf.lengths[.(ecfs.vs.ecfs.cov$query)]$length
ecfs.vs.ecfs.cov$slength <- ecf.lengths[.(ecfs.vs.ecfs.cov$subject)]$length
ecfs.vs.ecfs.cov[, cov.min := pmin(scov.min, qcov.min)]
ecfs.vs.ecfs.cov[, cov.max := pmax(scov.min, qcov.min)]

# filter ECF pairs with an extreme length difference
ecfs.vs.ecfs.cov <- ecfs.vs.ecfs.cov[!(ecfs.vs.ecfs.cov$qlength*ECF.LENGTH.DIFF.THR < ecfs.vs.ecfs.cov$slength | ecfs.vs.ecfs.cov$qlength/ECF.LENGTH.DIFF.THR > ecfs.vs.ecfs.cov$slength)]
# filter by min/max coverage
ecfs.vs.ecfs.cov <- ecfs.vs.ecfs.cov[ecfs.vs.ecfs.cov$cov.min >= ECF.MIN.COV.THR & ecfs.vs.ecfs.cov$cov.max >= ECF.MAX.COV.THR]

# create a MCL input
ecf.network <- data.table(from = c(ecfs.vs.ecfs.cov$query, all.ecfs),
                              to = c(ecfs.vs.ecfs.cov$subject, all.ecfs)
)

# create a temp mcl directory
mcl.run.dir <- file.path(data.tmp.dir, "mcl")
mcl.tmp.dir <- file.path(mcl.run.dir, "tmp")
make.dir(mcl.run.dir)
make.dir(mcl.tmp.dir)

# run MCL
input.filename <- file.path(mcl.run.dir, "mcl-input.txt")
fwrite(ecf.network, file = input.filename, sep = " ")
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
