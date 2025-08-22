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
ecf.table <- fread(ecf.table.output)
all.ecfs <- ecf.table$ecf.name
setkey(ecf.table, ecf.name)

# Calculate coverage
ecfs.vs.ecfs.cov <- fread(ecfs.vs.ecfs.cov.output)
ecfs.vs.ecfs.cov$qlength <- ecf.table[.(ecfs.vs.ecfs.cov$query)]$length
ecfs.vs.ecfs.cov$slength <- ecf.table[.(ecfs.vs.ecfs.cov$subject)]$length
ecfs.vs.ecfs.cov[, qcov := pmin(qcov.min, qcov.max)]
ecfs.vs.ecfs.cov[, scov := pmin(scov.min, scov.max)]
ecfs.vs.ecfs.cov[, cov.min := pmin(scov.min, qcov.min)]
ecfs.vs.ecfs.cov[, cov.max := pmax(scov.min, qcov.min)]

# calculate weights for ECF families
ecfs.vs.ecfs.cov[, weight := scov*qcov*prob/100]

ecf.nonself.hits <- ecfs.vs.ecfs.cov[qcov >= ECF.FAMILY.QCOV.THR & scov >= ECF.FAMILY.SCOV.THR] %>%
  select(from = query, to = subject, weight)
ecf.self.hits <- data.table(from = all.ecfs, to = all.ecfs, weight = 1)
ecf.mcl.table <- rbind(ecf.nonself.hits, ecf.self.hits)

# -- run MCL --
# create a temp mcl directory
mcl.run.dir <- file.path(data.tmp.dir, "mcl-ecf")
mcl.tmp.dir <- file.path(mcl.run.dir, "tmp")
make.dir(mcl.tmp.dir)

# run MCL
input.filename <- file.path(mcl.run.dir, "mcl-input.txt")
fwrite(ecf.mcl.table, file = input.filename, sep = " ")
output.filename <- file.path(mcl.run.dir, "mcl-output.txt")
run.mcl(input.filename, output.filename, inflation.factor = ECF.FAMILY.INFLATION.FACTOR, tmp.dir = mcl.tmp.dir)

# process MCL output
mcl.output <- readLines(output.filename)
Nf <- length(mcl.output)

mcl.output.dt <- data.table(index = 1:Nf, mcl_output = mcl.output)
mcl.output.dt <- mcl.output.dt[, .(ecf = unlist(tstrsplit(mcl_output, "\t", fixed = TRUE))), by = index]
mcl.output.dt[, ecf.family := sprintf("EF%.*d", nchar(Nf), index)]

# clear the directory
delete.dir(mcl.run.dir)

# save the output
ecf.family.table <- mcl.output.dt %>% select(ecf, family = ecf.family)
setkey(ecf.family.table, ecf)

fwrite(ecf.family.table, file = ecf.family.table.output)
