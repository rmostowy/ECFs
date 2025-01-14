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
ecfs.vs.ecfs.cov[, weight.fam := scov*qcov*prob/100]

# calculate weights for ECF superfamilies
ecfs.vs.ecfs.cov[, weight.sfam := cov.max*prob/100]

# filter tables as MCL inputs
ecf.mcl.fam.table <- ecfs.vs.ecfs.cov[qcov >= ECF.FAMILY.COV.THR & scov >= ECF.FAMILY.COV.THR] %>% 
  select(from = query, to = subject, weight = weight.fam)
ecf.mcl.sfam.table <- ecfs.vs.ecfs.cov[cov.max >= ECF.SUPERFAMILY.MAX.COV.THR] %>% 
  select(from = query, to = subject, weight = weight.sfam)

# add self hits
ecf.self.hits <- data.table(from = all.ecfs, to = all.ecfs, weight = 1)
ecf.mcl.fam.table <- rbind(ecf.mcl.fam.table, ecf.self.hits)
ecf.mcl.sfam.table <- rbind(ecf.mcl.sfam.table, ecf.self.hits)

# -- run MCL for families --
# create a temp mcl directory
mcl.run.dir <- file.path(data.tmp.dir, "mcl")
mcl.tmp.dir <- file.path(mcl.run.dir, "tmp-fam")
make.dir(mcl.run.dir)
make.dir(mcl.tmp.dir)

# run MCL
input.filename <- file.path(mcl.run.dir, "mcl-input-fam.txt")
fwrite(ecf.mcl.fam.table, file = input.filename, sep = " ")
output.filename <- file.path(mcl.run.dir, "mcl-output-fam.txt")
run.mcl(input.filename, output.filename, inflation.factor = 2, tmp.dir = mcl.tmp.dir)

# process MCL output
mcl.fam.output <- readLines(output.filename)
Nf <- length(mcl.fam.output)

mcl.fam.output.dt <- data.table(index = 1:Nf, mcl_output = mcl.fam.output)
mcl.fam.output.dt <- mcl.fam.output.dt[, .(ecf = unlist(tstrsplit(mcl_output, "\t", fixed = TRUE))), by = index]
mcl.fam.output.dt[, ecf.family := sprintf("FAM%.*d", nchar(Nf), index)]
# ecf.family.table <- mcl.vs.ecf[, .(ecf.family, ecf)]

# -- run MCL for superfamilies --
# create a temp mcl directory
mcl.tmp.dir <- file.path(mcl.run.dir, "tmp-sfam")
make.dir(mcl.tmp.dir)

# run MCL
input.filename <- file.path(mcl.run.dir, "mcl-input-sfam.txt")
fwrite(ecf.mcl.sfam.table, file = input.filename, sep = " ")
output.filename <- file.path(mcl.run.dir, "mcl-output-sfam.txt")
run.mcl(input.filename, output.filename, inflation.factor = 2, tmp.dir = mcl.tmp.dir)

# process MCL output
mcl.sfam.output <- readLines(output.filename)
Nsf <- length(mcl.sfam.output)

mcl.sfam.output.dt <- data.table(index = 1:Nsf, mcl_output = mcl.sfam.output)
mcl.sfam.output.dt <- mcl.sfam.output.dt[, .(ecf = unlist(tstrsplit(mcl_output, "\t", fixed = TRUE))), by = index]
mcl.sfam.output.dt[, ecf.superfamily := sprintf("SFAM%.*d", nchar(Nsf), index)]

# merge the tables
ecf.family.table <- mcl.fam.output.dt %>% select(ecf, fam = ecf.family)
ecf.superfamily.table <- mcl.sfam.output.dt %>% select(ecf, sfam = ecf.superfamily)

ecf.family.table.final <- merge(ecf.family.table, ecf.superfamily.table, by = "ecf")
fwrite(ecf.family.table.final, file = ecf.family.table.output)

