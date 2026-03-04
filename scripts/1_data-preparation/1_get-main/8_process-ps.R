### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(Biostrings)))
# suppressPackageStartupMessages(library(stringr))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

# load PhageScope metadata
ps.metadata <- fread(phagescope.protein.metadata.filename)
ps.seqs <- readAAStringSet(phagescope.protein.seqs.filename)
ps.seqs.n <- names(ps.seqs)

# load domain tables
ecf.table <- fread(ecf.table.output)
setkey(ecf.table, er.name)

# load domain hits for:

# a) ECFs
ps.vs.ecfs.all <- fread(ps.vs.ecf.filename)
ps.vs.ecfs.all[, hit.length := qend - qstart + 1]
ps.vs.ecfs.all[, qcov := round((qend - qstart + 1)/qlength, 4)]
ps.vs.ecfs.all[, scov := round((send - sstart + 1)/slength, 4)]
ps.vs.ecfs.all$sname <- ecf.table[.(ps.vs.ecfs.all$sname)]$ecf.name

ps.vs.ecfs <- ps.vs.ecfs.all %>%
  filter(prob >= MIN.PROB.FOR.HMM.FILTERING) %>%
  filter(eval <= MAX.EVAL.FOR.HMM.FILTERING) %>% 
  filter(scov >= MIN.SCOV.FOR.HMM.FILTERING) %>%
  filter(qcov >= MIN.QCOV.FOR.HMM.FILTERING) %>%
  filter(hit.length >= MIN.HIT.LENGTH.FOR.HMM.FILTERING)

# b) ECODs
ps.vs.ecod.all <- fread(ps.vs.ecod.filename)
ps.vs.ecod.all[, hit.length := qend - qstart + 1]
ps.vs.ecod.all[, qcov := round((qend - qstart + 1)/qlength, 4)]
ps.vs.ecod.all[, scov := round((send - sstart + 1)/slength, 4)]

ps.vs.ecod <- ps.vs.ecod.all %>%
  filter(prob >= MIN.PROB.FOR.HMM.FILTERING) %>%
  filter(eval <= MAX.EVAL.FOR.HMM.FILTERING) %>% 
  filter(scov >= MIN.SCOV.FOR.HMM.FILTERING) %>%
  filter(qcov >= MIN.QCOV.FOR.HMM.FILTERING) %>%
  filter(hit.length >= MIN.HIT.LENGTH.FOR.HMM.FILTERING)

# c) PFAMs
ps.vs.pfam.all <- fread(ps.vs.pfam.filename)
ps.vs.pfam.all[, hit.length := qend - qstart + 1]
ps.vs.pfam.all[, qcov := round((qend - qstart + 1)/qlength, 4)]
ps.vs.pfam.all[, scov := round((send - sstart + 1)/slength, 4)]
ps.vs.pfam.all[, accession := sub("\\..*", "", sname)]

ps.vs.pfam <- ps.vs.pfam.all %>%
  filter(prob >= MIN.PROB.FOR.HMM.FILTERING) %>%
  filter(eval <= MAX.EVAL.FOR.HMM.FILTERING) %>% 
  filter(scov >= MIN.SCOV.FOR.HMM.FILTERING) %>%
  filter(qcov >= MIN.QCOV.FOR.HMM.FILTERING) %>%
  filter(hit.length >= MIN.HIT.LENGTH.FOR.HMM.FILTERING)

# write data
fwrite(ps.metadata, file = phagescope.metadata.output)
writeXStringSet(ps.seqs, filepath = phagescope.sequences.output)

fwrite(ps.vs.ecfs, file = phagescope.vs.ecf)
fwrite(ps.vs.ecod, file = phagescope.vs.ecod)
fwrite(ps.vs.pfam, file = phagescope.vs.pfam)

