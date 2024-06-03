### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

#########################
# 1. ECF name conversion
#########################

# load ECF name table
ecf.name.table <- fread(ecf.name.table.filename, header = F)
setnames(ecf.name.table, c("er.name", "ecf.name"))
setkey(ecf.name.table, er.name)

#####################################
# 2. Process reprseq to ECF mapping
#####################################

# load reprseq to ecf mapping
reprseq.vs.ecfs.all <- fread(ecf.hits.filename)
reprseq.vs.ecfs.all[, hit.length := qend - qstart + 1]
reprseq.vs.ecfs.all[, qcov := round((qend - qstart + 1)/qlength, 4)]
reprseq.vs.ecfs.all[, scov := round((send - sstart + 1)/slength, 4)]

# filter ECF hits
reprseq.vs.ecfs <- reprseq.vs.ecfs.all %>%
  filter(eval <= MAX.EVAL.FOR.DOMAIN) %>% 
  filter(prob >= MIN.PROB.FOR.DOMAIN) %>%
  filter(hit.length >= MIN.DOMAIN.LENGTH)
reprseq.vs.ecfs$sname <- ecf.name.table[.(reprseq.vs.ecfs$sname)]$ecf.name

# write to file
fwrite(reprseq.vs.ecfs, file = reprseq.ecfs.output)

#####################################
# 2. Process ECF to ECF mapping
#####################################

ecfs.vs.ecfs.all <- fread(ecf.hits.vs.ecf.filename)
ecfs.vs.ecfs.all[, hit.length := qend - qstart + 1]
ecfs.vs.ecfs.all[, qcov := round((qend - qstart + 1)/qlength, 4)]
ecfs.vs.ecfs.all[, scov := round((send - sstart + 1)/slength, 4)]

# filter ECF-to-ECF hits
ecfs.vs.ecfs <- ecfs.vs.ecfs.all %>%
  filter(eval <= MAX.EVAL.FOR.ECF) %>% 
  filter(prob >= MIN.PROB.FOR.ECF) %>%
  filter(hit.length >= MIN.ECF.HIT.LENGTH)

# make ecf names consistent
# all(ecfs.vs.ecfs$sname %in% ecf.name.table$er.name)
ecfs.vs.ecfs$sname <- ecf.name.table[.(ecfs.vs.ecfs$sname)]$ecf.name
ecfs.vs.ecfs$qname <- gsub("_final", "", ecfs.vs.ecfs$qname)

# write to file
fwrite(ecfs.vs.ecfs, file = ecfs.vs.ecfs.output)

#####################################
# 3. Process ECF to ECOD mapping
#####################################

ecfs.vs.ecod.all <- fread(ecf.hits.vs.ecod.filename)
ecfs.vs.ecod.all[, hit.length := qend - qstart + 1]
ecfs.vs.ecod.all[, qcov := round((qend - qstart + 1)/qlength, 4)]
ecfs.vs.ecod.all[, scov := round((send - sstart + 1)/slength, 4)]

# filter ECF-to-ECF hits
ecfs.vs.ecod <- ecfs.vs.ecod.all %>%
  filter(eval <= MAX.EVAL.FOR.ECF) %>% 
  filter(prob >= MIN.PROB.FOR.ECF) %>%
  filter(hit.length >= MIN.ECF.HIT.LENGTH)

# make ecf names consistent
ecfs.vs.ecod$qname <- gsub("_final", "", ecfs.vs.ecod$qname)

# write to file
fwrite(ecfs.vs.ecod, file = ecfs.vs.ecod.output)

#####################################
# 4. Process ECF to PFAM mapping
#####################################

ecfs.vs.pfam.all <- fread(ecf.hits.vs.pfam.filename)
ecfs.vs.pfam.all[, hit.length := qend - qstart + 1]
ecfs.vs.pfam.all[, qcov := round((qend - qstart + 1)/qlength, 4)]
ecfs.vs.pfam.all[, scov := round((send - sstart + 1)/slength, 4)]

# filter ECF-to-ECF hits
ecfs.vs.pfam <- ecfs.vs.pfam.all %>%
  filter(eval <= MAX.EVAL.FOR.ECF) %>% 
  filter(prob >= MIN.PROB.FOR.ECF) %>%
  filter(hit.length >= MIN.ECF.HIT.LENGTH)

# make ecf names consistent
ecfs.vs.pfam$qname <- gsub("_final", "", ecfs.vs.pfam$qname)

# write to file
fwrite(ecfs.vs.pfam, file = ecfs.vs.pfam.output)


#####################################
# 5. Calculate ECF families (sECFs)
#####################################

# to be written

