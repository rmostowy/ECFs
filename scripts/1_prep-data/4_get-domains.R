### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

###########################################
# 1. Load and modify ECOD reference table
###########################################
ecod.domains <- fread(ecod.domains.filename, header = F, colClasses = list(character = "V1"))
setnames(ecod.domains, ecod.headers)

# Separate 'f_id' into 'x.level', 'h.level', 't.level', 'f.level'
ecod.domains[, c("x.level", "h.level", "t.level", "f.level") := tstrsplit(f_id, ".", fixed = TRUE)]

# Create hierarchical IDs
ecod.domains[, x_id := x.level]
ecod.domains[, h_id := paste(x.level, h.level, sep = ".")]
ecod.domains[, t_id := paste(x.level, h.level, t.level, sep = ".")]
ecod.domains[, f_id := paste(x.level, h.level, t.level, f.level, sep = ".")]

# propagate T names upwards if missing
missing.h <- ecod.domains$h_name == "NO_H_NAME"
missing.x <- ecod.domains$x_name == "NO_X_NAME"
ecod.domains[missing.h]$h_name <- ecod.domains[missing.h]$t_name
ecod.domains[missing.x]$x_name <- ecod.domains[missing.x]$t_name

ecod.domains$full.id <- sprintf("ECOD_%s_%s", ecod.domains$uid, ecod.domains$ecod_domain_id)

# write to a file
fwrite(ecod.domains, file = ecod.domains.output)

################################
# 2. Load and filter ECOD hits
################################
reprseq.domains.ecod.all <- fread(ecod.domain.hits.filename)
reprseq.domains.ecod.all[, hit.length := qend - qstart + 1]
reprseq.domains.ecod.all[, qcov := round((qend - qstart + 1)/qlength, 4)]
reprseq.domains.ecod.all[, scov := round((send - sstart + 1)/slength, 4)]

# filter ECOD hits
reprseq.domains.ecod <- reprseq.domains.ecod.all %>%
  filter(eval <= MAX.EVAL.FOR.DOMAIN) %>% 
  filter(prob >= MIN.PROB.FOR.DOMAIN) %>%
  filter(hit.length >= MIN.DOMAIN.LENGTH)

# write ECOD hits
fwrite(reprseq.domains.ecod, reprseq.domains.ecod.output)

################################
# 3. Load and filter PFAM hits
################################
reprseq.domains.pfam.all <- fread(pfam.domain.hits.filename)
reprseq.domains.pfam.all[, hit.length := qend - qstart + 1]
reprseq.domains.pfam.all[, qcov := round((qend - qstart + 1)/qlength, 4)]
reprseq.domains.pfam.all[, scov := round((send - sstart + 1)/slength, 4)]

# filter ECOD hits
reprseq.domains.pfam <- reprseq.domains.pfam.all %>%
  filter(eval <= MAX.EVAL.FOR.DOMAIN) %>% 
  filter(prob >= MIN.PROB.FOR.DOMAIN) %>%
  filter(hit.length >= MIN.DOMAIN.LENGTH)

# write ECOD hits
fwrite(reprseq.domains.pfam, reprseq.domains.pfam.output)
