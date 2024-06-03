### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

# load ECOD domain description table
ecod.domains <- fread(ecod.domains.filename, header = F)
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

# write to a file
fwrite(ecod.domains, file = ecod.domains.output)

