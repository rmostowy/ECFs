### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressPackageStartupMessages(library(stringr))

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

arch.dt <- data.table(arch = ecod.arch, arch.class = ecod.arch.class)
setkey(arch.dt, arch)

ecod.domains$arch.class <- arch.dt[.(ecod.domains$arch_name)]$arch.class

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

###########################################
# 3. Load and modify PFAM reference table
###########################################
# Read the contents of the SQL file
pfam.sql <- readLines(pfam.sql.filename)

# Extract column names from the SQL file
create_table_start <- grep("^CREATE TABLE `pfamA`", pfam.sql)
create_table_end <- grep("^\\) ENGINE=", pfam.sql)
create_table_statement <- pfam.sql[(create_table_start + 1):(create_table_end - 1)]
column_definitions <- create_table_statement[str_detect(create_table_statement, "^\\s*`")]
column_names <- str_extract(column_definitions, "`\\w+`")
column_names <- str_replace_all(column_names, "`", "")

# Extract relevant columns
pfam.domains.full <- fread(pfam.domains.filename) 
setnames(pfam.domains.full, column_names)
pfam.domains <- pfam.domains.full %>% select(accession = pfamA_acc, id = pfamA_id, description, length = average_length)

# Write to file
fwrite(pfam.domains, file = pfam.domains.output)

################################
# 4. Load and filter PFAM hits
################################
reprseq.domains.pfam.all <- fread(pfam.domain.hits.filename)
reprseq.domains.pfam.all[, hit.length := qend - qstart + 1]
reprseq.domains.pfam.all[, qcov := round((qend - qstart + 1)/qlength, 4)]
reprseq.domains.pfam.all[, scov := round((send - sstart + 1)/slength, 4)]
reprseq.domains.pfam.all[, accession := sub("\\..*", "", sname)]


# filter ECOD hits
reprseq.domains.pfam <- reprseq.domains.pfam.all %>%
  filter(eval <= MAX.EVAL.FOR.DOMAIN) %>% 
  filter(prob >= MIN.PROB.FOR.DOMAIN) %>%
  filter(hit.length >= MIN.DOMAIN.LENGTH)

# write ECOD hits
fwrite(reprseq.domains.pfam, reprseq.domains.pfam.output)
