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

# load reprseq tables
reprseq.table <- fread(reprseq.table.output)
all.reprseq <- unique(reprseq.table$reprseq)

###########################################
# 1. Load and modify ECOD reference table
###########################################
ecod.domains <- fread(ecod.domains.filename, colClasses = list(character = "#uid"))
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

# Calculation of total length of each ECOD
ecod.domains[, length := sum(
  sapply(strsplit(seqid_range, ","), function(range) {
    sum(sapply(range, function(r) {
      # Extract start and end positions using a robust regex pattern
      if (grepl(":", r) && grepl("-", r)) {
        positions <- as.numeric(unlist(strsplit(sub(".*:(\\d+-\\d+)", "\\1", r), "-")))
        positions[2] - positions[1] + 1  # Calculate length
      } else {
        0  # If pattern doesn't match, return 0
      }
    }))
  })
), by = seqid_range]

# write to a file
fwrite(ecod.domains, file = ecod.domains.output)

################################
# 2. Load and filter ECOD hits
################################

reprseq.domains.ecod.all <- fread(ecod.domain.hits.filename)
# reprseq.domains.ecod.all[, qname := sub("-hhsearch$", "", qname)]
reprseq.domains.ecod.all[, hit.length := qend - qstart + 1]
reprseq.domains.ecod.all[, qcov := round((qend - qstart + 1)/qlength, 4)]
reprseq.domains.ecod.all[, scov := round((send - sstart + 1)/slength, 4)]

# filter ECOD hits
reprseq.domains.ecod <- reprseq.domains.ecod.all %>%
  filter(eval <= MAX.EVAL.FOR.HMM.FILTERING) %>% 
  filter(prob >= MIN.PROB.FOR.HMM.FILTERING) %>%
  filter(hit.length >= MIN.HIT.LENGTH.FOR.HMM.FILTERING)

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
# Simplify and assign the length based on model_length (average_length is problematic)
pfam.domains <- pfam.domains.full %>% select(accession = pfamA_acc, id = pfamA_id, description, length = model_length)
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

# filter PFAM hits
reprseq.domains.pfam <- reprseq.domains.pfam.all %>%
  filter(prob >= MIN.PROB.FOR.HMM.FILTERING) %>%
  filter(eval <= MAX.EVAL.FOR.HMM.FILTERING) %>% 
  filter(scov >= MIN.SCOV.FOR.HMM.FILTERING) %>%
  filter(qcov >= MIN.QCOV.FOR.HMM.FILTERING) %>%
  filter(hit.length >= MIN.HIT.LENGTH.FOR.HMM.FILTERING)

# write PFAM hits
fwrite(reprseq.domains.pfam, reprseq.domains.pfam.output)
