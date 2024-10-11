### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

# reprseq.ecfs.output <- file.path(data.ecfs.dir, "reprseq-ecfs.txt") 
# ecfs.vs.ecfs.output <- file.path(data.ecfs.dir, "ecfs-vs-ecfs.txt") 
# ecfs.vs.ecod.output <- file.path(data.ecfs.dir, "ecfs-vs-ecod.txt") 
# ecfs.vs.pfam.output <- file.path(data.ecfs.dir, "ecfs-vs-pfam.txt") 

#########################
# 1. Load data
#########################

reprseq.table <- fread(reprseq.table.filename)
all.reprseq <- unique(reprseq.table$reprseq)

reprseq.ecfs.all <- fread(reprseq.ecfs.output)
reprseq.ecod.all <- fread(reprseq.domains.ecod.output)
reprseq.pfam.all <- fread(reprseq.domains.pfam.output)

#########################
# 2. Filter hits
#########################

reprseq.ecfs <- reprseq.ecfs.all %>%
  filter(scov >= MIN.SCOV.FOR.DOMAIN) %>%
  filter(qcov >= MIN.QCOV.FOR.DOMAIN)

setorder(reprseq.ecfs, qname, sname, -bitscore, eval)
reprseq.ecfs <- reprseq.ecfs[ , .SD[1], by = .(qname, sname)]

reprseq.ecod <- reprseq.ecod.all %>%
  filter(scov >= MIN.SCOV.FOR.DOMAIN) %>%
  filter(qcov >= MIN.QCOV.FOR.DOMAIN)

setorder(reprseq.ecod, qname, sname, -bitscore, eval)
reprseq.ecod <- reprseq.ecod[ , .SD[1], by = .(qname, sname)]

reprseq.pfam <- reprseq.pfam.all %>%
  filter(scov >= MIN.SCOV.FOR.DOMAIN) %>%
  filter(qcov >= MIN.QCOV.FOR.DOMAIN)

setorder(reprseq.pfam, qname, sname, -bitscore, eval)
reprseq.pfam <- reprseq.pfam[ , .SD[1], by = .(qname, sname)]


##################################################
# 3. Calculate overlaps between reprseq and ECOD
##################################################

reprseq.ecfs.overlap <- copy(reprseq.ecfs)
reprseq.ecods.overlap <- copy(reprseq.ecod)

# Rename columns for clarity
setnames(reprseq.ecfs.overlap, c("sname", "sstart", "send"), c("sname_db1", "sstart_db1", "send_db1"))
setnames(reprseq.ecods.overlap, c("sname", "sstart", "send", "qstart", "qend"), c("sname_db2", "sstart_db2", "send_db2", "qstart_db2", "qend_db2"))

# Perform a self-join to calculate coverage for all possible overlaps
combined.ecfs.ecods <- reprseq.ecfs.overlap[reprseq.ecods.overlap, on = .(qname), allow.cartesian = TRUE]

# Calculate overlap parameters
combined.ecfs.ecods[, c("overlap.start", "overlap.end", "overlap.length") := {
  res <- calculate_overlap_parms(qstart, qend, qstart_db2, qend_db2)
  list(res$overlap.start, res$overlap.end, res$overlap.length)
}]

# Remove rows with NA and 0 values
combined.ecfs.ecods <- na.omit(combined.ecfs.ecods, cols = "overlap.length")[overlap.length > 0]

# Order by qname, sname_db1, and decreasing coverage
setorder(combined.ecfs.ecods, qname, sname_db1, -overlap.length)

# Select the first row for each combination of qname, qstart, qend, sname_db1, sstart_db1, send_db1
result.ecfs.ecods <- combined.ecfs.ecods[, .SD[1], by = .(qname, sname_db1)]

reprseq.ecf.ecod.overlaps <- result.ecfs.ecods %>%
  filter(overlap.length >= MIN.OVERLAP.LENGTH) %>%
  select(qname, qlength, 
         qstart.ecf = qstart, qend.ecf = qend, hit.length.ecf = hit.length, sname.ecf = sname_db1, slength.ecf = slength,
         qstart.ecod = qstart_db2, qend.ecod = qend_db2, hit.length.ecod = i.hit.length, sname.ecod = sname_db2,
         overlap.start, overlap.end, overlap.length)

# write to file
fwrite(reprseq.ecf.ecod.overlaps, file = reprseq.ecf.vs.ecod.output)
  
##################################################
# 4. Calculate overlaps between reprseq and PFAM
##################################################

reprseq.pfam.overlap <- copy(reprseq.pfam)

# Rename columns for clarity
setnames(reprseq.pfam.overlap, c("sname", "sstart", "send", "qstart", "qend"), c("sname_db2", "sstart_db2", "send_db2", "qstart_db2", "qend_db2"))

# Perform a self-join to calculate coverage for all possible overlaps
combined.ecfs.pfam <- reprseq.ecfs.overlap[reprseq.pfam.overlap, on = .(qname), allow.cartesian = TRUE]

# Calculate overlap parameters
combined.ecfs.pfam[, c("overlap.start", "overlap.end", "overlap.length") := {
  res <- calculate_overlap_parms(qstart, qend, qstart_db2, qend_db2)
  list(res$overlap.start, res$overlap.end, res$overlap.length)
}]

# Remove rows with NA and 0 values
combined.ecfs.pfam <- na.omit(combined.ecfs.pfam, cols = "overlap.length")[overlap.length > 0]

# Order by qname, sname_db1, and decreasing coverage
setorder(combined.ecfs.ecods, qname, sname_db1, -overlap.length)

# Select the first row for each combination of qname, qstart, qend, sname_db1, sstart_db1, send_db1
result.ecfs.pfam <- combined.ecfs.pfam[, .SD[1], by = .(qname, sname_db1)]

reprseq.ecf.pfam.overlaps <- result.ecfs.pfam %>%
  filter(overlap.length >= MIN.OVERLAP.LENGTH) %>%
  select(qname, qlength, 
         qstart.ecf = qstart, qend.ecf = qend, hit.length.ecf = hit.length, sname.ecf = sname_db1, slength.ecf = slength,
         qstart.pfam = qstart_db2, qend.pfam = qend_db2, hit.length.pfam = i.hit.length, sname.pfam = sname_db2,
         overlap.start, overlap.end, overlap.length)

# Write to file
fwrite(reprseq.ecf.pfam.overlaps, file = reprseq.ecf.vs.pfam.output)

######################
# 5. Combine tables
######################

reprseq.ecf.ecod.overlaps.short <- reprseq.ecf.ecod.overlaps %>%
  select(qname, qlength, qstart.ecf, qend.ecf, hit.length.ecf, sname.ecf, slength.ecf,
         qstart.ecod, qend.ecod,
         overlap.ecod.start = overlap.start, overlap.ecod.end = overlap.end, overlap.ecod.length = overlap.length)

reprseq.pfam.ecod.overlaps.short <- reprseq.ecf.pfam.overlaps %>%
  select(qname, qlength, qstart.ecf, qend.ecf, hit.length.ecf, sname.ecf, slength.ecf,
         qstart.pfam, qend.pfam,
         overlap.pfam.start = overlap.start, overlap.pfam.end = overlap.end, overlap.pfam.length = overlap.length)

reprseq.ecf.ecod.pfam.overlaps.merged <- merge(reprseq.ecf.ecod.overlaps.short, reprseq.pfam.ecod.overlaps.short, by = c("qname", "qlength", "qstart.ecf", "qend.ecf", "hit.length.ecf", "sname.ecf", "slength.ecf"), all = TRUE)
fwrite(reprseq.ecf.ecod.pfam.overlaps.merged, file = reprseq.ecf.vs.ecod.and.pfam.output)

