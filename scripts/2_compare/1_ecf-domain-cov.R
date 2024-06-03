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


#########################
# 3. Calculate overlaps
#########################

this.reprseq <- "reprseq000005"

this.reprseq.ecfs <- reprseq.ecfs[reprseq.ecfs$qname == this.reprseq]
this.reprseq.ecods <- reprseq.ecod[reprseq.ecod$qname == this.reprseq]


#######################
# Proposed solution:
#######################

# Subset data tables based on the specific query
reprseq.ecfs.test <- copy(reprseq.ecfs)
reprseq.ecods.test <- copy(reprseq.ecod)

# Function to calculate overlap and coverage
calculate_overlap_and_coverage <- function(start1, end1, start2, end2) {
  overlap_length <- pmax(0, pmin(end1, end2) - pmax(start1, start2))
  coverage <- overlap_length / (end1 - start1)
  return(list(overlap_length = overlap_length, coverage = coverage))
}

# Rename columns for clarity
setnames(reprseq.ecfs.test, c("sname", "sstart", "send"), c("sname_db1", "sstart_db1", "send_db1"))
setnames(reprseq.ecods.test, c("sname", "sstart", "send", "qstart", "qend"), c("sname_db2", "sstart_db2", "send_db2", "qstart_db2", "qend_db2"))

# Perform a self-join to calculate coverage for all possible overlaps
combined <- reprseq.ecfs.test[reprseq.ecods.test, on = .(qname), allow.cartesian = TRUE]

# Calculate overlap and coverage
combined[, c("overlap_length", "cov") := {
  res <- calculate_overlap_and_coverage(qstart, qend, qstart_db2, qend_db2)
  list(res$overlap_length, res$coverage)
}]

# Remove rows with NA and 0 values
combined <- na.omit(combined, cols = "cov")[cov > 0]

# Find maximum coverage and corresponding protein name for each row in db1
# result <- combined[, .SD[which.max(cov)], by = .(qname, qstart, qend, sname_db1, sstart_db1, send_db1)]

# Order by qname, sname_db1, and decreasing coverage
setorder(combined, qname, sname_db1, -cov)

# Select the first row for each combination of qname, qstart, qend, sname_db1, sstart_db1, send_db1
result <- combined[, .SD[1], by = .(qname, qstart, qend, sname_db1, sstart_db1, send_db1)]

# Filter by hit length
result[result$overlap_length >= 20]




# 
# ###########################################################
# 
# # Function to calculate coverage
# calculate_coverage <- function(start1, end1, start2, end2) {
#   pmax(0, pmin(end1, end2) - pmax(start1, start2)) / (end1 - start1)
# }
# 
# # Rename columns for clarity
# setnames(this.reprseq.ecfs, c("sname", "sstart", "send"), c("sname_db1", "sstart_db1", "send_db1"))
# setnames(this.reprseq.ecods, c("sname", "sstart", "send"), c("sname_db2", "sstart_db2", "send_db2"))
# 
# # Perform a self-join to calculate coverage for all possible overlaps
# combined <- this.reprseq.ecfs[this.reprseq.ecods, on = .(qname), allow.cartesian = TRUE]
# 
# # Calculate coverage
# combined[, cov := calculate_coverage(qstart, qend, qstart.1, qend.1)]
# 
# # Find maximum coverage and corresponding protein name for each row in db1
# result <- combined[, .SD[which.max(cov)], by = .(qname, qstart, qend, sname_db1, sstart_db1, send_db1)]
# 
# 
# # Function to calculate coverage
# calculate_coverage <- function(start1, end1, start2, end2) {
#   pmax(0, pmin(end1, end2) - pmax(start1, start2)) / (end1 - start1)
# }
# 
# # Rename columns for clarity
# setnames(this.reprseq.ecfs, c("sname", "sstart", "send"), c("sname_db1", "sstart_db1", "send_db1"))
# setnames(this.reprseq.ecods, c("sname", "sstart", "send", "qstart", "qend"), c("sname_db2", "sstart_db2", "send_db2", "qstart_db2", "qend_db2"))
# 
# # Perform a self-join to calculate coverage for all possible overlaps
# combined <- this.reprseq.ecfs[this.reprseq.ecods, on = .(qname), allow.cartesian = TRUE]
# 
# # Calculate coverage
# combined[, cov := calculate_coverage(qstart, qend, qstart_db2, qend_db2)]
# 
# # Find maximum coverage and corresponding protein name for each row in this.reprseq.ecfs
# result <- combined[, .SD[which.max(cov)], by = .(qname, qstart, qend, sname_db1, sstart_db1, send_db1)]
# 
# # Select and rename columns for the final output
# result <- result[, .(qname, qstart, qend, sname = sname_db1, sstart = sstart_db1, send = send_db1, max_cov = cov, max_cov_name = sname_db2)]
# result[result$max_cov>0]
# 
# # Output result
# print(result)
# 
# 
# #########################
# 
# this.reprseq <- "reprseq000005"
# 
# # Subset data tables based on the specific query
# reprseq.ecfs.test <- copy(reprseq.ecfs)
# reprseq.ecods.test <- copy(reprseq.ecod)
# 
# # Function to calculate overlap and coverage
# calculate_overlap_and_coverage <- function(start1, end1, start2, end2) {
#   overlap_length <- pmax(0, pmin(end1, end2) - pmax(start1, start2))
#   coverage <- overlap_length / (end1 - start1)
#   return(list(overlap_length = overlap_length, coverage = coverage))
# }
# 
# # Rename columns for clarity
# setnames(reprseq.ecfs.test, c("sname", "sstart", "send"), c("sname_db1", "sstart_db1", "send_db1"))
# setnames(reprseq.ecods.test, c("sname", "sstart", "send", "qstart", "qend"), c("sname_db2", "sstart_db2", "send_db2", "qstart_db2", "qend_db2"))
# 
# # Perform a self-join to calculate coverage for all possible overlaps
# combined <- reprseq.ecfs.test[reprseq.ecods.test, on = .(qname), allow.cartesian = TRUE]
# 
# # Calculate overlap and coverage
# combined[, c("overlap_length", "cov") := calculate_overlap_and_coverage(qstart, qend, qstart_db2, qend_db2), by = .(qname, qstart, qend, sname_db1, sstart_db1, send_db1)]
# 
# # Find maximum coverage and corresponding protein name for each row in this.reprseq.ecfs
# result <- combined[, .SD[which.max(cov)], by = .(qname, qstart, qend, sname_db1, sstart_db1, send_db1)]
# 
# # Select and rename columns for the final output
# result <- result[, .(qname, qstart, qend, sname = sname_db1, sstart = sstart_db1, send = send_db1, overlap_length = overlap_length, max_cov = cov, max_cov_name = sname_db2)]
# result[result$overlap_length>=20]
# 
# # Output result
# print(result)
# 
# 
# ##########################
# 
# 
# # Perform a self-join to calculate coverage for all possible overlaps
# combined2 <- reprseq.ecfs.test[reprseq.ecods.test, on = .(qname), allow.cartesian = TRUE]
# 
# # Calculate overlap and coverage
# combined2[, c("overlap_length", "cov") := {
#   res <- calculate_overlap_and_coverage(qstart, qend, qstart_db2, qend_db2)
#   list(res$overlap_length, res$coverage)
# }]
# 
# # Find maximum coverage and corresponding protein name for each row in db1
# combined <- combined[order(-cov)]
# 
# result2 <- combined[, .SD[1], by = .(qname, qstart, qend, sname_db1, sstart_db1, send_db1)]
# 
# # Select and rename columns for the final output
# result2 <- result2[, .(qname, qstart, qend, sname = sname_db1, sstart = sstart_db1, send = send_db1, overlap_length, max_cov = cov, max_cov_name = sname_db2)]
# result2[result2$overlap_length>=20]
# 
# # Output result
# print(result)
# 
