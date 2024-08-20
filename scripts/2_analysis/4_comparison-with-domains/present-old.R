### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

# Load ECOD domain table
ecod.domains <- fread(ecod.domains.output)
setkey(ecod.domains, full.id)

# Load ECF table
ecfs.vs.qecfs <- fread(ecf.vs.qecf.output)
ecfs.vs.qecfs$reprseq <- substr(ecfs.vs.qecfs$qecf, 1, regexpr("_", ecfs.vs.qecfs$qecf)-1)
setkey(ecfs.vs.qecfs, ecf)
all.ecfs <- unique(ecfs.vs.qecfs$ecf)

# Load ECF vs ECOD table
ecfs.vs.ecod <- fread(ecfs.vs.ecod.output)
ecfs.vs.ecod[, score := round(qcov*scov*prob/100*pident/100, d = 6)]

# Look at the distribution of qcov/scov values for direct comparison

setorder(ecfs.vs.ecod, qname, -bitscore)
ecfs.best.ecod.hits <- ecfs.vs.ecod[, .SD[1], by = .(qname)]
ecfs.best.ecod.hits[, cov := qcov*scov]

ecfs.vs.ecod.raw <- fread(ecf.hits.vs.ecod.filename)
setorder(ecfs.vs.ecod.raw, qname, -prob)
ecfs.vs.ecod.max.prob <- ecfs.vs.ecod.raw[, .SD[1], by = .(qname)]

ggplot(ecfs.best.ecod.hits, aes(x = qcov, y = scov)) + geom_point(size = 1)
ggplot(ecfs.best.ecod.hits, aes(x = cov, y = pident)) + geom_point(size = 1)

# Closer look at subdomains
ecfs.subdomains.strong <- ecfs.best.ecod.hits %>% 
  filter(qcov >= 0.9 & scov < 0.6 & prob >= 95)
ecfs.subdomains.strong$t.name <- ecod.domains[.(ecfs.subdomains.strong$sname)]$t_name
# this set should roughly correspond to Bogna's Figure 3C

# 
ecfs.full.vs.ecod <- ecfs.vs.ecod %>%
  filter(qcov >= MIN.QCOV.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(prob >= MIN.PROB.FOR.DOMAIN.EQUIVALENCE)


ecf.index=1
this.ecf <- ecfs.subdomains.strong$qname[ecf.index]
ecfs.full.vs.ecod[ecfs.full.vs.ecod$qname == this.ecf]



# # look at individual subdomain ECFs
# 
# reprseq.ecfs.all <- fread(reprseq.ecfs.output)
# reprseq.ecfs <- reprseq.ecfs.all %>%
#   filter(scov >= MIN.SCOV.FOR.DOMAIN) %>%
#   filter(qcov >= MIN.QCOV.FOR.DOMAIN)
# 
# reprseq.ecod.all <- fread(reprseq.domains.ecod.output)
# reprseq.ecod <- reprseq.ecod.all %>%
#   filter(scov >= MIN.SCOV.FOR.DOMAIN) %>%
#   filter(qcov >= MIN.QCOV.FOR.DOMAIN)
# 
# reprseq.ecod$x.level <- ecod.domains[.(reprseq.ecod$sname)]$x_id
# reprseq.ecod.xlevels <- reprseq.ecod[, .(x.level = unique(x.level)), by = qname]
# 
# 
# 
# subdomains.in.reprseq <- reprseq.ecfs[reprseq.ecfs$sname %in% ecfs.subdomains.strong$qname] %>% select(sname, qname)
# 
# subdomains_with_levels <- merge(subdomains.in.reprseq, reprseq.ecod.xlevels, by = "qname", all.x = TRUE)
# subdomains_with_levels <- subdomains_with_levels[!is.na(subdomains_with_levels$x.level)]
# 
# subdomains_with_levels <- subdomains_with_levels[, .(x.levels = list(x.level)), by = .(sname, qname)]
# subdomains_with_levels[, dummy := 1]
# 
# # Define a function to find non-overlapping pairs
# find_non_overlapping_pairs <- function(dt) {
#   # Cross join by merging on a dummy column
#   pairs <- dt[dt, on = .(dummy), allow.cartesian = TRUE]
#   pairs <- pairs[qname != i.qname]  # Remove self-pairs
#   
#   # Function to check for non-overlapping x.levels
#   check_non_overlap <- function(x, y) {
#     intersect(unlist(x), unlist(y)) == integer(0)
#   }
#   
#   # Check for non-overlapping x.levels
#   pairs[, non_overlap := check_non_overlap(x.levels, i.x.levels), by = .(qname, i.qname)]
#   
#   # Filter pairs without shared x.level
#   non_shared_pairs <- pairs[non_overlap == TRUE, .(qname1 = qname, qname2 = i.qname)]
#   
#   return(non_shared_pairs)
# }
# 
# # Apply function for each sname
# results <- subdomains_with_levels[, find_non_overlapping_pairs(.SD), by = sname]
# 
# 
# # # Assuming 'this.ecf.in.reprseq' lists proteins with a given ECF and 'reprseq.ecod.xlevels' lists domain architectures
# # # First, join the data to filter only those qnames in this.ecf.in.reprseq
# # ecf_proteins_with_levels <- merge(this.ecf.in.reprseq, reprseq.ecod.xlevels, by = "qname", all.x = TRUE)
# # ecf_proteins_with_levels <- ecf_proteins_with_levels[!is.na(ecf_proteins_with_levels$x.level)]
# # 
# # # # Create all combinations using a cross join
# # # # Cross join by merging on a dummy column
# # # ecf_proteins_with_levels[, dummy := 1]
# # # pairs <- ecf_proteins_with_levels[ecf_proteins_with_levels, on = .(dummy), allow.cartesian = TRUE]
# # # pairs <- pairs[qname != i.qname]  # Remove self-pairs
# # # pairs[, shared := lengths(mapply(intersect, as.list(x.level), as.list(i.x.level))) > 0]
# # # non_shared_pairs <- pairs[shared == FALSE, .(qname1 = qname, qname2 = i.qname, x.level1 = x.level, x.level2 = i.x.level)]
# # 
# # 
# 
# 
# 
# 
# # Aggregate x.levels for each qname into a list
# 
# this.ecf <- "ecf_00119"
# this.ecf.in.reprseq <- reprseq.ecfs[reprseq.ecfs$sname == this.ecf] %>% select(sname, qname)
# ecf_proteins_with_levels <- merge(this.ecf.in.reprseq, reprseq.ecod.xlevels, by = "qname", all.x = TRUE)
# ecf_proteins_with_levels <- ecf_proteins_with_levels[!is.na(ecf_proteins_with_levels$x.level)]
# ecf_proteins_with_levels <- ecf_proteins_with_levels[, .(x.levels = list(unique(x.level))), by = qname]
# 
# ecf_proteins_with_levels[, dummy := 1]
# pairs <- ecf_proteins_with_levels[ecf_proteins_with_levels, on = .(dummy), allow.cartesian = TRUE]
# pairs <- pairs[qname != i.qname]  # Remove self-pairs
# # Function to check for non-overlapping x.levels
# check_non_overlap <- function(x, y) {
#   length(intersect(unlist(x), unlist(y))) == 0
# }
# pairs[, non_overlap := check_non_overlap(x.levels, i.x.levels), by = .(qname, i.qname)]
# non_shared_pairs <- pairs[non_overlap == TRUE, .(qname1 = qname, qname2 = i.qname, x.level1 = x.levels, x.level2 = i.x.levels)]
# non_shared_pairs
# 
