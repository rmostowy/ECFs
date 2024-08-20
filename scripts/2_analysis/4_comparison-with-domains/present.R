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

# Load ECF classification table
ecf.classification.table <- fread(ecf.classification.table.output)
ecf.subdomains <- ecf.classification.table[ecf.classification.table$type == "subdomain"]

#############
# 1. Identify all ECOD architechtures where subdomain ECFs are found via direct comparison (ECFs vs ECODs)
#############

ecfs.vs.ecod <- fread(ecfs.vs.ecod.output)
ecfs.full.vs.ecod <- ecfs.vs.ecod %>%
  # filter(qcov >= MIN.QCOV.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(qcov >= 0.7) %>%
  filter(prob >= MIN.PROB.FOR.DOMAIN.EQUIVALENCE)
setkey(ecfs.full.vs.ecod, qname)

x.levels.per.subdomain <- lapply(1:nrow(ecf.subdomains), function(subdomain.index){
  this.ecf <- ecf.subdomains$ecf[subdomain.index]
  this.ecf.hits <- ecfs.full.vs.ecod[.(this.ecf)]
  this.ecf.hits$x.level <- ecod.domains[.(this.ecf.hits$sname)]$x.level
  this.ecf.hits$x.name <- ecod.domains[.(this.ecf.hits$sname)]$x_name
  unique(this.ecf.hits$x.name)  
})
x.levels.per.subdomain.l <- sapply(x.levels.per.subdomain, length)
x.levels.per.subdomain[x.levels.per.subdomain.l>1]

# TODO: make a network of all X levels linked via subdomain ECFs

#############
# 2. Identify all ECOD architectures in reprseqs where subdomain ECFs are found 
#############

# A) First we check if putative subdomain ECFs, identidied via direct ECF-ECOD comparisons, are also exclusive subdomains
# when compared from the point of view of reprseqs

# Load reprseq vs ECFs
reprseq.ecfs.all <- fread(reprseq.ecfs.output)
reprseq.ecfs <- reprseq.ecfs.all %>%
  filter(scov >= MIN.SCOV.FOR.DOMAIN) %>%
  filter(qcov >= MIN.QCOV.FOR.DOMAIN)

# Load reprseq vs ECOD
reprseq.ecod.all <- fread(reprseq.domains.ecod.output)
reprseq.ecod <- reprseq.ecod.all %>%
  filter(scov >= MIN.SCOV.FOR.DOMAIN) %>%
  filter(qcov >= MIN.QCOV.FOR.DOMAIN)

# Calculate overlap between all ECFs and ECODs in all reprseq

reprseq.ecfs.overlap <- reprseq.ecfs %>% select(qname, qstart.ecf = qstart , qend.ecf = qend, qlength.ecf = qlength,
                                                sname.ecf = sname, sstart.ecf = sstart, send.ecf = send, slength.ecf = slength,
                                                pident.ecf = pident, hit.length.ecf = hit.length)

reprseq.ecods.overlap <- reprseq.ecod %>% select(qname, qstart.ecod = qstart , qend.ecod = qend, qlength.ecod = qlength,
                                                sname.ecod = sname, sstart.ecod = sstart, send.ecod = send, slength.ecod = slength,
                                                pident.ecod = pident, hit.length.ecod = hit.length, bitscore.ecod = bitscore)

combined.ecfs.ecods <- reprseq.ecfs.overlap[reprseq.ecods.overlap, on = .(qname), allow.cartesian = TRUE]

combined.ecfs.ecods[, c("overlap.start", "overlap.end", "overlap.length") := {
  res <- calculate_overlap_parms(qstart.ecf, qend.ecf, qstart.ecod, qend.ecod)
  list(res$overlap.start, res$overlap.end, res$overlap.length)
}]
combined.ecfs.ecods <- na.omit(combined.ecfs.ecods, cols = "overlap.length")[overlap.length > 0]

reprseq.ecf.ecod.overlaps <- combined.ecfs.ecods %>% filter(overlap.length >= MIN.OVERLAP.LENGTH)

# Filter overlaps only for subdomain ECFs
subdomain.ecods.overlaps <- reprseq.ecf.ecod.overlaps[reprseq.ecf.ecod.overlaps$sname.ecf %in% ecf.subdomains$ecf]
setkey(subdomain.ecods.overlaps, sname.ecf)

subdomain.ecods.overlaps$ecf.cov.overlap <- round(subdomain.ecods.overlaps$overlap.length/subdomain.ecods.overlaps$hit.length.ecf, d=4)
subdomain.ecods.overlaps$ecod.cov.overlap <- round(subdomain.ecods.overlaps$overlap.length/subdomain.ecods.overlaps$hit.length.ecod, d=4)
subdomain.ecods.overlaps[, scov.ecod := (send.ecod - sstart.ecod + 1)/slength.ecod]
subdomain.ecods.overlaps[, ecod.cov.by.ecf := round(scov.ecod * ecod.cov.overlap, d = 4)]

subdomain.ecods.overlaps$is.subdomain <- NA
subdomain.ecods.overlaps$is.subdomain[subdomain.ecods.overlaps$ecf.cov.overlap > 0.9 & subdomain.ecods.overlaps$ecod.cov.by.ecf >= 0.6] <- FALSE
subdomain.ecods.overlaps$is.subdomain[subdomain.ecods.overlaps$ecf.cov.overlap > 0.9 & subdomain.ecods.overlaps$ecod.cov.by.ecf < 0.6] <- TRUE

setkey(subdomain.ecods.overlaps, sname.ecf)

subdomain.reprseq.confirmed <- sapply(1:nrow(ecf.subdomains), function(index){
  this.subdomain <- ecf.subdomains$ecf[index]
  this.subdomain.overlaps <- subdomain.ecods.overlaps[.(this.subdomain)]
  this.subdomain.overlaps <- this.subdomain.overlaps[!is.na(this.subdomain.overlaps$is.subdomain)]
  confirmed.subdomain <- all(this.subdomain.overlaps$is.subdomain)
  confirmed.subdomain  
})
ecf.subdomains$reprseq.confirmed <- subdomain.reprseq.confirmed

ecf.subdomains.confirmed <- ecf.subdomains[ecf.subdomains$reprseq.confirmed]

# B) Then we check if 

x.levels.per.subdomain.reprseq <- lapply(1:nrow(ecf.subdomains.confirmed), function(confirmed.subdomain.index){
  this.conf.subdomain <- ecf.subdomains.confirmed$ecf[confirmed.subdomain.index]
  this.conf.subdomain.ecods.overlaps <- subdomain.ecods.overlaps[.(this.conf.subdomain)]
  this.conf.subdomain.ecods.overlaps$x.level <- ecod.domains[.(this.conf.subdomain.ecods.overlaps$sname.ecod)]$x.level
  setorder(this.conf.subdomain.ecods.overlaps, qname, sname.ecf, x.level, -bitscore.ecod)
  this.conf.subdomain.ecods.overlaps.filtered <- this.conf.subdomain.ecods.overlaps[, .SD[1], by = .(qname, sname.ecf, x.level)]
  unique(this.conf.subdomain.ecods.overlaps.filtered$x.level)  
})
x.levels.per.subdomain.reprseq.l <- sapply(x.levels.per.subdomain.reprseq, length)
x.levels.per.subdomain.reprseq[x.levels.per.subdomain.reprseq.l>1]



