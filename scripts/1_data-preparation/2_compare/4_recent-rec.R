### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

# Load refseq table filtered by NCBI
reprseq.table <- fread(reprseq.table.ncbi.filtered.filename)
all.reprseq <- unique(reprseq.table$reprseq)

# Load ECOD domain table
ecod.domains <- fread(ecod.domains.output)
setkey(ecod.domains, full.id)

# Load ECF families
ecf.family.table <- fread(ecf.family.table.output)
setkey(ecf.family.table, ecf)

# Load ECF table
ecf.class.table <- fread(ecf.classification.table.output)
ecf.class.table$family <- ecf.family.table[.(ecf.class.table$ecf)]$family
ecf2ecod <- ecf.class.table %>% select(ecf, family, t.id)
ecf2ecod$t.name <- sapply(1:nrow(ecf2ecod), function(k){
  out <- "unknown"
  this.t.id <- ecf2ecod$t.id[k]
  if(this.t.id != "unknown") out <- ecod.domains[t_id == this.t.id]$t_name[1]  
  return(out)
})
setorder(ecf2ecod, family)

# Map TID onto each unique family
ef2ecod <- ecf2ecod[
  , .(
    t.id = if (all(t.id == "unknown")) "unknown" else first(t.id[t.id != "unknown"])
  ), 
  by = family
]
setkey(ef2ecod, family)

# Load all-by-all reprseq
reprseq.vs.reprseq <- fread(reprseq.vs.reprseq.filename)
reprseq.vs.reprseq[, qcov := round((qend - qstart + 1)/qlength, d = 4)]
reprseq.vs.reprseq[, scov := round((send - sstart + 1)/slength, d = 4)]
reprseq.vs.reprseq <- reprseq.vs.reprseq[qname %chin% all.reprseq & sname %chin% all.reprseq]

# Load functional annotations
reprseq.functions.unique <- fread(reprseq.functions.unique.output)
setkey(reprseq.functions.unique, reprseq)

# Load reprseq vs ECFs
reprseq.ecfs.all <- fread(reprseq.ecfs.output)
reprseq.ecfs <- reprseq.ecfs.all %>%
  filter(prob >= MIN.PROB.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(scov >= MIN.SCOV.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE)
reprseq.ecfs$EF <- ecf.family.table[.(reprseq.ecfs$sname)]$family
setkey(reprseq.ecfs, qname)

# Load reprseq vs ECOD
reprseq.ecod.all <- fread(reprseq.domains.ecod.output)
reprseq.ecod <- reprseq.ecod.all %>%
  filter(prob >= MIN.PROB.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(scov >= MIN.SCOV.FOR.DOMAIN.EQUIVALENCE)
reprseq.ecod$t.id <- ecod.domains[.(reprseq.ecod$sname)]$t_id
setkey(reprseq.ecod, qname)

reprseq.cov.filename <- '/Users/rmostowy/MGG Dropbox/Projects/divRBP/phage-pp-workdir-refseq-hhblits/output/prot-families/families/dataset-full/repr-hits-pairwise-prob95.csv'
reprseq.cov.orig <- fread(reprseq.cov.filename)
reprseq.cov.orig[, max_cov := NULL]
reprseq.cov.orig[, min_cov := NULL]

# adjust pair names so that qname < sname
reprseq.cov.ab <- reprseq.cov.orig
reprseq.cov.ba <- reprseq.cov.orig %>% 
  select(qname = sname, sname = qname, prob,
         scov_min = qcov_min, scov_max = qcov_max,
         qcov_min = scov_min, qcov_max = scov_max, pident)
reprseq.cov.both <- rbind(reprseq.cov.ab, reprseq.cov.ba)
reprseq.cov <- reprseq.cov.both[qname < sname]
setkey(reprseq.cov, qname, sname)


# find hits with the highest pident for all reprseq pairs
reprseq.vs.reprseq.prefiltered <- reprseq.vs.reprseq[prob >= 95]
setorder(reprseq.vs.reprseq.prefiltered, qname, sname, -pident)
reprseq.vs.reprseq.filtered <- reprseq.vs.reprseq.prefiltered[, .SD[1], by = .(qname, sname)]
reprseq.vs.reprseq.filtered <- reprseq.vs.reprseq.filtered[qname < sname]
reprseq.vs.reprseq.filtered.qcov.scov <- reprseq.cov[.(reprseq.vs.reprseq.filtered[, .(qname, sname)])]
reprseq.vs.reprseq.filtered.qcov.scov[, qcov := pmin(qcov_min, qcov_max, na.rm = T)]
reprseq.vs.reprseq.filtered.qcov.scov[, scov := pmin(scov_min, scov_max, na.rm = T)]
rr <- cbind(reprseq.vs.reprseq.filtered, reprseq.vs.reprseq.filtered.qcov.scov %>% select(pident.full = pident, qcov.full = qcov, scov.full = scov))
rr[, pident.full := pident.full*100]
rr[, hit.length := qend - qstart + 1]

rr[, cov.max.full := pmax(qcov.full, scov.full, na.rm = T)]
rrr <- rr[pident > 90 & cov.max.full < 0.5 & hit.length >= 30]

