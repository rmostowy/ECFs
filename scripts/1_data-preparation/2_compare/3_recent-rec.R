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
reprseq.table <- fread(reprseq.table.filename)
all.reprseq <- unique(reprseq.table$reprseq)

# Load all-by-all reprseq
reprseq.vs.reprseq <- fread(reprseq.vs.reprseq.filename)
reprseq.vs.reprseq[, qcov := round((qend - qstart + 1)/qlength, d = 4)]
reprseq.vs.reprseq[, scov := round((send - sstart + 1)/slength, d = 4)]
reprseq.vs.reprseq <- reprseq.vs.reprseq[qname %chin% all.reprseq & sname %chin% all.reprseq]

# Load reprseq all-by-all cov table
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

# adjust qcov/scov table to match the high pident table
reprseq.vs.reprseq.filtered.qcov.scov <- reprseq.cov[.(reprseq.vs.reprseq.filtered[, .(qname, sname)])]
reprseq.vs.reprseq.filtered.qcov.scov[, qcov := pmax(qcov_min, qcov_max, na.rm = T)]
reprseq.vs.reprseq.filtered.qcov.scov[, scov := pmax(scov_min, scov_max, na.rm = T)]

reprseq.vs.reprseq.full <- cbind(reprseq.vs.reprseq.filtered, reprseq.vs.reprseq.filtered.qcov.scov %>% select(pident.full = pident, qcov.full = qcov, scov.full = scov))
reprseq.vs.reprseq.full[, pident.full := pident.full*100]
reprseq.vs.reprseq.full[, q.hit.length := qend - qstart + 1]
reprseq.vs.reprseq.full[, s.hit.length := send - sstart + 1]

recent.recombination.table <- reprseq.vs.reprseq.full %>% 
  filter(q.hit.length >= RECENT.REC.MIN.Q.HIT.LENGTH &
           s.hit.length >= RECENT.REC.MIN.S.HIT.LENGTH &
           pident >= RECENT.REC.MIN.PIDENT &
           qcov.full <= RECENT.REC.MAX.QCOV &
           scov.full <= RECENT.REC.MAX.SCOV)
recent.recombination.table.final <- recent.recombination.table %>% 
  select(qname, sname, qstart, qend, qlength, sstart, send, slength,
         pident, qcov, qcov.full, scov, scov.full, 
         hit.length.query = q.hit.length, hit.length.subject = s.hit.length)

# write to a file
fwrite(recent.recombination.table.final, file = reprseq.recent.rec.table.output)
         

