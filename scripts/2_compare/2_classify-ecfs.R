### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

ecfs.vs.qecfs <- fread(ecf.vs.qecf.output)
all.ecfs <- unique(ecfs.vs.qecfs$ecf)

ecfs.vs.ecod <- fread(ecfs.vs.ecod.output)
ecfs.with.ecod.hit <- unique(ecfs.vs.ecod$qname)

ecfs.equivalent <- ecfs.vs.ecod %>% 
  filter(prob >= MIN.PROB.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(qcov >= MIN.QCOV.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(scov >= MIN.SCOV.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE)

ecfs.subdomain <- ecfs.vs.ecod %>% 
  filter(prob >= MIN.PROB.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(qcov >= MIN.QCOV.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(scov <= MAX.SCOV.FOR.SUBDOMAIN) %>%
  filter(eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE)

ecfs.potentially.new <- setdiff(all.ecfs, ecfs.with.ecod.hit)


# ecfs.vs.ecod



# subdomain = ecfs.hits.vs.domain.described %>% 
#   filter(prob >= min.prob &  ecf.cov >= min.ecf.cov & domain.cov <= max.subdomain.cov & eval <= max.eval) %>%
#   distinct(ecf, t_id, t_name, h_id, h_name, ecf.length)



# # parameters for comparing ECFs to ECODs
# MIN.PROB.FOR.DOMAIN.EQUIVALENCE = 95
# MIN.QCOV.FOR.DOMAIN.EQUIVALENCE = 0.9
# MIN.SCOV.FOR.DOMAIN.EQUIVALENCE = 0.9
# MAX.EVAL.FOR.DOMAIN.EQUIVALENCE = 10^(-3)
# MAX.SCOV.FOR.SUBDOMAIN = 0.6


