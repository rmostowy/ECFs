### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(patchwork)))
suppressWarnings(suppressPackageStartupMessages(library(ggh4x)))
# suppressWarnings(suppressPackageStartupMessages(library(ComplexUpset)))
# suppressWarnings(suppressPackageStartupMessages(library(ggraph)))
# suppressWarnings(suppressPackageStartupMessages(library(broom)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

#################################
### Set up output directory
#################################
this.part.dir <- file.path(ecf.proj.dir, "part5_mosaic-struct")
make.dir(this.part.dir)

#################################
### Load all data:
### 1. ECF tables
### 2. reprseq lengths
### 3. reprseq annotations
### 4. reprseq vs ECF table
### 5. reprseq arch metrics
#################################

# ECF tables
ecf.table <- fread(ecf.table.output)
ecf.family.table <- fread(ecf.family.table.output)

ecod.domains <- fread(ecod.domains.output)

# reprseq lengths
reprseq.lengths <- fread(reprseq.lengths.filename)

# Load functional annotations
reprseq.functions.unique <- fread(reprseq.functions.unique.output)
reprseq.functions <- fread(reprseq.functions.all.output)

# Load reprseq vs ECFs
reprseq.vs.ecfs.all <- fread(reprseq.ecfs.output)

# Load reprseq vs ECOD
reprseq.vs.ecod.all <- fread(reprseq.domains.ecod.output)


#################################
### Process loaded tables
#################################

setkey(ecf.family.table, ecf)

setkey(ecod.domains, full.id)

setnames(reprseq.lengths, "name", "qname")
setkey(reprseq.lengths, qname)

# setkey(reprseq.functions, qname)
setkey(reprseq.functions.unique, reprseq)

reprseq.vs.ecfs <- reprseq.vs.ecfs.all %>% 
  filter(prob >= MIN.PROB.FOR.DOMAIN.DETECT.IN.RHMM,
         eval <= MAX.EVAL.FOR.DOMAIN.DETECT.IN.RHMM, 
         qcov >= MIN.RHMM.COV.FOR.DOMAIN.DETECT.IN.RHMM,
         scov >= MIN.DOMAIN.COV.FOR.DOMAIN.DETECT.IN.RHMM,
         hit.length >= MIN.HIT.LENGTH.FOR.DOMAIN.DETECT.IN.RHMM
  )
reprseq.vs.ecfs <- cbind(reprseq.vs.ecfs, family = ecf.family.table[.(reprseq.vs.ecfs$sname)]$family)
setkey(reprseq.vs.ecfs, qname)

reprseq.vs.ecod <- reprseq.vs.ecod.all %>% 
  filter(prob >= MIN.PROB.FOR.DOMAIN.DETECT.IN.RHMM,
         eval <= MAX.EVAL.FOR.DOMAIN.DETECT.IN.RHMM, 
         qcov >= MIN.RHMM.COV.FOR.DOMAIN.DETECT.IN.RHMM,
         scov >= MIN.DOMAIN.COV.FOR.DOMAIN.DETECT.IN.RHMM,
         hit.length >= MIN.HIT.LENGTH.FOR.DOMAIN.DETECT.IN.RHMM
  )
reprseq.vs.ecod <- cbind(reprseq.vs.ecod, t.id = ecod.domains[.(reprseq.vs.ecod$sname)]$t_id)



###########################################
###########################################

# 1. Identify EF which corresponds to depo

depo.t.id <- "207.2.1"

ecfs.vs.ecod <- fread(ecfs.vs.ecod.output)
ecfs.vs.ecod <- cbind(ecfs.vs.ecod, t.id = ecod.domains[.(ecfs.vs.ecod$sname)]$t_id)
ecfs.vs.ecod <- cbind(ecfs.vs.ecod, family = ecf.family.table[.(ecfs.vs.ecod$qname)]$family)
ecfs.vs.ecod.depo <- ecfs.vs.ecod[t.id == depo.t.id]

sgnh.t.id <- "2007.5.1"
ecfs.vs.ecod.sgn <- ecfs.vs.ecod[t.id == sgnh.t.id]
my.table(ecfs.vs.ecod.sgn$family)


# 2. Identify reprseqs with the depo fold
depo.reprseqs <- unique(reprseq.vs.ecod[t.id == depo.t.id]$qname)


ri=1

this.rhmm <- depo.reprseqs[ri]
this.rhmm.ecf <- reprseq.vs.ecfs[qname == this.rhmm]
this.rhmm.ecod.depo <- reprseq.vs.ecod[qname == this.rhmm & t.id == depo.t.id]
this.rhmm.depo.range <- reduce_ranges(this.rhmm.ecod.depo)
plot.ecf.distr.over.query(this.rhmm, this.rhmm.ecf)
ri <- ri + 1

# "reprseq057423"
# "reprseq087038"
# "reprseq035294"

this.ef <- "EF0255" #RBP N-terminal
all.rhmms.with.this.ef <- unique(reprseq.vs.ecfs[family == this.ef]$qname)
z <- reprseq.vs.ecfs[qname %in% all.rhmms.with.this.ef]
z[family == "EF0013"]

ri=1
this.rhmm <- all.rhmms.with.this.ef[ri]
this.rhmm.ecf <- reprseq.vs.ecfs[qname == this.rhmm]
this.rhmm.ecod.depo <- reprseq.vs.ecod[qname == this.rhmm & t.id == depo.t.id]
this.rhmm.depo.range <- reduce_ranges(this.rhmm.ecod.depo)
plot.ecf.distr.over.query(this.rhmm, this.rhmm.ecf)
ri <- ri + 1


this.rhmm <- "reprseq035294"
this.rhmm <- "reprseq040908"
this.rhmm.ecf <- reprseq.vs.ecfs[qname == this.rhmm]
this.rhmm.ecod.depo <- reprseq.vs.ecod[qname == this.rhmm & t.id == depo.t.id]
this.rhmm.depo.range <- reduce_ranges(this.rhmm.ecod.depo)
plot.ecf.distr.over.query(this.rhmm, this.rhmm.ecf)


# this.rhmm <- "reprseq057423"
# this.rhmm <- "reprseq087038"
# this.rhmm <- "reprseq035294"
reprseq.vs.ecod[qname == this.rhmm]

