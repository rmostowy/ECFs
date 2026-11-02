### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

# # Set up paths and directories
# Sys.setenv(PATH = paste("/opt/bio/hhsuite/bin", "/opt/bio/hhsuite/scripts", Sys.getenv("PATH"), sep = ":"),
#            HHLIB = "/opt/bio/hhsuite")

# db.a3m <- file.path(database.dir, "a3m")
# make.dir(db.a3m)


#################################
### Load all data:
### 1. ECF tables
### 3. domain tables
### 4. ECF vs domain tables
### 5. reprseq vs ecfs
### 6. reprseq vs pfam
### 7. reprseq functional annotations
#################################

# ECF tables
ecf.table         <- fread(ecf.table.output)
ecf.family.table  <- fread(ecf.family.table.output)

# domain tables
ecod.domains      <- fread(ecod.domains.output)
pfam.domains      <- fread(pfam.domains.output)

# ECF vs domain tables
ecfs.vs.ecod          <- fread(ecfs.vs.ecod.output)
ecfs.vs.pfam          <- fread(ecfs.vs.pfam.output)

# reprseq vs ecfs
reprseq.vs.ecfs.all   <- fread(reprseq.ecfs.output)

# reprseq functional annotations
reprseq.functions.unique <- fread(reprseq.functions.unique.output)
reprseq.functions.all    <- fread(reprseq.functions.all.output)


#################################
### Process loaded tables
#################################
## Keys and lookups used throughout. Both source files keyed `ecf.family.table`
## on `ecf` and `ecf.table` on `ecf.name` -- both are set explicitly here so
## every downstream `.( )` lookup is unambiguous regardless of analysis order.
setkey(ecf.table,                  ecf.name)
setkey(ecf.family.table,           ecf)
setkey(ecod.domains,               full.id)
setkey(pfam.domains,               accession)
setkey(reprseq.functions.unique,   reprseq)

all.ecfs <- unique(ecf.table$ecf.name)
Ne <- length(all.ecfs)

## Attach ECOD hierarchy levels and confidence-filter the ECF-vs-domain hit
## tables (former part4 process step; used in Steps 2 and 3 below).
ecfs.vs.ecod[, x.id    := ecod.domains[.(ecfs.vs.ecod$sname)]$x_id]
ecfs.vs.ecod[, h.id    := ecod.domains[.(ecfs.vs.ecod$sname)]$h_id]
ecfs.vs.ecod[, t.id    := ecod.domains[.(ecfs.vs.ecod$sname)]$t_id]
ecfs.vs.ecod[, x.level := ecod.domains[.(ecfs.vs.ecod$sname)]$x.level]
ecfs.vs.ecod[, h.level := ecod.domains[.(ecfs.vs.ecod$sname)]$h.level]
ecfs.vs.ecod[, t.level := ecod.domains[.(ecfs.vs.ecod$sname)]$t.level]

ecfs.vs.ecod.conf <- ecfs.vs.ecod %>%
  filter(eval <= MAX.EVALUE.FOR.ECF.DOMAIN.CONFIDENT &
           prob >= MIN.PROB.FOR.ECF.DOMAIN.CONFIDENT)
ecfs.vs.pfam.conf <- ecfs.vs.pfam %>%
  filter(eval <= MAX.EVALUE.FOR.ECF.DOMAIN.CONFIDENT &
           prob >= MIN.PROB.FOR.ECF.DOMAIN.CONFIDENT)

## Confidence-filter the rHMM-vs-domain hit tables (former part4 process step;
## used in Step 3 below).
reprseq.vs.ecfs <- reprseq.vs.ecfs.all %>%
  filter(prob >= MIN.PROB.FOR.DOMAIN.DETECT.IN.RHMM,
         eval <= MAX.EVAL.FOR.DOMAIN.DETECT.IN.RHMM,
         qcov >= MIN.RHMM.COV.FOR.DOMAIN.DETECT.IN.RHMM,
         scov >= MIN.DOMAIN.COV.FOR.DOMAIN.DETECT.IN.RHMM,
         hit.length >= MIN.HIT.LENGTH.FOR.DOMAIN.DETECT.IN.RHMM)

setkey(ecfs.vs.ecod.conf, qname)
setkey(ecfs.vs.pfam.conf, qname)
setkey(reprseq.vs.ecfs, sname)

ecf.metadata.list <- lapply(1:Ne, function(ecf.index){
  top.ecod.hit.tid <- top.ecod.hit.tname <- NA
  top.pfam.hit.accession <- top.pfam.hit.description <- NA
  top.functions <- NA 
  
  this.ecf <- all.ecfs[ecf.index]
  this.ecf.table.metadata <- ecf.table[.(this.ecf)]
  this.ef <- ecf.family.table[.(this.ecf)]$family
  
  
  # Get top ECOD hit
  this.ecf.ecods <- ecfs.vs.ecod.conf[.(this.ecf)]
  setorder(this.ecf.ecods, -bitscore, -prob)
  if(nrow(this.ecf.ecods)>0){
    top.ecod.hit.tid <- this.ecf.ecods$t.id[1]
    top.ecod.hit.tname <- ecod.domains[t_id == top.ecod.hit.tid]$t_name[1]
  }
  
  # Get top Pfam hit
  this.ecf.pfam <- ecfs.vs.pfam.conf[.(this.ecf)]
  setorder(this.ecf.pfam, -bitscore, -prob)
  if(nrow(this.ecf.pfam)>0){
    top.pfam.hit.accession <- this.ecf.pfam$accession[1]
    top.pfam.hit.description <- pfam.domains[accession == top.pfam.hit.accession]$description[1]
  }
  
  # Get top Functions
  reprseq.with.this.ecf <- reprseq.vs.ecfs[.(this.ecf)]$qname
  this.ecf.top.functions <- reprseq.functions.unique[.(reprseq.with.this.ecf)]
  this.ecf.top.functions <- this.ecf.top.functions[!is.na(score)]
  setorder(this.ecf.top.functions, -score)
  if(nrow(this.ecf.top.functions)>0){
    this.ecf.max.score <- this.ecf.top.functions$score[1]
    this.ecf.top.funct.hits <- unique(this.ecf.top.functions[score >= this.ecf.max.score]$funct)
    top.functions <- paste(this.ecf.top.funct.hits, collapse = "|")
  }
  
  data.table(ecf.name = this.ecf.table.metadata$ecf.name,
             EF = this.ef,
             ecf.length = this.ecf.table.metadata$length,
             ecf.size = this.ecf.table.metadata$size,
             top.ecod.hit.tid,
             top.ecod.hit.tname,
             top.pfam.hit.accession,
             top.pfam.hit.description,
             top.functions
  )
})

ecf.metadata <- rbindlist(ecf.metadata.list)
ecf.metadata.outfile <- file.path(database.dir, "ecf-metadata.csv")
fwrite(ecf.metadata, file = ecf.metadata.outfile)







