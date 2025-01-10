### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################
# Load reprseq tables
reprseq.table <- fread(reprseq.table.filename)
reprseq.table[, accession := sub("\\.[0-9]+:.*", "", cluster.name)]
all.reprseq <- unique(reprseq.table$reprseq)

# Load metadata table
metadata.table <- fread(metadata.table.output)

# Load ECOD domain table
ecod.domains <- fread(ecod.domains.output)
setkey(ecod.domains, full.id)

# Load PHROG annotations
reprseq.functions.unique <- fread(reprseq.functions.unique.output)
setkey(reprseq.functions.unique, reprseq)

# Load ECF table
ecf.name.table <- fread(ecf.name.table.output)
all.ecfs <- unique(ecf.name.table$ecf.name)
Ne <- length(all.ecfs)

# Load ECF vs REPRSEQ table
reprseq.ecfs <- fread(reprseq.ecfs.output)
reprseq.ecfs.conf <- reprseq.ecfs %>% filter(scov >= 0.95 & prob >= 0.95) 
ecf2reprseq <- reprseq.ecfs.conf %>% select(ecf = sname, reprseq = qname)
setkey(ecf2reprseq, ecf)

# Load ECF vs ECOD table
ecfs.vs.ecod <- fread(ecfs.vs.ecod.output)
ecfs.vs.ecod[, score := round(qcov*scov*prob/100*pident/100, d = 6)]

# Add T-level information
ecfs.vs.ecod$t.id <- ecod.domains[.(ecfs.vs.ecod$sname)]$t_id
ecfs.vs.ecod$arch.class <- ecod.domains[.(ecfs.vs.ecod$sname)]$arch.class

setkey(ecfs.vs.ecod, qname)



# ecfs.with.ecod.hit <- unique(ecfs.vs.ecod$qname)
# ecfs.potentially.new <- setdiff(all.ecfs, ecfs.with.ecod.hit)

ecf.classification.table.list <- lapply(1:Ne, function(index){
  # cat(index,"\n")
  # get linked PHROG annotations
  this.ecf <- all.ecfs[index]
  this.ecf.has.reprseqs <- this.ecf %in% ecf2reprseq$ecf
  
  this.ecf.phrogs.funct <- "unknown"
  this.ecf.phrogs.class <- "unknown"
  
  if(this.ecf.has.reprseqs){
    
    this.ecf.reprseqs <- ecf2reprseq[.(this.ecf)]$reprseq
    this.ecf.accessions <- unique(reprseq.table[reprseq %in% this.ecf.reprseqs]$accession)
    this.ecf.metadata <- metadata.table[accession %in% this.ecf.accessions]
    
    # ictv.classes <- unique(this.ecf.metadata$class.ictv)
    # ictv.classes <- ictv.classes[ictv.classes != ""]
    # ictv.families <- unique(this.ecf.metadata$family.ncbi)
    # ictv.classes <- ictv.classes[ictv.classes != ""]
    
    this.ecf.linked.annotations <- reprseq.functions.unique[.(this.ecf.reprseqs)]
    this.ecf.linked.annotations.known <- this.ecf.linked.annotations[this.ecf.linked.annotations$funct != ""]
    
    if(nrow(this.ecf.linked.annotations.known)>0){
      setorder(this.ecf.linked.annotations.known, -score)
      this.ecf.phrogs.funct <- this.ecf.linked.annotations.known$funct[1]
      this.ecf.phrogs.class <- this.ecf.linked.annotations.known$class[1]
    }
  }
  

  # get linked ECODs
  this.ecf.type <- "new"
  this.ecf.arch <- "unknown"
  this.ecf.t.level <- "unknown"
  if(this.ecf %in% ecfs.vs.ecod$qname){
    this.ecf.type <- "similar"
    this.ecfs.vs.ecod <- ecfs.vs.ecod[.(this.ecf)]
    
    this.ecfs.vs.ecod.equivalent <- this.ecfs.vs.ecod %>% 
      filter(prob >= MIN.PROB.FOR.DOMAIN.EQUIVALENCE) %>%
      filter(qcov >= MIN.QCOV.FOR.DOMAIN.EQUIVALENCE) %>%
      filter(scov >= MIN.SCOV.FOR.DOMAIN.EQUIVALENCE) %>%
      filter(eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE)
    is.equivalent <- nrow(this.ecfs.vs.ecod.equivalent)>0
    
    is.subdomain.cond1 <- all(this.ecfs.vs.ecod$scov < MAX.SCOV.FOR.SUBDOMAIN)
    this.ecfs.vs.ecod.subdomain <- this.ecfs.vs.ecod %>% 
      filter(prob >= MIN.PROB.FOR.SUBDOMAIN) %>%
      filter(qcov >= MIN.QCOV.FOR.SUBDOMAIN) %>%
      filter(eval <= MAX.EVAL.FOR.SUBDOMAIN)
    is.subdomain.cond2 <- nrow(this.ecfs.vs.ecod.subdomain)>0
    is.subdomain <- is.subdomain.cond1 & is.subdomain.cond2
    
    if(is.equivalent){
      this.ecfs.vs.ecod <- this.ecfs.vs.ecod.equivalent
      this.ecf.type <- "equivalent"
    } else if(is.subdomain){
      this.ecfs.vs.ecod <- this.ecfs.vs.ecod.subdomain
      this.ecf.type <- "subdomain"
    }
    
    setorder(this.ecfs.vs.ecod, -score)
    this.ecf.arch <- this.ecfs.vs.ecod$arch.class[1]
    this.ecf.t.level <- this.ecfs.vs.ecod$t.id[1]
    
  }
  data.table(ecf = this.ecf, type = this.ecf.type, 
             arch.class = this.ecf.arch, t.id = this.ecf.t.level,
             funct = this.ecf.phrogs.funct, class = this.ecf.phrogs.class)
  
})
ecf.classification.table <- rbindlist(ecf.classification.table.list)

has.functional.hit <- ecf.classification.table$funct != "unknown"
has.domain.hit <- ecf.classification.table$type != "new"

ecf.classification.table$hit.types <- "none"
ecf.classification.table$hit.types[has.functional.hit & !has.domain.hit] <- "F"
ecf.classification.table$hit.types[!has.functional.hit & has.domain.hit] <- "D"
ecf.classification.table$hit.types[has.functional.hit & has.domain.hit] <- "F+D"

fwrite(ecf.classification.table, file = ecf.classification.table.output, sep = "\t")
