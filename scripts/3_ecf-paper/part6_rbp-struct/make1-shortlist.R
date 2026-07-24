### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(patchwork)))
suppressWarnings(suppressPackageStartupMessages(library(ggh4x)))
suppressWarnings(suppressPackageStartupMessages(library(Biostrings)))
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
this.part.dir <- file.path(ecf.proj.dir, "part6_rbp-struct")
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

# Load ECFs vs ECOF
ecfs.vs.ecod <- fread(ecfs.vs.ecod.output)

# Protein arch metrics
reprseq.arch.metrics <- fread(reprseq.arch.metrics.output)


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

setkey(reprseq.arch.metrics, qname)

###########################################
###########################################

#############################################################
### Visualise ECF/EF architecture of many filtered RBPs
#############################################################

# sort ECOD hits
setorder(reprseq.vs.ecod, qname, t.id, -bitscore)
reprseq.vs.ecod.top <- reprseq.vs.ecod[, .SD[1], by = .(qname, t.id)]

# list all unique family combinations per protein
reprseqs.fam <- unique(reprseq.vs.ecfs[, .(qname, family)])
setkey(reprseqs.fam, qname)
reprseqs.fam.comb <- reprseqs.fam[, .(fams = list(sort(unique(family)))), by = qname]
reprseqs.fam.comb[, fam.key := vapply(fams, paste, collapse = "|", FUN.VALUE = character(1))]
reprseqs.fam.comb[, fam.first := sapply(fams, '[', 1)]
setkey(reprseqs.fam.comb, qname)

# Filtering criteria
MIN.N.FAMS=2
MIN.LENGTH=600
MIN.N.ECF.HITS=20

rbp.strings <- c("tail fiber", "tail spike")
reprseq.rbps <- reprseq.functions.unique[funct %chin% rbp.strings]$reprseq
rbps.vs.ecfs <- reprseq.vs.ecfs[qname %chin% reprseq.rbps]
all.rbps <- unique(rbps.vs.ecfs$qname)
rbps.arch.metrics <- reprseq.arch.metrics[qname %chin% all.rbps]
rbps.arch.metrics.filtered <- rbps.arch.metrics[n.fams >= MIN.N.FAMS & length >= MIN.LENGTH & n.ecf.hits >= MIN.N.ECF.HITS]
rbps.arch.metrics.filtered$fam.key <- reprseqs.fam.comb[.(rbps.arch.metrics.filtered$qname)]$fam.key
setorder(rbps.arch.metrics.filtered, fam.key, -n.ecf.hits)
rbps.arch.metrics.filtered.top.arch <- rbps.arch.metrics.filtered[, .SD[1], by = .(fam.key)]
repr.rbps <- rbps.arch.metrics.filtered.top.arch$qname


rbp.vis.dir <- file.path(this.part.dir, "vis-rbp-archs")
make.dir(rbp.vis.dir)

for(ri in 1:length(repr.rbps)){
  this.ri <-  repr.rbps[ri]
  this.ri.ecf.hits <- reprseq.vs.ecfs[qname == this.ri]
  this.ri.ecod.hits <- reprseq.vs.ecod.top[qname == this.ri]
  this.ri.ecod.hits.formatted <- this.ri.ecod.hits %>% select(qname, qstart, qend, qlength, sname, sstart, send, slength, pident, bitscore, eval, prob, pval, hit.length, qcov, scov, family = t.id)
  this.ri.all.hits <- rbind(this.ri.ecf.hits, this.ri.ecod.hits.formatted)
  this.ri.plot <- plot.ecf.distr.over.query(this.ri, this.ri.all.hits)
  this.ri.plot.outfile <- file.path(rbp.vis.dir, sprintf("%s.png", this.ri))
  ggsave(this.ri.plot.outfile, this.ri.plot, width = 7, height = 6)  
}


# shortlist2 <- rbps.arch.metrics.filtered.top.arch[grepl("EF0033", fam.key) & 
#                                       grepl("EF0080", fam.key) &
#                                       grepl("EF0079", fam.key)]$qname
# shortlist3 <- rbps.arch.metrics.filtered.top.arch[grepl("EF0212", fam.key) & 
#                                                    grepl("EF0255", fam.key)]$qname
# 
# sort(shortlist2)
# sort(shortlist3)


# this.shortlist.dir <- file.path(rbp.vis.dir, "af-model-shortlist2")
# make.dir(this.shortlist.dir)
# for(ri in 1:length(shortlist2)){
#   this.ri <-  shortlist2[ri]
#   this.ri.ecf.hits <- reprseq.vs.ecfs[qname == this.ri]
#   this.ri.ecod.hits <- reprseq.vs.ecod.top[qname == this.ri]
#   this.ri.ecod.hits.formatted <- this.ri.ecod.hits %>% select(qname, qstart, qend, qlength, sname, sstart, send, slength, pident, bitscore, eval, prob, pval, hit.length, qcov, scov, family = t.id)
#   this.ri.all.hits <- rbind(this.ri.ecf.hits, this.ri.ecod.hits.formatted)
#   this.ri.plot <- plot.ecf.distr.over.query(this.ri, this.ri.all.hits)
#   this.ri.plot.outfile <- file.path(this.shortlist.dir, sprintf("%s.png", this.ri))
#   ggsave(this.ri.plot.outfile, this.ri.plot, width = 7, height = 6)  
# }

af.model.set <- c("reprseq019028", "reprseq000317", "reprseq012341", "reprseq031110", "reprseq048029", "reprseq057378", "reprseq125186")
af.model.set.dir <- file.path(rbp.vis.dir, "af-model-set")
make.dir(af.model.set.dir)

reprseq.seqs <- readAAStringSet(reprseq.seqs.filename)
reprseq.seqs.n <- names(reprseq.seqs)
reprseq.seqs.af <- reprseq.seqs[reprseq.seqs.n %in% af.model.set]
writeXStringSet(reprseq.seqs.af, filepath = file.path(af.model.set.dir, "af-model-set.fasta"))


for(ri in 1:length(af.model.set)){
  this.ri <-  af.model.set[ri]
  this.ri.ecf.hits <- reprseq.vs.ecfs[qname == this.ri]
  this.ri.ecod.hits <- reprseq.vs.ecod.top[qname == this.ri]
  this.ri.ecod.hits.formatted <- this.ri.ecod.hits %>% select(qname, qstart, qend, qlength, sname, sstart, send, slength, pident, bitscore, eval, prob, pval, hit.length, qcov, scov, family = t.id)
  this.ri.all.hits <- rbind(this.ri.ecf.hits, this.ri.ecod.hits.formatted)
  this.ri.plot <- plot.ecf.distr.over.query(this.ri, this.ri.all.hits)
  this.ri.plot.outfile <- file.path(af.model.set.dir, sprintf("%s.png", this.ri))
  ggsave(this.ri.plot.outfile, this.ri.plot, width = 7, height = 6)  
}
