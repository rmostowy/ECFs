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
this.part.dir <- file.path(ecf.proj.dir, "part6_so-what")
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


shortlist <- rbps.arch.metrics.filtered.top.arch[grepl("EF0033", fam.key) & 
                                      grepl("EF0080", fam.key) &
                                      grepl("EF0079", fam.key)]$qname
sort(shortlist)

af.model.set <- c("reprseq019028", "reprseq000317", "reprseq012341", "reprseq031110", "reprseq048029", "reprseq057378", "reprseq125186")
af.model.set.dir <- file.path(rbp.vis.dir, "af-model-set")
make.dir(af.model.set.dir)

reprseq.seqs <- readAAStringSet(reprseq.seqs.filename)
reprseq.seqs.n <- names(reprseq.seqs)
reprseq.seqs.af <- reprseq.seqs[reprseq.seqs.n %in% af.model.set]
writeXStringSet(reprseq.seqs.af, filepath = file.path(af.model.set.dir, "af-model-shortlist1.fasta"))


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



# ###########################################
# ###########################################
# 
# # 1. Find all RBPs with PHROGs annotations
# # 2. Identify all EFs in those reprseqs
# # 3. Sort them by position
# # 4. Extract all EFs present in the first 
# #    100 AAs, but preserve ECF informaiton
# # 5. Optionally: exclude any known ECODs
# # 6. Make a list of all EFs which are
# #    putative N-termini
# # 7. Sort them by frequency and pick top Ntop
# # 8. For each N-terminal EF, identify
# #    unique EF combinations which contain
# #    multiple EFs and N-terminal as first
# # 9. For all reprseqs with the same EF combo,
# #    pick one with the highest number of ECFs
# # 10.Save unique sets of modular RBPs with
# #    representative proteins and plot the 
# #    ECF distribution in each, by including 
# #    ECODs as well
# 
# 
# # 1. Find all RBPs with PHROGs annotations
# 
# rbp.strings <- c("tail fiber", "tail spike")
# reprseq.rbps <- reprseq.functions.unique[funct %chin% rbp.strings]$reprseq
# 
# # 2. Identify all EFs in those reprseqs
# rbps.vs.ecfs <- reprseq.vs.ecfs[qname %chin% reprseq.rbps]
# 
# # 3. Sort them by position
# setorder(rbps.vs.ecfs, qname, qstart, family)
# 
# # 4. Extract all EFs present in the first 
# #    100 AAs, but preserve ECF informaiton
# rbps.vs.ecfs.overlap <- copy(rbps.vs.ecfs)
# rbps.vs.ecfs.overlap[, first.hundred.overlap := calculate_overlap_jaccard(qstart, qend, 1, 100)]
# rbps.vs.ecfs.first.hundred <- rbps.vs.ecfs.overlap[first.hundred.overlap>0]
# 
# # 5. Optionally: exclude any known ECODs
# rbps.vs.ecfs.first.hundred <- rbps.vs.ecfs.first.hundred[!sname %chin% ecfs.vs.ecod$qname]
# 
# # 6. Make a list of all EFs which are
# #    putative N-termini
# rbp.nterm.candidate.fam <- unique(rbps.vs.ecfs.first.hundred[, .(qname, family)])
# 
# # 7. Sort them by frequency and pick top Ntop
# Ntop <- 10
# rbp.nterm.candidate.fam.top <- rbp.nterm.candidate.fam[, .N, by = family][order(-N)][1:Ntop]
# 
# # 8. For each N-terminal EF, identify
# #    unique EF combinations which contain
# #    multiple EFs and N-terminal as first
# # 9. For all reprseqs with the same EF combo,
# #    pick one with the highest number of ECFs
# 
# # list all unique family combinations per protein
# reprseqs.fam <- unique(reprseq.vs.ecfs[, .(qname, family)])
# setkey(reprseqs.fam, qname)
# reprseqs.fam.comb <- reprseqs.fam[, .(fams = list(sort(unique(family)))), by = qname]
# reprseqs.fam.comb[, fam.key := vapply(fams, paste, collapse = "|", FUN.VALUE = character(1))]
# reprseqs.fam.comb[, fam.first := sapply(fams, '[', 1)]
# 
# 
# efi=4
# this.fam <- rbp.nterm.candidate.fam.top$family[efi]
# this.fam.reprseqs <- unique(rbps.vs.ecfs[family == this.fam]$qname)
# this.fam.reprseqs.fam.combs <- reprseqs.fam.comb[qname %chin% this.fam.reprseqs]
# this.fam.reprseqs.fam.combs$n.ecf.hits <- reprseq.arch.metrics[.(this.fam.reprseqs.fam.combs$qname)]$n.ecf.hits
# this.fam.reprseqs.fam.combs$n.fams <- reprseq.arch.metrics[.(this.fam.reprseqs.fam.combs$qname)]$n.fams
# this.fam.reprseqs.fam.combs.filtered <- this.fam.reprseqs.fam.combs[n.fams > 1]
# # this.fam.reprseqs.fam.combs.filtered <- this.fam.reprseqs.fam.combs[fam.first == this.fam & n.fams > 1]
# setorder(this.fam.reprseqs.fam.combs.filtered, fam.key, -n.ecf.hits)
# this.fam.repr.combs <- this.fam.reprseqs.fam.combs.filtered[, .SD[1], by = .(fam.key)]
# this.fam.repr.combs
# 
# 
# # 10.Save unique sets of modular RBPs with
# #    representative proteins and plot the 
# #    ECF distribution in each, by including 
# #    ECODs as well
# 
# ri=5
# this.ri <-  this.fam.repr.combs[ri]$qname
# this.ri.ecf.hits <- reprseq.vs.ecfs[qname == this.ri]
# this.ri.ecod.hits <- reprseq.vs.ecod[qname == this.ri]
# this.ri.ecod.hits


###########################################
###########################################
# 
# # 1. Identify EF which corresponds to depo
# 
# depo.t.id <- "207.2.1"
# 
# ecfs.vs.ecod <- fread(ecfs.vs.ecod.output)
# ecfs.vs.ecod <- cbind(ecfs.vs.ecod, t.id = ecod.domains[.(ecfs.vs.ecod$sname)]$t_id)
# ecfs.vs.ecod <- cbind(ecfs.vs.ecod, family = ecf.family.table[.(ecfs.vs.ecod$qname)]$family)
# ecfs.vs.ecod.depo <- ecfs.vs.ecod[t.id == depo.t.id]
# 
# sgnh.t.id <- "2007.5.1"
# ecfs.vs.ecod.sgn <- ecfs.vs.ecod[t.id == sgnh.t.id]
# my.table(ecfs.vs.ecod.sgn$family)
# 
# 
# # 2. Identify reprseqs with the depo fold
# depo.reprseqs <- unique(reprseq.vs.ecod[t.id == depo.t.id]$qname)
# 
# 
# ri=1
# 
# this.rhmm <- depo.reprseqs[ri]
# this.rhmm.ecf <- reprseq.vs.ecfs[qname == this.rhmm]
# this.rhmm.ecod.depo <- reprseq.vs.ecod[qname == this.rhmm & t.id == depo.t.id]
# this.rhmm.depo.range <- reduce_ranges(this.rhmm.ecod.depo)
# plot.ecf.distr.over.query(this.rhmm, this.rhmm.ecf)
# ri <- ri + 1
# 
# # "reprseq057423"
# # "reprseq087038"
# # "reprseq035294"
# 
# this.ef <- "EF0255" #RBP N-terminal
# all.rhmms.with.this.ef <- unique(reprseq.vs.ecfs[family == this.ef]$qname)
# z <- reprseq.vs.ecfs[qname %in% all.rhmms.with.this.ef]
# z[family == "EF0013"]
# 
# ri=1
# this.rhmm <- all.rhmms.with.this.ef[ri]
# this.rhmm.ecf <- reprseq.vs.ecfs[qname == this.rhmm]
# this.rhmm.ecod.depo <- reprseq.vs.ecod[qname == this.rhmm & t.id == depo.t.id]
# this.rhmm.depo.range <- reduce_ranges(this.rhmm.ecod.depo)
# plot.ecf.distr.over.query(this.rhmm, this.rhmm.ecf)
# ri <- ri + 1
# 
# 
# this.rhmm <- "reprseq035294"
# this.rhmm <- "reprseq040908"
# this.rhmm.ecf <- reprseq.vs.ecfs[qname == this.rhmm]
# this.rhmm.ecod.depo <- reprseq.vs.ecod[qname == this.rhmm & t.id == depo.t.id]
# this.rhmm.depo.range <- reduce_ranges(this.rhmm.ecod.depo)
# plot.ecf.distr.over.query(this.rhmm, this.rhmm.ecf)
# 
# 
# # this.rhmm <- "reprseq057423"
# # this.rhmm <- "reprseq087038"
# # this.rhmm <- "reprseq035294"
# reprseq.vs.ecod[qname == this.rhmm]
# 
