### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
# suppressWarnings(suppressPackageStartupMessages(library(seqinr)))
# suppressWarnings(suppressPackageStartupMessages(library(Biostrings)))
# suppressWarnings(suppressPackageStartupMessages(library(parallel)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

# Load refseq table
reprseq.table <- fread(reprseq.table.filename)
all.reprseq <- unique(reprseq.table$reprseq)

# Calculate NCBI overlaps per reprseq
max.overlap.per.reprseq <- reprseq.table[, .(reprseq, ncbi.overlap)]
setorder(max.overlap.per.reprseq, reprseq, -ncbi.overlap)
max.overlap.per.reprseq <- max.overlap.per.reprseq[, .SD[1], by = reprseq]
setkey(max.overlap.per.reprseq, reprseq)

# Load recent recombination table
recent.rec.table <- fread(reprseq.recent.rec.table.output)
recent.rec.table[, cov.max.full := pmax(qcov.full, scov.full, na.rm = T)]

# Load ECOD domain table
ecod.domains <- fread(ecod.domains.output)
setkey(ecod.domains, full.id)

# Load PFAM domain table
pfam.domains <- fread(pfam.domains.output)
setkey(pfam.domains, accession)

# Load ECF families
ecf.family.table <- fread(ecf.family.table.output)
setkey(ecf.family.table, ecf)

# Load ECF table
ecf.class.table <- fread(ecf.classification.table.output)
setkey(ecf.class.table, ecf)

ecf.class.table$family <- ecf.family.table[.(ecf.class.table$ecf)]$family
ecf2ecod <- ecf.class.table %>% select(ecf, family, t.id)
ecf2ecod$t.name <- sapply(1:nrow(ecf2ecod), function(k){
  out <- "unknown"
  this.t.id <- ecf2ecod$t.id[k]
  if(this.t.id != "unknown") out <- ecod.domains[t_id == this.t.id]$t_name[1]  
  return(out)
})
setorder(ecf2ecod, family)

# Load functional annotations
reprseq.functions <- fread(reprseq.functions.unique.output)
setkey(reprseq.functions, reprseq)

# Load recombination domain overlap table
rec.domain.overlaps <- fread("rec-domain-overlaps.txt")

get.accession <- function(string) sub("\\.\\d+$", "", string)


rec.domain.overlaps$prot1.ecod.id <- ecod.domains[.(rec.domain.overlaps$protein1.ecod.best)]$t_id
rec.domain.overlaps.ecod <- rec.domain.overlaps[which(protein1.ecod.score>0)]

most.frequent.ecods <- names(which(my.table(rec.domain.overlaps.ecod$prot1.ecod.id)>30))
most.frequent.ecods.names <- sapply(1:length(most.frequent.ecods), function(index) ecod.domains[t_id == most.frequent.ecods[index]]$t_name[1])

index=1
most.frequent.ecods.names[index]
hist(rec.domain.overlaps[prot1.ecod.id == most.frequent.ecods[index]]$protein1.ecod.score)

rbindlist(lapply(1:length(most.frequent.ecods), function(index){
  score.distr <- rec.domain.overlaps[prot1.ecod.id == most.frequent.ecods[index]]$protein1.ecod.score
  score.distr.quant <- quantile(score.distr)  
  data.table(most.frequent.ecods[index], score.med = score.distr.quant[3], score.low = score.distr.quant[2], score.high = score.distr.quant[4])
}))


# hs.ecod.domains <- rec.domain.overlaps[protein1.ecod.score > 0.8]$protein1.ecod.best
# ls.ecod.domains <- rec.domain.overlaps[protein1.ecod.score < 0.5 & protein1.ecod.score > 0]$protein1.ecod.best
# 
# my.table(ecod.domains[.(hs.ecod.domains)]$t_name)[1:50]
# my.table(ecod.domains[.(ls.ecod.domains)]$t_name)[1:50]
# 
# hs.pfam.domains <- rec.domain.overlaps[protein1.pfam.score > 0.8]$protein1.pfam.best
# hs.pfam.domains <- sapply(hs.pfam.domains, get.accession, USE.NAMES = F)
# ls.pfam.domains <- rec.domain.overlaps[protein1.pfam.score < 0.5 & protein1.pfam.score > 0]$protein1.pfam.best
# ls.pfam.domains <- sapply(ls.pfam.domains, get.accession, USE.NAMES = F)
# 
# my.table(pfam.domains[.(hs.pfam.domains)]$description)[1:50]
# my.table(pfam.domains[.(ls.pfam.domains)]$description)[1:50]


# x <- rec.domain.overlaps[protein1.name == "reprseq000005"]
# x.p <- sub("\\.\\d+$", "", x$protein1.pfam.best)
# pfam.domains[.(x.p)]
# x.e <- x$protein1.ecod.best
# ecod.domains[.(x.e)]
# ecf.family.table[.(x$protein1.ecf.best)]



####################################################################
####################################################################
####################################################################




ggplot(rec.domain.overlaps, aes(x = protein1.ecf.score, y = protein2.ecf.score)) + geom_point()


reprseq.ecfs.all <- fread(reprseq.ecfs.output)
setkey(reprseq.ecfs.all, qname)

# rec.domain.overlaps[protein1.ecf.score>0.9 & protein2.ecf.score == 0]
weird.index <- which(rec.domain.overlaps$protein1.ecf.score>0.9 & rec.domain.overlaps$protein2.ecf.score == 0)
# recent.rec.table[ri]

weird.index.can.be.explained <- sapply(1:length(weird.index), function(index){
  this.rr <- recent.rec.table[weird.index[index]]
  this.rr.overlaps <- rec.domain.overlaps[weird.index[index]]
  relevant.ecf <- this.rr.overlaps$protein1.ecf.best
  all.ecfs.other.protein <- reprseq.ecfs.all[.(this.rr$sname)]$sname
  relevant.ecf %in% all.ecfs.other.protein  
})

very.weird.index <- weird.index[!weird.index.can.be.explained]
recent.rec.table[very.weird.index]






#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################








# # Load reprseq vs ECFs
# reprseq.ecfs.all <- fread(reprseq.ecfs.output)
# reprseq.ecfs <- reprseq.ecfs.all %>%
#   filter(prob >= MIN.PROB.FOR.DOMAIN.EQUIVALENCE) %>%
#   filter(scov >= MIN.SCOV.FOR.DOMAIN.EQUIVALENCE) %>%
#   filter(eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE)
# reprseq.ecfs$EF <- ecf.family.table[.(reprseq.ecfs$sname)]$family
# setkey(reprseq.ecfs, qname)
# 
# # Load reprseq vs ECOD
# reprseq.ecod.all <- fread(reprseq.domains.ecod.output)
# reprseq.ecod <- reprseq.ecod.all %>%
#   filter(prob >= MIN.PROB.FOR.DOMAIN.EQUIVALENCE) %>%
#   filter(eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE)
# reprseq.ecod$t.id <- ecod.domains[.(reprseq.ecod$sname)]$t_id
# setkey(reprseq.ecod, qname)
# 
# # Load reprseq vs Pfam
# reprseq.pfam.all <- fread(reprseq.domains.pfam.output)
# reprseq.pfam <- reprseq.pfam.all %>%
#   filter(prob >= MIN.PROB.FOR.DOMAIN.EQUIVALENCE) %>%
#   filter(eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE)
# setkey(reprseq.pfam, qname)
# 
# 
# get.domain.overlap.score <- function(rec.start, rec.end, domain.start, domain.end, domain.cov){
#   domain.overlap.parms <- calculate_overlap_parms(rec.start, rec.end, domain.start, domain.end)
#   domain.hit.length <- domain.end - domain.start + 1
#   prop.of.domain.hit.that.overlaps.with.rec <- domain.overlap.parms$overlap.length/domain.hit.length
#   domain.overlap.score <- prop.of.domain.hit.that.overlaps.with.rec * domain.cov
#   return(domain.overlap.score)
# }
# 
# 
# 
# rec.domain.overlap.table.list <- lapply(1:nrow(recent.rec.table), function(pair.index){
#   this.rr <- recent.rec.table[pair.index]
#   
#   # overlaps with protein1
#   protein1.name <- this.rr$qname
#   protein1.length <- this.rr$qlength
#   aln1.start <- this.rr$qstart
#   aln1.end <- this.rr$qend
#   
#   # ECOD
#   protein1.ecod.hits <- reprseq.ecod[.(protein1.name)]
#   protein1.ecod.hits[, rr.overlap.jaccard := calculate_overlap_jaccard(aln1.start, aln1.end, qstart, qend)]
#   setorder(protein1.ecod.hits, -rr.overlap.jaccard)
#   protein1.ecod.best.overlap <- protein1.ecod.hits[1]
#   protein1.ecod.overlap.score <- get.domain.overlap.score(rec.start = aln1.start, rec.end = aln1.end,
#                                                           domain.start = protein1.ecod.best.overlap$qstart, domain.end = protein1.ecod.best.overlap$qend,
#                                                           domain.cov = protein1.ecod.best.overlap$scov)
#   protein1.ecod.best.hit <- protein1.ecod.best.overlap$sname
#   
#   # PFAM
#   protein1.pfam.hits <- reprseq.pfam[.(protein1.name)]
#   protein1.pfam.hits[, rr.overlap.jaccard := calculate_overlap_jaccard(aln1.start, aln1.end, qstart, qend)]
#   setorder(protein1.pfam.hits, -rr.overlap.jaccard)
#   protein1.pfam.best.overlap <- protein1.pfam.hits[1]
#   protein1.pfam.overlap.score <- get.domain.overlap.score(rec.start = aln1.start, rec.end = aln1.end,
#                                                           domain.start = protein1.pfam.best.overlap$qstart, domain.end = protein1.pfam.best.overlap$qend,
#                                                           domain.cov = protein1.pfam.best.overlap$scov)
#   protein1.pfam.best.hit <- protein1.pfam.best.overlap$sname
#   
#   # ECF
#   protein1.ecf.hits <- reprseq.ecfs[.(protein1.name)]
#   protein1.ecf.hits[, rr.overlap.jaccard := calculate_overlap_jaccard(aln1.start, aln1.end, qstart, qend)]
#   setorder(protein1.ecf.hits, -rr.overlap.jaccard)
#   protein1.ecf.best.overlap <- protein1.ecf.hits[1]
#   protein1.ecf.overlap.score <- get.domain.overlap.score(rec.start = aln1.start, rec.end = aln1.end,
#                                                          domain.start = protein1.ecf.best.overlap$qstart, domain.end = protein1.ecf.best.overlap$qend,
#                                                          domain.cov = protein1.ecf.best.overlap$scov)
#   protein1.ecf.best.hit <- protein1.ecf.best.overlap$sname
#   
#   
#   # overlaps with protein2
#   protein2.name <- this.rr$sname
#   protein2.length <- this.rr$slength
#   aln2.start <- this.rr$sstart
#   aln2.end <- this.rr$send
#   
#   # ECOD
#   protein2.ecod.hits <- reprseq.ecod[.(protein2.name)]
#   protein2.ecod.hits[, rr.overlap.jaccard := calculate_overlap_jaccard(aln1.start, aln1.end, qstart, qend)]
#   setorder(protein2.ecod.hits, -rr.overlap.jaccard)
#   protein2.ecod.best.overlap <- protein2.ecod.hits[1]
#   protein2.ecod.overlap.score <- get.domain.overlap.score(rec.start = aln1.start, rec.end = aln1.end,
#                                                           domain.start = protein2.ecod.best.overlap$qstart, domain.end = protein2.ecod.best.overlap$qend,
#                                                           domain.cov = protein2.ecod.best.overlap$scov)
#   protein2.ecod.best.hit <- protein2.ecod.best.overlap$sname
#   
#   # PFAM
#   protein2.pfam.hits <- reprseq.pfam[.(protein2.name)]
#   protein2.pfam.hits[, rr.overlap.jaccard := calculate_overlap_jaccard(aln1.start, aln1.end, qstart, qend)]
#   setorder(protein2.pfam.hits, -rr.overlap.jaccard)
#   protein2.pfam.best.overlap <- protein2.pfam.hits[1]
#   protein2.pfam.overlap.score <- get.domain.overlap.score(rec.start = aln1.start, rec.end = aln1.end,
#                                                           domain.start = protein2.pfam.best.overlap$qstart, domain.end = protein2.pfam.best.overlap$qend,
#                                                           domain.cov = protein2.pfam.best.overlap$scov)
#   protein2.pfam.best.hit <- protein2.pfam.best.overlap$sname
#   
#   # ECF
#   protein2.ecf.hits <- reprseq.ecfs[.(protein2.name)]
#   protein2.ecf.hits[, rr.overlap.jaccard := calculate_overlap_jaccard(aln1.start, aln1.end, qstart, qend)]
#   setorder(protein2.ecf.hits, -rr.overlap.jaccard)
#   protein2.ecf.best.overlap <- protein2.ecf.hits[1]
#   protein2.ecf.overlap.score <- get.domain.overlap.score(rec.start = aln1.start, rec.end = aln1.end,
#                                                          domain.start = protein2.ecf.best.overlap$qstart, domain.end = protein2.ecf.best.overlap$qend,
#                                                          domain.cov = protein2.ecf.best.overlap$scov)
#   protein2.ecf.best.hit <- protein2.ecf.best.overlap$sname
#   
#   this.rr.out <- data.table(protein1.name, 
#                             protein1.ecod.best = protein1.ecod.best.hit, protein1.ecod.score = round(protein1.ecod.overlap.score, d=2),
#                             protein1.pfam.best = protein1.pfam.best.hit, protein1.pfam.score = round(protein1.pfam.overlap.score, d=2),
#                             protein1.ecf.best = protein1.ecf.best.hit, protein1.ecf.score = round(protein1.ecf.overlap.score, d=2),
#                             protein2.name, 
#                             protein2.ecod.best = protein2.ecod.best.hit, protein2.ecod.score = round(protein2.ecod.overlap.score, d=2),
#                             protein2.pfam.best = protein2.pfam.best.hit, protein2.pfam.score = round(protein2.pfam.overlap.score, d=2),
#                             protein2.ecf.best = protein2.ecf.best.hit, protein2.ecf.score = round(protein2.ecf.overlap.score, d=2)
#   )
#   return(this.rr.out)  
# })
# 
# rec.domain.overlap.table <- rbindlist(rec.domain.overlap.table.list)
# fwrite(rec.domain.overlap.table, file = "rec-domain-overlaps.txt")
# 
