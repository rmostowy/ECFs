### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(igraph)))

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

# Filter hits 
ecfs.vs.ecod <- fread(ecfs.vs.ecod.output)
ecfs.full.vs.ecod <- ecfs.vs.ecod %>%
  filter(qcov >= MIN.QCOV.FOR.SUBDOMAIN) %>%
  filter(prob >= MIN.PROB.FOR.SUBDOMAIN) %>%
  filter(eval <= MAX.EVAL.FOR.SUBDOMAIN)
setkey(ecfs.full.vs.ecod, qname)

#################
# 1. Visualise the similarities between ECFs and ECODs
#################

ecfs.vs.ecod.filtered.p90 <- ecfs.vs.ecod %>% filter(prob >= 90)
ecfs.vs.ecod.filtered.p93 <- ecfs.vs.ecod %>% filter(prob >= 93)
ecfs.vs.ecod.filtered.p95 <- ecfs.vs.ecod %>% filter(prob >= 95)
ecfs.vs.ecod.filtered.p98 <- ecfs.vs.ecod %>% filter(prob >= 98)
ecfs.vs.ecod.filtered.p99 <- ecfs.vs.ecod %>% filter(prob >= 99)
ggplot(ecfs.vs.ecod, aes(x = qcov, y = scov)) + geom_bin2d(bins = 300)
ggplot(ecfs.vs.ecod.filtered.p90, aes(x = qcov, y = scov)) + geom_bin2d(bins = 300)
ggplot(ecfs.vs.ecod.filtered.p93, aes(x = qcov, y = scov)) + geom_bin2d(bins = 300)

ggplot(ecfs.vs.ecod.filtered.p98, aes(x = eval, y = prob)) + scale_x_log10() + geom_bin2d(bins = 300)

ecfs.vs.ecod.fq <- ecfs.vs.ecod %>% filter(qcov >= 0.99)
ggplot(ecfs.vs.ecod.fq, aes(x = eval, y = prob)) + scale_x_log10() + geom_bin2d(bins = 300)


#################
# 2. List all x levels linked via ECF confirmed as subdomains
#################

x.levels.per.subdomain <- lapply(1:nrow(ecf.subdomains), function(subdomain.index){
  this.ecf <- ecf.subdomains$ecf[subdomain.index]
  this.ecf.hits <- ecfs.full.vs.ecod[.(this.ecf)]
  this.ecf.hits$x.level <- ecod.domains[.(this.ecf.hits$sname)]$x.level
  this.ecf.hits$x.name <- ecod.domains[.(this.ecf.hits$sname)]$x_name
  unique(this.ecf.hits$x.level)  
})
x.levels.per.subdomain.l <- sapply(x.levels.per.subdomain, length)

ecf.subdomains$multiple.x <- x.levels.per.subdomain.l>1
ecf.subdomains.subset <- ecf.subdomains[ecf.subdomains$multiple.x]

linked.architectures.list <- lapply(1:nrow(ecf.subdomains.subset), function(this.ecf.index){
  this.ecf.name <- ecf.subdomains.subset$ecf[this.ecf.index]
  this.ecf.hits <- ecfs.full.vs.ecod[.(this.ecf.name)]
  this.ecf.hits$x.level <- ecod.domains[.(this.ecf.hits$sname)]$x.level
  unique.x.levels <- unique(this.ecf.hits$x.level)
  all.x.pairs <- as.data.table(expand.grid(x1 = unique.x.levels, x2 = unique.x.levels))
  all.x.pairs <- all.x.pairs[all.x.pairs$x1 < all.x.pairs$x2]
  all.x.pairs <- cbind(all.x.pairs, ecf = this.ecf.name)
  all.x.pairs  
})
linked.architectures <- rbindlist(linked.architectures.list)

ecfs.per.architectures <- linked.architectures[, .N, by = .(x1, x2)]

graph <- graph_from_data_frame(ecfs.per.architectures, directed = FALSE)
v.names.x.level <- as.numeric(V(graph)$name)
n.nodes <- length(v.names.x.level)
v.names.x.name <- sapply(1:n.nodes, function(index) ecod.domains$x_name[which(ecod.domains$x.level == v.names.x.level[index])][1])

pdf("plot1_subdomain-linked-x.pdf", width = 12, height = 10)
plot(graph, 
     layout = layout_with_fr(graph),  # Using the Fruchterman-Reingold layout
     vertex.size = 10,
     vertex.label.cex = 0.8,  # Adjust label size
     vertex.label = v.names.x.name,  # Use node names as labels
     edge.width = ecfs.per.architectures$N,  # Use 'N' as the width of the edges
     edge.label = ecfs.per.architectures$N,  # Optionally, add edge labels showing 'N'
     edge.arrow.size = 0.5)
dev.off()

#################
# 2. List all x levels linked via any ECFs
#################

x.levels.per.ecf <- lapply(1:nrow(ecf.classification.table), function(ecf.index){
  this.ecf <- ecf.classification.table$ecf[ecf.index]
  this.ecf.hits <- ecfs.full.vs.ecod[.(this.ecf)]
  this.ecf.hits$x.level <- ecod.domains[.(this.ecf.hits$sname)]$x.level
  unique(this.ecf.hits$x.level)  
})
x.levels.per.ecf.l <- sapply(x.levels.per.ecf, length)

ecf.classification.table$multiple.x <- x.levels.per.ecf.l>1
ecf.classification.table.subset <- ecf.classification.table[ecf.classification.table$multiple.x]


linked.architectures.all.list <- lapply(1:nrow(ecf.classification.table.subset), function(this.ecf.index){
  this.ecf.name <- ecf.classification.table.subset$ecf[this.ecf.index]
  this.ecf.hits <- ecfs.full.vs.ecod[.(this.ecf.name)]
  this.ecf.hits$x.level <- ecod.domains[.(this.ecf.hits$sname)]$x.level
  unique.x.levels <- unique(this.ecf.hits$x.level)
  all.x.pairs <- as.data.table(expand.grid(x1 = unique.x.levels, x2 = unique.x.levels))
  all.x.pairs <- all.x.pairs[all.x.pairs$x1 < all.x.pairs$x2]
  all.x.pairs <- cbind(all.x.pairs, ecf = this.ecf.name)
  all.x.pairs  
})
linked.architectures.all <- rbindlist(linked.architectures.all.list)

ecfs.per.architectures.all <- linked.architectures.all[, .N, by = .(x1, x2)]

graph <- graph_from_data_frame(ecfs.per.architectures.all, directed = FALSE)
v.names.x.level <- as.numeric(V(graph)$name)
n.nodes <- length(v.names.x.level)
v.names.x.name <- sapply(1:n.nodes, function(index) ecod.domains$x_name[which(ecod.domains$x.level == v.names.x.level[index])][1])

pdf("plot2_ecf-linked-x.pdf", width = 12, height = 10)
plot(graph, 
     layout = layout_with_fr(graph),  # Using the Fruchterman-Reingold layout
     vertex.size = 10,
     vertex.label.cex = 0.8,  # Adjust label size
     vertex.label = v.names.x.name,  # Use node names as labels
     edge.width = ecfs.per.architectures.all$N,  # Use 'N' as the width of the edges
     edge.label = ecfs.per.architectures.all$N,  # Optionally, add edge labels showing 'N'
     edge.arrow.size = 0.5)
dev.off()

# TODO: make a network of all X levels linked via subdomain ECFs

# #############
# # 2. Identify all ECOD architectures in reprseqs where subdomain ECFs are found 
# #############
# 
# # A) First we check if putative subdomain ECFs, identidied via direct ECF-ECOD comparisons, are also exclusive subdomains
# # when compared from the point of view of reprseqs
# 
# # Load reprseq vs ECFs
# reprseq.ecfs.all <- fread(reprseq.ecfs.output)
# reprseq.ecfs <- reprseq.ecfs.all %>%
#   filter(scov >= MIN.SCOV.FOR.DOMAIN) %>%
#   filter(qcov >= MIN.QCOV.FOR.DOMAIN)
# 
# # Load reprseq vs ECOD
# reprseq.ecod.all <- fread(reprseq.domains.ecod.output)
# reprseq.ecod <- reprseq.ecod.all %>%
#   filter(scov >= MIN.SCOV.FOR.DOMAIN) %>%
#   filter(qcov >= MIN.QCOV.FOR.DOMAIN)
# 
# # Calculate overlap between all ECFs and ECODs in all reprseq
# 
# reprseq.ecfs.overlap <- reprseq.ecfs %>% select(qname, qstart.ecf = qstart , qend.ecf = qend, qlength.ecf = qlength,
#                                                 sname.ecf = sname, sstart.ecf = sstart, send.ecf = send, slength.ecf = slength,
#                                                 pident.ecf = pident, hit.length.ecf = hit.length)
# 
# reprseq.ecods.overlap <- reprseq.ecod %>% select(qname, qstart.ecod = qstart , qend.ecod = qend, qlength.ecod = qlength,
#                                                 sname.ecod = sname, sstart.ecod = sstart, send.ecod = send, slength.ecod = slength,
#                                                 pident.ecod = pident, hit.length.ecod = hit.length, bitscore.ecod = bitscore)
# 
# combined.ecfs.ecods <- reprseq.ecfs.overlap[reprseq.ecods.overlap, on = .(qname), allow.cartesian = TRUE]
# 
# combined.ecfs.ecods[, c("overlap.start", "overlap.end", "overlap.length") := {
#   res <- calculate_overlap_parms(qstart.ecf, qend.ecf, qstart.ecod, qend.ecod)
#   list(res$overlap.start, res$overlap.end, res$overlap.length)
# }]
# combined.ecfs.ecods <- na.omit(combined.ecfs.ecods, cols = "overlap.length")[overlap.length > 0]
# 
# reprseq.ecf.ecod.overlaps <- combined.ecfs.ecods %>% filter(overlap.length >= MIN.OVERLAP.LENGTH)
# 
# # Filter overlaps only for subdomain ECFs
# subdomain.ecods.overlaps <- reprseq.ecf.ecod.overlaps[reprseq.ecf.ecod.overlaps$sname.ecf %in% ecf.subdomains$ecf]
# setkey(subdomain.ecods.overlaps, sname.ecf)
# 
# # 'ecf.cov.overlap' is the proportion of the ECF hit that is covered by the overlap 
# # (if high, then most of the ECF is covered by the ECF-ECOD)
# subdomain.ecods.overlaps$ecf.cov.overlap <- round(subdomain.ecods.overlaps$overlap.length/subdomain.ecods.overlaps$hit.length.ecf, d=4)
# 
# # 'ecod.cov.overlap' is the proportion of the ECOD hit (not entire ECOD) covered by the overlap
# # (again, if ~1, then the ECOD hit is equivalent to the overlap)
# subdomain.ecods.overlaps$ecod.cov.overlap <- round(subdomain.ecods.overlaps$overlap.length/subdomain.ecods.overlaps$hit.length.ecod, d=4)
# 
# # 'scov.ecod' is the proportion of the ECOD length that hits the query
# subdomain.ecods.overlaps[, scov.ecod := round((send.ecod - sstart.ecod + 1)/slength.ecod, d = 4)]
# 
# # 'ecod.cov.by.ecf' is the proportion of the entire ECOD covered by the overlap,
# # and it is calculated by weighing the proportion of full ECOD's length that hits the query
# # by the proportion of that hit length that overlaps with the ECF.
# subdomain.ecods.overlaps[, ecod.cov.by.ecf := round(scov.ecod * ecod.cov.overlap, d = 4)]
# 
# subdomain.ecods.overlaps$is.subdomain <- NA
# subdomain.ecods.overlaps$is.subdomain[subdomain.ecods.overlaps$ecf.cov.overlap > 0.9 & subdomain.ecods.overlaps$ecod.cov.by.ecf >= 0.6] <- FALSE
# subdomain.ecods.overlaps$is.subdomain[subdomain.ecods.overlaps$ecf.cov.overlap > 0.9 & subdomain.ecods.overlaps$ecod.cov.by.ecf < 0.6] <- TRUE
# 
# setkey(subdomain.ecods.overlaps, sname.ecf)
# 
# subdomain.reprseq.confirmed <- sapply(1:nrow(ecf.subdomains), function(index){
#   this.subdomain <- ecf.subdomains$ecf[index]
#   this.subdomain.overlaps <- subdomain.ecods.overlaps[.(this.subdomain)]
#   this.subdomain.overlaps <- this.subdomain.overlaps[!is.na(this.subdomain.overlaps$is.subdomain)]
#   confirmed.subdomain <- all(this.subdomain.overlaps$is.subdomain)
#   confirmed.subdomain  
# })
# ecf.subdomains$reprseq.confirmed <- subdomain.reprseq.confirmed
# 
# ecf.subdomains.confirmed <- ecf.subdomains[ecf.subdomains$reprseq.confirmed]
# 
# # B) Then we check if 
# 
# x.levels.per.subdomain.reprseq <- lapply(1:nrow(ecf.subdomains.confirmed), function(confirmed.subdomain.index){
#   this.conf.subdomain <- ecf.subdomains.confirmed$ecf[confirmed.subdomain.index]
#   this.conf.subdomain.ecods.overlaps <- subdomain.ecods.overlaps[.(this.conf.subdomain)]
#   this.conf.subdomain.ecods.overlaps$x.level <- ecod.domains[.(this.conf.subdomain.ecods.overlaps$sname.ecod)]$x.level
#   setorder(this.conf.subdomain.ecods.overlaps, qname, sname.ecf, x.level, -bitscore.ecod)
#   this.conf.subdomain.ecods.overlaps.filtered <- this.conf.subdomain.ecods.overlaps[, .SD[1], by = .(qname, sname.ecf, x.level)]
#   unique(this.conf.subdomain.ecods.overlaps.filtered$x.level)  
# })
# x.levels.per.subdomain.reprseq.l <- sapply(x.levels.per.subdomain.reprseq, length)
# x.levels.per.subdomain.reprseq[x.levels.per.subdomain.reprseq.l>1]
# 
# #############
# # 3. Sequence similarity of shared (subdomain) fragments
# #############
# 
# # Load sequence data
# reprseq.seqs <- read.fasta(reprseq.seqs.filename)
# reprseq.seqs.s <- sapply(reprseq.seqs, c2s)
# reprseq.seqs.n <- names(reprseq.seqs)
# 
# # Find all cases of shared subdomain ECFs between different Xs
# 
# confirmed.subdomain.index=5
# this.conf.subdomain <- ecf.subdomains.confirmed$ecf[confirmed.subdomain.index]
# this.conf.subdomain.ecods.overlaps <- subdomain.ecods.overlaps[.(this.conf.subdomain)]
# this.conf.subdomain.ecods.overlaps$x.level <- ecod.domains[.(this.conf.subdomain.ecods.overlaps$sname.ecod)]$x.level
# setorder(this.conf.subdomain.ecods.overlaps, qname, sname.ecf, x.level, -bitscore.ecod)
# this.conf.subdomain.ecods.overlaps.filtered <- this.conf.subdomain.ecods.overlaps[, .SD[1], by = .(qname, sname.ecf, x.level)]
# 
# View(this.conf.subdomain.ecods.overlaps.filtered)
# 
# this.conf.subdomain.ecods.overlaps.filtered[c(7, 151), ]
# 
# 
# reprseq.i.rn <- 1
# reprseq.j.rn <- 41
# reprseq.i.table <- this.conf.subdomain.ecods.overlaps.filtered[reprseq.i.rn,] 
# reprseq.j.table <- this.conf.subdomain.ecods.overlaps.filtered[reprseq.j.rn,] 
# 
# reprseq.name.i <- reprseq.i.table$qname
# reprseq.name.i.ecf.start <- reprseq.i.table$qstart.ecf
# reprseq.name.i.ecf.end <- reprseq.i.table$qend.ecf
# reprseq.seq.i <- reprseq.seqs.s[reprseq.seqs.n == reprseq.name.i]
# reprseq.seq.i.frag <- substr(reprseq.seq.i, reprseq.name.i.ecf.start, reprseq.name.i.ecf.end)
# 
# reprseq.name.j <- reprseq.j.table$qname
# reprseq.name.j.ecf.start <- reprseq.j.table$qstart.ecf
# reprseq.name.j.ecf.end <- reprseq.j.table$qend.ecf
# reprseq.seq.j <- reprseq.seqs.s[reprseq.seqs.n == reprseq.name.j]
# reprseq.seq.j.frag <- substr(reprseq.seq.j, reprseq.name.j.ecf.start, reprseq.name.j.ecf.end)
# 
# seq1 <- AAString(reprseq.seq.i.frag)
# seq2 <- AAString(reprseq.seq.j.frag)
# alignment <- pairwiseAlignment(pattern = seq1, subject = seq2, substitutionMatrix = "BLOSUM62", gapOpening = 10, gapExtension = 1, type = "global")
# aligned_pattern <- as.character(alignment@pattern)
# aligned_subject <- as.character(alignment@subject)
# # Function to calculate percentage identity
# calculate_percentage_identity <- function(aligned_seq1, aligned_seq2) {
#   matches <- sum(aligned_seq1 == aligned_seq2)  # Count matches
#   alignment_length <- length(aligned_seq1)     # Total positions in the alignment
#   percent_identity <- (matches / alignment_length) * 100  # Calculate percentage
#   return(percent_identity)
# }
# percentage_identity <- calculate_percentage_identity(s2c(aligned_pattern), s2c(aligned_subject))
# 
# 
# 
# 
# 
# 
