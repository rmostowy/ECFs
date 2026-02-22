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

# Load ECOD domain table
ecod.domains <- fread(ecod.domains.output)
setkey(ecod.domains, full.id)

# Load all-by-all reprseq
reprseq.vs.reprseq <- fread(reprseq.vs.reprseq.filename)
reprseq.vs.reprseq[, qcov := round((qend - qstart + 1)/qlength, d = 4)]
reprseq.vs.reprseq[, scov := round((send - sstart + 1)/slength, d = 4)]
reprseq.vs.reprseq <- reprseq.vs.reprseq[qname %chin% all.reprseq & sname %chin% all.reprseq]

# Load functional annotations
reprseq.functions.unique <- fread(reprseq.functions.unique.output)
setkey(reprseq.functions.unique, reprseq)

# Load taxonomy data
reprseq2taxonomy <- unique(reprseq.table[, c("reprseq", "cluster.name")])
reprseq2taxonomy[, accession := sub("\\..*", "", cluster.name)]
reprseq2taxonomy[, cluster.name := NULL]
metadata.table <- fread(metadata.table.output)
setkey(metadata.table, accession)

reprseq2taxonomy$family.ncbi <- metadata.table[.(reprseq2taxonomy$accession)]$family.ncbi
reprseq2taxonomy$genus.ncbi <- metadata.table[.(reprseq2taxonomy$accession)]$genus.ncbi
setkey(reprseq2taxonomy, reprseq)


# ############
# 1. Find all overlapping ECFs and ECODs in reprseqs and calculate overlaps
# ############

# Load reprseq vs ECFs
reprseq.ecfs.all <- fread(reprseq.ecfs.output)
reprseq.ecfs <- reprseq.ecfs.all %>%
  filter(prob >= MIN.PROB.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(scov >= MIN.SCOV.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE)

# Load reprseq vs ECOD
reprseq.ecod.all <- fread(reprseq.domains.ecod.output)
reprseq.ecod <- reprseq.ecod.all %>%
  filter(prob >= MIN.PROB.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(scov >= MIN.SCOV.FOR.DOMAIN.EQUIVALENCE)

# Calculate overlap between all ECFs and ECODs in all reprseq
reprseq.ecfs.overlap <- reprseq.ecfs %>% select(qname, qlength, qstart.ecf = qstart , qend.ecf = qend,
                                                sname.ecf = sname, sstart.ecf = sstart, send.ecf = send, slength.ecf = slength,
                                                pident.ecf = pident, hit.length.ecf = hit.length)
reprseq.ecods.overlap <- reprseq.ecod %>% select(qname, qstart.ecod = qstart , qend.ecod = qend,
                                                 sname.ecod = sname, sstart.ecod = sstart, send.ecod = send, slength.ecod = slength,
                                                 pident.ecod = pident, hit.length.ecod = hit.length, bitscore.ecod = bitscore)
combined.ecfs.ecods <- reprseq.ecfs.overlap[reprseq.ecods.overlap, on = .(qname), allow.cartesian = TRUE]
combined.ecfs.ecods[, c("overlap.start", "overlap.end", "overlap.length") := {
  res <- calculate_overlap_parms(qstart.ecf, qend.ecf, qstart.ecod, qend.ecod)
  list(res$overlap.start, res$overlap.end, res$overlap.length)
}]
combined.ecfs.ecods <- na.omit(combined.ecfs.ecods, cols = "overlap.length")[overlap.length > 0]

# Filter by minimum length
reprseq.ecf.ecod.overlaps <- combined.ecfs.ecods %>% filter(overlap.length >= MIN.OVERLAP.LENGTH)
setkey(reprseq.ecf.ecod.overlaps, sname.ecf)

# Calculate overlaps between ECFs and ECODs:

# a) 'ecf.cov.overlap' is the proportion of the ECF hit (from the query's point of view) 
# that is covered by the overlap (if high, then most of the ECF is covered by the ECF-ECOD)
reprseq.ecf.ecod.overlaps[, ecf.cov.overlap := round(overlap.length/hit.length.ecf, d = 4)]

# b) 'ecod.cov.overlap' is the proportion of the ECOD hit (not entire ECOD) covered by the overlap
# (again, if ~1, then the ECOD hit is equivalent to the overlap)
reprseq.ecf.ecod.overlaps[, ecod.cov.overlap := round(overlap.length/hit.length.ecod, d = 4)]

# c) 'scov.ecod' is the proportion of the ECOD length that hits the query
reprseq.ecf.ecod.overlaps[, scov.ecod := round((send.ecod - sstart.ecod + 1)/slength.ecod, d = 4)]

# d) 'ecod.cov.by.ecf' is the proportion of the entire ECOD covered by the overlap,
# and it is calculated by weighing the proportion of full ECOD's length that hits the query
# by the proportion of that hit length that overlaps with the ECF.
reprseq.ecf.ecod.overlaps[, ecod.cov.by.ecf := round(scov.ecod * ecod.cov.overlap, d = 4)]

# Consider only ECFs which are nearly covered by ECODs
reprseq.ecf.ecod.overlaps.shortlist <- reprseq.ecf.ecod.overlaps[reprseq.ecf.ecod.overlaps$ecf.cov.overlap > MIN.QCOV.FOR.DOMAIN.EQUIVALENCE]

# ############
# 2. Identify all pairs of ECFs which connect ECOF
# ############

# Identify all ECFs which overlap with at least 2 different H level ECODs
ecf.ecod.overlapping <- reprseq.ecf.ecod.overlaps.shortlist %>% select(ecf = sname.ecf, ecod = sname.ecod)
ecf.ecod.overlapping$x_id <- ecod.domains[.(ecf.ecod.overlapping$ecod)]$x_id
ecf.ecod.overlapping$h_id <- ecod.domains[.(ecf.ecod.overlapping$ecod)]$h_id
ecf.ecod.overlapping$t_id <- ecod.domains[.(ecf.ecod.overlapping$ecod)]$t_id
no.t.per.ecf <- ecf.ecod.overlapping[, .(no.t.levels = uniqueN(t_id)), by = .(ecf)]
ecfs.with.multiple.t.levels <- sort(no.t.per.ecf[no.t.levels>1]$ecf)

# Of all reprseq-ecf pairs, consider only a single ECOD hit amongst those with the same H level (with the highest bitscore)
reprseq.ecf.ecod.overlaps.filtered1 <- reprseq.ecf.ecod.overlaps.shortlist[sname.ecf %in% ecfs.with.multiple.t.levels]
reprseq.ecf.ecod.overlaps.filtered1$x.id <- ecod.domains[.(reprseq.ecf.ecod.overlaps.filtered1$sname.ecod)]$x_id
reprseq.ecf.ecod.overlaps.filtered1$h.id <- ecod.domains[.(reprseq.ecf.ecod.overlaps.filtered1$sname.ecod)]$h_id
reprseq.ecf.ecod.overlaps.filtered1$t.id <- ecod.domains[.(reprseq.ecf.ecod.overlaps.filtered1$sname.ecod)]$t_id
setorder(reprseq.ecf.ecod.overlaps.filtered1, qname, sname.ecf, t.id, -bitscore.ecod)
reprseq.ecf.ecod.overlaps.filtered2 <- reprseq.ecf.ecod.overlaps.filtered1[, .SD[1], by = .(qname, sname.ecf, t.id)]

# Inner-join the table to list of pairs of reprseqs
left.join.table <- reprseq.ecf.ecod.overlaps.filtered2 %>% select(query.a = qname, qlength.a = qlength, shared.ecf = sname.ecf, ecf.length = slength.ecf,
                                         ecf.start.a = qstart.ecf, ecf.end.a = qend.ecf, 
                                         ecod.a = sname.ecod, x.id.a = x.id, h.id.a = h.id, t.id.a = t.id, ecod.start.a = qstart.ecod, ecod.end.a = qend.ecod, ecod.cov.by.ecf.a = ecod.cov.by.ecf)

right.join.table <- reprseq.ecf.ecod.overlaps.filtered2 %>% select(query.b = qname, qlength.b = qlength, shared.ecf = sname.ecf,
                                          ecf.start.b = qstart.ecf, ecf.end.b = qend.ecf, 
                                          ecod.b = sname.ecod, x.id.b = x.id, h.id.b = h.id, t.id.b = t.id, ecod.start.b = qstart.ecod, ecod.end.b = qend.ecod, ecod.cov.by.ecf.b = ecod.cov.by.ecf)

putative.recomb.pairs.temp <- left.join.table[right.join.table, on = .(shared.ecf), allow.cartesian = TRUE]
putative.recomb.pairs.temp <- putative.recomb.pairs.temp[query.a < query.b]
# putative.recomb.pairs.temp <- putative.recomb.pairs.temp[query.a < query.b & t.id.a != t.id.b]

putative.recomb.pairs <- putative.recomb.pairs.temp %>% 
  select(query.a, query.b, shared.ecf, qlength.a, qlength.b,
         ecf.start.a, ecf.end.a, ecf.start.b, ecf.end.b,
         ecod.a, ecod.b, x.id.a, h.id.a, x.id.b, h.id.b, t.id.a, t.id.b,
         ecod.start.a, ecod.end.a, ecod.start.b, ecod.end.b,
         ecod.cov.by.ecf.a, ecod.cov.by.ecf.b)
# putative.recomb.pairs[ , different.x := x.id.a != x.id.b]


# ############
# 3. Merge-in the information about reprseq alignment
# ############

# Filter 'reprseq.vs.reprseq' to leave only the strongest hit between each pair
# setorder(reprseq.vs.reprseq, qname, sname, -bitscore)
setorder(reprseq.vs.reprseq, qname, sname, -pident)
reprseq.vs.reprseq.filtered <- reprseq.vs.reprseq[, .SD[1], by = .(qname, sname)] %>%
  select(qname, sname, aln.start.a = qstart, aln.end.a = qend, aln.start.b = sstart, aln.end.b = send,
         cov.a = qcov, cov.b = scov, aln.prob = prob, aln.eval = eval, aln.pident = pident)

# Merge full reprseq alignments to add to the 'putative.recomb.pairs' table
putative.recomb.pairs.aln <- merge(putative.recomb.pairs, reprseq.vs.reprseq.filtered,
                                      by.x = c("query.a", "query.b"),
                                      by.y = c("qname", "sname"),
                                      all.x = TRUE, all.y = FALSE)
putative.recomb.pairs.aln <- putative.recomb.pairs.aln[!is.na(aln.pident)]
putative.recomb.pairs.aln[, max.cov := pmax(cov.a, cov.b, na.rm = T)]

ggplot(putative.recomb.pairs.aln[aln.pident>20], aes(x = max.cov, y = aln.pident)) + geom_bin2d()

hist(putative.recomb.pairs.aln$aln.pident)
sh <- putative.recomb.pairs.aln[max.cov < 0.5 & aln.pident > 70]



setorder(reprseq.vs.reprseq, qname, sname, -pident)
rr <- reprseq.vs.reprseq[, .SD[1], by = .(qname, sname)]
rrr <- rr[prob>95 & pident > 95 & max.cov < 0.3]


# # ggplot(putative.recomb.pairs.aln[aln.pident>20], aes(x = cov.a, y = ecod.cov.by.ecf.a)) + geom_bin2d()
# 
# # ############
# # 4. Align ECF-containing fragments 
# # ############
# 
# # Load sequence data
# reprseq.seqs <- read.fasta(reprseq.seqs.filename)
# reprseq.seqs.s <- sapply(reprseq.seqs, c2s)
# reprseq.seqs.n <- names(reprseq.seqs)
# reprseq.seqs.table <- data.table(name = reprseq.seqs.n, string = toupper(reprseq.seqs.s))
# setkey(reprseq.seqs.table, name)
# 
# # Define an object to work on
# # putative.recomb.pairs.ecf.coords <- putative.recomb.pairs %>% select(query.a, query.b, qlength.a, ecf.start.a, ecf.end.a, qlength.b, ecf.start.b, ecf.end.b)
# putative.recomb.pairs.ecf.coords <- putative.recomb.pairs.aln
# 
# # Add full sequence information
# putative.recomb.pairs.ecf.coords$seq.a <- reprseq.seqs.table[.(putative.recomb.pairs.ecf.coords$query.a)]$string
# putative.recomb.pairs.ecf.coords$seq.b <- reprseq.seqs.table[.(putative.recomb.pairs.ecf.coords$query.b)]$string
# 
# # Remove reprseqs for which there's no sequence
# # These are repseqs that should have been removed at the inital filtering due to Xs in their sequence
# # TODO!!! WE NEED TO CORRECT THE DATA PROCESSING TO AVOID THIS
# 
# # seqs.are.ok <- !is.na(putative.recomb.pairs.ecf.coords$seq.a) & !is.na(putative.recomb.pairs.ecf.coords$seq.b)
# # putative.recomb.pairs.aln <- putative.recomb.pairs.aln[seqs.are.ok]
# # putative.recomb.pairs.ecf.coords <- putative.recomb.pairs.ecf.coords[seqs.are.ok]
# 
# # Extract subsequences for query.a using the correct string column
# putative.recomb.pairs.ecf.coords[, ecf.seq.a := substr(seq.a, ecf.start.a, ecf.end.a)]
# # Extract subsequences for query.b using the correct string column
# putative.recomb.pairs.ecf.coords[, ecf.seq.b := substr(seq.b, ecf.start.b, ecf.end.b)]
# # Drop full sequences
# putative.recomb.pairs.ecf.coords[, c("seq.a", "seq.b") := NULL]  # Removing the full sequence columns to save memory
# 
# calculate_percentage_identity <- function(aligned_seq1, aligned_seq2) {
#   matches <- sum(aligned_seq1 == aligned_seq2)  # Count matches
#   alignment_length <- length(aligned_seq1)     # Total positions in the alignment
#   percent_identity <- (matches / alignment_length) * 100  # Calculate percentage
#   return(percent_identity)
# }
# get.pident <- function(alignment){
#   percentage_identity <- NA
#   alignment.is.null <- is.null(alignment)
#   if(!alignment.is.null){
#     aligned_pattern <- as.character(alignment@pattern)
#     aligned_subject <- as.character(alignment@subject)
#     percentage_identity <- calculate_percentage_identity(s2c(aligned_pattern), s2c(aligned_subject))
#   }
#   return(percentage_identity)
# }
# 
# # Align all ECF fragments
# 
# # #|
# # set.seed(1234)
# # sample.size <- 0.1
# # nrow.pairs <- nrow(putative.recomb.pairs.ecf.coords)
# # nrows.sample <- ceiling(nrow.pairs*sample.size)
# # rows2sample <- sort(sample(nrow.pairs, nrows.sample))
# # putative.recomb.pairs.ecf.coords <- putative.recomb.pairs.ecf.coords[rows2sample]
# # #|
# 
# # N <- nrow(putative.recomb.pairs.ecf.coords)
# # 
# # 
# # NO.CORES=8
# # all.alignments.list <- mclapply(1:N, function(pair.index){
# #   seq1.raw <- putative.recomb.pairs.ecf.coords$ecf.seq.a[pair.index]
# #   seq2.raw <- putative.recomb.pairs.ecf.coords$ecf.seq.b[pair.index]
# #   seq1 <- AAString(toupper(seq1.raw))
# #   seq2 <- AAString(toupper(seq2.raw))
# #   alignment <- pairwiseAlignment(pattern = seq1, subject = seq2, substitutionMatrix = "BLOSUM62", type = "global")
# #   # aligned_pattern <- as.character(alignment@pattern)
# #   # aligned_subject <- as.character(alignment@subject)
# #   # # Function to calculate percentage identity
# #   # percentage_identity <- calculate_percentage_identity(s2c(aligned_pattern), s2c(aligned_subject))
# #   return(alignment)
# # }, mc.cores = NO.CORES)
# # 
# # all.alignments.pident <- sapply(all.alignments.list, get.pident)
# # putative.recomb.pairs.ecf.coords$ecf.pident <- all.alignments.pident
# # 
# # 
# 
# ############
# 
# putative.recomb.pairs.ecf.coords.subset <- putative.recomb.pairs.ecf.coords
# NO.CORES=10
# all.alignments.list <- mclapply(1:nrow(putative.recomb.pairs.ecf.coords.subset), function(pair.index){
# # all.alignments.list <- lapply(1:nrow(putative.recomb.pairs.ecf.coords.subset), function(pair.index){
#   cat(pair.index,"\n")
#   seq1.raw <- putative.recomb.pairs.ecf.coords.subset$ecf.seq.a[pair.index]
#   seq2.raw <- putative.recomb.pairs.ecf.coords.subset$ecf.seq.b[pair.index]
#   alignment <- NULL
#   check <- !is.na(seq1.raw) & !is.na(seq2.raw)
#   if(check){
#     seq1 <- AAString(toupper(seq1.raw))
#     seq2 <- AAString(toupper(seq2.raw))
#     alignment <- pairwiseAlignment(pattern = seq1, subject = seq2, substitutionMatrix = "BLOSUM62", type = "global")
#     
#   }
#   return(alignment)
# }, mc.cores = NO.CORES)
# # })
# all.alignments.pident <- sapply(all.alignments.list, get.pident)
# putative.recomb.pairs.ecf.coords$ecf.pident <- all.alignments.pident
# 
# fwrite(putative.recomb.pairs.ecf.coords, file = ecf.sharing.pairs.output)
# 
# 
# putative.recomb.pairs.ecf.coords$genus.a <- name.table.unique[.(putative.recomb.pairs.ecf.coords$query.a)]$genus.ncbi
# putative.recomb.pairs.ecf.coords$genus.b <- name.table.unique[.(putative.recomb.pairs.ecf.coords$query.b)]$genus.ncbi
# putative.recomb.pairs.ecf.coords$funct.a <- reprseq.functions.unique[.(putative.recomb.pairs.ecf.coords$query.a)]$funct
# putative.recomb.pairs.ecf.coords$class.a <- reprseq.functions.unique[.(putative.recomb.pairs.ecf.coords$query.a)]$class
# putative.recomb.pairs.ecf.coords$funct.b <- reprseq.functions.unique[.(putative.recomb.pairs.ecf.coords$query.b)]$funct
# putative.recomb.pairs.ecf.coords$class.b <- reprseq.functions.unique[.(putative.recomb.pairs.ecf.coords$query.b)]$class
# putative.recomb.pairs.ecf.coords$t.name.a <- ecod.domains[.(putative.recomb.pairs.ecf.coords$ecod.a)]$t_name
# putative.recomb.pairs.ecf.coords$t.id.a <- ecod.domains[.(putative.recomb.pairs.ecf.coords$ecod.a)]$t_id
# putative.recomb.pairs.ecf.coords$t.name.b <- ecod.domains[.(putative.recomb.pairs.ecf.coords$ecod.b)]$t_name
# putative.recomb.pairs.ecf.coords$t.id.b <- ecod.domains[.(putative.recomb.pairs.ecf.coords$ecod.b)]$t_id
# 
# 
# View(putative.recomb.pairs.ecf.coords[ecf.pident > 50 & cov.a < 0.5 & cov.b < 0.5])
# 
# source("examples.R")
# 
# subdomain.shortlist <- putative.recomb.pairs.ecf.coords[ecf.pident > 50 & cov.a < 0.5 & cov.b < 0.5]
# 
# row.index <- 143
# # row.index <- 76
# subdomain.shortlist[row.index]
# 
# protein1_name <- subdomain.shortlist$query.a[row.index]
# protein2_name <- subdomain.shortlist$query.b[row.index]
# 
# protein1_length <- subdomain.shortlist$qlength.a[row.index]
# protein2_length <- subdomain.shortlist$qlength.b[row.index]
# aln1_start <- subdomain.shortlist$aln.start.a[row.index]
# aln1_end <- subdomain.shortlist$aln.end.a[row.index]
# aln2_start <- subdomain.shortlist$aln.start.b[row.index]
# aln2_end <- subdomain.shortlist$aln.end.b[row.index]
# 
# ecf1_start <- subdomain.shortlist$ecf.start.a[row.index]
# ecf1_end <- subdomain.shortlist$ecf.end.a[row.index]
# ecf2_start <- subdomain.shortlist$ecf.start.b[row.index]
# ecf2_end <- subdomain.shortlist$ecf.end.b[row.index]
# 
# aln_identity <- subdomain.shortlist$aln.pident[row.index]/100
# ecf_identity <- subdomain.shortlist$ecf.pident[row.index]/100
# 
# ecod1_start <- subdomain.shortlist$ecod.start.a[row.index]
# ecod1_end <- subdomain.shortlist$ecod.end.a[row.index]
# ecod2_start <- subdomain.shortlist$ecod.start.b[row.index]
# ecod2_end <- subdomain.shortlist$ecod.end.b[row.index]
# 
# protein1_funct <- subdomain.shortlist$funct.a[row.index]
# protein1_class <- subdomain.shortlist$class.a[row.index]
# protein1_tname <- subdomain.shortlist$t.name.a[row.index]
# protein1_tid <- subdomain.shortlist$t.id.a[row.index]
# protein1_description <- sprintf("%s | %s\n%s (%s)", protein1_funct, protein1_class, protein1_tname, protein1_tid)
# 
# protein2_funct <- subdomain.shortlist$funct.b[row.index]
# protein2_class <- subdomain.shortlist$class.b[row.index]
# protein2_tname <- subdomain.shortlist$t.name.b[row.index]
# protein2_tid <- subdomain.shortlist$t.id.b[row.index]
# protein2_description <- sprintf("%s | %s\n%s (%s)", protein2_funct, protein2_class, protein2_tname, protein2_tid)
# 
# # source("examples.R")
# plot_protein_alignment(protein1_name, protein2_name,
#                        protein1_length, protein2_length, 
#                        aln1_start, aln1_end, aln2_start, aln2_end, 
#                        ecf1_start, ecf1_end, ecf2_start, ecf2_end, 
#                        aln_identity, ecf_identity, 
#                        ecod1_start, ecod1_end, ecod2_start, ecod2_end,
#                        protein1_description, protein2_description)
# 
# 
# ############################
# 
# 
# # putative.recomb.pairs.aln$ecf.seq.a <- putative.recomb.pairs.ecf.coords$ecf.seq.a
# # putative.recomb.pairs.aln$ecf.seq.b <- putative.recomb.pairs.ecf.coords$ecf.seq.b
# # putative.recomb.pairs.aln$ecf.pident <- putative.recomb.pairs.ecf.coords$ecf.pident
# # 
# # fwrite(putative.recomb.pairs.aln, file = ecf.sharing.pairs.output)
# 
# # z <- putative.recomb.pairs.aln[putative.recomb.pairs.aln$ecf.pident>70 & putative.recomb.pairs.aln$different.x] %>% select(ecod.a, ecod.b, x.id.a, x.id.b)
# # z$x.name.a <- ecod.domains[.(z$ecod.a)]$x_name
# # z$x.name.b <- ecod.domains[.(z$ecod.b)]$x_name
# 
# 
# # # TEST ON A SUBSET OF 10%
# # set.seed(1234)
# # sample.size <- 0.1
# # nrow.pairs <- nrow(putative.recomb.pairs.ecf.coords)
# # nrows.sample <- ceiling(nrow.pairs*sample.size)
# # rows2sample <- sort(sample(nrow.pairs, nrows.sample))
# # 
# # putative.recomb.pairs.ecf.coords.subsample <- putative.recomb.pairs.ecf.coords[rows2sample]
# # N <- nrow(putative.recomb.pairs.ecf.coords.subsample)
# # 
# # NO.CORES=12
# # all.alignments.list <- mclapply(1:N, function(pair.index){
# #   seq1.raw <- putative.recomb.pairs.ecf.coords.subsample$ecf.seq.a[pair.index]
# #   seq2.raw <- putative.recomb.pairs.ecf.coords.subsample$ecf.seq.b[pair.index]
# #   seq1 <- AAString(toupper(seq1.raw))
# #   seq2 <- AAString(toupper(seq2.raw))
# #   alignment <- pairwiseAlignment(pattern = seq1, subject = seq2, substitutionMatrix = "BLOSUM62", type = "global")
# #   # aligned_pattern <- as.character(alignment@pattern)
# #   # aligned_subject <- as.character(alignment@subject)
# #   # # Function to calculate percentage identity
# #   # percentage_identity <- calculate_percentage_identity(s2c(aligned_pattern), s2c(aligned_subject))
# #   return(alignment)
# # }, mc.cores = NO.CORES)
# # 
# # all.alignments.pident <- sapply(all.alignments.list, get.pident)
# # 
# # putative.recomb.pairs.ecf.coords.subsample$ecf.pident <- all.alignments.pident
# 
# 
# 
# 
# ####################################################################
# ## |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# ####################################################################
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # #############
# # # 1. Identify all ECOD architectures in reprseqs where subdomain ECFs are found
# # #############
# # 
# # # A) First we check if putative subdomain ECFs, identidied via direct ECF-ECOD comparisons, are also exclusive subdomains
# # # when compared from the point of view of reprseqs
# # 
# # # Load reprseq vs ECFs
# # reprseq.ecfs.all <- fread(reprseq.ecfs.output)
# # reprseq.ecfs <- reprseq.ecfs.all %>%
# #   filter(prob >= MIN.PROB.FOR.DOMAIN.EQUIVALENCE) %>%
# #   filter(eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE) %>%
# #   filter(scov >= MIN.SCOV.FOR.DOMAIN) %>%
# #   filter(qcov >= MIN.QCOV.FOR.DOMAIN)
# # 
# # # Load reprseq vs ECOD
# # reprseq.ecod.all <- fread(reprseq.domains.ecod.output)
# # reprseq.ecod <- reprseq.ecod.all %>%
# #   filter(prob >= MIN.PROB.FOR.DOMAIN.EQUIVALENCE) %>%
# #   filter(eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE) %>%
# #   filter(scov >= MIN.SCOV.FOR.DOMAIN) %>%
# #   filter(qcov >= MIN.QCOV.FOR.DOMAIN)
# # 
# # # Calculate overlap between all ECFs and ECODs in all reprseq
# # 
# # reprseq.ecfs.overlap <- reprseq.ecfs %>% select(qname, qlength, qstart.ecf = qstart , qend.ecf = qend,
# #                                                 sname.ecf = sname, sstart.ecf = sstart, send.ecf = send, slength.ecf = slength,
# #                                                 pident.ecf = pident, hit.length.ecf = hit.length)
# # 
# # reprseq.ecods.overlap <- reprseq.ecod %>% select(qname, qstart.ecod = qstart , qend.ecod = qend,
# #                                                 sname.ecod = sname, sstart.ecod = sstart, send.ecod = send, slength.ecod = slength,
# #                                                 pident.ecod = pident, hit.length.ecod = hit.length, bitscore.ecod = bitscore)
# # 
# # combined.ecfs.ecods <- reprseq.ecfs.overlap[reprseq.ecods.overlap, on = .(qname), allow.cartesian = TRUE]
# # 
# # combined.ecfs.ecods[, c("overlap.start", "overlap.end", "overlap.length") := {
# #   res <- calculate_overlap_parms(qstart.ecf, qend.ecf, qstart.ecod, qend.ecod)
# #   list(res$overlap.start, res$overlap.end, res$overlap.length)
# # }]
# # combined.ecfs.ecods <- na.omit(combined.ecfs.ecods, cols = "overlap.length")[overlap.length > 0]
# # 
# # reprseq.ecf.ecod.overlaps <- combined.ecfs.ecods %>% filter(overlap.length >= MIN.OVERLAP.LENGTH)
# # 
# # # Filter overlaps only for subdomain ECFs
# # subdomain.ecods.overlaps <- reprseq.ecf.ecod.overlaps[reprseq.ecf.ecod.overlaps$sname.ecf %in% ecf.subdomains$ecf]
# # setkey(subdomain.ecods.overlaps, sname.ecf)
# # 
# # # 'ecf.cov.overlap' is the proportion of the ECF hit that is covered by the overlap 
# # # (if high, then most of the ECF is covered by the ECF-ECOD)
# # subdomain.ecods.overlaps$ecf.cov.overlap <- round(subdomain.ecods.overlaps$overlap.length/subdomain.ecods.overlaps$hit.length.ecf, d=4)
# # 
# # # 'ecod.cov.overlap' is the proportion of the ECOD hit (not entire ECOD) covered by the overlap
# # # (again, if ~1, then the ECOD hit is equivalent to the overlap)
# # subdomain.ecods.overlaps$ecod.cov.overlap <- round(subdomain.ecods.overlaps$overlap.length/subdomain.ecods.overlaps$hit.length.ecod, d=4)
# # 
# # # 'scov.ecod' is the proportion of the ECOD length that hits the query
# # subdomain.ecods.overlaps[, scov.ecod := round((send.ecod - sstart.ecod + 1)/slength.ecod, d = 4)]
# # 
# # # 'ecod.cov.by.ecf' is the proportion of the entire ECOD covered by the overlap,
# # # and it is calculated by weighing the proportion of full ECOD's length that hits the query
# # # by the proportion of that hit length that overlaps with the ECF.
# # subdomain.ecods.overlaps[, ecod.cov.by.ecf := round(scov.ecod * ecod.cov.overlap, d = 4)]
# # 
# # subdomain.ecods.overlaps$is.subdomain <- NA
# # subdomain.ecods.overlaps$is.subdomain[subdomain.ecods.overlaps$ecf.cov.overlap > MIN.QCOV.FOR.DOMAIN.EQUIVALENCE & subdomain.ecods.overlaps$ecod.cov.by.ecf >= MAX.SCOV.FOR.SUBDOMAIN] <- FALSE
# # subdomain.ecods.overlaps$is.subdomain[subdomain.ecods.overlaps$ecf.cov.overlap > MIN.QCOV.FOR.DOMAIN.EQUIVALENCE & subdomain.ecods.overlaps$ecod.cov.by.ecf < MAX.SCOV.FOR.SUBDOMAIN] <- TRUE
# # 
# # setkey(subdomain.ecods.overlaps, sname.ecf)
# # 
# # subdomain.reprseq.confirmed <- sapply(1:nrow(ecf.subdomains), function(index){
# #   this.subdomain <- ecf.subdomains$ecf[index]
# #   this.subdomain.overlaps <- subdomain.ecods.overlaps[.(this.subdomain)]
# #   this.subdomain.overlaps <- this.subdomain.overlaps[!is.na(this.subdomain.overlaps$is.subdomain)]
# #   confirmed.subdomain <- all(this.subdomain.overlaps$is.subdomain)
# #   confirmed.subdomain  
# # })
# # ecf.subdomains$reprseq.confirmed <- subdomain.reprseq.confirmed
# # 
# # ecf.subdomains.confirmed <- ecf.subdomains[ecf.subdomains$reprseq.confirmed]
# # 
# # # B) Then we check if 
# # 
# # x.levels.per.subdomain.reprseq <- lapply(1:nrow(ecf.subdomains.confirmed), function(confirmed.subdomain.index){
# #   this.conf.subdomain <- ecf.subdomains.confirmed$ecf[confirmed.subdomain.index]
# #   this.conf.subdomain.ecods.overlaps <- subdomain.ecods.overlaps[.(this.conf.subdomain)]
# #   this.conf.subdomain.ecods.overlaps$x.level <- ecod.domains[.(this.conf.subdomain.ecods.overlaps$sname.ecod)]$x.level
# #   setorder(this.conf.subdomain.ecods.overlaps, qname, sname.ecf, x.level, -bitscore.ecod)
# #   this.conf.subdomain.ecods.overlaps.filtered <- this.conf.subdomain.ecods.overlaps[, .SD[1], by = .(qname, sname.ecf, x.level)]
# #   unique(this.conf.subdomain.ecods.overlaps.filtered$x.level)  
# # })
# # x.levels.per.subdomain.reprseq.l <- sapply(x.levels.per.subdomain.reprseq, length)
# # x.levels.per.subdomain.reprseq[x.levels.per.subdomain.reprseq.l>1]
# # 
# # #############
# # # 3. Sequence similarity of shared (subdomain) fragments
# # #############
# # 
# # # Load sequence data
# # reprseq.seqs <- read.fasta(reprseq.seqs.filename)
# # reprseq.seqs.s <- sapply(reprseq.seqs, c2s)
# # reprseq.seqs.n <- names(reprseq.seqs)
# # reprseq.seqs.table <- data.table(name = reprseq.seqs.n, string = reprseq.seqs.s)
# # setkey(reprseq.seqs.table, name)
# # 
# # # Find all cases of shared subdomain ECFs between different Xs
# # 
# # confirmed.subdomain.index=1
# # # this.conf.subdomain <- ecf.subdomains.confirmed$ecf[confirmed.subdomain.index]
# # # this.conf.subdomain <- ecf.subdomains.confirmed$ecf[c(1:5)]
# # # this.conf.subdomain.ecods.overlaps <- subdomain.ecods.overlaps[.(this.conf.subdomain)]
# # this.conf.subdomain.ecods.overlaps <- subdomain.ecods.overlaps
# # this.conf.subdomain.ecods.overlaps$x.level <- ecod.domains[.(this.conf.subdomain.ecods.overlaps$sname.ecod)]$x.level
# # setorder(this.conf.subdomain.ecods.overlaps, qname, sname.ecf, x.level, -bitscore.ecod)
# # this.conf.subdomain.ecods.overlaps.filtered <- this.conf.subdomain.ecods.overlaps[, .SD[1], by = .(qname, sname.ecf, x.level)]
# # 
# # # Prototype ECF-recombination table
# # 
# # this.input <- this.conf.subdomain.ecods.overlaps.filtered
# # 
# # left.join.table <- this.input %>% select(query.a = qname, qlength.a = qlength, shared.ecf = sname.ecf, ecf.length = slength.ecf,
# #                                          ecf.start.a = qstart.ecf, ecf.end.a = qend.ecf, 
# #                                          ecod.a = sname.ecod, x.level.a = x.level, ecod.start.a = qstart.ecod, ecod.end.a = qend.ecod)
# # 
# # right.join.table <- this.input %>% select(query.b = qname, qlength.b = qlength, shared.ecf = sname.ecf,
# #                                           ecf.start.b = qstart.ecf, ecf.end.b = qend.ecf, 
# #                                           ecod.b = sname.ecod, x.level.b = x.level, ecod.start.b = qstart.ecod, ecod.end.b = qend.ecod)
# # 
# # putative.recomb.pairs.temp <- left.join.table[right.join.table, on = .(shared.ecf), allow.cartesian = TRUE]
# # putative.recomb.pairs.temp <- putative.recomb.pairs.temp[query.a < query.b & x.level.a != x.level.b]
# # 
# # putative.recomb.pairs <- putative.recomb.pairs.temp %>% 
# #   select(query.a, query.b, shared.ecf, qlength.a, qlength.b,
# #          ecf.start.a, ecf.end.a, ecf.start.b, ecf.end.b,
# #          ecod.a, ecod.b, x.level.a, x.level.b,
# #          ecod.start.a, ecod.end.a, ecod.start.b, ecod.end.b)
# # 
# # # Filter 'reprseq.vs.reprseq' to leave only the strongest hit between each pair
# # setorder(reprseq.vs.reprseq, qname, sname, -bitscore)
# # reprseq.vs.reprseq.filtered <- reprseq.vs.reprseq[, .SD[1], by = .(qname, sname)] %>%
# #   select(qname, sname, aln.start.a = qstart, aln.end.a = qend, aln.start.b = sstart, aln.end.b = send,
# #          cov.a = qcov, cov.b = scov, aln.prob = prob, aln.eval = eval, aln.pident = pident)
# # 
# # # Merge full reprseq alignments to add to the 'putative.recomb.pairs' table
# # putative.recomb.pairs.merged <- merge(putative.recomb.pairs, reprseq.vs.reprseq.filtered,
# #                                by.x = c("query.a", "query.b"),
# #                                by.y = c("qname", "sname"),
# #                                all.x = TRUE, all.y = FALSE)
# # 
# # putative.recomb.pairs.filtered <- putative.recomb.pairs.merged[!is.na(putative.recomb.pairs.merged$cov.a)]
# # 
# # # Exclude rows where full sequences do not align
# # # putative.recomb.pairs.filtered <- putative.recomb.pairs.merged[!is.na(putative.recomb.pairs.merged$cov.a)] %>%
# # #   select(query.a, query.b, shared.ecf, ecf.start.a, ecf.end.a, ecf.start.b, ecf.end.b, 
# # #          aln.start.a, aln.end.a, aln.start.b, aln.end.b, aln.pident)
# # 
# # putative.recomb.pairs.filtered[, aln.length.a := aln.end.a - aln.start.a + 1]
# # putative.recomb.pairs.filtered[, ecf.length.a := ecf.end.a - ecf.start.a + 1]
# # putative.recomb.pairs.filtered[, aln.length.b := aln.end.b - aln.start.b + 1]
# # putative.recomb.pairs.filtered[, ecf.length.b := ecf.end.b - ecf.start.b + 1]
# # putative.recomb.pairs.filtered[, overlap.length.a := {
# #   res <- calculate_overlap_parms(ecf.start.a, ecf.end.a, aln.start.a, aln.end.a)
# #   res$overlap.length
# # }]
# # putative.recomb.pairs.filtered[, overlap.length.b := {
# #   res <- calculate_overlap_parms(ecf.start.b, ecf.end.b, aln.start.b, aln.end.b)
# #   res$overlap.length
# # }]
# # putative.recomb.pairs.filtered[, prop.overlap2aln.a := overlap.length.a/aln.length.a]
# # putative.recomb.pairs.filtered[, prop.overlap2ecf.a := overlap.length.a/ecf.length.a]
# # putative.recomb.pairs.filtered[, prop.overlap2aln.b := overlap.length.b/aln.length.b]
# # putative.recomb.pairs.filtered[, prop.overlap2ecf.b := overlap.length.b/ecf.length.b]
# # 
# # OVERLAP.PARM <- 0.8 # move to config later!
# # 
# # putative.recomb.pairs.filtered[, ecf.aln.overlap := (prop.overlap2aln.a > OVERLAP.PARM) & 
# #                                  (prop.overlap2ecf.a > OVERLAP.PARM) &
# #                                  (prop.overlap2aln.b > OVERLAP.PARM) &
# #                                  (prop.overlap2ecf.b > OVERLAP.PARM)]
# # 
# # putative.recomb.pairs.shortlist <- putative.recomb.pairs.filtered[putative.recomb.pairs.filtered$ecf.aln.overlap]
# # putative.recomb.pairs.shortlist[putative.recomb.pairs.shortlist$aln.pident > 50]
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # putative.recomb.pairs.filtered2 <- putative.recomb.pairs.merged[interesting.shortlist]
# # 
# # 
# # # putative.recomb.pairs.filtered2 <- putative.recomb.pairs.filtered %>%
# # #   filter(prop.overlap2aln.a > 0.8,
# # #          prop.overlap2ecf.a > 0.8,
# # #          prop.overlap2aln.b > 0.8,
# # #          prop.overlap2ecf.b > 0.8
# # #   )
# # 
# # 
# # putative.recomb.pairs.merged[1]
# # 
# # x <- putative.recomb.pairs.merged[1]
# # 
# # # aln.length.a <- x$aln.end.a - x$aln.start.a + 1
# # # ecf.length.a <- x$ecf.end.a - x$ecf.start.a + 1
# # # overlap.length.a <- calculate_overlap_parms(x$ecf.start.a, x$ecf.end.a, x$aln.start.a, x$aln.end.a)$overlap.length
# # # prop.overlap2aln.a <- overlap.length.a/aln.length.a
# # # prop.overlap2ecf.a <- overlap.length.a/ecf.length.a
# # # 
# # # aln.length.b <- x$aln.end.b - x$aln.start.b + 1
# # # ecf.length.b <- x$ecf.end.b - x$ecf.start.b + 1
# # # overlap.length.b <- calculate_overlap_parms(x$ecf.start.b, x$ecf.end.b, x$aln.start.b, x$aln.end.b)$overlap.length
# # # prop.overlap2aln.b <- overlap.length.b/aln.length.b
# # # prop.overlap2ecf.b <- overlap.length.b/ecf.length.b
# # 
# # 
# # 
# # # putative.recomb.pairs.merged[, calculate_overlap_parms()
# # 
# # # # extract sequences of all reprseq from the 'putative.recomb.pairs' table
# # #  
# # # # Merge for query.a
# # # putative.recomb.pairs <- merge(putative.recomb.pairs, reprseq.seqs.table, by.x = "query.a", by.y = "name", all.x = TRUE, suffixes = c("", ".a"))
# # # setnames(putative.recomb.pairs, "string", "seq.a")
# # # # Merge for query.b
# # # putative.recomb.pairs <- merge(putative.recomb.pairs, reprseq.seqs.table, by.x = "query.b", by.y = "name", all.x = TRUE, suffixes = c("", ".b"))
# # # setnames(putative.recomb.pairs, "string", "seq.b")
# # # 
# # # # Extract subsequences for query.a using the correct string column
# # # putative.recomb.pairs[, ecf.seq.a := substr(seq.a, ecf.start.a, ecf.end.a)]
# # # # Extract subsequences for query.b using the correct string column
# # # putative.recomb.pairs[, ecf.seq.b := substr(seq.b, ecf.start.b, ecf.end.b)]
# # # # Drop full sequences
# # # putative.recomb.pairs[, c("seq.a", "seq.b") := NULL]  # Removing the full sequence columns to save memory
# # # 
# # # # set column order
# # # all.col.names <- colnames(putative.recomb.pairs)
# # # all.col.names.reorderd <- c(all.col.names[2], all.col.names[1], all.col.names[-(1:2)])
# # # setcolorder(putative.recomb.pairs, all.col.names.reorderd)
# # # 
# # # 
# # # 
# # # 
# # # 
# # # 
# # # 
# # # # reprseq.ecfs.overlap <- reprseq.ecfs %>% select(qname, qlength, qstart.ecf = qstart , qend.ecf = qend,
# # # #                                                 sname.ecf = sname, sstart.ecf = sstart, send.ecf = send, slength.ecf = slength,
# # # #                                                 pident.ecf = pident, hit.length.ecf = hit.length)
# # # 
# # # 
# # # 
# # # 
# # # 
# # # 
# # # # View(this.conf.subdomain.ecods.overlaps.filtered)
# # # # this.conf.subdomain.ecods.overlaps.filtered[c(7, 151), ]
# # # # 
# # # # 
# # # # reprseq.i.rn <- 1
# # # # reprseq.j.rn <- 41
# # # # reprseq.i.table <- this.conf.subdomain.ecods.overlaps.filtered[reprseq.i.rn,] 
# # # # reprseq.j.table <- this.conf.subdomain.ecods.overlaps.filtered[reprseq.j.rn,] 
# # # # 
# # # # reprseq.name.i <- reprseq.i.table$qname
# # # # reprseq.name.i.ecf.start <- reprseq.i.table$qstart.ecf
# # # # reprseq.name.i.ecf.end <- reprseq.i.table$qend.ecf
# # # # reprseq.seq.i <- reprseq.seqs.s[reprseq.seqs.n == reprseq.name.i]
# # # # reprseq.seq.i.frag <- substr(reprseq.seq.i, reprseq.name.i.ecf.start, reprseq.name.i.ecf.end)
# # # # 
# # # # reprseq.name.j <- reprseq.j.table$qname
# # # # reprseq.name.j.ecf.start <- reprseq.j.table$qstart.ecf
# # # # reprseq.name.j.ecf.end <- reprseq.j.table$qend.ecf
# # # # reprseq.seq.j <- reprseq.seqs.s[reprseq.seqs.n == reprseq.name.j]
# # # # reprseq.seq.j.frag <- substr(reprseq.seq.j, reprseq.name.j.ecf.start, reprseq.name.j.ecf.end)
# # # # 
# # # # seq1 <- AAString(reprseq.seq.i.frag)
# # # # seq2 <- AAString(reprseq.seq.j.frag)
# # # # alignment <- pairwiseAlignment(pattern = seq1, subject = seq2, substitutionMatrix = "BLOSUM62", gapOpening = 10, gapExtension = 1, type = "global")
# # # # aligned_pattern <- as.character(alignment@pattern)
# # # # aligned_subject <- as.character(alignment@subject)
# # # # # Function to calculate percentage identity
# # # # calculate_percentage_identity <- function(aligned_seq1, aligned_seq2) {
# # # #   matches <- sum(aligned_seq1 == aligned_seq2)  # Count matches
# # # #   alignment_length <- length(aligned_seq1)     # Total positions in the alignment
# # # #   percent_identity <- (matches / alignment_length) * 100  # Calculate percentage
# # # #   return(percent_identity)
# # # # }
# # # # percentage_identity <- calculate_percentage_identity(s2c(aligned_pattern), s2c(aligned_subject))
# # # 
