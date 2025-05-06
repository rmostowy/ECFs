### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(Biostrings)))
# suppressWarnings(suppressPackageStartupMessages(library(reshape2)))
# suppressWarnings(suppressPackageStartupMessages(library(ggh4x)))
# suppressWarnings(suppressPackageStartupMessages(library(tidyr)))
# suppressWarnings(suppressPackageStartupMessages(library(seqinr)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

# Load all datasets
ecf.table <- fread(ecf.table.output)
ecf.family.table <- fread(ecf.family.table.output)

reprseq.vs.ecfs <- fread(reprseq.ecfs.output)

ecod.domains <- fread(ecod.domains.output)
pfam.domains <- fread(pfam.domains.output)

reprseq.vs.ecod <- fread(reprseq.domains.ecod.output)
reprseq.vs.ecod <- reprseq.vs.ecod %>% filter(qcov >= MIN.QCOV.FOR.DOMAIN, scov >= MIN.SCOV.FOR.DOMAIN)

reprseq.vs.pfam <- fread(reprseq.domains.pfam.output)
reprseq.vs.pfam <- reprseq.vs.pfam %>% filter(qcov >= MIN.QCOV.FOR.DOMAIN, scov >= MIN.SCOV.FOR.DOMAIN)

reprseq.seqs <- readAAStringSet(reprseq.seqs.filename)
reprseq.seqs.n <- names(reprseq.seqs)
reprseq.seqs.l <- width(reprseq.seqs)

reprseqs.with.domains <- union(unique(reprseq.vs.ecfs$qname), unique(reprseq.vs.ecod$qname))
reprseqs.with.domains <- union(reprseqs.with.domains, unique(reprseq.vs.pfam$qname))


###############################################################
### 1. Plot the distribution of the number of reprseqs per ECF
###############################################################
this.results.dir <- "analysis1_distr"
make.dir(this.results.dir)

reprseq.per.ecf <- reprseq.vs.ecfs[, .(no.reprseqs = uniqueN(qname)), by = sname]

p.repr.per.ecf <- ggplot(reprseq.per.ecf, aes(x = no.reprseqs)) +
  geom_histogram(bins = 100, fill = "blue", color = "black", alpha = 0.7) + scale_x_continuous(limits = c(1,200))

p.repr.per.ecf.logx <- ggplot(reprseq.per.ecf, aes(x = no.reprseqs)) +
  geom_histogram(bins = 50, fill = "blue", color = "black", alpha = 0.7) + scale_x_log10()

p.repr.per.ecf.fn <- file.path(this.results.dir, "plot1a_repr-per-ecf.png")
p.repr.per.ecf.logx.fn <- file.path(this.results.dir, "plot1b_repr-per-ecf-log.png")
ggsave(p.repr.per.ecf.fn, p.repr.per.ecf, width = 6, height = 5)
ggsave(p.repr.per.ecf.logx.fn, p.repr.per.ecf.logx, width = 5, height = 5)

###############################################################
### 2. Plot the distribution of the number of ECFs per EF
###############################################################

# Rank the families by size
ecfs.per.family <- ecf.family.table[, .(no.ecfs = uniqueN(ecf)), by = family]
# Define bins for the number of ECFs
ecfs.per.family[, bin := cut(
  no.ecfs,
  breaks = c(-Inf, 1, 5, 10, 50, 100, 500, Inf),  # Define bin edges
  labels = c("1", "2–5", "6–10", "11–50", "51-100", "101-500", ">500"),
  right = TRUE
)]

# Count the number of families in each bin
family_counts <- ecfs.per.family[, .N, by = bin]
setorder(family_counts, bin)  # Sort by bin order

# Plot the binned data
p.no.fam.per.no.ecfs <- ggplot(family_counts, aes(x = bin, y = N)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = N), vjust = -0.5, size = 4) +
  labs(
    x = "ECF count",
    y = "Number of Families with this ECF count"
  ) +
  theme_bw(base_size = 14)

p.no.fam.per.no.ecfs.fn <- file.path(this.results.dir, "plot2a_no-fams-per-no-ecfs.png")
ggsave(p.no.fam.per.no.ecfs, file = p.no.fam.per.no.ecfs.fn, width = 6, height = 5)


###############################################################
### 3. Plot the distribution of ECF lengths (as well as for ECODs and Pfams)
###############################################################

reprseq.lengths <- data.table(name = reprseq.seqs.n, length = reprseq.seqs.l)
reprseq.lengths[, db := "rHMMs"]
reprseq.lengths.filtered <- reprseq.lengths[name %in% reprseqs.with.domains]

ecf.lengths <- ecf.table %>% select(name = ecf.name, length)
ecf.lengths[, db := "ECF"]
ecod.lengths <- ecod.domains %>% select(name = full.id, length)
ecod.lengths[, db := "ECOD"]
pfam.lengths <- pfam.domains %>% select(name = accession, length)
pfam.lengths[, length := round(length, d = 0)]
pfam.lengths[, db := "PFAM"]

all.lengths <- rbind(reprseq.lengths.filtered, ecf.lengths, ecod.lengths, pfam.lengths)

# Calculate medians for each database
medians <- all.lengths[, .(median_length = median(length)), by = db]
medians$db <- factor(medians$db, levels = unique(all.lengths$db))

means <- all.lengths[, .(mean_length = mean(length)), by = db]
means$db <- factor(means$db, levels = unique(all.lengths$db))

all.lengths$db <- factor(all.lengths$db, levels = unique(all.lengths$db))
# Create the plot
p.length.distr <- ggplot(all.lengths, aes(x = length, fill = db)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.6, bins = 80) +  # Histogram with transparency
  # geom_vline(data = medians, aes(xintercept = median_length, color = db), 
  #            linetype = "dashed", size = 1) +  # Median line
  geom_vline(data = means, aes(xintercept = mean_length, color = db), 
             linetype = "dashed", size = 1) +  # Median line
  facet_grid(db ~ ., scales = "free_x") +  # Align distributions vertically
  scale_fill_manual(values = c("rHMMs" = "red", "ECF" = "#1f77b4", "PFAM" = "#ff7f0e", "ECOD" = "#2ca02c")) +  # Custom colors
  scale_color_manual(values = c("rHMMs" = "red", "ECF" = "#1f77b4", "PFAM" = "#ff7f0e", "ECOD" = "#2ca02c")) +  # Matching colors
  scale_x_log10() +
  labs(
    x = "Length",
    y = "Density",
    fill = "Database"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = NULL,
    strip.text = element_text(face = "bold")
  )

p.length.distr.fn <- file.path(this.results.dir, "plot3-length-distr.png")
ggsave(p.length.distr.fn, p.length.distr, width = 6, height = 5)


