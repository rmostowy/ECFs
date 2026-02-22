### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(igraph)))
suppressWarnings(suppressPackageStartupMessages(library(patchwork)))
# suppressWarnings(suppressPackageStartupMessages(library(seqinr)))
# suppressWarnings(suppressPackageStartupMessages(library(Biostrings)))
# suppressWarnings(suppressPackageStartupMessages(library(parallel)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

# Define output folder
table.dir <- "tables"
output.dir <- file.path(visual.dir, "part5_rem", "rr-network")
make.dir(output.dir)

# Load reprseq table
reprseq.table <- fread(reprseq.table.output)
all.reprseq <- unique(reprseq.table$reprseq)

# Load recent recombination table
rem.table <- fread(rem.table.output)
rem.table[, cov.max.full := pmax(qcov.full, scov.full, na.rm = T)]

# Load functional annotations
reprseq.functions <- fread(reprseq.functions.unique.output)
setkey(reprseq.functions, reprseq)

# Load ECOD domain table
ecod.domains <- fread(ecod.domains.output)
setkey(ecod.domains, full.id)

# Load reprseq vs ECOD
reprseq.ecod.all <- fread(reprseq.domains.ecod.output)
reprseq.ecod <- reprseq.ecod.all %>%
  filter(prob >= MIN.PROB.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(scov >= MIN.SCOV.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE)
setkey(reprseq.ecod, qname)

reprseq2ecod <- reprseq.ecod %>% select(qname, sname, bitscore)
setorder(reprseq2ecod, qname, -bitscore)
reprseq2ecod <- reprseq2ecod[, .SD[1], by = .(qname)]
reprseq2ecod$t.id <- ecod.domains[.(reprseq2ecod$sname)]$t_id
reprseq2ecod$t.name <- ecod.domains[.(reprseq2ecod$sname)]$t_name
reprseq2ecod.final <- reprseq2ecod %>% select(reprseq = qname, t.id, t.name)

reprseq.metadata <- reprseq.functions %>% select(reprseq, phrogs.score = score, funct, class, annotation)
reprseq.metadata <- merge(reprseq.metadata, reprseq2ecod.final, by = 'reprseq', all.x = T) 
reprseq.metadata.outfile <- file.path(output.dir, "reprseq-table.txt")
fwrite(reprseq.metadata, file = reprseq.metadata.outfile, sep = "\t")
# setkey(reprseq.metadata, reprseq)

rr.network <- rem.table %>% select(qname, sname, qcov, scov, pident, hit.length = hit.length.query)
rr.network.outfile <- file.path(output.dir, "rem-network.txt")
fwrite(rr.network, file = rr.network.outfile, sep = "\t")

########################################
## Plot most frequent PHROGs and ECODs 
# across the full dataset and in the 
# recent-rec dataset
########################################

rr.reprseq <- sort(unique(rem.table$qname, rem.table$sname))

# Function to plot top N categories
plot_top_n <- function(data, column, n = 20, title) {
  data[, (column) := fifelse(is.na(get(column)) | get(column) == "", "Unknown", get(column))]  # Handle missing values
  top_n <- data[, .N, by = get(column)][order(-N)][1:n]  # Get top N most frequent categories
  
  ggplot(top_n, aes(x = reorder(get, -N), y = N, fill = get)) +
    geom_bar(stat = "identity") +
    coord_flip() +  # Horizontal bars for readability
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = title, x = column, y = "Frequency")
  
}

funct.plot.all <- plot_top_n(reprseq.metadata[funct != "Unknown"], "funct", 20, "Top 20 Functions (all proteins)")
funct.plot.rr <- plot_top_n(reprseq.metadata[funct != "Unknown" & reprseq %chin% rr.reprseq], "funct", 20, "Top 20 Functions (REM proteins)")
funct.plots <- funct.plot.all / funct.plot.rr  # Stack function plots
funct.plots.outname <- file.path(output.dir, "plot1_funct-plot.png")
ggsave(funct.plots.outname, funct.plots, width = 6, height = 12)

ecod.plot.all <- plot_top_n(reprseq.metadata[t.name != "Unknown"], "t.name", 20, "Top 20 ECODs (all proteins)")
ecod.plot.rr <- plot_top_n(reprseq.metadata[t.name != "Unknown" & reprseq %chin% rr.reprseq], "t.name", 20, "Top 20 ECODs (REM proteins)")
ecod.plots <- ecod.plot.all / ecod.plot.rr  # Stack function plots
ecod.plots.outname <- file.path(output.dir, "plot2_ecod-plot.png")
ggsave(ecod.plots.outname, ecod.plots, width = 8, height = 12)

# nc.pairs <- relative.classes[(protein1.class == 'N' & protein2.class == "C") | (protein1.class == "C" & protein2.class == "N")]
# nc.pairs.reprseqs <- unique(nc.pairs$protein1.name, nc.pairs$protein2.name)
# plot_top_n(reprseq.metadata[t.name != "Unknown" & reprseq %chin% nc.pairs.reprseqs], "t.name", 10, "Top 20 Functions (all proteins)")
