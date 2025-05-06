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

# Load data
ecod.domains <- fread(ecod.domains.output)
setkey(ecod.domains, full.id)
ecf.table <- fread(ecf.table.output)

# Filter hits
ecfs.vs.ecod <- fread(ecfs.vs.ecod.output)
ecfs.full.vs.ecod <- ecfs.vs.ecod %>%
  filter(qcov >= MIN.QCOV.FOR.SUBDOMAIN) %>%
  filter(prob >= MIN.PROB.FOR.SUBDOMAIN) %>%
  filter(eval <= MAX.EVAL.FOR.SUBDOMAIN)
setkey(ecfs.full.vs.ecod, qname)


# index=1
# this.ecf <- ecf.name.table$ecf.name[index]
# this.ecf.hits <- ecfs.full.vs.ecod[.(this.ecf)]
# this.ecf.hits$t.id <- ecod.domains[.(this.ecf.hits$sname)]$t_id

# ecfs.ecod1.overlap <- this.ecf.hits %>% select(qname, qlength,
#                                                 sname.ecod1 = sname, sstart.ecod1 = sstart, send.ecod1 = send, slength.ecod1 = slength,
#                                                 pident.ecod1 = pident, hit.length.ecod1 = hit.length, t.id.ecod1 = t.id)
# ecfs.ecod2.overlap <- this.ecf.hits %>% select(qname,
#                                                sname.ecod2 = sname, sstart.ecod2 = sstart, send.ecod2 = send, slength.ecod2 = slength,
#                                                pident.ecod2 = pident, hit.length.ecod2 = hit.length, t.id.ecod2 = t.id)

# ecf.subset <- ecf.name.table$ecf.name[1:1000]
# ecfs.full.vs.ecod.sub <- ecfs.full.vs.ecod[qname %in% ecf.subset]
# ecfs.full.vs.ecod.sub$t.id <- ecod.domains[.(ecfs.full.vs.ecod.sub$sname)]$t_id

ecfs.full.vs.ecod$t.id <- ecod.domains[.(ecfs.full.vs.ecod$sname)]$t_id

ecfs.ecod1.overlap <- ecfs.full.vs.ecod %>% select(qname, qlength, sname.ecod1 = sname, t.id.ecod1 = t.id,
                                                       sstart.ecod1 = sstart, send.ecod1 = send, slength.ecod1 = slength,
                                                       scov.ecod1 = scov, pident.ecod1 = pident)
ecfs.ecod2.overlap <- ecfs.full.vs.ecod %>% select(qname, sname.ecod2 = sname, t.id.ecod2 = t.id,
                                                       sstart.ecod2 = sstart, send.ecod2 = send, slength.ecod2 = slength,
                                                       scov.ecod2 = scov, pident.ecod2 = pident)
shared.ecf.table <- ecfs.ecod1.overlap[ecfs.ecod2.overlap, on = .(qname), allow.cartesian = TRUE]
shared.ecf.table <- shared.ecf.table[sname.ecod1 < sname.ecod2]

###
#||

linked.ecod.pairs.all <- shared.ecf.table %>%
  select(ecod1.name = sname.ecod1, ecod1.t.id = t.id.ecod1,
         ecod2.name = sname.ecod2, ecod2.t.id = t.id.ecod2,
         ecf = qname, ecf.length = qlength,
         ecod1.start = sstart.ecod1, ecod1.end = send.ecod1, ecod1.length = slength.ecod1, 
         ecod2.start = sstart.ecod2, ecod2.end = send.ecod2, ecod2.length = slength.ecod2,
         ecod1.cov = scov.ecod1, ecod1.pident = pident.ecod1,
         ecod2.cov = scov.ecod2, ecod2.pident = pident.ecod2)

# Remove redundant pairs
linked.ecod.pairs.all[, `:=`(
  ecod1_sym = pmin(ecod1.name, ecod2.name),
  ecod2_sym = pmax(ecod1.name, ecod2.name)
)]
unique_pairs <- unique(linked.ecod.pairs.all[, .(ecod1_sym, ecod2_sym)])
unique_pairs[, pair.index := .I]
linked.ecod.pairs.all.index <- merge(
  linked.ecod.pairs.all,
  unique_pairs[, .(ecod1_sym, ecod2_sym, pair.index)],
  by = c("ecod1_sym", "ecod2_sym"),
  all.x = TRUE
)
setorder(linked.ecod.pairs.all.index, pair.index, ecf)
# linked.ecod.pairs.all.filtered <- linked.ecod.pairs.all.index[, .SD[1], by = .(pair.index, ecf)] %>% select(c(pair.index, output.colnames))
linked.ecod.pairs.all.filtered <- linked.ecod.pairs.all.index[, .SD[1], by = .(pair.index, ecf)]
linked.ecod.pairs.all.filtered[, cov := round(ecod1.cov * ecod2.cov, d = 4)]

#||
###


linked.ecod.pairs <- shared.ecf.table[t.id.ecod1 != t.id.ecod2] %>%
  select(ecod1.name = sname.ecod1, ecod1.t.id = t.id.ecod1,
         ecod2.name = sname.ecod2, ecod2.t.id = t.id.ecod2,
         ecf = qname, ecf.length = qlength,
         ecod1.start = sstart.ecod1, ecod1.end = send.ecod1, ecod1.length = slength.ecod1, 
         ecod2.start = sstart.ecod2, ecod2.end = send.ecod2, ecod2.length = slength.ecod2,
         ecod1.cov = scov.ecod1, ecod1.pident = pident.ecod1,
         ecod2.cov = scov.ecod2, ecod2.pident = pident.ecod2)
output.colnames <- colnames(linked.ecod.pairs)

# Remove redundant pairs
linked.ecod.pairs[, `:=`(
  ecod1_sym = pmin(ecod1.name, ecod2.name),
  ecod2_sym = pmax(ecod1.name, ecod2.name)
)]
unique_pairs <- unique(linked.ecod.pairs[, .(ecod1_sym, ecod2_sym)])
unique_pairs[, pair.index := .I]
linked.ecod.pairs.index <- merge(
  linked.ecod.pairs,
  unique_pairs[, .(ecod1_sym, ecod2_sym, pair.index)],
  by = c("ecod1_sym", "ecod2_sym"),
  all.x = TRUE
)
setorder(linked.ecod.pairs.index, pair.index, ecf)
linked.ecod.pairs.filtered <- linked.ecod.pairs.index[, .SD[1], by = .(pair.index, ecf)]
# linked.ecod.pairs.filtered <- linked.ecod.pairs.index[, .SD[1], by = .(pair.index, ecf)] %>% select(c(pair.index, output.colnames))

# Add ECOD comparisons
# Define a function to extract hierarchical levels
extract_ecod_levels <- function(ecod_id, levels = 1) {
  # Split the ECOD ID by "."
  parts <- strsplit(ecod_id, "\\.", fixed = FALSE)
  # Extract the desired number of levels and collapse into a string
  sapply(parts, function(x) paste(x[1:levels], collapse = "."))
}
linked.ecod.pairs.filtered[, `:=`(
  same.x.level = extract_ecod_levels(ecod1.t.id, levels = 1) ==  extract_ecod_levels(ecod2.t.id, levels = 1),
  same.h.level = extract_ecod_levels(ecod1.t.id, levels = 2) ==  extract_ecod_levels(ecod2.t.id, levels = 2)
)]

# calculate bidirectional coverage
linked.ecod.pairs.filtered[, cov := round(ecod1.cov * ecod2.cov, d = 4)]

# assign ECOD comparison categories
linked.ecod.pairs.filtered$ecod.comp <- "same.H"
linked.ecod.pairs.filtered[same.x.level & !same.h.level]$ecod.comp <- "diff.H"
linked.ecod.pairs.filtered[same.x.level == F]$ecod.comp <- "diff.X"


setorder(linked.ecod.pairs.filtered, pair.index, -cov)
linked.ecod.pairs.filtered2 <- linked.ecod.pairs.filtered[, .SD[1], by = .(pair.index)]

setorder(linked.ecod.pairs.filtered, ecf, -cov)
linked.ecod.pairs.filtered3 <- linked.ecod.pairs.filtered[, .SD[1], by = .(ecf)]

#######################
# Visualise
# Step 1: Create bins for the `cov` column
linked.ecod.pairs.filtered2[, cov_bin := cut(cov, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)]

# Step 2: Aggregate data to count occurrences of `ecod.comp` classes in each bin
agg_data <- linked.ecod.pairs.filtered2[, .N, by = .(cov_bin, ecod.comp)]
agg_data$ecod.comp <- factor(agg_data$ecod.comp, levels = c("same.H", "diff.H", "diff.X"))

# Step 3: Create the bar plot
p.ecf.vs.ecod.prop <- ggplot(agg_data, aes(x = cov_bin, y = N, fill = ecod.comp)) +
  geom_bar(stat = "identity", position = "fill") +  # Position = "fill" for proportions
  scale_fill_brewer(palette = "Set3") +  # Customise fill colours
  theme_minimal() +
  labs(
    title = "Distribution of ECOD Comp Classes Across Cov Ranges",
    x = "Cov Range",
    y = "Proportion",
    fill = "ECOD Comp Class"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

# Step 3: Create the bar plot with actual counts
p.ecf.vs.ecod.abs <- ggplot(agg_data, aes(x = cov_bin, y = N, fill = ecod.comp)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bar to show actual counts
  scale_fill_brewer(palette = "Set3") +  # Customise fill colours
  theme_minimal() +
  labs(
    title = "Actual Counts of ECOD Comp Classes Across Cov Ranges",
    x = "Cov Range",
    y = "Count",
    fill = "ECOD Comp Class"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +   # Rotate x-axis labels
  scale_y_sqrt()

ggsave("plo1_ecf-vs-ecod-absolute.pdf", p.ecf.vs.ecod.abs, width = 7, height = 5)
ggsave("plo2_ecf-vs-ecod-prop.pdf", p.ecf.vs.ecod.prop, width = 7, height = 5)


