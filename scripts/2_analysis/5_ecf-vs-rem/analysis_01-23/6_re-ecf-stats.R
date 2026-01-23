### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(stringr)))
suppressWarnings(suppressPackageStartupMessages(library(ggalluvial)))
suppressWarnings(suppressPackageStartupMessages(library(Biostrings)))
# suppressWarnings(suppressPackageStartupMessages(library(parallel)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

source("functions.R")

# Define output folder
table.dir <- "tables"
output.dir <- file.path(visual.dir, "analysis_02-21", "rec-equiv-ecfs")
make.dir(output.dir)

# Load refseq table
reprseq.table <- fread(reprseq.table.filename)
setkey(reprseq.table, reprseq)
all.reprseq <- unique(reprseq.table$reprseq)

# Calculate NCBI overlaps per reprseq
max.overlap.per.reprseq <- reprseq.table[, .(reprseq, ncbi.overlap)]
setorder(max.overlap.per.reprseq, reprseq, -ncbi.overlap)
max.overlap.per.reprseq <- max.overlap.per.reprseq[, .SD[1], by = reprseq]
setkey(max.overlap.per.reprseq, reprseq)

# Load recent recombination table
recent.rec.table <- fread(reprseq.recent.rec.table.output)
recent.rec.table[, cov.max.full := pmax(qcov.full, scov.full, na.rm = T)]

# Load recombination domain overlap table
rec.domain.overlap.table.outfile <- file.path(table.dir, "rec-domain-overlaps.txt")
rec.domain.overlaps <- fread(rec.domain.overlap.table.outfile)
rec.domain.overlaps[, row.id := .I]

# Load ECF table
ecf.table <- fread(ecf.table.output)
setkey(ecf.table, ecf.name)

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

# Load reprseq vs ECFs
reprseq.ecfs.all <- fread(reprseq.ecfs.output)
reprseq.ecfs <- reprseq.ecfs.all %>%
  filter(prob >= MIN.PROB.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(scov >= MIN.SCOV.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE)
reprseq.ecfs$EF <- ecf.family.table[.(reprseq.ecfs$sname)]$family
setkey(reprseq.ecfs, sname)

# Load reprseq vs ECFs relaxed
reprseq.ecfs.relaxed <- fread(reprseq.ecfs.all.output)

# Load reprseq vs ECOD
reprseq.ecod.all <- fread(reprseq.domains.ecod.output)
reprseq.ecod <- reprseq.ecod.all %>%
  filter(prob >= MIN.PROB.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(scov >= MIN.SCOV.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE)
reprseq.ecod$t.id <- ecod.domains[.(reprseq.ecod$sname)]$t_id
setkey(reprseq.ecod, qname)

# Load reprseq vs Pfam
reprseq.pfam.all <- fread(reprseq.domains.pfam.output)
reprseq.pfam <- reprseq.pfam.all %>%
  filter(prob >= MIN.PROB.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(scov >= MIN.SCOV.FOR.DOMAIN.EQUIVALENCE) %>%
  filter(eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE)
setkey(reprseq.pfam, qname)

# Load ECFs vs ECOD
ecfs.vs.ecod <- fread(ecfs.vs.ecod.output)
setkey(ecfs.vs.ecod, qname)

# Load ECFs vs PFAM
ecfs.vs.pfam <- fread(ecfs.vs.pfam.output)
setkey(ecfs.vs.pfam, qname)


match2type <- function(qcov.value, scov.value){
  final.type <- "new"
  both.values.non.na <- !is.na(qcov.value) & !is.na(scov.value)
  if(both.values.non.na){
    if(qcov.value >= MIN.QCOV.FOR.DOMAIN.EQUIVALENCE & scov.value >= MIN.SCOV.FOR.DOMAIN.EQUIVALENCE){
      final.type <- "equivalent"
    } else if(qcov.value >= MIN.QCOV.FOR.DOMAIN.EQUIVALENCE & scov.value < MIN.SCOV.FOR.DOMAIN.EQUIVALENCE){
      final.type <- "subdomain" 
    } else if(qcov.value < MIN.QCOV.FOR.DOMAIN.EQUIVALENCE & scov.value >= MIN.SCOV.FOR.DOMAIN.EQUIVALENCE){
      final.type <- "domain+"
    } else{
      final.type <- "similar"
    }
  }
  return(final.type)
}


re.ecf.common.outfile <- file.path(output.dir, "re-ecf-common.txt")
re.ecf.common <- fread(re.ecf.common.outfile)

# Compare ECFs to domains
ecf.best.match.parms.list <- lapply(1:nrow(re.ecf.common), function(index){
  # cat(index,"\n")
  this.ecf <- re.ecf.common$ecf[index]
  this.family <- re.ecf.common$family[index]
  
  this.ecf.reprseq <- reprseq.ecfs[.(this.ecf)]
  
  this.ecf.freq <- length(unique(this.ecf.reprseq$qname))
  this.ecf.has.ecod.hits <- this.ecf %in% ecfs.vs.ecod$qname
  best.ecod.match.type <-  "new"
  best.ecod.match.name <- this.ecf.most.freq.function <- "unknown"
  best.ecod.match.sname <- best.ecod.match.id <- NA
  if(this.ecf.has.ecod.hits){
    this.ecf.vs.ecod <- ecfs.vs.ecod[.(this.ecf)]  
    this.ecf.vs.ecod.max.scov <- this.ecf.vs.ecod[which.max(scov)]  
    best.ecod.match.type <- match2type(this.ecf.vs.ecod.max.scov$qcov, this.ecf.vs.ecod.max.scov$scov)
    best.ecod.match.sname <- this.ecf.vs.ecod.max.scov$sname
    best.ecod.match.id <- ecod.domains[.(best.ecod.match.sname)]$t_id
    best.ecod.match.name <- ecod.domains[.(best.ecod.match.sname)]$t_name
  }
  
  this.ecf.has.pfam.hits <- this.ecf %in% ecfs.vs.pfam$qname
  best.pfam.match.type <- "new"
  best.pfam.match.name <- "unknown"
  best.pfam.match.sname <- best.pfam.match.id <- NA
  if(this.ecf.has.pfam.hits){
    this.ecf.vs.pfam <- ecfs.vs.pfam[.(this.ecf)]  
    this.ecf.vs.pfam.max.scov <- this.ecf.vs.pfam[which.max(scov)]  
    best.pfam.match.type <- match2type(this.ecf.vs.pfam.max.scov$qcov, this.ecf.vs.pfam.max.scov$scov)
    best.pfam.match.sname <- get.pfam.accession(this.ecf.vs.pfam.max.scov$sname)
    best.pfam.match.id <- pfam.domains[.(best.pfam.match.sname)]$id
    best.pfam.match.name <- pfam.domains[.(best.pfam.match.sname)]$description
  }
  
  this.ecf.functions <- reprseq.functions[.(this.ecf.reprseq$qname)]$funct
  this.ecf.functions <- this.ecf.functions[this.ecf.functions != ""]
  if(length(this.ecf.functions) > 0) this.ecf.most.freq.function <- names(my.table(this.ecf.functions)[1])

  data.table(ecf = this.ecf, family = this.family, freq = this.ecf.freq,
             ecod.type = best.ecod.match.type, ecod.id = best.ecod.match.id, ecod.name = best.ecod.match.name,
             pfam.type = best.pfam.match.type, pfam.id = best.pfam.match.id, pfam.name = best.pfam.match.name,
             funct = this.ecf.most.freq.function
             )
})
ecf.best.match.parms <- rbindlist(ecf.best.match.parms.list)
setorder(ecf.best.match.parms, family, ecf)

#############################################
## Plot: visualise best match table:
## Category vs ECOD
#############################################

# Define number of top ECFs to retain
top_n_ecfs <- 30  # Change as needed

# Filter for the top N most frequent ECFs
top_ecfs <- ecf.best.match.parms %>%
  arrange(desc(freq)) %>%
  head(top_n_ecfs)

# Create a heatmap-friendly dataframe
heatmap_data <- top_ecfs %>%
  select(ecf, ecod.type, ecod.name) %>%
  mutate(ecod.name = ifelse(is.na(ecod.name), "unknown", ecod.name))

# Define color palette for known ecod.name values
ecod_colors <- RColorBrewer::brewer.pal(n = length(unique(heatmap_data$ecod.name)) - 1, "Paired")
fill_colors <- c(setNames(ecod_colors, unique(heatmap_data$ecod.name)[unique(heatmap_data$ecod.name) != "unknown"]), "unknown" = "grey")

p.most.freq.re.ecfs <- ggplot(heatmap_data, aes(x = ecf, y = factor(ecod.type, levels = c("new", "equivalent", "subdomain", "domain+", "similar")), fill = ecod.name)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = fill_colors, na.value = "grey") +
  theme_bw() +
  labs(
    title = "Heatmap of Top ECFs by ECOD Type",
    x = "ECF",
    y = "Category",
    fill = "ECOD (T-Name)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p.most.freq.re.ecfs.outfile <- file.path(output.dir, "plot1_most-freq-RE-ECFs.png")
ggsave(p.most.freq.re.ecfs.outfile, p.most.freq.re.ecfs, width = 15, height = 4)


#############################################
## Plot: visualise unknown ECFs
## PFAM and PHROG annotations
#############################################

# Subset data: Focus only on ECFs where ECOD name is "unknown" (i.e., no ECOD hit)
ecf_no_ecod <- ecf.best.match.parms %>%
  filter(ecod.name == "unknown", !is.na(pfam.name), !is.na(funct))

# Count occurrences for Sankey plot
sankey_data <- ecf_no_ecod %>%
  count(pfam.name, funct)

# Create Sankey plot
plot.unknown.sankey <- ggplot(sankey_data, aes(axis1 = pfam.name, axis2 = funct, y = n)) +
  geom_alluvium(aes(fill = pfam.name), alpha = 0.8) +
  geom_stratum(fill = "gray", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("PFAM Name", "Function (PHROGs)")) +
  theme_bw() +
  labs(title = "Sankey Diagram: Mapping PFAM Names to Functions for ECFs with No ECOD Hit",
       x = NULL,
       y = "Count",
       fill = "PFAM Name") +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        legend.position = "none")

p.unknown.outfile <- file.path(output.dir, "plot2_unknown-RE-ECFs.png")
ggsave(p.unknown.outfile, plot.unknown.sankey, width = 10, height = 6)


##############################################################################
##############################################################################
##############################################################################
##############################################################################


# # # Compare ECFs to PHROGs
# # ecf2funct.list <- lapply(1:nrow(ecf.best.match.parms), function(index){
# #   this.ecf <- ecf.best.match.parms$ecf[index]
# #   this.family <- ecf.best.match.parms$family[index]
# #   this.ecf.ecod.id <- ecf.best.match.parms$ecod.id[index]
# #   this.ecf.ecod.name <- ecf.best.match.parms$ecod.name[index]
# #   this.ecf.reprseqs <- reprseq.ecfs[.(this.ecf)]$qname
# #   data.table(ecf = this.ecf, family = this.family, ecod.id = this.ecf.ecod.id, ecod.name = this.ecf.ecod.name, reprseq = this.ecf.reprseqs)  
# # })
# # ecf2funct <- rbindlist(ecf2funct.list)
# # ecf2funct$funct <- reprseq.functions[.(ecf2funct$reprseq)]$funct
# # ecf2funct$class <- reprseq.functions[.(ecf2funct$reprseq)]$class
# # ecf2funct[funct == "", funct := "unknown"]
# # ecf2funct[class == "", class := "unknown"]
# # setorder(ecf2funct, reprseq)
# 
# 
# #############################################
# ## Plot: compare ECOD types with PFAM types
# #############################################
# 
# # Filter out NA values
# sankey_data <- ecf.best.match.parms %>%
#   filter(!is.na(ecod.type), !is.na(pfam.type)) %>%
#   count(ecod.type, pfam.type)
# 
# # Sankey plot
# plot.type.assignment.ecod.pfam <- ggplot(sankey_data, aes(axis1 = ecod.type, axis2 = pfam.type, y = n)) +
#   geom_alluvium(aes(fill = ecod.type), alpha = 0.8) +
#   geom_stratum(fill = "grey") +
#   geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 4) +
#   scale_x_discrete(limits = c("ECOD Type", "PFAM Type")) +
#   theme_minimal() +
#   labs(title = "Sankey Diagram: Mapping ECOD Types to PFAM Types",
#        x = NULL,
#        y = "Count",
#        fill = "ECOD Type") +
#   theme(axis.text.x = element_text(size = 12, face = "bold"))
# 
# #############################################
# ## Plot: domain types with their functions
# #############################################
# 
# # Aggregate count
# heatmap_data_ecod <- ecf.best.match.parms %>%
#   count(ecod.type, ecod.name)
# heatmap_data_pfam <- ecf.best.match.parms %>%
#   count(pfam.type, pfam.name)
# 
# # Plot
# p.heatmap.ecod <- ggplot(heatmap_data, aes(x = ecod.type, y = ecod.name, fill = n)) +
#   geom_tile() +
#   scale_fill_viridis_c() +
#   theme_minimal() +
#   labs(title = "Heatmap of ECF Type Frequency per Domain Function",
#        x = "ECF Type",
#        y = "Domain Function",
#        fill = "Count")
# 
# p.heatmap.pfam <- ggplot(heatmap_data_pfam, aes(x = pfam.type, y = pfam.name, fill = n)) +
#   geom_tile() +
#   scale_fill_viridis_c() +
#   theme_minimal() +
#   labs(title = "Heatmap of ECF Type Frequency per Domain Function",
#        x = "ECF Type",
#        y = "Domain Function",
#        fill = "Count")
# 
# #############################################
# ## Plot: relationship between domains and functions of RE-ECFs
# #############################################
# 
# 
# # ###########################################################################
# # ###########################################################################
# # ###########################################################################
# # ###########################################################################
# # 
# # # ecf.list <- re.ecf.common[type != 'new']$ecf
# # 
# # 
# # make.qcov.scov.plot <- function(table){
# #   plot <- ggplot(table, aes(x = qcov, y = scov)) + geom_point(size = 4) 
# #   plot <- plot + scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0, 1))
# #   plot <- plot + theme_bw()
# #   return(plot)
# # }
# # 
# # index=5
# # # this.ecf <- ecf.list[index]
# # 
# # 
# # # 1. Find all ECODs linked to this ECF
# # ecf.best.match.parms.list <- lapply(1:nrow(re.ecf.common), function(index){
# #   this.ecf <- re.ecf.common$ecf[index]
# #   this.family <- re.ecf.common$family[index]
# #   this.ecf.has.ecod.hits <- this.ecf %in% ecfs.vs.ecod$qname
# #   best.match.scov <- best.match.qcov <- best.match.ecod <- NA
# #   if(this.ecf.has.ecod.hits){
# #     this.ecf.vs.ecod <- ecfs.vs.ecod[.(this.ecf)]  
# #     this.ecf.vs.ecod.max.scov <- this.ecf.vs.ecod[which.max(scov)]  
# #     best.match.qcov <- this.ecf.vs.ecod.max.scov$qcov
# #     best.match.scov <- this.ecf.vs.ecod.max.scov$scov
# #     best.match.ecod <- this.ecf.vs.ecod.max.scov$sname
# #   }
# #   data.table(ecf = this.ecf, family = this.family, best.match.qcov, best.match.scov, best.match.ecod)  
# # })
# # ecf.best.match.parms <- rbindlist(ecf.best.match.parms.list)
# # ecf.best.match.parms$id <- ecod.domains[.(ecf.best.match.parms$best.match.ecod)]$t_id
# # ecf.best.match.parms$name <- ecod.domains[.(ecf.best.match.parms$best.match.ecod)]$t_name
# # setorder(ecf.best.match.parms, family)
# # 
# # 
# # # Plot best matches with the corresponding scov/qcov values
# # ecf.best.match.parms[, shape := ifelse(best.match.scov >= 0.8, 22, 21)]  # Square = 22, Circle = 21
# # ecf.best.match.parms[, ecf := factor(ecf, levels = unique(ecf))]
# # 
# # # Generate distinct colors for all unique ECOD IDs
# # set.seed(42)  # Ensures consistent color assignment
# # unique_ecod_ids <- unique(na.omit(ecf.best.match.parms$id))  # Remove NA IDs
# # ecod_colors <- setNames(colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(unique_ecod_ids)), unique_ecod_ids)
# # 
# # # Plot
# # ggplot(ecf.best.match.parms, aes(x = ecf, y = best.match.qcov, shape = factor(shape), fill = id)) +
# #   geom_point(size = 4, alpha = 0.6, na.rm = FALSE) + 
# #   scale_shape_manual(values = c("21" = 21, "22" = 22), labels = c("S < 0.8", "S ≥ 0.8")) +
# #   scale_fill_manual(values = ecod_colors, na.translate = FALSE) +  # Assign colors to IDs
# #   scale_y_continuous(limits = c(0,1)) +
# #   theme_minimal() +
# #   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
# #   labs(title = "Best Q-Coverage Matches for Each ECF",
# #        x = "ECF",
# #        y = "Best Match Q-Coverage",
# #        shape = "S-Coverage Threshold",
# #        fill = "ECOD ID") +
# #   guides(fill = guide_legend(override.aes = list(size = 4)))  # Adjust legend size for readability
# # 
