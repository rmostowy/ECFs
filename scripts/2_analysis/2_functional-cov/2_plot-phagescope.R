### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
# suppressWarnings(suppressPackageStartupMessages(library(Biostrings)))
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

# Define output directories
this.results.dir <- file.path(visual.dir, "part2_funct-cov", "analysis2-phagescope")
make.dir(this.results.dir)

# Load data
phagescope.metadata <- fread(phagescope.metadata.output)

ps.vs.ecf <- fread(phagescope.vs.ecf)
ps.vs.ecod <- fread(phagescope.vs.ecod)
ps.vs.pfam <- fread(phagescope.vs.pfam)

##########################

phagescope.metadata[, ecf.hit := MGG_Protein_ID %in% ps.vs.ecf$qname]
phagescope.metadata[, ecod.hit := MGG_Protein_ID %in% ps.vs.ecod$qname]
phagescope.metadata[, pfam.hit := MGG_Protein_ID %in% ps.vs.pfam$qname]

protein.vs.domains <- phagescope.metadata[, c("MGG_Protein_ID", "ecf.hit", "ecod.hit", "pfam.hit")]

protein.vs.category <- phagescope.metadata[
  , .(category = unlist(strsplit(Protein_classification, ";", fixed = TRUE))),
  by = .(MGG_Protein_ID)
]





phagescope.metadata$Protein_classification






# 
# # Load all datasets
# ecf.table <- fread(ecf.table.output)
# all.ecfs <- unique(ecf.table$ecf.name)
# Ne <- length(all.ecfs)
# 
# reprseq.ecfs <- fread(reprseq.ecfs.output)
# reprseq.ecfs <- reprseq.ecfs %>% filter(qcov >= MIN.QCOV.FOR.DOMAIN, scov >= MIN.SCOV.FOR.DOMAIN)
# 
# reprseq.domains.ecod <- fread(reprseq.domains.ecod.output)
# reprseq.ecod <- reprseq.domains.ecod %>% filter(qcov >= MIN.QCOV.FOR.DOMAIN, scov >= MIN.SCOV.FOR.DOMAIN)
# 
# reprseq.domains.pfam <- fread(reprseq.domains.pfam.output)
# reprseq.pfam <- reprseq.domains.pfam %>% filter(qcov >= MIN.QCOV.FOR.DOMAIN, scov >= MIN.SCOV.FOR.DOMAIN)
# 
# ecf.family.table <- fread(ecf.family.table.output)
# 
# reprseq.with.ecf.hits <- unique(reprseq.ecfs$qname)
# # all.ecfs <- unique(reprseq.ecfs$sname)
# 
# reprseq.functions.all <- fread(reprseq.functions.all.output)
# # reprseq.functions.all <- reprseq.functions.all[reprseq.functions.all$qname %in% all.reprseq]
# setkey(reprseq.functions.all, annotation)
# 
# top.annotations <- reprseq.functions.all[reprseq.functions.all$top.annotation] %>% select(funct, class, annotation) %>% distinct()
# 
# # ecf.family.table <- ecf.family.table[ecf.family.table$ecf %in% all.ecfs]
# no.families <- length(unique(ecf.family.table$family))
# 
# reprseq.table <- fread(reprseq.table.output)
# # reprseq.table[, accession := tstrsplit(cluster.name, ":")[[1]]]
# # reprseq.table[, accession := tstrsplit(accession, ".", fixed = T)[[1]]]
# 
# metadata.table <- fread(metadata.table.output)
# setkey(metadata.table, accession)
# 
# ecod.domains <- fread(ecod.domains.output)
# pfam.domains <- fread(pfam.domains.output)
# 
# ###############################################################
# ### 4. Plot the domain coverage (ECFs vs Pfam/ECODs) per function:
# ###    a) number of proteins with any hit
# ###    b) proportion of aa space with a domain hit
# ###############################################################
# 
# 
# 
# funct.cov.list <- lapply(1:nrow(top.annotations), function(annotation.index){
#   this.annotation <- top.annotations$annotation[annotation.index]
#   this.funct <- top.annotations$funct[annotation.index]
#   this.class <- top.annotations$class[annotation.index]
#   this.reprseq.functions.all <- reprseq.functions.all[.(this.annotation)]
#   all.reprseqs.with.this.annotation <- unique(this.reprseq.functions.all$qname)
#   all.reprseqs.with.this.annotation.has.ecf <- all.reprseqs.with.this.annotation %in% reprseq.with.ecf.hits
#   all.reprseqs.with.this.annotation.has.ecod <- all.reprseqs.with.this.annotation %in% reprseq.ecod$qname
#   all.reprseqs.with.this.annotation.has.pfam <- all.reprseqs.with.this.annotation %in% reprseq.pfam$qname
#   
#   data.table(funct = this.funct, class = this.class, annotation = this.annotation,
#              no.reprseq = length(all.reprseqs.with.this.annotation),
#              no.reprseq.with.ecf = lw(all.reprseqs.with.this.annotation.has.ecf),
#              no.reprseq.with.ecod = lw(all.reprseqs.with.this.annotation.has.ecod),
#              no.reprseq.with.pfam = lw(all.reprseqs.with.this.annotation.has.pfam)
#   )  
# })
# 
# funct.cov <- rbindlist(funct.cov.list)
# 
# setDT(funct.cov)
# 
# # Calculate domain coverage proportions using set to avoid copying
# funct.cov[, `:=`(ecf_coverage = no.reprseq.with.ecf / no.reprseq,
#                  ecod_coverage = no.reprseq.with.ecod / no.reprseq,
#                  pfam_coverage = no.reprseq.with.pfam / no.reprseq)]
# 
# funct.cov <- funct.cov %>% filter(no.reprseq >= MIN.NUM.PROT.PER.ANNOT)
# 
# 
# 
# # Melt the data for ggplot
# melted_data <- melt(funct.cov, id.vars = c("funct", "class"), 
#                     measure.vars = c("ecf_coverage", "ecod_coverage", "pfam_coverage"),
#                     variable.name = "database", value.name = "coverage")
# setDT(melted_data)
# # Define the colours for the domain coverages
# coverage_colours <- c("ecf_coverage" = "#08306B", "ecod_coverage" = "#006D2C", "pfam_coverage" = "#800026")
# # Set the order of the facets
# melted_data$class <- factor(melted_data$class, levels = phrogs.class.order)
# # Count the number of functions per class
# n_funct <- melted_data[, .N, by = class]
# # Define strip backgrounds using ggh4x
# # strip <- strip_themed(
# #   background_y = elem_list_rect(fill = phrogs.col[match(levels(melted_data$class), phrogs.class)])
# # )
# strip <- strip_themed(
#   background_x = elem_list_rect(fill = phrogs.col[match(levels(melted_data$class), phrogs.class)])
# )
# 
# # Create the ggplot object with customized strips and varying plot heights
# p.domain.cov <- ggplot(melted_data) +
#   geom_col(aes(x = funct, y = coverage, fill = database), 
#            position = "dodge", width = 0.7, linewidth = 0.2) +
#   scale_fill_manual(values = coverage_colours) +
#   theme(strip.text.y = element_text(face = "bold", size = 7, angle = 360)) +
#   labs(x = "Function", y = "Domain Coverage") +  
#   facet_grid2(cols = vars(class), scales = "free_x",  
#               labeller = labeller(class = phrog.class.labels),
#               strip = strip, axes = "x", space = "free") +  # Apply the `strip` object
#   guides(color = "none") + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),  
#         strip.placement = "outside",
#         panel.spacing = unit(0.1, "lines"))
# 
# # Draw the plot
# p.domain.cov.fn <- file.path(this.results.dir, "plot1_domain-cov.pdf")
# ggsave(file = p.domain.cov.fn, p.domain.cov, width = 15, height = 5)
# 
# # Calculate and plot domain coverage in reprseq of known vs. unknown function
# all.reprseq <- unique(reprseq.table$reprseq)
# has.assigned.function <- all.reprseq %in% reprseq.functions.all$qname
# has.ecf.hit <- all.reprseq %in% reprseq.with.ecf.hits
# has.ecod.hit <- all.reprseq %in% reprseq.ecod$qname
# has.pfam.hit <- all.reprseq %in% reprseq.pfam$qname
# reprseq.domain.coverage <- data.table(reprseq = all.reprseq, has.assigned.function, has.ecf.hit, has.ecod.hit, has.pfam.hit)
# 
# # Transforming the data to a long format
# long_data <- melt(reprseq.domain.coverage, id.vars = "has.assigned.function", 
#                   measure.vars = c("has.ecf.hit", "has.ecod.hit", "has.pfam.hit"),
#                   variable.name = "domain", value.name = "hit")
# 
# setDT(long_data)
# 
# # Calculate the proportions
# proportions <- long_data[, .(proportion = mean(hit, na.rm = TRUE)), by = .(has.assigned.function, domain)]
# 
# # Plot
# p.known.unknown <- ggplot(proportions, aes(x = domain, y = proportion, fill = domain)) +
#   geom_col(position = position_dodge()) +
#   facet_wrap(~ has.assigned.function, scales = "free_y", labeller = labeller(has.assigned.function = c(`TRUE` = "With Function", `FALSE` = "Without Function"))) +
#   labs(title = "Domain Coverage in Proteins", 
#        x = "Domain", 
#        y = "Proportion of Proteins with Hit",
#        fill = "Domain") +
#   scale_y_continuous(limits = c(0,1))
#   
# p.known.unknown.fn <- file.path(this.results.dir, "plot2_domain-cov-unknown.pdf")
# ggsave(p.known.unknown.fn, p.known.unknown, width = 8, height = 5)
