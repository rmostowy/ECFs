### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(reshape2)))
suppressWarnings(suppressPackageStartupMessages(library(ggh4x)))
suppressWarnings(suppressPackageStartupMessages(library(tidyr)))
suppressWarnings(suppressPackageStartupMessages(library(seqinr)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

# Load all datasets
ecfs.vs.qecfs <- fread(ecf.vs.qecf.output)
all.ecfs <- unique(ecfs.vs.qecfs$ecf)
Ne <- length(all.ecfs)

reprseq.ecfs <- fread(reprseq.ecfs.output)
reprseq.ecfs <- reprseq.ecfs %>% filter(qcov >= MIN.QCOV.FOR.DOMAIN, scov >= MIN.SCOV.FOR.DOMAIN)

reprseq.domains.ecod <- fread(reprseq.domains.ecod.output)
reprseq.ecod <- reprseq.domains.ecod %>% filter(qcov >= MIN.QCOV.FOR.DOMAIN, scov >= MIN.SCOV.FOR.DOMAIN)

reprseq.domains.pfam <- fread(reprseq.domains.pfam.output)
reprseq.pfam <- reprseq.domains.pfam %>% filter(qcov >= MIN.QCOV.FOR.DOMAIN, scov >= MIN.SCOV.FOR.DOMAIN)

ecf.family.table <- fread(ecf.family.table.output)

reprseq.with.ecf.hits <- unique(reprseq.ecfs$qname)
# all.ecfs <- unique(reprseq.ecfs$sname)

reprseq.functions.all <- fread(reprseq.functions.all.output)
# reprseq.functions.all <- reprseq.functions.all[reprseq.functions.all$qname %in% all.reprseq]
setkey(reprseq.functions.all, annotation)

top.annotations <- reprseq.functions.all[reprseq.functions.all$top.annotation] %>% select(funct, class, annotation) %>% distinct()

# ecf.family.table <- ecf.family.table[ecf.family.table$ecf %in% all.ecfs]
no.families <- length(unique(ecf.family.table$ecf.family))

reprseq.table <- fread(reprseq.table.filename)
reprseq.table[, accession := tstrsplit(cluster.name, ":")[[1]]]
reprseq.table[, accession := tstrsplit(accession, ".", fixed = T)[[1]]]

metadata.table <- fread(metadata.table.output)
setkey(metadata.table, accession)

ecod.domains <- fread(ecod.domains.output)
pfam.domains <- fread(pfam.domains.output)


###############################################################
### 1. Plot the distribution of the number of reprseqs per ECF
###############################################################
this.results.dir <- "analysis1_distr"
make.dir(this.results.dir)

reprseq.per.ecf <- reprseq.ecfs[, .(no.reprseqs = uniqueN(qname)), by = sname]

p.repr.per.ecf <- ggplot(reprseq.per.ecf, aes(x = no.reprseqs)) +
  geom_histogram(bins = 100, fill = "blue", color = "black", alpha = 0.7) + scale_x_continuous(limits = c(1,200))

p.repr.per.ecf.logx <- ggplot(reprseq.per.ecf, aes(x = no.reprseqs)) +
  geom_histogram(bins = 20, fill = "blue", color = "black", alpha = 0.7) + scale_x_log10()

p.repr.per.ecf.fn <- file.path(this.results.dir, "plot1a_repr-per-ecf.png")
p.repr.per.ecf.logx.fn <- file.path(this.results.dir, "plot1b_repr-per-ecf-log.png")
ggsave(p.repr.per.ecf.fn, p.repr.per.ecf, width = 6, height = 5)
ggsave(p.repr.per.ecf.logx.fn, p.repr.per.ecf.logx, width = 5, height = 5)

###############################################################
### 2. Plot the distribution of the number of ECFs per EF
###############################################################

ecfs.per.family <- ecf.family.table[, .(no.ecfs = uniqueN(ecf)), by = ecf.family]

p.ecfs.per.fam.logx <- ggplot(ecfs.per.family, aes(x = no.ecfs)) +
  geom_histogram(bins = 15, fill = "darkblue", color = "black", alpha = 0.7) + scale_x_log10() + scale_y_log10()

p.ecfs.per.fam.logx.fn <- file.path(this.results.dir, "plot2a_ecfs-per-fam.png")
ggsave(p.ecfs.per.fam.logx.fn, p.ecfs.per.fam.logx, width = 6, height = 5)

###############################################################
### 3. Plot the distribution of ECF lengths (as well as for ECODs and Pfams)
###############################################################

ecf.lengths <- reprseq.ecfs %>% select(ecf = sname, length = slength) %>% distinct()

p.ecf.lengths <- ggplot(ecf.lengths, aes(x = length)) +
  geom_histogram(bins = 30, fill = "darkred", color = "black", alpha = 0.7) + scale_x_log10()
p.ecf.lengths.fn <- file.path(this.results.dir, "plot3a_ecf-lengths.png")
ggsave(p.ecf.lengths.fn, p.ecf.lengths, width = 6, height = 5)

ecod.domains.seqs <- read.fasta(ecod.domains.fasta.filename)
ecod.domains.lengths <- sapply(ecod.domains.seqs, length)
# hist(ecod.domains.lengths, breaks = 30)
ecod.lengths <- data.table(ecod = names(ecod.domains.seqs), length = ecod.domains.lengths)

p.ecod.lengths <- ggplot(ecod.lengths, aes(x = length)) +
  geom_histogram(bins = 30, fill = "darkred", color = "black", alpha = 0.7) + scale_x_log10()
p.ecod.lengths.fn <- file.path(this.results.dir, "plot3b_ecod-lengths.png")
ggsave(p.ecod.lengths.fn, p.ecod.lengths, width = 6, height = 5)

p.pfam.lengths <- ggplot(pfam.domains, aes(x = length)) +
  geom_histogram(bins = 30, fill = "darkred", color = "black", alpha = 0.7) + scale_x_log10()
p.pfam.lengths.fn <- file.path(this.results.dir, "plot3c_pfam-lengths.png")
ggsave(p.pfam.lengths.fn, p.pfam.lengths, width = 6, height = 5)

###############################################################
### 4. Plot the domain coverage (ECFs vs Pfam/ECODs) per function:
###    a) number of proteins with any hit
###    b) proportion of aa space with a domain hit
###############################################################

this.results.dir <- "analysis2_funct-cov"
make.dir(this.results.dir)

funct.cov.list <- lapply(1:nrow(top.annotations), function(annotation.index){
  this.annotation <- top.annotations$annotation[annotation.index]
  this.funct <- top.annotations$funct[annotation.index]
  this.class <- top.annotations$class[annotation.index]
  this.reprseq.functions.all <- reprseq.functions.all[.(this.annotation)]
  all.reprseqs.with.this.annotation <- unique(this.reprseq.functions.all$qname)
  all.reprseqs.with.this.annotation.has.ecf <- all.reprseqs.with.this.annotation %in% reprseq.with.ecf.hits
  all.reprseqs.with.this.annotation.has.ecod <- all.reprseqs.with.this.annotation %in% reprseq.ecod$qname
  all.reprseqs.with.this.annotation.has.pfam <- all.reprseqs.with.this.annotation %in% reprseq.pfam$qname
  
  data.table(funct = this.funct, class = this.class, annotation = this.annotation,
             no.reprseq = length(all.reprseqs.with.this.annotation),
             no.reprseq.with.ecf = lw(all.reprseqs.with.this.annotation.has.ecf),
             no.reprseq.with.ecod = lw(all.reprseqs.with.this.annotation.has.ecod),
             no.reprseq.with.pfam = lw(all.reprseqs.with.this.annotation.has.pfam)
  )  
})

funct.cov <- rbindlist(funct.cov.list)

setDT(funct.cov)

# Calculate domain coverage proportions using set to avoid copying
funct.cov[, `:=`(ecf_coverage = no.reprseq.with.ecf / no.reprseq,
                 ecod_coverage = no.reprseq.with.ecod / no.reprseq,
                 pfam_coverage = no.reprseq.with.pfam / no.reprseq)]

funct.cov <- funct.cov %>% filter(no.reprseq >= MIN.NUM.PROT.PER.ANNOT)



# Melt the data for ggplot
melted_data <- melt(funct.cov, id.vars = c("funct", "class"), 
                    measure.vars = c("ecf_coverage", "ecod_coverage", "pfam_coverage"),
                    variable.name = "database", value.name = "coverage")
setDT(melted_data)
# Define the colours for the domain coverages
coverage_colours <- c("ecf_coverage" = "#08306B", "ecod_coverage" = "#006D2C", "pfam_coverage" = "#800026")
# Set the order of the facets
melted_data$class <- factor(melted_data$class, levels = phrogs.class.order)
# Count the number of functions per class
n_funct <- melted_data[, .N, by = class]
# Define strip backgrounds using ggh4x
strip <- strip_themed(
  background_y = elem_list_rect(fill = phrogs.col[match(levels(melted_data$class), phrogs.class)])
)

# Create the ggplot object with customized strips and varying plot heights
p.domain.cov <- ggplot(melted_data) +
  geom_col(aes(y = funct, x = coverage, fill = database), 
           position = "dodge", width = 0.7, size = 0.2) +
  scale_fill_manual(values = coverage_colours) +
  # theme(strip.text = element_text(face = "bold", size = 7, angle = 90)) + 
  theme(strip.text.y = element_text(face = "bold", size = 7, angle = 360)) +
  labs(y = "Function", x = "Domain Coverage") +
  facet_grid2(rows = vars(class), scales = "free_y", 
              labeller = labeller(class = phrog.class.labels),
              strip = strip, axes = "y", space = "free") +
  guides(color = "none") +
  theme(axis.text.y = element_text(angle = 0, hjust = 1),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.spacing = unit(0.1, "lines"))

# Draw the plot
p.domain.cov.fn <- file.path(this.results.dir, "plot1_domain-cov.pdf")
ggsave(file = p.domain.cov.fn, p.domain.cov, width = 7, height = 10)

# Calculate and plot domain coverage in reprseq of known vs. unknown function
all.reprseq <- unique(reprseq.table$reprseq)
has.assigned.function <- all.reprseq %in% reprseq.functions.all$qname
has.ecf.hit <- all.reprseq %in% reprseq.with.ecf.hits
has.ecod.hit <- all.reprseq %in% reprseq.ecod$qname
has.pfam.hit <- all.reprseq %in% reprseq.pfam$qname
reprseq.domain.coverage <- data.table(reprseq = all.reprseq, has.assigned.function, has.ecf.hit, has.ecod.hit, has.pfam.hit)

# Transforming the data to a long format
long_data <- melt(reprseq.domain.coverage, id.vars = "has.assigned.function", 
                  measure.vars = c("has.ecf.hit", "has.ecod.hit", "has.pfam.hit"),
                  variable.name = "domain", value.name = "hit")

setDT(long_data)

# Calculate the proportions
proportions <- long_data[, .(proportion = mean(hit, na.rm = TRUE)), by = .(has.assigned.function, domain)]

# Plot
p.known.unknown <- ggplot(proportions, aes(x = domain, y = proportion, fill = domain)) +
  geom_col(position = position_dodge()) +
  facet_wrap(~ has.assigned.function, scales = "free_y", labeller = labeller(has.assigned.function = c(`TRUE` = "With Function", `FALSE` = "Without Function"))) +
  labs(title = "Domain Coverage in Proteins", 
       x = "Domain", 
       y = "Proportion of Proteins with Hit",
       fill = "Domain") +
  scale_y_continuous(limits = c(0,1))
  
p.known.unknown.fn <- file.path(this.results.dir, "plot2_domain-cov-unknown.pdf")
ggsave(p.known.unknown.fn, p.known.unknown, width = 8, height = 5)

#######################################
### 5. Calculate reprseq frequency
#######################################

frequency.thresholds <- seq(0,1,0.2)
no.ranges <- length(frequency.thresholds)-1

# Calculate frequency based on reprseq cluster size
reprseq.cluster.size <- reprseq.table[, .N, by = .(reprseq, family)]
max.size <- max(reprseq.cluster.size$N)
freq.per.thr <- sapply(1:max.size, function(size.thr) sum(reprseq.cluster.size[reprseq.cluster.size$N==size.thr]$N)/sum(reprseq.cluster.size$N))
cum.freq.per.thr <- cumsum(freq.per.thr)

reprseq.per.freq.quintile.list <- lapply(1:no.ranges, function(range.index){
  this.range.lower <- frequency.thresholds[range.index]
  if(range.index == no.ranges){
    this.range.upper <- 1
  } else{
    this.range.upper <- frequency.thresholds[range.index+1]  
  }
  sizes.within.this.range <- which((cum.freq.per.thr > this.range.lower) & (cum.freq.per.thr <= this.range.upper))
  reprseqs.within.this.range <- reprseq.cluster.size[reprseq.cluster.size$N %in% sizes.within.this.range]$reprseq
  data.table(reprseq = reprseqs.within.this.range, frequency.quintile = range.index)  
})
reprseq.per.freq.cluster.quintile <- rbindlist(reprseq.per.freq.quintile.list)

# Calculate frequency based on taxonomic genus
genus.freq <- metadata.table[, .N, by = genus.ncbi]
genus.freq <- genus.freq[!genus.freq$genus.ncbi %in% c("Unclassified", "")]

ncbi.genus <- metadata.table[.(reprseq.table$accession)]$genus.ncbi
ncbi.genus[ncbi.genus == "Unclassified" | ncbi.genus == ""] <- NA
reprseq.table$genus <- ncbi.genus

max.genus.size <- max(genus.freq$N)
genus.freq.per.thr <- sapply(1:max.genus.size, function(size.thr) sum(genus.freq[genus.freq$N==size.thr]$N)/sum(genus.freq$N))
cum.genus.freq.per.thr <- cumsum(genus.freq.per.thr)

reprseq.per.freq.quintile.list <- lapply(1:no.ranges, function(range.index){
  this.range.lower <- frequency.thresholds[range.index]
  if(range.index == no.ranges){
    this.range.upper <- 1
  } else{
    this.range.upper <- frequency.thresholds[range.index+1]  
  }
  sizes.within.this.range <- which((cum.genus.freq.per.thr > this.range.lower) & (cum.genus.freq.per.thr <= this.range.upper))
  genera.within.this.range <- genus.freq[genus.freq$N %in% sizes.within.this.range]$genus.ncbi
  reprseq.within.this.range <- unique(reprseq.table[reprseq.table$genus %in% genera.within.this.range]$reprseq)
  data.table(reprseq = reprseq.within.this.range, frequency.quintile = range.index)
})
reprseq.per.freq.genus.quintile <- rbindlist(reprseq.per.freq.quintile.list)


##########################################
### 6. Visualise domain coverage per freq
##########################################

cov.per.freq.list <- lapply(1:no.ranges, function(this.quantile){# get domain coverage within this quantile for cluster frequency
  reprseq.per.freq.cluster.within <- reprseq.per.freq.cluster.quintile[reprseq.per.freq.cluster.quintile$frequency.quintile == this.quantile]$reprseq
  reprseq.per.freq.cluster.within.with.ecf <- reprseq.per.freq.cluster.within[reprseq.per.freq.cluster.within %in% reprseq.with.ecf.hits]
  reprseq.per.freq.cluster.within.with.ecod <- reprseq.per.freq.cluster.within[reprseq.per.freq.cluster.within %in% reprseq.ecod$qname]
  reprseq.per.freq.cluster.within.with.pfam <- reprseq.per.freq.cluster.within[reprseq.per.freq.cluster.within %in% reprseq.pfam$qname]
  
  out.cluster <- data.table(quantile = this.quantile, freq.type = "cluster",
                            no.reprseq = length(reprseq.per.freq.cluster.within),
                            no.reprseq.with.ecf = length(reprseq.per.freq.cluster.within.with.ecf),
                            no.reprseq.with.ecod = length(reprseq.per.freq.cluster.within.with.ecod),
                            no.reprseq.with.pfam = length(reprseq.per.freq.cluster.within.with.pfam)
  )
  
  # get domain coverage within this quantile for genus frequency
  reprseq.per.freq.genus.within <- reprseq.per.freq.genus.quintile[reprseq.per.freq.genus.quintile$frequency.quintile == this.quantile]$reprseq
  reprseq.per.freq.genus.within.with.ecf <- reprseq.per.freq.genus.within[reprseq.per.freq.genus.within %in% reprseq.with.ecf.hits]
  reprseq.per.freq.genus.within.with.ecod <- reprseq.per.freq.genus.within[reprseq.per.freq.genus.within %in% reprseq.ecod$qname]
  reprseq.per.freq.genus.within.with.pfam <- reprseq.per.freq.genus.within[reprseq.per.freq.genus.within %in% reprseq.pfam$qname]
  
  out.genus <- data.table(quantile = this.quantile, freq.type = "genus",
                          no.reprseq = length(reprseq.per.freq.genus.within),
                          no.reprseq.with.ecf = length(reprseq.per.freq.genus.within.with.ecf),
                          no.reprseq.with.ecod = length(reprseq.per.freq.genus.within.with.ecod),
                          no.reprseq.with.pfam = length(reprseq.per.freq.genus.within.with.pfam)
  )
  
  rbind(out.cluster, out.genus)
})
cov.per.freq <- rbindlist(cov.per.freq.list)

# Calculate proportions
cov.per.freq[, `:=`(ecf.cov = no.reprseq.with.ecf / no.reprseq,
                    ecod.cov = no.reprseq.with.ecod / no.reprseq,
                    pfam.cov = no.reprseq.with.pfam / no.reprseq)]


# Convert the data from wide to long format
long_data.abs <- pivot_longer(cov.per.freq, 
                          cols = c("no.reprseq", "no.reprseq.with.ecf", "no.reprseq.with.ecod", "no.reprseq.with.pfam"),
                          names_to = "category",
                          values_to = "count")

# Plot
p.abs <- ggplot(long_data.abs, aes(x = as.factor(quantile), y = count, fill = category)) +
  geom_col(position = position_dodge()) +
  facet_wrap(~freq.type, scales = "free_y") +  # free_y allows each facet to have its own y scale
  labs(title = "Coverage per Frequency Type across Quintiles",
       x = "Quintile",
       y = "Count",
       fill = "Category") +
  theme_minimal()

# Convert the data from wide to long format
long_data.prop <- pivot_longer(cov.per.freq, 
                          cols = c("ecf.cov", "ecod.cov", "pfam.cov"),
                          names_to = "category",
                          values_to = "proportion")

# Plot
p.prop <- ggplot(long_data.prop, aes(x = as.factor(quantile), y = proportion, fill = category)) +
  geom_col(position = position_dodge()) +
  facet_wrap(~freq.type, scales = "free_y") +  # Allows each facet to have its own y scale
  labs(title = "Proportional Coverage per Frequency Type across Quintiles",
       x = "Quintile",
       y = "Proportion",
       fill = "Coverage Type") +
  theme_minimal() + scale_y_continuous(limits = c(0,1))

p.abs.fn <- file.path(this.results.dir, "plot3a_domain-cov-freq-abs.pdf")
p.prop.fn <- file.path(this.results.dir, "plot3b_domain-cov-freq-rel.pdf")
ggsave(file = p.abs.fn, p.abs, width = 6, height = 4)
ggsave(file = p.prop.fn, p.prop, width = 6, height = 4)


  
