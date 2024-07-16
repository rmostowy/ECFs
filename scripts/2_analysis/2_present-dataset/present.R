### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(reshape2)))
suppressWarnings(suppressPackageStartupMessages(library(ggh4x)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

# Load all datasets
reprseq.ecfs <- fread(reprseq.ecfs.output)
reprseq.ecfs <- reprseq.ecfs %>% filter(qcov >= MIN.QCOV.FOR.DOMAIN, scov >= MIN.SCOV.FOR.DOMAIN)

reprseq.domains.ecod <- fread(reprseq.domains.ecod.output)
reprseq.ecod <- reprseq.domains.ecod %>% filter(qcov >= MIN.QCOV.FOR.DOMAIN, scov >= MIN.SCOV.FOR.DOMAIN)

reprseq.domains.pfam <- fread(reprseq.domains.pfam.output)
reprseq.pfam <- reprseq.domains.pfam %>% filter(qcov >= MIN.QCOV.FOR.DOMAIN, scov >= MIN.SCOV.FOR.DOMAIN)

ecf.family.table <- fread(ecf.family.table.output)

all.reprseq <- unique(reprseq.ecfs$qname)
all.ecfs <- unique(reprseq.ecfs$sname)

reprseq.functions.all <- fread(reprseq.functions.all.output)
# reprseq.functions.all <- reprseq.functions.all[reprseq.functions.all$qname %in% all.reprseq]
setkey(reprseq.functions.all, annotation)

top.annotations <- reprseq.functions.all[reprseq.functions.all$top.annotation] %>% select(funct, class, annotation) %>% distinct()

ecf.family.table <- ecf.family.table[ecf.family.table$ecf %in% all.ecfs]
no.families <- length(unique(ecf.family.table$ecf.family))

# 1. Plot the distribution of the number of reprseqs per ECF
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


# 2. Plot the distribution of the number of ECFs per EF
ecfs.per.family <- ecf.family.table[, .(no.ecfs = uniqueN(ecf)), by = ecf.family]

p.ecfs.per.fam.logx <- ggplot(ecfs.per.family, aes(x = no.ecfs)) +
  geom_histogram(bins = 15, fill = "darkblue", color = "black", alpha = 0.7) + scale_x_log10() + scale_y_log10()

p.ecfs.per.fam.logx.fn <- file.path(this.results.dir, "plot2a_ecfs-per-fam.png")
ggsave(p.ecfs.per.fam.logx.fn, p.ecfs.per.fam.logx, width = 6, height = 5)


# 3. Plot the distribution of ECF lengths (as well as for ECODs and Pfams)

ecf.lengths <- reprseq.ecfs %>% select(ecf = sname, length = slength) %>% distinct()

p.ecf.lengths <- ggplot(ecf.lengths, aes(x = length)) +
  geom_histogram(bins = 30, fill = "darkred", color = "black", alpha = 0.7) + scale_x_log10()
p.ecf.lengths.fn <- file.path(this.results.dir, "plot3a_ecf-lengths.png")
ggsave(p.ecf.lengths.fn, p.ecf.lengths, width = 6, height = 5)


# 4. Plot the domain coverage (ECFs vs Pfam/ECODs) per function:
#    a) number of proteins with any hit
#    b) proportion of aa space with a domain hit

this.results.dir <- "analysis2_funct-cov"
make.dir(this.results.dir)

funct.cov.list <- lapply(1:nrow(top.annotations), function(annotation.index){
  this.annotation <- top.annotations$annotation[annotation.index]
  this.funct <- top.annotations$funct[annotation.index]
  this.class <- top.annotations$class[annotation.index]
  this.reprseq.functions.all <- reprseq.functions.all[.(this.annotation)]
  all.reprseqs.with.this.annotation <- unique(this.reprseq.functions.all$qname)
  all.reprseqs.with.this.annotation.has.ecf <- all.reprseqs.with.this.annotation %in% all.reprseq
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

# 5. Plot the domain coverage per frequency






















