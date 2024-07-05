### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

# this.output.dir <- 

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

lw(all.reprseqs.with.this.annotation.has.ecf)/length(all.reprseqs.with.this.annotation)



# 5. Plot the domain coverage per frequency




