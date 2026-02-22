### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(ComplexUpset)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

# Define output directories
this.visual.dir <- file.path(visual.dir, "part3_novel-ecfs")
this.results.dir <- this.visual.dir
make.dir(this.results.dir)

# Load reprseq tables
reprseq.table <- fread(reprseq.table.output)
reprseq.table[, accession := sub(":.*$", "", ncbi.id)]
proteins.per.reprseq <- reprseq.table[, .N, by = reprseq]

all.reprseq <- unique(reprseq.table$reprseq)

# Load metadata table
metadata.table <- fread(metadata.table.output)

# Load ECOD domain table
ecod.domains <- fread(ecod.domains.output)
setkey(ecod.domains, full.id)

# Load PFAM domain table
pfam.domains <- fread(pfam.domains.output)
setkey(pfam.domains, accession)

# Load PHROG annotations
reprseq.functions.unique <- fread(reprseq.functions.unique.output)
setkey(reprseq.functions.unique, reprseq)

# Load ECF table
ecf.table <- fread(ecf.table.output)
all.ecfs <- unique(ecf.table$ecf.name)
ecf.lengths <- ecf.table$length
Ne <- length(all.ecfs)

# Load ECF vs REPRSEQ table
reprseq.ecfs <- fread(reprseq.ecfs.output)
reprseq.ecfs.conf <- reprseq.ecfs %>% filter(scov >= 0.95 & prob >= 0.95) 
ecf2reprseq <- reprseq.ecfs.conf %>% select(ecf = sname, reprseq = qname)
setkey(ecf2reprseq, ecf)

# Load ECF vs ECOD table
ecfs.vs.ecod <- fread(ecfs.vs.ecod.output)
ecfs.vs.ecod[, score := round(qcov*scov*prob/100*pident/100, d = 6)]

# Load ECF vs PFAM table
ecfs.vs.pfam <- fread(ecfs.vs.pfam.output)
ecfs.vs.pfam[, score := round(qcov*scov*prob/100*pident/100, d = 6)]

setkey(ecfs.vs.ecod, qname)
setkey(ecfs.vs.pfam, qname)


ecf.domain.table.list <- lapply(1:Ne, function(index){
  cat(index,"\n")
  this.ecf <- all.ecfs[index]
  this.ecf.length <- ecf.lengths[index]
  this.ecf.has.reprseqs <- this.ecf %in% ecf2reprseq$ecf
  
  this.ecf.funct <- "unknown"
  this.ecf.taxa <- "unknown"
  this.ecf.no.reprseq <- this.ecf.no.proteins <- this.ecf.no.genera <- 0
  
  # get linked functions, taxa, diversity
  if(this.ecf.has.reprseqs){
    # get linked functions
    this.ecf.reprseqs <- ecf2reprseq[.(this.ecf)]$reprseq
    this.ecf.linked.annotations <- reprseq.functions.unique[.(this.ecf.reprseqs)]
    this.ecf.linked.annotations.known <- this.ecf.linked.annotations[this.ecf.linked.annotations$funct != ""]
    
    if(nrow(this.ecf.linked.annotations.known)>0){
      setorder(this.ecf.linked.annotations.known, -score)
      this.ecf.funct <- this.ecf.linked.annotations.known$funct[1]
    }
    
    this.ecf.accessions <- unique(reprseq.table[reprseq %in% this.ecf.reprseqs]$accession)
    this.ecf.metadata <- metadata.table[accession %in% this.ecf.accessions]
    this.ecf.genera <- unique(this.ecf.metadata$genus.ncbi)
    this.ecf.genera.known <- this.ecf.genera[this.ecf.genera != "Unclassified"]
    if(any(this.ecf.genera != "Unclassified")) this.ecf.taxa <- "known"
    
    this.ecf.no.reprseq <- length(this.ecf.reprseqs)
    this.ecf.no.proteins <- nrow(reprseq.table[reprseq %in% this.ecf.reprseqs])
    this.ecf.no.genera <- length(unique(this.ecf.genera.known))
  }
  
  # get linked domain hits
  this.ecf.ecod.hit <- "no"
  if(this.ecf %in% ecfs.vs.ecod$qname) this.ecf.ecod.hit <- "yes"
  this.ecf.pfam.hit <- "no"
  if(this.ecf %in% ecfs.vs.pfam$qname) this.ecf.pfam.hit <- "yes"
  
  data.table(name = this.ecf, length = this.ecf.length, ecod.hit = this.ecf.ecod.hit, pfam.hit = this.ecf.pfam.hit, 
             funct = this.ecf.funct, taxa = this.ecf.taxa, 
             no.reprseq = this.ecf.no.reprseq, no.proteins = this.ecf.no.proteins, no.genera = this.ecf.no.genera)
})
ecf.domain.table <- rbindlist(ecf.domain.table.list)

fwrite(ecf.domain.table, file = "ecf-domain-table.txt")

#####################
###### EXPLORE ######
#####################

ecf.domain.table[, domain.hit := ecod.hit == "yes" | pfam.hit == "yes"]
ecf.domain.table[, funct.known := funct != "unknown"]

p.prop.domain.hits <- ggplot(ecf.domain.table, aes(x = domain.hit, fill = domain.hit)) +
  geom_bar() +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "lightgrey")) +
  theme_minimal() +
  labs(
    title = "Proportion of Domain Hits",
    x = "Domain Hit",
    y = "Count",
    fill = "Domain Hit"
  )
p.prop.domain.hits.outfile <- file.path(this.results.dir, "plot1-prop-domhits.pdf")
ggsave(p.prop.domain.hits.outfile, p.prop.domain.hits, width = 4, height = 4)

p.domain.hit.by.length <- ggplot(ecf.domain.table, aes(x = domain.hit, y = length, fill = domain.hit)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "lightgrey")) +
  theme_bw(base_size = 8) +
  labs(
    title = "Distribution of Length by Domain Hit Status",
    x = "Domain Hit",
    y = "Length"
  ) + scale_y_log10()
p.domain.hit.by.length.outfile <- file.path(this.results.dir, "plot2-domhit-by-length.pdf")
ggsave(p.domain.hit.by.length.outfile, p.domain.hit.by.length, width = 4, height = 4)

p.domain.hit.by.div <- ggplot(ecf.domain.table, aes(x = domain.hit, y = no.proteins, fill = domain.hit)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "lightgrey")) +
  theme_bw(base_size = 8) +
  labs(
    title = "No proteins by Domain Hit Status",
    x = "Domain Hit",
    y = "No proteins"
  )  + scale_y_log10()
p.domain.hit.by.div.outfile <- file.path(this.results.dir, "plot3a-domhit-by-div1.pdf")
ggsave(p.domain.hit.by.div.outfile, p.domain.hit.by.div, width = 4, height = 4)

p.domain.hit.by.div2 <- ggplot(ecf.domain.table, aes(x = domain.hit, y = no.reprseq, fill = domain.hit)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "lightgrey")) +
  theme_bw(base_size = 8) +
  labs(
    title = "No rHMMs by Domain Hit Status",
    x = "Domain Hit",
    y = "No reprseq"
  )  + scale_y_log10()
p.domain.hit.by.div2.outfile <- file.path(this.results.dir, "plot3b-domhit-by-div2.pdf")
ggsave(p.domain.hit.by.div2.outfile, p.domain.hit.by.div2, width = 4, height = 4)

p.domain.hit.by.div3 <- ggplot(ecf.domain.table, aes(x = domain.hit, y = no.genera, fill = domain.hit)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "lightgrey")) +
  theme_bw(base_size = 8) +
  labs(
    title = "No genera by Domain Hit Status",
    x = "Domain Hit",
    y = "No genera"
  )  + scale_y_log10()
p.domain.hit.by.div3.outfile <- file.path(this.results.dir, "plot3c-domhit-by-div3.pdf")
ggsave(p.domain.hit.by.div3.outfile, p.domain.hit.by.div3, width = 4, height = 4)

p.domain.hit.by.funct <- ggplot(ecf.domain.table, aes(x = funct.known, fill = domain.hit)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "lightgrey")) +
  theme_bw(base_size = 8) +
  labs(
    title = "Functional Category by Domain Hit Status",
    x = "Functional Category",
    y = "Proportion",
    fill = "Domain Hit"
  )
p.domain.hit.by.funct.outfile <- file.path(this.results.dir, "plot4-domhit-by-funct.pdf")
ggsave(p.domain.hit.by.funct.outfile, p.domain.hit.by.funct, width = 4, height = 4)

p.domain.hit.by.taxonomy <- ggplot(ecf.domain.table, aes(x = taxa, fill = domain.hit)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "lightgrey")) +
  theme_bw(base_size = 8) +
  labs(
    title = "Taxonomic assignmnet by Domain Hit Status",
    x = "Genus assignments",
    y = "Proportion",
    fill = "Domain Hit"
  )
p.domain.hit.by.funct.outfile <- file.path(this.results.dir, "plot5-domhit-by-taxonomy.pdf")
ggsave(p.domain.hit.by.fįunct.outfile, p.domain.hit.by.taxonomy, width = 4, height = 4)


# Prepare data for the UpSet plot
ecf.domain.table[, ECOD := ecod.hit == "yes"]
ecf.domain.table[, PFAM := pfam.hit == "yes"]

# Select relevant columns for the upset plot
upset_data <- ecf.domain.table[, .(ECOD, PFAM)]

# Create the UpSet plot
p.upset.domain <- upset(
  upset_data,
  intersect = c("ECOD", "PFAM"),
  mode = "exclusive_intersection",  # Correct mode
  themes = list(
    'overall_sizes' = theme(
      axis.title.x = element_text(),
      axis.title.y = element_text()
    ),
    'intersections_matrix' = theme(
      axis.title.x = element_text(),
      axis.title.y = element_text()
    ),
    'intersections_sizes' = theme(
      axis.title.x = element_text(),
      axis.title.y = element_text()
    )
  )
)
p.upset.domain.outfile <- file.path(this.results.dir, "plot6-upset-domain.pdf")
ggsave(p.upset.domain.outfile, p.upset.domain, width = 7, height = 4)
