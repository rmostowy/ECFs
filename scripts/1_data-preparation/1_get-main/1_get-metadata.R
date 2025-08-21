### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(readxl)))
# suppressWarnings(suppressPackageStartupMessages(library(rentrez)))
# suppressWarnings(suppressPackageStartupMessages(library(rtracklayer)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

# Load metadata table
metadata.table <- fread(metadata.table.filename)

# Load ICTV
ictv.release <- as.data.table(read_excel(ictv.release.filename, 2))
setnames(ictv.release, old = "Virus REFSEQ accession", new = "refseq")
setkey(ictv.release, refseq)

ictv.release.bac <- ictv.release[ictv.release$`Host source` == "bacteria"]

metadata2ictv <- ictv.release[.(metadata.table$Accession)]
metadata.table.ictv <- metadata2ictv %>% select(realm.ictv = "Realm", kingdom.ictv = "Kingdom", 
                         phylum.ictv = "Phylum", class.ictv = "Class", 
                         order.ictv  = "Order", family.ictv = "Family",
                         genus.ictv = "Genus", species.ictv = "Species", 
                         genome.composition.ictv = "Genome composition")


metadata.table.filtered <- metadata.table %>% select(accession = "Accession", description = "Description", 
                          genome.length = "Genome Length (bp)", 
                          is.jumbo = "Jumbophage",
                          mol.gc = "molGC (%)", 
                          molecule = "Molecule", 
                          no.cds = "Number CDS", 
                          perc.pos.strand = "Positive Strand (%)",
                          perc.neg.strand = "Negative Strand (%)",
                          coding.capacity = "Coding Capacity (%)",
                          no.trnas = "tRNAs",
                          genbank.division = "Genbank Division",
                          # bacphlip.score = "bacphlip_virulent_score",
                          family.ncbi = "Family", 
                          genus.ncbi = "Genus", 
                          host = "Host")

metadata.table.final <- cbind(metadata.table.filtered, metadata.table.ictv)
fwrite(metadata.table.final, file = metadata.table.output)
