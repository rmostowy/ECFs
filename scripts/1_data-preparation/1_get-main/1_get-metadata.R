### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(readxl)))
suppressWarnings(suppressPackageStartupMessages(library(rentrez)))
suppressWarnings(suppressPackageStartupMessages(library(rtracklayer)))

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
                          bacphlip.score = "bacphlip_virulent_score",
                          family.ncbi = "Family_ncbi_2023", 
                          genus.ncbi = "Genus_ncbi_2023", 
                          host = "Host")

metadata.table.final <- cbind(metadata.table.filtered, metadata.table.ictv)

fwrite(metadata.table.final, file = metadata.table.output)

# Extract all accessions in which reprseqs were found 
# (some of these are missing from the metadata table)
reprseq.table <- fread(reprseq.table.filename)
reprseq.table[, accession := gsub("\\.(\\d+):.*", "", ncbi.id)]
all.accessions <- unique(reprseq.table$accession)


# Define the batch size and delay
batch_size <- 300
delay <- 1  # 1 second between batch requests

# Split accessions into batches
accession_batches <- split(all.accessions, ceiling(seq_along(all.accessions) / batch_size))

# Initialise an empty list to store results
results <- list()

# Loop over each batch
for (i in seq_along(accession_batches)) {
  cat("Processing batch", i, "of", length(accession_batches), "...\n")
  
  # Fetch GFF3 annotations for the current batch
  batch <- accession_batches[[i]]
  
  # Fetch GFF3 annotations in text format
  gff3_text <- tryCatch(
    entrez_fetch(db = "nuccore", id = batch, rettype = "gff3", retmode = "text"),
    error = function(e) {
      message("Error fetching batch ", i, ": ", e$message)
      return(NULL)
    }
  )
  
  if (!is.null(gff3_text)) {
    # Write GFF3 text to a temporary file
    temp_gff3 <- tempfile(fileext = ".gff3")
    writeLines(gff3_text, temp_gff3)
    
    # Import the GFF3 file
    gff3 <- tryCatch(
      import(temp_gff3),
      error = function(e) {
        message("Error importing GFF3 for batch ", i, ": ", e$message)
        return(NULL)
      }
    )
    
    # Filter for CDS features and convert to data.table
    if (!is.null(gff3)) {
      cds <- gff3[gff3$type == "CDS"]
      cds_dt <- as.data.table(cds)
      cds_dt[, accession := gsub("\\.\\d+$", "", seqnames)]  # Remove version number from accession
      
      # Append to results list
      results[[i]] <- cds_dt[, .(accession, start, end, width, strand)]
    }
  }
  
  # Delay between requests
  Sys.sleep(delay)
}

accession.table <- rbindlist(results, use.names = TRUE, fill = TRUE)
fwrite(accession.table, file = accessions.cds.output)

