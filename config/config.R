root <- rprojroot::find_rstudio_root_file()

#################
## Paths
#################

main.path <- "/Users/rmostowy"
main.data.dir <- file.path(main.path, "Data")
main.data.tmp.dir <- file.path(main.data.dir, "tmp-ecfs")


#################
## Labels
#################

phrogs.class <- c("integration and excision", "head and packaging", 
                  "transcription regulation", "connector", 
                  "tail", "lysis", 
                  "DNA, RNA and nucleotide metabolism", "other",
                  "moron, auxiliary metabolic gene and host takeover")
phrogs.col <- c("#FF7F00", "#C7E9C0", "#984EA3", "#41AB5D", "#006D2C", "#E41A1C", "#377EB8", "#999999", "#F781BF")

phrogs.class.order <- c("head and packaging", "connector", "tail", 
                        "DNA, RNA and nucleotide metabolism", "integration and excision",
                        "moron, auxiliary metabolic gene and host takeover", "transcription regulation",
                        "lysis", "other")

phrog.class.labels <- c(
  'integration and excision' = "integration\nand excision",
  'head and packaging' = "head and\npackaging",
  "connector" = "connector",
  "tail" = "tail",
  "lysis" = "lysis",
  'transcription regulation' = "transcription\nregulation",
  'DNA, RNA and nucleotide metabolism' = "DNA, RNA and\nnucleotide\nmetabolism",
  "other" = "other",
  'moron, auxiliary metabolic gene and host takeover' = "moron, auxiliary\nmetabolic gene\nand host takeover"
)
