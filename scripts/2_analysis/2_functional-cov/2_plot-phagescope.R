### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(ComplexUpset)))
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

protein.vs.category[category %in% c("unsorted", "hypothetical")]$category <- "unknown"
unique.categories <- unique(protein.vs.category$category)


################################################
# 1. Plot the PhageScope protein dataset 
# coverage by the three domain databases
################################################

my_themes <- list(
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

df_upset <- protein.vs.domains[, .(
  ECF  = as.integer(as.logical(ecf.hit)),
  ECOD = as.integer(as.logical(ecod.hit)),
  Pfam = as.integer(as.logical(pfam.hit))
)] |> as.data.frame()


this.plot.filename <- file.path(this.results.dir, "plot1_relcov-upset.pdf")

domain.upset.plot <- upset(
  df_upset,
  intersect = colnames(df_upset),
  themes = my_themes
)

ggsave(plot = domain.upset.plot, filename = this.plot.filename, width = 10, height = 4.5)


################################################
# 2. Plot relative coverage per functional category
################################################

functional.cov.domains.table.list <- lapply(1:length(unique.categories), function(index){
  this.category <- unique.categories[index]
  this.category.proteins <- protein.vs.category[category %chin% this.category]$MGG_Protein_ID
  this.category.proteins.no <- length(this.category.proteins)
  this.category.proteins.vs.domains <- protein.vs.domains[MGG_Protein_ID %chin% this.category.proteins]
  this.category.ecf.coverage <- length(which(this.category.proteins.vs.domains$ecf.hit))
  this.category.ecod.coverage <- length(which(this.category.proteins.vs.domains$ecod.hit))
  this.category.pfam.coverage <- length(which(this.category.proteins.vs.domains$pfam.hit))
  this.category.out <- data.table(category = this.category,
                              ecf.coverage = this.category.ecf.coverage/this.category.proteins.no,
                              ecod.coverage = this.category.ecod.coverage/this.category.proteins.no,
                              pfam.coverage = this.category.pfam.coverage/this.category.proteins.no
  )
  return(this.category.out)
})
functional.cov.domains.table <- rbindlist(functional.cov.domains.table.list)

dt_long <- melt(
  functional.cov.domains.table,
  id.vars = "category",
  variable.name = "database",
  value.name = "coverage"
)

coverage_colours <- c(
  ecf.coverage  = "#08306B",
  ecod.coverage = "#006D2C",
  pfam.coverage = "#800026"
)


phagescope.cov.outname <- file.path(this.results.dir, "plot2_funct-cov.pdf")
phagescope.cov <- ggplot(dt_long, aes(x = category, y = coverage, fill = database)) +
  geom_col(
    position = position_dodge(width = 0.7),
    width = 0.6,
    linewidth = 0.2
  ) +
  scale_fill_manual(values = coverage_colours) +
  scale_y_continuous(limits = c(0,1)) +
  labs(
    x = "Functional category",
    y = "Domain coverage",
    fill = "Database"
  ) +
  # theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.placement = "outside",
    panel.spacing = unit(0.1, "lines")
  )

ggsave(plot = phagescope.cov, filename = phagescope.cov.outname, width = 7, height = 4.5)



