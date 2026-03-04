### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(igraph)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(ggraph)))
suppressWarnings(suppressPackageStartupMessages(library(patchwork)))
suppressWarnings(suppressPackageStartupMessages(library(broom)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

#################################
### Set up output directory
#################################
this.part.dir <- file.path(ecf.proj.dir, "part2_containment")
make.dir(this.part.dir)

# figures.dir <- file.path(visual.dir, "figures")
# suppl.dir <- file.path(visual.dir, "figures")
# make.dir(figures.dir)

#################################
### Load all data:
### 1. ECF tables
### 2. ECF comparsion tables
#################################

# ECF tables
ecf.table <- fread(ecf.table.output)
ecf.family.table <- fread(ecf.family.table.output)

# ECF comparison tables
ecfs.vs.ecfs <- fread(ecfs.vs.ecfs.output) # raw
ecfs.vs.ecfs.cov <- fread(ecfs.vs.ecfs.cov.output) # processed

#################################
### Process loaded tables
#################################

all.ecfs <- unique(ecf.table$ecf.name)
setkey(ecf.family.table, ecf)
setkey(ecf.table, ecf.name)

ecfs.vs.ecfs[, max.cov := pmax(qcov, scov)]
setkey(ecfs.vs.ecfs, qname, sname)

ecfs.vs.ecfs.cov$qlength <- ecf.table[.(ecfs.vs.ecfs.cov$query)]$length
ecfs.vs.ecfs.cov$slength <- ecf.table[.(ecfs.vs.ecfs.cov$subject)]$length
ecfs.vs.ecfs.cov[, qcov := pmin(qcov.min, qcov.max)]
ecfs.vs.ecfs.cov[, scov := pmin(scov.min, scov.max)]
ecfs.vs.ecfs.cov[, cov.min := pmin(scov.min, qcov.min)]
ecfs.vs.ecfs.cov[, cov.max := pmax(scov.min, qcov.min)]
eve.conf <- ecfs.vs.ecfs.cov[prob >= CONTAINMENT.PROB.THR]

# generate network table
ecf.network.ab <- eve.conf %>% select(query, subject, prob, qcov, scov)
ecf.network.ba <- eve.conf %>% select(query = subject, subject = query, prob, qcov = scov, scov = qcov)
ecf.network <- rbind(ecf.network.ab, ecf.network.ba)

# Define base theme for plotting
base_theme <- theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0)
  )

#################################
### 1. Generate per-ECF 
### summary table
#################################

dt <- copy(ecf.network)
dt[, subject.ef := ecf.family.table[.(subject)]$family]
ecf.no.fams <- dt[, .(no.fam = uniqueN(subject.ef)), by = query]

ecf.connectivity.stats <- dt[, .(
  mean.qcov = mean(qcov),
  mean.scov = mean(scov),
  med.qcov = median(qcov),
  med.scov = median(scov),
  degree = .N
  ), by = query]
ecf.connectivity.stats[, family := ecf.family.table[.(ecf.connectivity.stats$query)]$family]
ecf.connectivity.stats[, length := ecf.table[.(ecf.connectivity.stats$query)]$length]
ecf.connectivity.stats[, asym_mean := mean.qcov - mean.scov]
ecf.connectivity.stats[, asym_med := med.qcov - med.scov]
ecf.connectivity.stats[, strength_mean := (mean.qcov + mean.scov)/2]
ecf.connectivity.stats[, strength_med := (med.qcov + med.scov)/2]

ecf.connectivity.stats <- merge(ecf.connectivity.stats, ecf.no.fams, by = "query")

#################################
### 2. Show qcov/scov heatmap
#################################

# heatmap_fill <- scale_fill_gradient(
#   low = "#F7FBFF",
#   high = "#2C7FB8",
#   name = "Count"
# )

heatmap_fill <- scale_fill_viridis_c(
  option = "cividis",
  trans = "sqrt",
  name = "Count"
)

p.cov.heatmap <- ggplot(eve.conf, aes(x = qcov, y = scov)) +
  geom_bin2d(bins = 60) +
  heatmap_fill +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    x = "Query coverage",
    y = "Subject coverage"
  ) +
  base_theme

p.cov.heatmap.outfile <- file.path(this.part.dir, "plot1_cov-heatmap.png")
ggsave(p.cov.heatmap.outfile, p.cov.heatmap, width = 6, height = 5)

#################################
### 3. Show qcov/scov heatmap
### for per-ECF averages
#################################

p.connect.cov.heatmap <- ggplot(
  ecf.connectivity.stats,
  aes(x = mean.qcov, y = mean.scov)
) +
  geom_bin2d(bins = 60) +
  heatmap_fill +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    x = "Mean query coverage",
    y = "Mean subject coverage"
  ) +
  base_theme

p.connect.cov.heatmap.outfile <- file.path(this.part.dir, "plot2_connect-cov-heatmap.png")
ggsave(p.connect.cov.heatmap.outfile, p.connect.cov.heatmap, width = 6, height = 5)

#################################
### 4. Asymmetry distribution
#################################

p.asym.dist <- ggplot(ecf.connectivity.stats, aes(x = asym_med)) +
  geom_histogram(bins = 60, fill = "#2c7fb8", color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 14) +
  labs(
    x = "Asymmetry (median qcov − median scov)",
    y = "Number of ECFs",
    title = "Distribution of hierarchical asymmetry"
  )

ggsave(file.path(this.part.dir, "plot3_asym-dist.png"),
       p.asym.dist, width = 6, height = 4)

#################################
### 4. Asymmetry vs. Length
#################################

p.asym.len <- ggplot(
  ecf.connectivity.stats,
  aes(x = asym_mean, y = length)
) +
  geom_bin2d(bins = 60) +
  heatmap_fill +
  scale_y_continuous(trans = "log10") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6) +
  labs(
    x = "Asymmetry (mean qcov − mean scov)",
    y = "ECF length (aa, log10 scale)"
  ) +
  base_theme

ggsave(file.path(this.part.dir, "plot4_asym-vs-length.png"),
       p.asym.len, width = 6, height = 5)

#################################
### 5. Asymmetry vs. Strength
#################################

p.asym.str <- ggplot(ecf.connectivity.stats,
                     aes(x = asym_med, y = strength_med)) +
  geom_bin2d(bins = 70) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 14) +
  labs(
    x = "Asymmetry (median qcov − median scov)",
    y = "Median strength ((qcov + scov)/2)",
    title = "Hierarchical role vs match strength"
  )

ggsave(file.path(this.part.dir, "plot5_asym-vs-strength.png"),
       p.asym.str, width = 6, height = 5)

#################################
### 6. Median qcov vs. scov
### clean diagonal view
#################################

p.qs.node <- ggplot(ecf.connectivity.stats,
                    aes(x = med.qcov, y = med.scov)) +
  geom_bin2d(bins = 70) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_bw(base_size = 14) +
  labs(
    x = "Median qcov",
    y = "Median scov",
    title = "Per-ECF similarity geometry"
  )

ggsave(file.path(this.part.dir, "plot6_node-qcov-scov.png"),
       p.qs.node, width = 6, height = 5)


#################################
### 7. Asymmetry vs. Degree
#################################

p.asym.deg <- ggplot(ecf.connectivity.stats,
                     aes(x = asym_med, y = degree)) +
  geom_bin2d(bins = 70) +
  scale_y_continuous(trans = "log10") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 14) +
  labs(
    x = "Asymmetry (median)",
    y = "Degree (log10 scale)",
    title = "Hierarchy vs connectivity"
  )

ggsave(file.path(this.part.dir, "plot7_asym-vs-degree.png"),
       p.asym.deg, width = 6, height = 5)

#################################
### 8. Asymmetry vs. Neighbour
### family count
#################################

p.asym.fam <- ggplot(ecf.connectivity.stats,
                     aes(x = asym_med, y = no.fam)) +
  geom_bin2d(bins = 60) +
  scale_y_continuous() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 14) +
  labs(
    x = "Asymmetry (median)",
    y = "Number of distinct neighbor families (log10 scale)",
    title = "Hierarchy vs family diversity"
  )

ggsave(file.path(this.part.dir, "plot8_asym-vs-family-div.png"),
       p.asym.fam, width = 6, height = 5)

#################################
### 9. Strength vs. Length
#################################

p.str.len <- ggplot(ecf.connectivity.stats,
                    aes(x = strength_med, y = length)) +
  geom_bin2d(bins = 70) +
  scale_y_continuous(trans = "log10") +
  theme_bw(base_size = 14) +
  labs(
    x = "Median strength",
    y = "ECF length (log10 scale)",
    title = "Match strength vs fragment size"
  )

ggsave(file.path(this.part.dir, "plot9_strength-vs-length.png"),
       p.str.len, width = 6, height = 5)

#################################
### 10. Statistical models to 
### resolve relationship between
### no.fam, asymm, length and degree
#################################

asym.fam.model <- copy(ecf.connectivity.stats)
asym.fam.model[, log_degree := log1p(degree)]      # log(1+degree) avoids log(0)
asym.fam.model[, log_length := log10(length)]

m1 <- glm(no.fam ~ log_degree,
          data = asym.fam.model,
          family = poisson())
m2 <- glm(no.fam ~ log_degree + asym_med,
          data = asym.fam.model,
          family = poisson())
m3 <- glm(no.fam ~ log_degree + asym_med + log_length,
          data = asym.fam.model,
          family = poisson())
summary(m3)

coef_df <- tidy(m3, conf.int = TRUE, exponentiate = TRUE)
coef_df <- coef_df[coef_df$term != "(Intercept)", ]
coef_df$term <- c(
  "log(Degree)",
  "Asymmetry",
  "log(Length)"
)
p.glm.model <- ggplot(coef_df, aes(x = estimate, y = term)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), width = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(
    x = "Incidence rate ratio (exp(coef))",
    y = NULL
  ) +
  theme_classic(base_size = 14)



#################################
### 10. Containment network
#################################

containment.network <- ecf.network[qcov >= DOMAIN.CONTAINMENT.COV.THR]
# convert into igraph object
containment.graph <- graph_from_data_frame(containment.network, directed = TRUE)
containment.graph.nodes <- V(containment.graph)$name
no.nodes.cg <- length(containment.graph.nodes)

# 1. Find all pairs of ECFs that contain an ECF contained in both of them
# 2. Calculate cov.max between all of them

contain.neighbour.cov.max.list <- lapply(1:no.nodes.cg, function(node.index){
  # cat(node.index,"\n")
  out <- NULL
  this.node <- containment.graph.nodes[node.index]
  # this.node.transitivity <- full.network.tr[which(names(full.network.tr)==this.node)]
  # 2. Identify all possible pairwise connections between neighbours in the containment network
  neighbours <- neighbors(containment.graph, this.node, mode = "out")
  no.neighbours <- length(neighbours)
  if(no.neighbours >= 2){
    neighbour.pairs <- combn(neighbours$name, 2, simplify = FALSE)
    neighbour.pairs.dt <- as.data.table(do.call(rbind, neighbour.pairs))
    no.pairs <- nrow(neighbour.pairs.dt)
    neighbour.pair.weigths <- ecfs.vs.ecfs[.(neighbour.pairs.dt)]$max.cov
    neighbour.pair.weigths <- fcoalesce(neighbour.pair.weigths, 0)
    
    # 5. Return 
    out <- data.table(ecf = this.node, neighbour.cov.max = neighbour.pair.weigths)  
  }
  return(out)
})

contain.neighbour.cov.max <- rbindlist(contain.neighbour.cov.max.list)
contain.neighbour.cov.max.avg <- contain.neighbour.cov.max[
  , .(avg.neighbour.cov.max = mean(neighbour.cov.max, na.rm = TRUE)),
  by = ecf
]
contain.neighbour.cov.max.avg <- merge(contain.neighbour.cov.max.avg, contain.neighbour.cov.max[, .N, by = ecf], by = "ecf")

data2plot <- rbind(data.table(type = "full", value = contain.neighbour.cov.max$neighbour.cov.max),
                   data.table(type = "average", value = contain.neighbour.cov.max.avg$avg.neighbour.cov.max))

# Create the histogram with facets and independent Y axes
p.qij <- ggplot(data2plot, aes(x = value)) +
  geom_histogram(bins = 30, fill = "#69b3a2", alpha = 0.7, color = "black") +  # Histogram
  facet_wrap(~ type, scales = "free_y") +  # Independent Y axis per panel
  labs(
    x = "Value",
    y = "Frequency"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank()
  )

p.qij.fn <- file.path(this.part.dir, "plot10_neighbour-cov.png")
ggsave(p.qij.fn, p.qij, width = 7, height = 3)

# rgb(51, 91, 170, maxColorValue = 255)
"#335BAA"

# Create the histogram with facets and independent Y axes
p.cont.cont <- ggplot(data2plot[type == "full"], aes(x = value)) +
  geom_histogram(
    bins = 40,
    fill = "#335BAA",
    color = "white",
    linewidth = 0.2
  ) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.6) +
  scale_x_continuous(
    breaks = seq(0, 1, 0.2),
    expand = c(0.03, 0.03)
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    x = "Container–container similarity (max coverage)",
    y = "Number of container pairs"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    plot.margin = margin(10, 10, 10, 10)
  )
p.cont.cont


#################################
### Generate Figure 3
#################################

# p.cov.heatmap # panel A
# p.connect.cov.heatmap # panel B
# p.asym.len # panel C
# p.cont.cont # panel D

fig3 <- (
  (p.cov.heatmap | p.connect.cov.heatmap) /
    (p.asym.len    | p.cont.cont)
) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0.02, 0.98)
  )
# fig3

fig3.outfile <- file.path(ecf.figures.dir, "figure3.png")
ggsave(fig3.outfile, fig3, width = 9, height = 7)

#################################
### Generate Supplementary Figs
#################################

# p.asym.dist
# p.asym.fam
# p.asym.deg
# p.glm.model

p.asym.dist.outfile <- file.path(ecf.suppl.figures.dir, "suppl-3p1_asymm-dist.png")
ggsave(p.asym.dist.outfile, p.asym.dist, width = 5, height = 4)

p.no.fam.asymm <- (
  (p.asym.fam / p.asym.deg / p.glm.model)
) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0.02, 0.98)
  )
p.no.fam.asymm.outfile <- file.path(ecf.suppl.figures.dir, "suppl-3p2_asymm-nofam.png")
ggsave(p.no.fam.asymm.outfile, p.no.fam.asymm, width = 5, height = 9)

# polish the GLM table
tab <- tidy(m3)  # non-exponentiated
## Remove intercept (optional — usually removed in papers)
tab <- tab[tab$term != "(Intercept)", ]
## Compute IRR (exp(beta)) and 95% CI
tab$IRR  <- exp(tab$estimate)
tab$CIlo <- exp(tab$estimate - 1.96 * tab$std.error)
tab$CIhi <- exp(tab$estimate + 1.96 * tab$std.error)
## Format numbers to 2 significant digits
fmt2 <- function(x) signif(x, 2)
tab$estimate  <- fmt2(tab$estimate)
tab$std.error <- fmt2(tab$std.error)
tab$statistic <- fmt2(tab$statistic)
tab$IRR       <- fmt2(tab$IRR)
tab$CIlo      <- fmt2(tab$CIlo)
tab$CIhi      <- fmt2(tab$CIhi)
## Format p-values nicely
tab$p.value <- ifelse(tab$p.value < 0.001,
                      "<0.001",
                      formatC(tab$p.value, digits = 2, format = "f"))
## Rename terms for readability
tab$term <- dplyr::recode(tab$term,
                          log_degree = "log1p(degree)",
                          asym_med   = "asymmetry",
                          log_length = "log10(length)")
## Select and reorder columns
tab_final <- tab[, c("term", "estimate", "std.error",
                     "IRR", "CIlo", "CIhi", "statistic", "p.value")]
colnames(tab_final) <- c(
  "Predictor",
  "β",
  "SE",
  "IRR",
  "CI 2.5%",
  "CI 97.5%",
  "z",
  "p"
)
# tab_final

tab.outfile <- file.path(ecf.suppl.figures.dir, "suppl-tab-3p1_GLM-results.csv")
write.csv(tab_final, tab.outfile, row.names = FALSE)


## Function to convert tab_final into LaTeX table
make_latex_table <- function(tab_final,
                             file,
                             caption = "Poisson GLM results for neighbour-family diversity.",
                             label = "tab:glm_asymmetry") {
  
  sink(file)
  
  cat("\\begin{table}[ht]\n")
  cat("\\centering\n")
  cat(paste0("\\caption{", caption, "}\n"))
  cat(paste0("\\label{", label, "}\n"))
  cat("\\begin{tabular}{lrrrrrrr}\n")
  cat("\\hline\n")
  
  cat(paste(colnames(tab_final), collapse = " & "))
  cat(" \\\\\n")
  cat("\\hline\n")
  
  for (i in 1:nrow(tab_final)) {
    cat(paste(tab_final[i, ], collapse = " & "))
    cat(" \\\\\n")
  }
  
  cat("\\hline\n")
  cat("\\end{tabular}\n")
  cat("\\end{table}\n")
  
  sink()
}

tab.latex.outfile <- file.path(ecf.suppl.figures.dir, "suppl-tab-3p1_GLM-results-latex.txt")
make_latex_table(tab_final, file = tab.latex.outfile)

