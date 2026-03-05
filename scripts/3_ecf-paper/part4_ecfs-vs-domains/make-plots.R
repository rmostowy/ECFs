### load libraries
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(patchwork)))
suppressWarnings(suppressPackageStartupMessages(library(ggh4x)))
# suppressWarnings(suppressPackageStartupMessages(library(ComplexUpset)))
# suppressWarnings(suppressPackageStartupMessages(library(ggraph)))
# suppressWarnings(suppressPackageStartupMessages(library(broom)))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.R")
source(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

#################################
### Set up output directory
#################################
this.part.dir <- file.path(ecf.proj.dir, "part4_ecfs-vs-domains")
make.dir(this.part.dir)

#################################
### Load all data:
### 1. ECF tables
### 2. domain tables
### 3. ecf vs domain tables
### 4. reprseq vs ecfs
### 5. 
### 6. 
### 7. 
#################################

# ECF tables
ecf.table <- fread(ecf.table.output)
ecf.family.table <- fread(ecf.family.table.output)

# domain tables
ecod.domains <- fread(ecod.domains.output)
pfam.domains <- fread(pfam.domains.output)

# ecf vs domains tables
ecfs.vs.ecod <- fread(ecfs.vs.ecod.output)
ecfs.vs.pfam <- fread(ecfs.vs.pfam.output)
ecfs.vs.ecod.hhsenser <- fread(ecfs.vs.ecod.hhsenser.output)
ecfs.vs.pfam.hhsenser <- fread(ecfs.vs.pfam.hhsenser.output)

# reprseq vs ecfs
reprseq.vs.ecfs.all <- fread(reprseq.ecfs.output)

# functional annotations
reprseq.functions.unique <- fread(reprseq.functions.unique.output)
reprseq.functions.all <- fread(reprseq.functions.all.output)

#################################
### Process loaded tables
#################################

setkey(ecf.family.table, ecf)

all.ecfs <- unique(ecf.table$ecf.name)
Ne <- length(all.ecfs)

setkey(ecod.domains, full.id)
setkey(pfam.domains, accession)

ecfs.vs.ecod[, x.id := ecod.domains[.(ecfs.vs.ecod$sname)]$x_id]
ecfs.vs.ecod[, h.id := ecod.domains[.(ecfs.vs.ecod$sname)]$h_id]
ecfs.vs.ecod[, t.id := ecod.domains[.(ecfs.vs.ecod$sname)]$t_id]
ecfs.vs.ecod[, x.level := ecod.domains[.(ecfs.vs.ecod$sname)]$x.level]
ecfs.vs.ecod[, h.level := ecod.domains[.(ecfs.vs.ecod$sname)]$h.level]
ecfs.vs.ecod[, t.level := ecod.domains[.(ecfs.vs.ecod$sname)]$t.level]

ecfs.vs.ecod.conf <- ecfs.vs.ecod %>%
  filter(eval <= MAX.EVALUE.FOR.ECF.DOMAIN.CONFIDENT & prob >= MIN.PROB.FOR.ECF.DOMAIN.CONFIDENT)
ecfs.vs.pfam.conf <- ecfs.vs.pfam %>%
  filter(eval <= MAX.EVALUE.FOR.ECF.DOMAIN.CONFIDENT & prob >= MIN.PROB.FOR.ECF.DOMAIN.CONFIDENT)

reprseq.vs.ecfs <- reprseq.vs.ecfs.all %>% 
  filter(prob >= MIN.PROB.FOR.DOMAIN.DETECT.IN.RHMM,
         eval <= MAX.EVAL.FOR.DOMAIN.DETECT.IN.RHMM, 
         qcov >= MIN.RHMM.COV.FOR.DOMAIN.DETECT.IN.RHMM,
         scov >= MIN.DOMAIN.COV.FOR.DOMAIN.DETECT.IN.RHMM,
         hit.length >= MIN.HIT.LENGTH.FOR.DOMAIN.DETECT.IN.RHMM
  )

############################################
### 1. Compare how many ECFs hit
### domain databases
############################################

ecf.domain.coverage <- data.table(ecf = all.ecfs)

# ECOD
ecf.domain.coverage[, ecod.eval0p1 := all.ecfs %in% ecfs.vs.ecod$qname]
ecf.domain.coverage[, ecod.eval0p001 := all.ecfs %in% ecfs.vs.ecod.conf$qname]
ecf.domain.coverage[, ecod.eval10 := all.ecfs %in% ecfs.vs.ecod.hhsenser$qname]
# PFAM
ecf.domain.coverage[, pfam.eval0p1 := all.ecfs %in% ecfs.vs.pfam$qname]
ecf.domain.coverage[, pfam.eval0p001 := all.ecfs %in% ecfs.vs.pfam.conf$qname]
ecf.domain.coverage[, pfam.eval10 := all.ecfs %in% ecfs.vs.pfam.hhsenser$qname]

# dt <- copy(ecf.domain.coverage)

# Helper: assign the 4-category label for a given pair of booleans
assign_cat <- function(has_ecod, has_pfam) {
  fifelse(
    has_ecod & has_pfam, "Both (ECOD + Pfam)",
    fifelse(
      has_ecod & !has_pfam, "ECOD only",
      fifelse(!has_ecod & has_pfam, "Pfam only", "Neither")
    )
  )
}

# Build a long table: one row per ECF × threshold
# (thresholds ordered: strict -> lenient)
long_list <- list(
  data.table(ecf = ecf.domain.coverage$ecf,
             threshold = "1e-3",
             has_ecod = ecf.domain.coverage$ecod.eval0p001,
             has_pfam = ecf.domain.coverage$pfam.eval0p001),
  data.table(ecf = ecf.domain.coverage$ecf,
             threshold = "1e-1",
             has_ecod = ecf.domain.coverage$ecod.eval0p1,
             has_pfam = ecf.domain.coverage$pfam.eval0p1),
  data.table(ecf = ecf.domain.coverage$ecf,
             threshold = "10",
             has_ecod = ecf.domain.coverage$ecod.eval10,
             has_pfam = ecf.domain.coverage$pfam.eval10)
)
long_dt <- rbindlist(long_list)

# Category + counts
long_dt[, category := assign_cat(has_ecod, has_pfam)]

counts_dt <- long_dt[, .N, by = .(threshold, category)]

# Factor ordering for nicer plot
counts_dt[, threshold := factor(threshold, levels = c("1e-3", "1e-1", "10"))]
counts_dt[, category := factor(
  category,
  levels = c("Both (ECOD + Pfam)", "ECOD only", "Pfam only", "Neither")
)]

# (Optional) a clean, colorblind-friendly palette
cat_cols <- c(
  "Both (ECOD + Pfam)" = "#009E73",
  "ECOD only"          = "#0072B2",
  "Pfam only"          = "#D55E00",
  "Neither"            = "grey70"
)

# Plot
p.ecf.domain.coverage <- ggplot(counts_dt, aes(x = threshold, y = N, fill = category)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.25) +
  scale_fill_manual(values = cat_cols, name = NULL) +
  labs(
    x = "E-value threshold (strict → lenient)",
    y = "Number of ECFs"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    legend.text  = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.key.width  = unit(0.4, "cm"),
    legend.key.height = unit(0.4, "cm"),
    legend.margin = margin(t = 0),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )
p.ecf.domain.coverage.outfile <- file.path(this.part.dir, "plot1_ecf-domain-cov.png")
ggsave(p.ecf.domain.coverage.outfile, p.ecf.domain.coverage, width = 6, height = 5)


############################################
### 2. Classify ECFs based on domain hits and
### analyse the class distribution based
### on ECF length
############################################

# A) ECOD comparisons

mt <- copy(ecfs.vs.ecod)
mt[, conf := (prob >= MIN.PROB.FOR.DOMAIN.RELATEDNESS)]
mt[, is.equiv := conf & (prob >= MIN.PROB.FOR.DOMAIN.EQUIVALENCE & qcov >= MIN.QCOV.FOR.DOMAIN.EQUIVALENCE &
                    scov >= MIN.SCOV.FOR.DOMAIN.EQUIVALENCE & eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE)]
mt[, is.subdom := conf & (prob >= MIN.PROB.FOR.SUBDOMAIN & qcov >= MIN.QCOV.FOR.SUBDOMAIN &
                     scov <= MAX.SCOV.FOR.SUBDOMAIN & eval <= MAX.EVAL.FOR.SUBDOMAIN)]
mt[, is.superdom := conf & (prob >= MIN.PROB.FOR.SUPERDOMAIN & qcov <= MAX.QCOV.FOR.SUPERDOMAIN &
                       scov >= MIN.SCOV.FOR.SUPERDOMAIN & eval <= MAX.EVAL.FOR.SUPERDOMAIN)]
ecf.domain.assignment <- mt[, {
  
  is.conf <- any(conf)
  is.equiv.class <- any(is.equiv)
  is.subdom.class <- any(is.subdom)
  is.superdom.class <- any(is.superdom)
  class.assignment <- "distant"
  
  if(is.conf){
    class.assignment <- "related"
    if(is.equiv.class){
      class.assignment <- "equivalent"
    } else if(is.subdom.class & !is.superdom.class){
      class.assignment <- "subdomain"
    } else if(!is.subdom.class & is.superdom.class){
      class.assignment <- "superdomain"
    } else if(is.subdom.class & is.superdom.class){
      ambiguous.rows <- which(is.subdom | is.superdom)
      ambiguous.row.with.longest.hit <- ambiguous.rows[which.max(hit.length[ambiguous.rows])]
      this.is.subdomain <- is.subdom[ambiguous.row.with.longest.hit]
      class.assignment <- if(this.is.subdomain) "subdomain" else "superdomain"
    }
  }

 .(class.assignment = class.assignment) 
  
}, by = qname]

ecf.ecod.domain.assignment <- merge(ecf.domain.assignment, ecf.table[, .(qname = ecf.name, length, size)], by = "qname")
ecf.ecod.domain.assignment[, len_bin := cut(length,
                                       breaks = quantile(length, probs = seq(0, 1, by = 0.1), na.rm = TRUE),
                                       include.lowest = TRUE)]

class.plotting.order <- c("equivalent", "subdomain", "superdomain", "related", "distant")
prop_dt <- ecf.ecod.domain.assignment[, .N, by = .(len_bin, class.assignment)]
prop_dt[, prop := N / sum(N), by = len_bin]
prop_dt[, num := N, by = len_bin]
prop_dt[, class.assignment := factor(class.assignment, levels = rev(class.plotting.order))]

p.ecod.domain.cat.by.length <- ggplot(prop_dt, aes(x = len_bin, y = prop, fill = class.assignment)) +
  geom_col(width = 0.9) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "ECF length (quantile bins)", y = "Proportion", fill = "Class") +
   theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    legend.text  = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.key.width  = unit(0.4, "cm"),
    legend.key.height = unit(0.4, "cm"),
    legend.margin = margin(t = 0),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) 
p.ecod.domain.cat.by.length.outfile <- file.path(this.part.dir, "plot2_dom-class-ecod.png")
ggsave(p.ecod.domain.cat.by.length.outfile, p.ecod.domain.cat.by.length, width = 6, height = 5)


# B) PFAM comparisons

mt <- copy(ecfs.vs.pfam)
mt[, conf := (prob >= MIN.PROB.FOR.DOMAIN.RELATEDNESS)]
mt[, is.equiv := conf & (prob >= MIN.PROB.FOR.DOMAIN.EQUIVALENCE & qcov >= MIN.QCOV.FOR.DOMAIN.EQUIVALENCE &
                           scov >= MIN.SCOV.FOR.DOMAIN.EQUIVALENCE & eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE)]
mt[, is.subdom := conf & (prob >= MIN.PROB.FOR.SUBDOMAIN & qcov >= MIN.QCOV.FOR.SUBDOMAIN &
                            scov <= MAX.SCOV.FOR.SUBDOMAIN & eval <= MAX.EVAL.FOR.SUBDOMAIN)]
mt[, is.superdom := conf & (prob >= MIN.PROB.FOR.SUPERDOMAIN & qcov <= MAX.QCOV.FOR.SUPERDOMAIN &
                              scov >= MIN.SCOV.FOR.SUPERDOMAIN & eval <= MAX.EVAL.FOR.SUPERDOMAIN)]
ecf.domain.assignment <- mt[, {
  
  is.conf <- any(conf)
  is.equiv.class <- any(is.equiv)
  is.subdom.class <- any(is.subdom)
  is.superdom.class <- any(is.superdom)
  class.assignment <- "distant"
  
  if(is.conf){
    class.assignment <- "related"
    if(is.equiv.class){
      class.assignment <- "equivalent"
    } else if(is.subdom.class & !is.superdom.class){
      class.assignment <- "subdomain"
    } else if(!is.subdom.class & is.superdom.class){
      class.assignment <- "superdomain"
    } else if(is.subdom.class & is.superdom.class){
      ambiguous.rows <- which(is.subdom | is.superdom)
      ambiguous.row.with.longest.hit <- ambiguous.rows[which.max(hit.length[ambiguous.rows])]
      this.is.subdomain <- is.subdom[ambiguous.row.with.longest.hit]
      class.assignment <- if(this.is.subdomain) "subdomain" else "superdomain"
    }
  }
  
  .(class.assignment = class.assignment) 
  
}, by = qname]

ecf.pfam.domain.assignment <- merge(ecf.domain.assignment, ecf.table[, .(qname = ecf.name, length, size)], by = "qname")
ecf.pfam.domain.assignment[, len_bin := cut(length,
                                       breaks = quantile(length, probs = seq(0, 1, by = 0.1), na.rm = TRUE),
                                       include.lowest = TRUE)]

class.plotting.order <- c("equivalent", "subdomain", "superdomain", "related", "distant")
prop_dt <- ecf.pfam.domain.assignment[, .N, by = .(len_bin, class.assignment)]
prop_dt[, prop := N / sum(N), by = len_bin]
prop_dt[, num := N, by = len_bin]
prop_dt[, class.assignment := factor(class.assignment, levels = rev(class.plotting.order))]

p.pfam.domain.cat.by.length <- ggplot(prop_dt, aes(x = len_bin, y = prop, fill = class.assignment)) +
  geom_col(width = 0.9) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "ECF length (quantile bins)", y = "Proportion", fill = "Class") +
  theme_bw(base_size = 12) +
  theme(
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) 
p.pfam.domain.cat.by.length.outfile <- file.path(this.part.dir, "plot3_dom-class-pfam.png")
ggsave(p.pfam.domain.cat.by.length.outfile, p.pfam.domain.cat.by.length, width = 6, height = 5)

# C) Compare the two sets

ecod <- ecf.ecod.domain.assignment[, .(qname,
                 ecod_class = class.assignment,
                 length,
                 size)]
pfam <- ecf.pfam.domain.assignment[, .(qname,
                 pfam_class = class.assignment,
                 length,
                 size)]

# Union of ECF ids
all_ecf <- data.table(qname = sort(unique(c(ecod$qname, pfam$qname))))
all_ecf <- merge(all_ecf, ecf.table[, .(qname = ecf.name, length, size)], by = "qname")

# Merge
dt <- merge(all_ecf, ecod, by = c("qname", "length", "size"), all.x = TRUE)
dt <- merge(dt, pfam, by = c("qname", "length", "size"), all.x = TRUE)

# Add explicit "absent" class
dt[is.na(ecod_class), ecod_class := "absent"]
dt[is.na(pfam_class), pfam_class := "absent"]

# harmonize factor levels (make plots stable)
class.plotting.order2 <- c(class.plotting.order, "absent")
dt[, ecod_class := factor(ecod_class, levels = class.plotting.order2)]
dt[, pfam_class := factor(pfam_class, levels = class.plotting.order2)]

# Quick agreement stats
dt[, agree := (as.character(ecod_class) == as.character(pfam_class))]


cm <- dt[, .N, by = .(ecod_class, pfam_class)]

# Option 1: counts
p_cm_counts <- ggplot(cm, aes(x = pfam_class, y = ecod_class, fill = N)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(trans = "sqrt") +
  labs(x = "Pfam assignment class", y = "ECOD assignment class", fill = "Count") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Option 2: row-normalized proportions (within ECOD class)
cm[, prop_row := N / sum(N), by = ecod_class]
p_cm_rowprop <- ggplot(cm, aes(x = pfam_class, y = ecod_class, fill = prop_row)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Pfam assignment class", y = "ECOD assignment class", fill = "Row %") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p.ecod.pfam.agreement <- (p_cm_counts / p_cm_rowprop) + plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 24),
    plot.tag.position = c(0.02, 0.98)
  )
p.ecod.pfam.agreement.outfile <- file.path(this.part.dir, "plot4_ecod-pfam-agreement.png")
ggsave(p.ecod.pfam.agreement.outfile, p.ecod.pfam.agreement, width = 6, height = 10)

# Run a statistical chi-squared test 
tab <- table(dt$ecod_class, dt$pfam_class)
N <- sum(tab)
r <- nrow(tab); c <- ncol(tab)
chi <- suppressWarnings(chisq.test(tab))
V <- sqrt(unname(chi$statistic) / (N * min(r - 1, c - 1)))
diag_agree <- sum(diag(tab)) / N
# 4) Format p-value nicely
pval_txt <- if (is.na(chi$p.value)) {
  "NA"
} else if (chi$p.value < 2.2e-16) {
  "< 2.2e-16"
} else {
  format(chi$p.value, scientific = TRUE, digits = 2)
}

# Manuscript-ready sentence
sent <- sprintf(
  paste0(
    "ECOD- and Pfam-based assignment classes show a clear correspondence across ECFs ",
    "(Pearson's chi-squared test: χ²(%d) = %.0f, p %s; Cramér's V = %.3f). ",
    "Agreement is strongest in the domain-like categories (equivalent/related/subdomain), ",
    "while most mismatches shift into 'related', consistent with differences in reference coverage and boundaries."
  ),
  unname(chi$parameter),
  unname(chi$statistic),
  pval_txt,
  V
)
sent.outfile <- file.path(this.part.dir, "plot4_stats.txt")
write(sent, file = sent.outfile, ncolumns = 1)

# D) Look at class distribution within families
ecf.ecod.domain.assignment[, len_bin := NULL]
ecf.ecod.domain.assignment$family <- ecf.family.table[.(ecf.ecod.domain.assignment$qname)]$family
ecf.ecod.domain.assignment[, w := 1L]        # unweighted
fam_size <- ecf.ecod.domain.assignment[, .N, by = family]
fam_p <- ecf.ecod.domain.assignment[, .(w = sum(w)), by = .(family, class.assignment)]
fam_p[, p := w / sum(w), by = family]
setorder(fam_p, family, -p, -w, class.assignment)
dominant_dt <- fam_p[, .(
  dominant_class = class.assignment[1]
), by = family]
fam_div <- fam_p[, .(
  n_classes = .N,
  eff_shannon = exp(-sum(p * log(p)))
), by = family]
# merge dominant class info
fam_div <- merge(fam_div, dominant_dt, by = "family", all.x = TRUE)
fam_final <- merge(fam_div, fam_size, by = "family", all.x = TRUE)

# A) histogram: number of classes
pA <- ggplot(fam_final, aes(n_classes)) +
  geom_bar() +
  theme_classic(base_size = 12) +
  labs(x = "Number of unique domain classes", y = "No. ECF families")

# B) histogram: effective Shannon
pB <- ggplot(fam_final, aes(eff_shannon)) +
  geom_histogram(bins = 30) +
  theme_classic(base_size = 12) +
  labs(x = "Effective Shannon (exp(H))", y = "No. ECF families")

# C) eff_shannon vs N
pC <- ggplot(fam_final, aes(N, eff_shannon)) +
  geom_bin_2d(bins = 50) +
  # geom_point(alpha = 0.99, size = 1) +
  scale_x_continuous(trans = "log10") +
  theme_classic(base_size = 12) +
  labs(x = "Family size N (log10)", y = "Effective Shannon (exp(H))")

# D) dominant-class counts for low-diversity families
low <- fam_final[eff_shannon < 1.5]
pDdat <- rbind(
  low[, .N, by = dominant_class][, group := "All (exp(H)<1.5)"],
  low[N > 1, .N, by = dominant_class][, group := "N>1 (exp(H)<1.5)"]
)
pDdat[, dominant_class := factor(dominant_class, levels = rev(class.plotting.order))]
pDdat[, group := factor(group, levels = c("All (exp(H)<1.5)", "N>1 (exp(H)<1.5)"))]
pD <- ggplot(pDdat, aes(dominant_class, N, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75, color = "white", linewidth = 0.2) +
  coord_flip() +
  theme_classic(base_size = 12) +
  labs(x = NULL, y = "No. ECF families") +
  theme(legend.position = "top", legend.title = element_blank())

ecf.ecod.domain.assignment[, w := NULL]

# Combine into 2x2 with tags
p.domain.class.per.fam <- (pA | pB) / (pC | pD) + plot_annotation(tag_levels = "A")
p.domain.class.per.fam.outfile <- file.path(this.part.dir, "plot5_domain-class-per-fam.png")
ggsave(p.domain.class.per.fam.outfile, p.domain.class.per.fam, width = 8, height = 8)

# Old part of the code, as a sanity check and for educational purposes
# 
# get.class <- function(dt.all){
#   class.assignment <- "distant"
#   dt <- dt.all[prob >= MIN.PROB.FOR.DOMAIN.RELATEDNESS]
#   any.hits.conf <- length(dt$eval)>0
#   if(any.hits.conf){
#     class.assignment <- "related"
#     is.equiv <- dt$prob >= MIN.PROB.FOR.DOMAIN.EQUIVALENCE &
#       dt$qcov >= MIN.QCOV.FOR.DOMAIN.EQUIVALENCE &
#       dt$scov >= MIN.SCOV.FOR.DOMAIN.EQUIVALENCE &
#       dt$eval <= MAX.EVAL.FOR.DOMAIN.EQUIVALENCE
#     is.subdom <- dt$prob >= MIN.PROB.FOR.SUBDOMAIN &
#       dt$qcov >= MIN.QCOV.FOR.SUBDOMAIN &
#       dt$scov <= MAX.SCOV.FOR.SUBDOMAIN &
#       dt$eval <= MAX.EVAL.FOR.SUBDOMAIN
#     is.superdom <- dt$prob >= MIN.PROB.FOR.SUPERDOMAIN &
#       dt$qcov <= MAX.QCOV.FOR.SUPERDOMAIN &
#       dt$scov >= MIN.SCOV.FOR.SUPERDOMAIN &
#       dt$eval <= MAX.EVAL.FOR.SUPERDOMAIN
#     
#     if(any(is.equiv)){
#       class.assignment <- "equivalent"
#     } else if(any(is.subdom) & !any(is.superdom)){
#       class.assignment <- "subdomain"
#     } else if(!any(is.subdom) & any(is.superdom)){
#       class.assignment <- "superdomain"
#     } else if(any(is.subdom) & any(is.superdom)){
#       dt.subsup <- cbind(dt, data.table(is.equiv, is.subdom, is.superdom))
#       dt.subsup <- dt.subsup[is.subdom == T | is.superdom == T]
#       dt.subsup <- dt.subsup[which.max(hit.length)]
#       if(dt.subsup$is.subdom==T){
#         class.assignment <- "subdomain"
#       } else{
#         class.assignment <- "superdomain"
#       }
#     } else{
#       class.assignment <- "related"
#     }
#   }
#   return(class.assignment)
# }
# z <- ecfs.vs.ecod[, .(dom.class = get.class(.SD)), by = qname]
# my.table(z$dom.class)

############################################
### 3. Investigate ECFs which link ECODs
### and quantify how they link different 
### ECOD architectures
############################################

ecfs.vs.ecod.link <- ecfs.vs.ecod[prob >= MIN.PROB.FOR.ECF.DOMAIN.CONFIDENT]
# ecfs.vs.ecod.link[qname == "ecf_00001"]

ecf.ecod.linking.stats <- ecfs.vs.ecod.link[, .(
  no.unique.x = uniqueN(x.id),
  no.unique.h = uniqueN(h.id),
  no.unique.t = uniqueN(t.id)  
), by = qname]
ecf.ecod.linking.stats <- merge(ecf.ecod.linking.stats, ecf.ecod.domain.assignment, by = "qname", all.x = T)


plot_dt <- ecf.ecod.linking.stats[
  class.assignment %chin% c("equivalent","related","subdomain","superdomain") &
    !is.na(no.unique.x)
]
plot_dt[, class.assignment := factor(
  class.assignment,
  levels = class.plotting.order
)]

p.ecod.links <- ggplot(plot_dt, aes(x = class.assignment, y = no.unique.x, fill = class.assignment)) +
  geom_violin(width = .95, trim = TRUE, scale = "width",
              colour = "grey25", linewidth = .35) +
  geom_point(position = position_jitter(width = .14, height = 0),
             alpha = .2, size = 1, colour = "grey10") +
  scale_y_continuous(breaks = sort(unique(plot_dt$no.unique.x))) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(5, "Set2")[c(5,4,3,2)]) +
  labs(x = "ECF–ECOD relationship", y = "ECOD X-groups linked") +
  theme_classic(base_size = 12) +
  guides(fill = "none") +
  theme(
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    axis.text.x  = element_text(angle = 0, vjust = 1),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = .25),
    panel.grid.minor.y = element_blank()
  )
p.ecod.links.outfile <- file.path(this.part.dir, "plot6_ecod-links.png")
ggsave(p.ecod.links.outfile, p.ecod.links, width = 4, height = 4)


############################################
### 4. Look at the putative functions
### of reprseq which contain "new" ECFs
############################################

# new.ecfs <- ecf.domain.coverage[ecod.eval0p1 == F]$ecf
# reprseq.vs.new.ecfs <- reprseq.vs.ecfs[sname %chin% new.ecfs]
# new.ecfs.linked.rhmms <- unique(reprseq.vs.new.ecfs$qname)

new_ecfs <- ecf.domain.coverage[ecod.eval0p1 == FALSE, unique(ecf)]
reprseq_new <- unique(reprseq.vs.ecfs[sname %chin% new_ecfs, qname])

# new.ecfs.vs.pfam <- ecfs.vs.pfam.conf[qcov >= MIN.QCOV.FOR.DOMAIN.EQUIVALENCE & qname %chin% new_ecfs]
# new.ecfs.linked.pfams <- unique(new.ecfs.vs.pfam$accession)

# -------------------------
# 1) Build analysis table (one row per reprseq)
# -------------------------
dt <- copy(reprseq.functions.unique)
setnames(dt, "reprseq", "qname")
dt[, is_new := qname %chin% reprseq_new]
dt[, funct2 := fifelse(is.na(funct) | funct == "", "unknown", funct)]
dt[, class2 := fifelse(is.na(class) | class == "", "unknown", class)]

# -------------------------
# 2) Enrichment per function (Fisher), include "unknown" in dt but exclude from tests
# -------------------------
dt_test <- dt[funct2 != "unknown"]
enr <- dt_test[
  , {
    a <- sum(is_new)                               # NEW & this funct
    b <- .N - a                                    # not-NEW & this funct
    c <- sum(dt_test$is_new) - a                   # NEW & other funct
    d <- (nrow(dt_test) - sum(dt_test$is_new)) - b # not-NEW & other funct
    
    ft <- fisher.test(matrix(c(a,b,c,d), nrow = 2))
    .(n_new = a,
      n_total = .N,
      odds = unname(ft$estimate),
      p = ft$p.value)
  },
  by = .(funct2, class2)
]

enr[, `:=`(
  fdr = p.adjust(p, method = "BH"),
  log2OR = log2(odds),
  frac_new = n_new / n_total
)]

# -------------------------
# 3) Choose what to plot: significant if any, else top N by FDR then effect size
# -------------------------
FDR_CUT <- 0.05
SHOW_N  <- 40L
NTOP <- 6

plot_dt <- enr[order(fdr)][fdr <= FDR_CUT, head(.SD, NTOP), by = class2]

# plot_dt <- {
#   sig <- enr[fdr <= FDR_CUT]
#   if (nrow(sig) > 0L) sig else enr
# }[order(fdr, -abs(log2OR))][1:min(.N, SHOW_N)]

# order classes in facets
plot_dt[, class2 := factor(class2, levels = phrogs.class.order)]
# make a facet-aware y key and order it by log2OR (globally, but unique per facet)
plot_dt[, ykey := paste(funct2, class2, sep = "___")]
plot_dt[, ykey := reorder(ykey, log2OR)]


# -------------------------
# 4) Plot
# -------------------------
strip <- ggh4x::strip_themed(
  background_y = ggh4x::elem_list_rect(
    fill = grDevices::adjustcolor(
      phrogs.col[match(phrogs.class.order, phrogs.class)],
      alpha.f = 0.25
    ),
    colour = "black",
    linewidth = 0.4
  )
)

p.new.ecf.enrichment <- ggplot(plot_dt, aes(x = log2OR, y = ykey)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
  geom_point(aes(size = n_total, colour = -log10(fdr)), alpha = 0.9) +
  ggh4x::facet_grid2(
    rows = vars(class2),
    scales = "free_y",
    space  = "free_y",
    labeller = labeller(class2 = phrog.class.labels),
    strip = strip
  ) +
  scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
  scale_colour_gradient(
    low = "#D7301F", high = "darkblue",
    trans = "sqrt",
    limits = c(0, max(-log10(plot_dt$fdr), na.rm = TRUE)),
    name = expression(-log[10]("FDR"))
  ) +
  scale_size_continuous(name = "N reprseqs", range = c(1.5, 5)) +
  labs(
    x = "Enrichment among NEW-ECF-associated reprseqs (log2 odds ratio)",
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    strip.placement = "outside",
    strip.text.y = element_text(angle = 0, hjust = 0, face = "bold", size = 9),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 9),
    panel.spacing = unit(0.35, "lines"),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = .25)
  )

p.new.ecf.enrichment.outfile <- file.path(this.part.dir, "plot7_new-ecf-enrichment.png")
ggsave(p.new.ecf.enrichment.outfile, p.new.ecf.enrichment, width = 8, height = 10)


############################################
### Make Figures
############################################

fig5 <- (
  (p.ecf.domain.coverage / p.ecod.domain.cat.by.length / p.ecod.links) |
    plot_spacer() |
    p.new.ecf.enrichment
) +
  plot_layout(widths = c(1.3, 0.08, 1)) +   # middle number controls spacing
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 24),
    plot.tag.position = c(0.02, 0.98)
  )
fig5.outfile <- file.path(ecf.figures.dir, "figure5.png")
ggsave(fig5.outfile, fig5, width = 12, height = 10)

##########################################################
##########################################################
##########################################################
##########################################################



# out <- ecfs.vs.ecod[qname == "ecf_00077"]
# 
# mean.pident <- function(dt) mean(dt$pident)
# ecfs.vs.ecod[, .(mean_pident = mean.pident(.SD)), by = qname]
# 
# mean.vals <- function(dt){
#   list(
#     mean_pident = mean(dt$pident),
#     min_pident = min(dt$pident),
#     max_pident = max(dt$pident)
#   )
# }
# mean.val <- function(dt){
#   list(
#     mean_pident = mean(dt$pident)
#   )
# }
# ecfs.vs.ecod[, mean.vals(.SD), by = qname]
# ecfs.vs.ecod[, mean.val(.SD), by = qname]
# 
# 
# mean.vals.fast <- function(pident) {
#   list(
#     mean_pident = mean(pident),
#     min_pident = min(pident),
#     max_pident = max(pident)
#   )
# }
# ecfs.vs.ecod[, mean.vals.fast(pident), by = qname]