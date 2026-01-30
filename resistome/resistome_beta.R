#### OBJETIVO DE ESTE SCRIPT:
# Generar una gráfica de beta diversidad del resistoma al gusto de las demás y calcular si hay diferencias entre grupos
library(vegan)
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(GUniFrac)
#library(DESeq2)
library(MicEco)
library(dplyr)
#library(MiscMetabar)
library(patchwork)
library(tibble)
#library(ecole)
library(pairwiseAdonis)
library(mia)

# directorios útiles
script_dir = getwd()
setwd("../data/r_data")
data_dir = getwd()
setwd("../../results/beta_res")
res_dir = getwd()

setwd(data_dir)
meta_tab = read.csv("metadata_resistome.csv", sep = ",", row.names = 1)
meta_tab[order(row.names(meta_tab)), ] #por consistencia con luego las abundancias
META = sample_data(meta_tab)
arg_freqs_tab = read.csv("count_arg.csv", sep = ",", header = TRUE, row.names = 1)
mge_freqs_tab = read.csv("count_mge.csv", sep = ",", header = TRUE, row.names = 1)
ARG = otu_table(arg_freqs_tab, taxa_are_rows = TRUE)
MGE = otu_table(mge_freqs_tab, taxa_are_rows = TRUE)
resis = phyloseq(ARG, META)
transp = phyloseq(MGE, META)

# Once I have the phyloseq objects, I want to fix the naming of various things
resis@sam_data$sampling_point[resis@sam_data$sampling_point == "uru"] = "Lakes"
resis@sam_data$sampling_point[resis@sam_data$sampling_point == "ion"] = "Lakes"
resis@sam_data$sampling_point[resis@sam_data$sampling_point == "ard"] = "Ardley"

transp@sam_data$sampling_point[transp@sam_data$sampling_point == "uru"] = "Lakes"
transp@sam_data$sampling_point[transp@sam_data$sampling_point == "ion"] = "Lakes"
transp@sam_data$sampling_point[transp@sam_data$sampling_point == "ard"] = "Ardley"

resis@sam_data$type_g[resis@sam_data$type_g == "plastic"] = "Plastic"
resis@sam_data$type_g[resis@sam_data$type_g == "control"] = "Surrounding env."

transp@sam_data$type_g[transp@sam_data$type_g == "plastic"] = "Plastic"
transp@sam_data$type_g[transp@sam_data$type_g == "control"] = "Surrounding env"

##### ALL AT ONCE ####
setwd(res_dir)

#and obtain the equivalent distances
res_pcoa = ordinate(resis, method = "PCoA", distance = "bray")
res_nmds = ordinate(resis, method = "NMDS", distance = "bray")
transp_pcoa = ordinate(transp, method = "PCoA", distance = "bray")
transp_nmds = ordinate(transp, method = "NMDS", distance = "bray")

## Quick aside to calculate stress values. It works well enough, but should be improved upon  
# Get the distance matrices
bc_res_matrix = as.matrix(phyloseq::distance(resis, method = "bray"))
bc_transp_matrix = as.matrix(phyloseq::distance(transp, method = "bray"))
# Get the distance matrices of each ordination. Requieres mia and this is the hacky part
# The NMDS plots have the information available from the get-go, so I just need to calculate the PCoAs
resis_tse = convertFromPhyloseq(resis)
transp_tse = convertFromPhyloseq(transp)

resis_tse <- addMDS(resis_tse, FUN = getDissimilarity,
                    name = "PCoA", method = "bray", assay.type = "counts")
resis_pcoa_matrix <- as.matrix(dist(reducedDim(resis_tse, "PCoA")))
transp_tse <- addMDS(transp_tse, FUN = getDissimilarity,
                    name = "PCoA", method = "bray", assay.type = "counts")
transp_pcoa_matrix <- as.matrix(dist(reducedDim(transp_tse, "PCoA")))

stress_res = sum((resis_pcoa_matrix - bc_res_matrix)^2) / sum(bc_res_matrix^2)
stress_res
stress_transp = sum((transp_pcoa_matrix - bc_transp_matrix)^2) / sum(bc_transp_matrix^2)
stress_transp



# Define the palette
palette = c("#B22C2C", "#85B22C", "#B26F2C", "#2C85B2")

# Bray PCOA
png("bc_pcoa_res.png", width = 6, height = 5, units = "in", res = 300)
plot_ordination(resis, res_pcoa, "samples", color = "type_f", shape = "sampling_point", 
                title = "RESISTOME BRAY-CURTIS PCoA") +
  geom_point(size=6, alpha=0.55) +
  scale_color_manual(values = palette) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1)
invisible(dev.off())

# I'll do each plot twice and then stitch them together with patchwork. This version has less info because
# the intention is to make a 2x2 grid with all the info (such as the legend) in one of the squares
pcoa_res = plot_ordination(resis, res_pcoa, "samples", color = "type_f", shape = "sampling_point", title = "") +
  geom_point(size=6, alpha=0.55) +
  scale_color_manual(values = palette) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1, legend.position = "none",
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 12))


png("bc_pcoa_transp.png", width = 6, height = 5, units = "in", res = 300)
plot_ordination(transp, transp_pcoa, "samples", color = "type_f", shape = "sampling_point", 
                title = "MGEs BRAY-CURTIS PCoA") +
  scale_color_manual(values = palette) + 
  geom_point(size=6, alpha=0.55) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1)
invisible(dev.off())

pcoa_transp = plot_ordination(transp, transp_pcoa, "samples", color = "type_f", shape = "sampling_point", title = "") +
  geom_point(size=6, alpha=0.55) +
  scale_color_manual(values = palette) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1, legend.position = "none",
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 12))


# NMDS now
png("bc_nmds_res.png", width = 6, height = 5, units = "in", res = 300)
plot_ordination(resis, res_nmds, "samples", color = "type_f", shape = "sampling_point", 
                title = "RESISTOME BRAY-CURTIS NMDS") +
  scale_color_manual(values = palette) + 
  geom_point(size=6, alpha=0.55) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1)
invisible(dev.off())
nmds_resis = plot_ordination(resis, res_nmds, "samples", color = "type_f", shape = "sampling_point", title = "") +
  geom_point(size=6, alpha=0.55) +
  scale_color_manual(values = palette) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1, legend.position = "none",
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 12))

png("bc_nmds_transp.png", width = 6, height = 5, units = "in", res = 300)
plot_ordination(transp, transp_nmds, "samples", color = "type_f", shape = "sampling_point", 
                title = "MGEs BRAY-CURTIS NMDS") +
  scale_color_manual(values = palette) + 
  geom_point(size=6, alpha=0.55) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1)
invisible(dev.off())
nmds_transp = plot_ordination(transp, transp_nmds, "samples", color = "type_f", shape = "sampling_point", title = "") +
  geom_point(size=6, alpha=0.55) +
  scale_color_manual(values = palette) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1, legend.position = "none",
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 12))


# Now, betadispersion and PERMANOVAS. But, first, we need to generate distance objects
resis_dist = phyloseq::distance(resis, "bray")
transp_dist = phyloseq::distance(transp, "bray")

# betadisp
beta_resis = betadisper(resis_dist, group = c(data.frame(sample_data(resis))$tf.sp))
write.csv(permutest(beta_resis, permutations = 10000)$tab, "betadisp_resis.csv")

beta_transp = betadisper(transp_dist, group = c(data.frame(sample_data(transp))$tf.sp))
write.csv(permutest(beta_transp, permutations = 10000)$tab, "betadisp_transp.csv")


# PERMANOVA
perma_resis = adonis3(resis_dist ~ type_g + sampling_point + type_g:sampling_point, data = data.frame(sample_data(resis)), permutations = 1000)
write.csv(perma_resis$aov.tab,"perma_resis.csv")
pp_resis = pairwise.adonis(resis_dist, factors = resis@sam_data$tf.sp, p.adjust.m = "fdr", perm = 10000)
write.csv(pp_resis, "pp_resis_all.csv")
pp_resis = pairwise.adonis(resis_dist, factors = resis@sam_data$type_g, p.adjust.m = "fdr", perm = 10000)
write.csv(pp_resis, "pp_resis_tg.csv")
pp_resis = pairwise.adonis(resis_dist, factors = resis@sam_data$type_f, p.adjust.m = "fdr", perm = 10000)
write.csv(pp_resis, "pp_resis_tf.csv")
pp_resis = pairwise.adonis(resis_dist, factors = resis@sam_data$sampling_point, p.adjust.m = "fdr", perm = 10000)
write.csv(pp_resis, "pp_resis_area.csv")

perma_transp = adonis3(transp_dist ~ type_g + sampling_point + type_g:sampling_point, data = data.frame(sample_data(transp)), permutations = 1000)
write.csv(perma_transp$aov.tab,"perma_transp.csv")
pp_transp = pairwise.adonis(transp_dist, factors = transp@sam_data$tf.sp, p.adjust.m = "fdr", perm = 10000)
write.csv(pp_transp, "pp_transp_all.csv")
pp_transp = pairwise.adonis(transp_dist, factors = transp@sam_data$type_g, p.adjust.m = "fdr", perm = 10000)
write.csv(pp_transp, "pp_transp_tg.csv")
pp_transp = pairwise.adonis(transp_dist, factors = transp@sam_data$type_f, p.adjust.m = "fdr", perm = 10000)
write.csv(pp_transp, "pp_transp_tf.csv")
pp_transp = pairwise.adonis(transp_dist, factors = transp@sam_data$sampling_point, p.adjust.m = "fdr", perm = 10000)
write.csv(pp_transp, "pp_transp_area.csv")

##### ONLY ARDLEY #####
# As there are quite significant differences explainable by sampling site (ardley has nothing to do with ards), 
# maybe we're losing information by studying everything together and not by separate. For that reason, let's
# analyze just ardley samples first and then just ard samples second
ard_resis = subset_samples(resis, sampling_point == "Ardley")
ard_transp = subset_samples(transp, sampling_point == "Ardley")

# Just for tidyness sake
setwd("Ard")
#distances
res_pcoa = ordinate(ard_resis, method = "PCoA", distance = "bray")
res_nmds = ordinate(ard_resis, method = "NMDS", distance = "bray")
transp_pcoa = ordinate(ard_transp, method = "PCoA", distance = "bray")
transp_nmds = ordinate(ard_transp, method = "NMDS", distance = "bray")

# Now the whole stress thing
bc_res_matrix = as.matrix(phyloseq::distance(ard_resis, method = "bray"))
bc_transp_matrix = as.matrix(phyloseq::distance(ard_transp, method = "bray"))
resis_tse = convertFromPhyloseq(ard_resis)
transp_tse = convertFromPhyloseq(ard_transp)

resis_tse <- addMDS(resis_tse, FUN = getDissimilarity,
                    name = "PCoA", method = "bray", assay.type = "counts")
resis_pcoa_matrix <- as.matrix(dist(reducedDim(resis_tse, "PCoA")))
transp_tse <- addMDS(transp_tse, FUN = getDissimilarity,
                     name = "PCoA", method = "bray", assay.type = "counts")
transp_pcoa_matrix <- as.matrix(dist(reducedDim(transp_tse, "PCoA")))

stress_res = sum((resis_pcoa_matrix - bc_res_matrix)^2) / sum(bc_res_matrix^2)
stress_res
stress_transp = sum((transp_pcoa_matrix - bc_transp_matrix)^2) / sum(bc_transp_matrix^2)
stress_transp



# Bray PCOA
#tab10_cols = c("#ff7f0e", "#1f77b4")

png("bc_pcoa_res.png", width = 6, height = 5, units = "in", res = 300)
plot_ordination(ard_resis, res_pcoa, "samples", color = "type_f", shape = "sampling_point", 
                title = "RESISTOME BRAY-CURTIS PCoA") +
  scale_color_manual(values = palette) + 
  geom_point(size=6, alpha=0.55) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1)
invisible(dev.off())
ard_pcoa_res = plot_ordination(ard_resis, res_pcoa, "samples", color = "type_f", shape = "sampling_point", title = "") +
  geom_point(size=6, alpha=0.55) +
  scale_color_manual(values = palette) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1, legend.position = "none",
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 12))


png("bc_pcoa_transp.png", width = 6, height = 5, units = "in", res = 300)
plot_ordination(ard_transp, transp_pcoa, "samples", color = "type_f", shape = "sampling_point", 
                title = "MGEs BRAY-CURTIS PCoA") +
  scale_color_manual(values = palette) + 
  geom_point(size=6, alpha=0.55) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1)
invisible(dev.off())
ard_pcoa_transp = plot_ordination(ard_transp, transp_pcoa, "samples", color = "type_f", shape = "sampling_point", title = "") +
  geom_point(size=6, alpha=0.55) +
  scale_color_manual(values = palette) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1, legend.position = "none",
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 12))

# NMDS now
png("bc_nmds_res.png", width = 6, height = 5, units = "in", res = 300)
plot_ordination(ard_resis, res_nmds, "samples", color = "type_f", shape = "sampling_point", 
                title = "RESISTOME BRAY-CURTIS NMDS") +
  scale_color_manual(values = palette) + 
  geom_point(size=6, alpha=0.55) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1)
invisible(dev.off())
ard_nmds_res = plot_ordination(ard_resis, res_nmds, "samples", color = "type_f", shape = "sampling_point", title = "") +
  geom_point(size=6, alpha=0.55) +
  scale_color_manual(values = palette) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1, legend.position = "none",
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 12))


png("bc_nmds_transp.png", width = 6, height = 5, units = "in", res = 300)
plot_ordination(ard_transp, transp_nmds, "samples", color = "type_f", shape = "sampling_point", 
                title = "MGEs BRAY-CURTIS NMDS") +
  scale_color_manual(values = palette) + 
  geom_point(size=6, alpha=0.55) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1)
invisible(dev.off())
ard_nmds_transp = plot_ordination(ard_transp, transp_nmds, "samples", color = "type_f", shape = "sampling_point", title = "") +
  geom_point(size=6, alpha=0.55) +
  scale_color_manual(values = palette) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1, legend.position = "none",
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 12))



# betadispersion and PERMANOVAS. First we need to generate distance objects
resis_dist = phyloseq::distance(ard_resis, "bray")
transp_dist = phyloseq::distance(ard_transp, "bray")

# betadisp
beta_resis = betadisper(resis_dist, group = c(data.frame(sample_data(ard_resis))$type_f))
write.csv(permutest(beta_resis, permutations = 10000)$tab, "betadisp_resis.csv")

beta_transp = betadisper(transp_dist, group = c(data.frame(sample_data(ard_transp))$type_f))
write.csv(permutest(beta_transp, permutations = 10000)$tab, "betadisp_transp.csv")


# PERMANOVA
perma_resis = adonis3(resis_dist ~ type_f, data = data.frame(sample_data(ard_resis)), permutations = 1000)
write.csv(perma_resis$aov.tab,"perma_resis.csv")
pp_resis = pairwise.adonis(resis_dist, factors = ard_resis@sam_data$type_g, p.adjust.m = "fdr", perm = 10000)
write.csv(pp_resis, "pp_resis_tg.csv")
pp_resis = pairwise.adonis(resis_dist, factors = ard_resis@sam_data$type_f, p.adjust.m = "fdr", perm = 10000)
write.csv(pp_resis, "pp_resis_tf.csv")
perma_transp = adonis3(transp_dist ~ type_f, data = data.frame(sample_data(ard_transp)), permutations = 1000)
write.csv(perma_transp$aov.tab,"perma_transp.csv")
pp_transp = pairwise.adonis(transp_dist, factors = ard_transp@sam_data$type_g, p.adjust.m = "fdr", perm = 10000)
write.csv(pp_transp, "pp_transp_tg.csv")
pp_transp = pairwise.adonis(transp_dist, factors = ard_transp@sam_data$type_f, p.adjust.m = "fdr", perm = 10000)
write.csv(pp_transp, "pp_transp_tf.csv")

##### ONLY LAKES ####
lake_resis = subset_samples(resis, sampling_point == "Lakes")
lake_transp = subset_samples(transp, sampling_point == "Lakes")

# Just for tidyness sake
setwd("../Lakes")
#distances
res_pcoa = ordinate(lake_resis, method = "PCoA", distance = "bray")
res_nmds = ordinate(lake_resis, method = "NMDS", distance = "bray")
transp_pcoa = ordinate(lake_transp, method = "PCoA", distance = "bray")
transp_nmds = ordinate(lake_transp, method = "NMDS", distance = "bray")

# Now the whole stress thing
bc_res_matrix = as.matrix(phyloseq::distance(lake_resis, method = "bray"))
bc_transp_matrix = as.matrix(phyloseq::distance(lake_transp, method = "bray"))
resis_tse = convertFromPhyloseq(lake_resis)
transp_tse = convertFromPhyloseq(lake_transp)

resis_tse <- addMDS(resis_tse, FUN = getDissimilarity,
                    name = "PCoA", method = "bray", assay.type = "counts")
resis_pcoa_matrix <- as.matrix(dist(reducedDim(resis_tse, "PCoA")))
transp_tse <- addMDS(transp_tse, FUN = getDissimilarity,
                     name = "PCoA", method = "bray", assay.type = "counts")
transp_pcoa_matrix <- as.matrix(dist(reducedDim(transp_tse, "PCoA")))

stress_res = sum((resis_pcoa_matrix - bc_res_matrix)^2) / sum(bc_res_matrix^2)
stress_res
stress_transp = sum((transp_pcoa_matrix - bc_transp_matrix)^2) / sum(bc_transp_matrix^2)
stress_transp



# Bray PCOA
#tab10_cols = c("#ff7f0e", "#1f77b4")

png("bc_pcoa_res.png", width = 6, height = 5, units = "in", res = 300)
plot_ordination(lake_resis, res_pcoa, "samples", color = "type_f", shape = "sampling_point", 
                title = "RESISTOME BRAY-CURTIS PCoA") +
  scale_color_manual(values = palette) + 
  geom_point(size=6, alpha=0.55) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1)
invisible(dev.off())
lake_pcoa_res = plot_ordination(lake_resis, res_pcoa, "samples", color = "type_f", shape = "sampling_point", title = "") +
  geom_point(size=6, alpha=0.55) +
  scale_color_manual(values = palette) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1, legend.position = "none",
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 12))

png("bc_pcoa_transp.png", width = 6, height = 5, units = "in", res = 300)
plot_ordination(lake_transp, transp_pcoa, "samples", color = "type_f", shape = "sampling_point", 
                title = "MGEs BRAY-CURTIS PCoA") +
  scale_color_manual(values = palette) + 
  geom_point(size=6, alpha=0.55) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1)
invisible(dev.off())
lake_pcoa_transp = plot_ordination(lake_transp, transp_pcoa, "samples", color = "type_f", shape = "sampling_point", title = "") +
  geom_point(size=6, alpha=0.55) +
  scale_color_manual(values = palette) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1, legend.position = "none",
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 12))

# NMDS now
png("bc_nmds_res.png", width = 6, height = 5, units = "in", res = 300)
plot_ordination(lake_resis, res_nmds, "samples", color = "type_f", shape = "sampling_point", 
                title = "RESISTOME BRAY-CURTIS NMDS") +
  scale_color_manual(values = palette) + 
  geom_point(size=6, alpha=0.55) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1)
invisible(dev.off())
lake_nmds_res = plot_ordination(lake_resis, res_nmds, "samples", color = "type_f", shape = "sampling_point", title = "") +
  geom_point(size=6, alpha=0.55) +
  scale_color_manual(values = palette) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1, legend.position = "none",
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 12))


png("bc_nmds_transp.png", width = 6, height = 5, units = "in", res = 300)
plot_ordination(lake_transp, transp_nmds, "samples", color = "type_f", shape = "sampling_point", 
                title = "MGEs BRAY-CURTIS NMDS") +
  scale_color_manual(values = palette) + 
  geom_point(size=6, alpha=0.55) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1)
invisible(dev.off())
lake_nmds_transp = plot_ordination(lake_transp, transp_nmds, "samples", color = "type_f", shape = "sampling_point", title = "") +
  geom_point(size=6, alpha=0.55) +
  scale_color_manual(values = palette) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1, legend.position = "none",
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 12))

# betadispersion and PERMANOVAS
resis_dist = phyloseq::distance(lake_resis, "bray")
transp_dist = phyloseq::distance(lake_transp, "bray")

# betadisp
beta_resis = betadisper(resis_dist, group = c(data.frame(sample_data(lake_resis))$type_f))
write.csv(permutest(beta_resis, permutations = 10000)$tab, "betadisp_resis.csv")

beta_transp = betadisper(transp_dist, group = c(data.frame(sample_data(lake_transp))$type_f))
write.csv(permutest(beta_transp, permutations = 10000)$tab, "betadisp_transp.csv")

# PERMANOVA
perma_resis = adonis3(resis_dist ~ type_f, data = data.frame(sample_data(lake_resis)), permutations = 1000)
write.csv(perma_resis$aov.tab,"perma_resis.csv")
pp_resis = pairwise.adonis(resis_dist, factors = lake_resis@sam_data$type_g, p.adjust.m = "fdr", perm = 10000)
write.csv(pp_resis, "pp_resis_tg.csv")
pp_resis = pairwise.adonis(resis_dist, factors = lake_resis@sam_data$type_f, p.adjust.m = "fdr", perm = 10000)
write.csv(pp_resis, "pp_resis_tf.csv")
perma_transp = adonis3(transp_dist ~ type_f, data = data.frame(sample_data(lake_transp)), permutations = 1000)
write.csv(perma_transp$aov.tab,"perma_transp.csv")
pp_transp = pairwise.adonis(transp_dist, factors = lake_transp@sam_data$type_g, p.adjust.m = "fdr", perm = 10000)
write.csv(pp_transp, "pp_transp_tg.csv")
pp_transp = pairwise.adonis(transp_dist, factors = lake_transp@sam_data$type_f, p.adjust.m = "fdr", perm = 10000)
write.csv(pp_transp, "pp_transp_tf.csv")

### And now we use a whole new framework to make the dbRDA and get its weights. As it is the last of the plots we will use
## this is when we stitch everything together
setwd(res_dir)
library("miaViz")

resis@sam_data[["Sample..type"]] = resis@sam_data$type_g
resis@sam_data[["Sampling..area"]] = resis@sam_data$sampling_point

resis_tse = makeTreeSEFromPhyloseq(resis)
resis_tse = transformAssay(resis_tse, # Bray Curtis, relativas
                           assay.type = "counts",
                           method = "relabundance",
                           pseudocount = TRUE)
resis_tse = runRDA(resis_tse,
                   assay.type = "relabundance",
                   formula = assay ~ Sample..type + Sampling..area, # + type_g:sam_area,
                   distance = "bray",
                   na.action = na.exclude)

# Now, the whole stress fiasco
bc_res_matrix = as.matrix(phyloseq::distance(resis, method = "bray"))
resis_dbrda_matrix <- as.matrix(dist(reducedDim(resis_tse, "RDA")))
stress_res = sum((resis_dbrda_matrix - bc_res_matrix)^2) / sum(bc_res_matrix^2)
stress_res

# dbRDA plot
png("RDA.png", width = 6, height = 4, units = "in", res = 300)
plotRDA(resis_tse, "RDA", colour_by = "type_f", 
        shape_by = "Sampling..area", point_size = 5, add.ellipse = "colour",
        arrow.size = 0.50, vec.size = 0.5, label.size = 3, vec.text = FALSE) +
  scale_color_manual(values = palette) +
  theme(aspect.ratio = 1)
invisible(dev.off())
dbrda_resis = plotRDA(resis_tse, "RDA", colour_by = "type_f", 
               shape_by = "Sampling..area", point_size = 6, add.ellipse = "colour", ellipse.linewidth = 1,
               add.vectors = FALSE) +
  scale_color_manual(values = palette) +
  theme_bw() +
  theme(aspect.ratio = 1, legend.position = "none",
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 12)) +
  labs(x = "dbRDA 1", y = "dbRDA 2") 

# Now with only three
three_joined = (pcoa_res | nmds_resis)/(dbrda_resis | dbrda_resis)
three_joined
ggsave(filename = "three_joined_arg.tif", plot = three_joined, device = "tiff",
       width = 10, height = 8, dpi = 300)


transp@sam_data[["Sample..type"]] = transp@sam_data$type_g
transp@sam_data[["Sampling..area"]] = transp@sam_data$sampling_point

transp_tse = makeTreeSEFromPhyloseq(transp)
transp_tse = transformAssay(transp_tse, # Bray Curtis, relativas
                           assay.type = "counts",
                           method = "relabundance",
                           pseudocount = TRUE)
transp_tse = runRDA(transp_tse,
                   assay.type = "relabundance",
                   formula = assay ~ Sample..type + Sampling..area, # + type_g:sam_area,
                   distance = "bray",
                   na.action = na.exclude)

# Stress
bc_transp_matrix = as.matrix(phyloseq::distance(transp, method = "bray"))
transp_dbrda_matrix <- as.matrix(dist(reducedDim(transp_tse, "RDA")))
stress_res = sum((transp_dbrda_matrix - bc_res_matrix)^2) / sum(bc_res_matrix^2)
stress_res


png("RDA_mge.png", width = 6, height = 4, units = "in", res = 300)
plotRDA(transp_tse, "RDA", colour_by = "type_f", 
        shape_by = "Sampling..area", point_size = 5, add.ellipse = "colour",
        arrow.size = 0.50, vec.size = 0.5, label.size = 3, vec.text = FALSE) +
  scale_color_manual(values = palette) +
  theme(aspect.ratio = 1)
invisible(dev.off())
dbrda_transp = plotRDA(transp_tse, "RDA", colour_by = "type_f", 
                      shape_by = "Sampling..area", point_size = 6, add.ellipse = "colour", ellipse.linewidth = 1,
                      add.vectors = FALSE) +
  scale_color_manual(values = palette) +
  theme_bw() +
  theme(aspect.ratio = 1, legend.position = "none",
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 12)) +
  labs(x = "dbRDA 1", y = "dbRDA 2") 

three_joined = (pcoa_transp | nmds_transp)/(dbrda_transp | dbrda_transp)
three_joined
ggsave(filename = "three_joined_mge.tif", plot = three_joined, device = "tiff",
       width = 10, height = 8, dpi = 300)
transp_nmds

# Now ardley
resis_ard = subset_samples(resis, sampling_point == "Ardley")
ardR_tse = makeTreeSEFromPhyloseq(resis_ard)
ardR_tse = transformAssay(ardR_tse, # Bray Curtis, relativas
                           assay.type = "counts",
                           method = "relabundance",
                           pseudocount = TRUE)
ardR_tse = runRDA(ardR_tse,
                   assay.type = "relabundance",
                   formula = assay ~ type_f,
                   distance = "bray",
                   na.action = na.exclude)

bc_res_matrix = as.matrix(phyloseq::distance(resis_ard, method = "bray"))
resis_dbrda_matrix <- as.matrix(dist(reducedDim(ardR_tse, "RDA")))
stress_res = sum((resis_dbrda_matrix - bc_res_matrix)^2) / sum(bc_res_matrix^2)
stress_res


setwd("Ard")
png("dbRDA.png", width = 6, height = 4, units = "in", res = 300)
plotRDA(ardR_tse, "RDA", colour_by = "type_f",
        point_size = 5, add.ellipse = "colour",
        arrow.size = 0.50, vec.size = 0.5, label.size = 3, vec.text = FALSE) +
  scale_color_manual(values = palette) +
  theme(aspect.ratio = 1)
invisible(dev.off())
dbrda_ardR = plotRDA(ardR_tse, "RDA", colour_by = "type_f", 
                       shape_by = "Sampling..area", point_size = 6, add.ellipse = "colour", ellipse.linewidth = 1,
                       add.vectors = FALSE) +
  scale_color_manual(values = palette) +
  theme_bw() +
  theme(aspect.ratio = 1, legend.position = "none",
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 12)) +
  labs(x = "dbRDA 1", y = "dbRDA 2") 

all_together_ardR = (ard_pcoa_res | ard_nmds_res)/(dbrda_ardR | dbrda_ardR)
print(all_together_ardR)
ggsave("all_joined_arg.png", all_together_ardR)

transp_ard = subset_samples(transp, sampling_point == "Ardley")

ardT_tse = makeTreeSEFromPhyloseq(transp_ard)
ardT_tse = transformAssay(ardT_tse, # Bray Curtis, relativas
                            assay.type = "counts",
                            method = "relabundance",
                            pseudocount = TRUE)
ardT_tse = runRDA(ardT_tse,
                    assay.type = "relabundance",
                    formula = assay ~ type_f, 
                    distance = "bray",
                    na.action = na.exclude)

bc_transp_matrix = as.matrix(phyloseq::distance(transp_ard, method = "bray"))
transp_dbrda_matrix <- as.matrix(dist(reducedDim(ardT_tse, "RDA")))
stress_transp = sum((transp_dbrda_matrix - bc_transp_matrix)^2) / sum(bc_transp_matrix^2)
stress_transp


png("dbRDA_mge.png", width = 6, height = 4, units = "in", res = 300)
plotRDA(ardT_tse, "RDA", colour_by = "type_f",
        point_size = 5, add.ellipse = "colour",
        arrow.size = 0.50, vec.size = 0.5, label.size = 3, vec.text = FALSE) +
  scale_color_manual(values = palette) +
  theme(aspect.ratio = 1)
invisible(dev.off())

dbrda_ardT = plotRDA(ardT_tse, "RDA", colour_by = "type_f", 
                     point_size = 6, add.ellipse = "colour", ellipse.linewidth = 1,
                     add.vectors = FALSE) +
  scale_color_manual(values = palette) +
  theme_bw() +
  theme(aspect.ratio = 1, legend.position = "none",
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 12)) +
  labs(x = "dbRDA 1", y = "dbRDA 2") 

all_together_ardT = (ard_pcoa_transp | ard_nmds_transp)/ (dbrda_ardT | dbrda_ardT)
print(all_together_ardT)
ggsave("all_joined_MGE.png", all_together_ardT)

# Now Lakes
resis_lake = subset_samples(resis, sampling_point == "Lakes")

lakeR_tse = makeTreeSEFromPhyloseq(resis_lake)
lakeR_tse = transformAssay(lakeR_tse, # Bray Curtis, relativas
                          assay.type = "counts",
                          method = "relabundance",
                          pseudocount = TRUE)
lakeR_tse = runRDA(lakeR_tse,
                  assay.type = "relabundance",
                  formula = assay ~ Sample..type + type_f,
                  distance = "bray",
                  na.action = na.exclude)
bc_res_matrix = as.matrix(phyloseq::distance(resis_lake, method = "bray"))
resis_dbrda_matrix <- as.matrix(dist(reducedDim(lakeR_tse, "RDA")))
stress_res = sum((resis_dbrda_matrix - bc_res_matrix)^2) / sum(bc_res_matrix^2)
stress_res


setwd("../Lakes")

png("RDA.png", width = 6, height = 4, units = "in", res = 300)
plotRDA(lakeR_tse, "RDA", colour_by = "type_f", 
        shape_by = "Sampling..area", point_size = 5, add.ellipse = "colour",
        arrow.size = 0.50, vec.size = 0.5, label.size = 3, vec.text = FALSE) +
  scale_color_manual(values = palette) +
  theme(aspect.ratio = 1)
invisible(dev.off())

dbrda_lakeR = plotRDA(lakeR_tse, "RDA", colour_by = "type_f", 
                     point_size = 6, add.ellipse = "colour", ellipse.linewidth = 1,
                     add.vectors = FALSE) +
  scale_color_manual(values = palette) +
  theme_bw() +
  theme(aspect.ratio = 1, legend.position = "none",
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 12)) +
  labs(x = "dbRDA 1", y = "dbRDA 2") 

all_together_lakeR = (lake_pcoa_res | lake_nmds_res)/ (dbrda_lakeR | dbrda_lakeR)
print(all_together_lakeR)
ggsave("all_joined_ARG.png", all_together_lakeR)


transp_lake = subset_samples(transp, sampling_point == "Lakes")

lakeT_tse = makeTreeSEFromPhyloseq(transp_lake)
lakeT_tse = transformAssay(lakeT_tse, # Bray Curtis, relativas
                          assay.type = "counts",
                          method = "relabundance",
                          pseudocount = TRUE)
lakeT_tse = runRDA(lakeT_tse,
                  assay.type = "relabundance",
                  formula = assay ~ Sample..type + type_f, # + type_g:sam_area,
                  distance = "bray",
                  na.action = na.exclude)

bc_transp_matrix = as.matrix(phyloseq::distance(transp_lake, method = "bray"))
transp_dbrda_matrix <- as.matrix(dist(reducedDim(lakeT_tse, "RDA")))
stress_transp = sum((transp_dbrda_matrix - bc_transp_matrix)^2) / sum(bc_transp_matrix^2)
stress_transp

png("RDA_mge.png", width = 6, height = 4, units = "in", res = 300)
plotRDA(lakeT_tse, "RDA", colour_by = "type_f", 
        shape_by = "Sampling..area", point_size = 5, add.ellipse = "colour",
        arrow.size = 0.50, vec.size = 0.5, label.size = 3, vec.text = FALSE) +
  scale_color_manual(values = palette) +
  theme(aspect.ratio = 1)
invisible(dev.off())
dbrda_lakeT = plotRDA(lakeT_tse, "RDA", colour_by = "type_f", 
                     point_size = 6, add.ellipse = "colour", ellipse.linewidth = 1,
                     add.vectors = FALSE) +
  scale_color_manual(values = palette) +
  theme_bw() +
  theme(aspect.ratio = 1, legend.position = "none",
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 12)) +
  labs(x = "dbRDA 1", y = "dbRDA 2") 
all_together_lakeT = (lake_pcoa_transp | lake_nmds_transp)/ (dbrda_lakeT | dbrda_lakeT)
print(all_together_lakeT)
ggsave("all_joined_MGE.png", all_together_lakeT)
