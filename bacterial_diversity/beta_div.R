#### BETA DIVERSITY ####

#### LIBRARIES ####
# frameworks
library("mia")
library("phyloseq")
library("dplyr")
library("tibble")
# plotting 
library("ggplot2")
library("patchwork")
library("miaViz")
library("ggrepel")
# stats
library("vegan")
library("dendextend")
library("GUniFrac")
library("pairwiseAdonis")
#library("ecole")

#### DIRECTORIOS ####
script_dir = getwd()
setwd("../data/r_backup")
data_dir = getwd()
setwd("../../results/beta")
beta_dir = getwd()

#### RAREFACTAR ####
setwd(data_dir)
ph = readRDS("phystuff_1")

ph@sam_data$sampling_point[ph@sam_data$sampling_point == "Uru"] = "Uruguay"
ph@sam_data$sampling_point[ph@sam_data$sampling_point == "Ion"] = "Ionosférico"
ph@sam_data$sampling_point[ph@sam_data$sampling_point == "Ard"] = "Ardley"

ph@sam_data$type_g[ph@sam_data$type_g == "plastic"] = "Plastic"
ph@sam_data$type_g[ph@sam_data$type_g == "control"] = "Surrounding env."

ph@sam_data$type_f[ph@sam_data$type_f == "Suelo"] = "Soil"
ph@sam_data$type_f[ph@sam_data$type_f == "Agua"] = "Water"

ph@sam_data$tf.sp = paste(ph@sam_data$sampling_point, ph@sam_data$type_f, sep = " - ")  

# I need to make a new variable to group lakes together
ph@sam_data$spf = ph@sam_data$sampling_point
ph@sam_data$spf[ph@sam_data$sampling_point == "Uruguay"] = "Lakes"
ph@sam_data$spf[ph@sam_data$sampling_point == "Ionosférico"] = "Lakes"

# And recreate the tf.sp variable in it's new form
ph@sam_data$tf.spf = paste(ph@sam_data$spf, ph@sam_data$type_f, sep = " - ")
ph@sam_data$tg.spf = paste(ph@sam_data$spf, ph@sam_data$type_g, sep = " - ")

phtse = convertFromPhyloseq(ph)

tab10_cols = c("#ff7f0e", "#1f77b4")
chosen4_cols = c("#B22C2C", "#85B22C", "#B26F2C", "#2C85B2")

# I'm going to rarefact the counts to account for uneven sequencing depth
rarefaction_depth = assay(phtse, "counts") %>% 
                    colSums() %>% 
                    min() 

phtse = rarefyAssay(phtse, 
                     assay.type = "counts",
                     sample = rarefaction_depth,
                     replace = FALSE,
                     name = "rarefied")

# For curiosity's sake, we can check which zOTUs went away
print(setdiff(row.names(assay(phtse, "counts")), 
              row.names(assay(phtse, "rarefied"))
              )
      ) # No zOTUs were removed by rarefaction

# Change the counts in the original phyloseq object
otu_table(ph) = otu_table(assay(phtse, "rarefied"), taxa_are_rows = TRUE)


#### OBTAIN DISTANCES ####
# First the phylogenetic distances
ph@phy_tree = phangorn::midpoint(ph@phy_tree) # the tree needs to be rooted
gunifracs = GUniFrac(t(otu_table(ph)), 
                     phy_tree(ph), 
                     alpha = c(0, 0.5, 1),
                     verbose = TRUE)
dist_uni = as.dist(gunifracs$unifracs[, , "d_UW"])       # Unweighted UniFrac
dist_wuni = as.dist(gunifracs$unifracs[, , "d_1"])        # Weighted UniFrac
dist_guni = as.dist(gunifracs$unifracs[, , "d_0.5"])      # GUniFrac with alpha 0.5
dist_bc = distance(ph, "bray")     # Bray-Curtis

phtse = addDissimilarity(phtse,
                          method = "robust.aitchison",
                          assay.type = "rarefied",
                          niter = 100)
dist_ait = metadata(phtse)[["robust.aitchison"]]
colnames(dist_ait) = colnames(otu_table(ph))
rownames(dist_ait) = colnames(dist_ait)
dist_ait = as.dist(dist_ait)
dist_ait


#### BETADISPERSION ####
setwd(beta_dir)
# Unifrac
beta_disp = betadisper(dist_uni, # distancia a emplear
                       group = data.frame(sample_data(ph))$type_g) # cómo agrupar las muestras
write.csv(permutest(beta_disp, permutations = 10000)$tab, "betadis_uni_tg.csv")
beta_disp = betadisper(dist_uni, # distancia a emplear
                       group = data.frame(sample_data(ph))$type_f) # cómo agrupar las muestras
write.csv(permutest(beta_disp, permutations = 10000)$tab, "betadis_uni_tf.csv")
# GUnifrac
beta_disp = betadisper(dist_guni, 
                       group = data.frame(sample_data(ph))$type_g)
write.csv(permutest(beta_disp, permutations = 10000)$tab, "betadis_guni_tg.csv")
beta_disp = betadisper(dist_guni, 
                       group = data.frame(sample_data(ph))$type_f)
write.csv(permutest(beta_disp, permutations = 10000)$tab, "betadis_guni_tf.csv")
# WUnifrac
write.csv(permutest(betadisper(dist_wuni, group = data.frame(sample_data(ph))$type_g),
                    permutations = 10000)$tab,"betadis_wuni_tg.csv")
write.csv(permutest(betadisper(dist_wuni, group = data.frame(sample_data(ph))$type_f),
                    permutations = 10000)$tab,"betadis_wuni_tf.csv")
# Bray-Curtis
beta_disp = betadisper(dist_bc, group = data.frame(sample_data(ph))$type_g) %>%
            permutest(, permutations = 10000)
write.csv(beta_disp$tab, "betadis_bc_tg.csv")
beta_disp = betadisper(dist_bc, group = data.frame(sample_data(ph))$type_f) %>%
  permutest(, permutations = 10000)
write.csv(beta_disp$tab, "betadis_bc_tf.csv")
# aitchison
beta_disp = betadisper(dist_ait, 
                       group = data.frame(sample_data(ph))$type_g)
write.csv(permutest(beta_disp, permutations = 10000)$tab, "betadis_ait_tg.csv")
beta_disp = betadisper(dist_uni, 
                       group = data.frame(sample_data(ph))$type_f)
write.csv(permutest(beta_disp, permutations = 10000)$tab, "betadis_ait_tf.csv")

# The results are:
# - BC is good in all cases
# - uni only in type_g
# - Guni in all
# - Uni only in tf
# - Wuni in all

#### PERMANOVA ####
# I'm using the global permanova from GUnifrac to run all at once and not lose statistical power
dists = list(dist_ait, dist_bc, dist_uni, dist_wuni, dist_guni)
perma = PermanovaG2(dists ~ type_g, data = data.frame(sample_data(ph)))
# PermanovaG2 da resultados de cada matriz de distancia por separado y luego uno conjunto, global. Nos interesa quedarnoslo todo
write.csv(perma$aov.tab.list[[1]], "perma_ait_tg.csv")
write.csv(perma$aov.tab.list[[2]], "perma_bc_tg.csv")
write.csv(perma$aov.tab.list[[3]], "perma_uni_tg.csv")
write.csv(perma$aov.tab.list[[4]], "perma_wuni_tg.csv")
write.csv(perma$aov.tab.list[[5]], "perma_guni_tg.csv")
write.csv(perma$p.tab, "perma_all_tg.csv")

# Now, type_f
perma = PermanovaG2(dists ~ type_f, data = data.frame(sample_data(ph)))
write.csv(perma$aov.tab.list[[1]], "perma_ait_tf.csv")
write.csv(perma$aov.tab.list[[2]], "perma_bc_tf.csv")
write.csv(perma$aov.tab.list[[3]], "perma_uni_tf.csv")
write.csv(perma$aov.tab.list[[4]], "perma_wuni_tf.csv")
write.csv(perma$aov.tab.list[[5]], "perma_guni_tf.csv")
write.csv(perma$p.tab, "perma_all_tf.csv")

# Let's check only sampling_point
perma = PermanovaG2(dists ~ spf, data = data.frame(sample_data(ph)))
write.csv(perma$aov.tab.list[[1]], "perma_ait_spf.csv")
write.csv(perma$aov.tab.list[[2]], "perma_bc_spf.csv")
write.csv(perma$aov.tab.list[[3]], "perma_uni_spf.csv")
write.csv(perma$aov.tab.list[[4]], "perma_wuni_spf.csv")
write.csv(perma$aov.tab.list[[5]], "perma_guni_spf.csv")
write.csv(perma$p.tab, "perma_all_spf.csv")

# And let's combine both factors
perma = PermanovaG2(dists ~ tg.spf, data = data.frame(sample_data(ph))) # at tg level
write.csv(perma$aov.tab.list[[1]], "perma_ait_tg-spf.csv")
write.csv(perma$aov.tab.list[[2]], "perma_bc_tg-spf.csv")
write.csv(perma$aov.tab.list[[3]], "perma_uni_tg-spf.csv")
write.csv(perma$aov.tab.list[[4]], "perma_wuni_tg-spf.csv")
write.csv(perma$aov.tab.list[[5]], "perma_guni_tg-spf.csv")
write.csv(perma$p.tab, "perma_all_tg-spf.csv")

perma = PermanovaG2(dists ~ tf.spf, data = data.frame(sample_data(ph))) # at tf level
write.csv(perma$aov.tab.list[[1]], "perma_ait_tf-spf.csv")
write.csv(perma$aov.tab.list[[2]], "perma_bc_tf-spf.csv")
write.csv(perma$aov.tab.list[[3]], "perma_uni_tf-spf.csv")
write.csv(perma$aov.tab.list[[4]], "perma_wuni_tf-spf.csv")
write.csv(perma$aov.tab.list[[5]], "perma_guni_tf-spf.csv")
write.csv(perma$p.tab, "perma_all_tf-spf.csv")

#### PERMANOVA PAIRWISE ####
# Once we have differences, we need to weed out where those differences are
# Bray-Curtis
pp_transp = pairwise.adonis(dist_bc, factors = ph@sam_data$type_f, p.adjust.m = "fdr", perm = 10000)
write.csv(pp_transp, "pp_tf_bc.csv")
pp_transp = pairwise.adonis(dist_bc, factors = ph@sam_data$tg.spf, p.adjust.m = "fdr", perm = 10000)
write.csv(pp_transp, "pp_tg-spf_bc.csv")
pp_transp = pairwise.adonis(dist_bc, factors = ph@sam_data$tf.spf, p.adjust.m = "fdr", perm = 10000)
write.csv(pp_transp, "pp_tf-spf_bc.csv")

# aitchison
pp_transp = pairwise.adonis(dist_ait, factors = ph@sam_data$type_f, p.adjust.m = "fdr", perm = 10000)
write.csv(pp_transp, "pp_tf_ait.csv")
pp_transp = pairwise.adonis(dist_ait, factors = ph@sam_data$tg.spf, p.adjust.m = "fdr", perm = 10000)
write.csv(pp_transp, "pp_tg-spf_ait.csv")
pp_transp = pairwise.adonis(dist_ait, factors = ph@sam_data$tf.spf, p.adjust.m = "fdr", perm = 10000)
write.csv(pp_transp, "pp_tf-spf_ait.csv")

# Unifrac
pp_transp = pairwise.adonis(dist_uni, factors = ph@sam_data$type_f, p.adjust.m = "fdr", perm = 10000)
write.csv(pp_transp, "pp_tf_uni.csv")
pp_transp = pairwise.adonis(dist_uni, factors = ph@sam_data$tg.spf, p.adjust.m = "fdr", perm = 10000)
write.csv(pp_transp, "pp_tg-spf_uni.csv")
pp_transp = pairwise.adonis(dist_uni, factors = ph@sam_data$tf.spf, p.adjust.m = "fdr", perm = 10000)
write.csv(pp_transp, "pp_tf-spf_uni.csv")

# WUnifrac
pp_transp = pairwise.adonis(dist_wuni, factors = ph@sam_data$type_f, p.adjust.m = "fdr", perm = 10000)
write.csv(pp_transp, "pp_tf_wuni.csv")
pp_transp = pairwise.adonis(dist_wuni, factors = ph@sam_data$tg.spf, p.adjust.m = "fdr", perm = 10000)
write.csv(pp_transp, "pp_tg-spf_wuni.csv")
pp_transp = pairwise.adonis(dist_wuni, factors = ph@sam_data$tf.spf, p.adjust.m = "fdr", perm = 10000)
write.csv(pp_transp, "pp_tf-spf_wuni.csv")

# GUniFrac
pp_transp = pairwise.adonis(dist_guni, factors = ph@sam_data$type_f, p.adjust.m = "fdr", perm = 10000)
write.csv(pp_transp, "pp_tf_guni.csv")
pp_transp = pairwise.adonis(dist_guni, factors = ph@sam_data$tg.spf, p.adjust.m = "fdr", perm = 10000)
write.csv(pp_transp, "pp_tg-spf_guni.csv")
pp_transp = pairwise.adonis(dist_guni, factors = ph@sam_data$tf.spf, p.adjust.m = "fdr", perm = 10000)
write.csv(pp_transp, "pp_tf-spf_guni.csv")

#### ORDINATE BRAY-CURTIS ####
bc_nmds = ordinate(ph, method = "NMDS", distance = dist_bc)
bc_nmds$stress
bc_nmds_plot_L = plot_ordination(ph, bc_nmds, type = "samples", color = "type_f", shape = "sampling_point", 
                               title = "BRAY-CURTIS NMDS") +
  scale_color_manual(values = chosen4_cols) + 
  geom_point(size=6, alpha=0.55) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1) # This one has a legend
bc_nmds_plot = plot_ordination(ph, bc_nmds, type = "samples", color = "type_f", shape = "sampling_point", 
                                 title = "Bray-Curtis NMDS") +
  scale_color_manual(values = chosen4_cols) + 
  geom_point(size=6, alpha=0.55) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1, legend.position = "none") #this one does not



bc_pcoa = ordinate(ph, method = "PCoA", distance = dist_bc)
resis_pcoa_matrix <- as.matrix(dist(bc_pcoa$vectors))
stress_res = sum((resis_pcoa_matrix - as.matrix(dist_bc)^2) / sum(as.matrix(dist_bc))^2)
stress_res
bc_pcoa_plot = plot_ordination(ph, bc_pcoa, type = "samples", color = "type_f", shape = "sampling_point", 
                title = "Bray-Curtis PCoA") +
  scale_color_manual(values = chosen4_cols) + 
  geom_point(size=6, alpha=0.55) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1, legend.position = "none")


metadata = as.data.frame(data.matrix(ph@sam_data)) # Necesario hacer esto de antemano para que vegan reconozca el objeto
bc_dbrda = dbrda(dist_bc ~ type_f + spf + type_f:spf, data = metadata) #calcular el dbrda
scores_dbrda = scores(bc_dbrda) #sacar los valores de los puntos
bip = scores_dbrda$biplot #y de los vectores
sites_dbrda = as.data.frame(scores_dbrda$sites) # sacar los valores de los sitios

gg_tbl = tibble(sites_dbrda$dbRDA1, sites_dbrda$dbRDA2, 
                ph@sam_data$type_f, ph@sam_data$spf, 
                1:nrow(ph@sam_data), row.names(ph@sam_data))
colnames(gg_tbl) = c("X", "Y", "colour_by", "shape_by", "order_by", "rows")
gg_tbl = column_to_rownames(gg_tbl, var = "rows")

gg_tbl = rownames_to_column(gg_tbl, var = "label") %>% # cambiamos los rownames a columna para no perder la info de qué muestra es cuál
  add_column(plot_type = "sites", .before = "label") %>% # añadimos una columna para saber qué datos se tienen que usar en qué plot
  add_row(label = row.names(bip)[1], plot_type = "biplot", #añadimos la fila de la flecha de type_g
          X = bip[row.names(bip)[1], colnames(bip)[1]], Y = bip[row.names(bip)[1], colnames(bip)[2]]) %>%
  add_row(label = row.names(bip)[2], plot_type = "biplot", #añadimos la fila de la flecha de Year
          X = bip[row.names(bip)[2], colnames(bip)[1]], Y = bip[row.names(bip)[2], colnames(bip)[2]]) %>%
  add_row(label = row.names(bip)[3], plot_type = "biplot", #añadimos la fila de la flecha de Year
          X = bip[row.names(bip)[3], colnames(bip)[1]], Y = bip[row.names(bip)[2], colnames(bip)[2]])

point_mask = gg_tbl$plot_type == "sites"
biplot_mask = gg_tbl$plot_type == "biplot"

perma_model = as.data.frame(anova.cca(bc_dbrda, by = NULL)) 
perma_margin = as.data.frame(anova.cca(bc_dbrda, by = "terms", model = "direct")) 
perma_all = rbind(perma_model[1, ], perma_margin)
perma_all[ , "Total variance"] <- perma_all["Model", 2] + perma_all["Residual", 2]
perma_all[ , "Explained variance"] <- perma_all[ , 2] / perma_all[ , "Total variance"]

for (var_name in row.names(perma_all)) {
  if (var_name %in% gg_tbl$label){
    var_data = perma_all[var_name, ]
    explained_var = round(var_data[, "Explained variance"] * 100, digits = 1) # paso a pct y redondeo por estética
    p_val = round(var_data[, "Pr(>F)"], digits = 3) # redondeo por estética
    final_label_text = paste0(var_name, " (", explained_var, ", p = ", p_val, ")")
    gg_tbl = gg_tbl %>%
      mutate(label = replace(label, label == var_name, final_label_text))
  }
} 

bc_dbrda_plot = ggplot() + 
  geom_point(data = gg_tbl[point_mask, ],
             aes(x = X, y = Y, color = colour_by, shape = shape_by, size=8, alpha=0.55)) +
  geom_point(data = gg_tbl[point_mask, ],
             aes(x = X, y = Y, color = colour_by, shape = shape_by, size=1)) +
  #geom_segment(data = gg_tbl[biplot_mask, ], 
  #             aes(x = 0, y = 0, xend = X, yend = Y), 
  #             arrow = arrow(length = unit(0.5,"cm"))) +
  #geom_label_repel(data = gg_tbl[biplot_mask, ], # cambia el comando a geom_label_repel()
  #                 aes(x = X, y = Y, label = label),
  #                 fill = "white", # rellena la caja de blanco
  #                 xlim = c(-Inf, Inf), ylim = c(-Inf, Inf)) +
  scale_color_manual(values = chosen4_cols) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1, legend.position = "none") +
  labs(title ="Bray-Curtis dbRDA",
       x = "dbRDA1", y = "dbRDA2") #+
  #scale_x_continuous(expand = expansion(mult = 0.5))
perma_margin #check what is and isn't explained by the model

# Stitch everything together
three_joined = (bc_nmds_plot | bc_pcoa_plot)/(bc_dbrda_plot | bc_nmds_plot_L)
three_joined
ggsave(filename = "BC_joined.tif", plot = three_joined, device = "tiff",
       width = 10, height = 8, dpi = 300)

#### ORDINATE AITCHISON ####
ait_nmds = ordinate(ph, method = "NMDS", distance = dist_ait)
ait_nmds$stress
ait_nmds_plot_L = plot_ordination(ph, ait_nmds, type = "samples", color = "type_f", shape = "sampling_point", 
                title = "AITCHISON NMDS") +
  scale_color_manual(values = chosen4_cols) + 
  geom_point(size=6, alpha=0.55) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1)
ait_nmds_plot= plot_ordination(ph, ait_nmds, type = "samples", color = "type_f", shape = "sampling_point", 
                             title = "Aitchison NMDS") +
  scale_color_manual(values = chosen4_cols) + 
  geom_point(size=6, alpha=0.55) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1, legend.position = "none")


ait_pcoa = ordinate(ph, method = "PCoA", distance = dist_ait)
resis_pcoa_matrix <- as.matrix(dist(ait_pcoa$vectors))
stress_res = sum((resis_pcoa_matrix - as.matrix(dist_ait)^2) / sum(as.matrix(dist_ait))^2)
stress_res
ait_pcoa_plot = plot_ordination(ph, ait_pcoa, type = "samples", color = "type_f", shape = "sampling_point", 
                title = "Aitchison PCoA") +
  scale_color_manual(values = chosen4_cols) + 
  geom_point(size=6, alpha=0.55) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1, legend.position = "none")


metadata = as.data.frame(data.matrix(ph@sam_data)) # Necesario hacer esto de antemano para que vegan reconozca el objeto
ait_dbrda = dbrda(dist_ait ~ type_f + spf + type_f:spf, data = metadata) #calcular el dbrda
scores_dbrda = scores(ait_dbrda) #sacar los valores de los puntos
bip = scores_dbrda$biplot #y de los vectores
sites_dbrda = as.data.frame(scores_dbrda$sites) # sacar los valores de los sitios

gg_tbl = tibble(sites_dbrda$dbRDA1, sites_dbrda$dbRDA2, 
                ph@sam_data$type_f, ph@sam_data$spf, 
                1:nrow(ph@sam_data), row.names(ph@sam_data))
colnames(gg_tbl) = c("X", "Y", "colour_by", "shape_by", "order_by", "rows")
gg_tbl = column_to_rownames(gg_tbl, var = "rows")

gg_tbl = rownames_to_column(gg_tbl, var = "label") %>% # cambiamos los rownames a columna para no perder la info de qué muestra es cuál
  add_column(plot_type = "sites", .before = "label") %>% # añadimos una columna para saber qué datos se tienen que usar en qué plot
  add_row(label = row.names(bip)[1], plot_type = "biplot", #añadimos la fila de la flecha de type_g
          X = bip[row.names(bip)[1], colnames(bip)[1]], Y = bip[row.names(bip)[1], colnames(bip)[2]]) %>%
  add_row(label = row.names(bip)[2], plot_type = "biplot", #añadimos la fila de la flecha de Year
          X = bip[row.names(bip)[2], colnames(bip)[1]], Y = bip[row.names(bip)[2], colnames(bip)[2]]) %>%
  add_row(label = row.names(bip)[3], plot_type = "biplot", #añadimos la fila de la flecha de Year
          X = bip[row.names(bip)[3], colnames(bip)[1]], Y = bip[row.names(bip)[2], colnames(bip)[2]])

point_mask = gg_tbl$plot_type == "sites"
biplot_mask = gg_tbl$plot_type == "biplot"

perma_model = as.data.frame(anova.cca(ait_dbrda, by = NULL)) 
perma_margin = as.data.frame(anova.cca(ait_dbrda, by = "terms", model = "direct")) 
perma_all = rbind(perma_model[1, ], perma_margin)
perma_all[ , "Total variance"] <- perma_all["Model", 2] + perma_all["Residual", 2]
perma_all[ , "Explained variance"] <- perma_all[ , 2] / perma_all[ , "Total variance"]

for (var_name in row.names(perma_all)) {
  if (var_name %in% gg_tbl$label){
    var_data = perma_all[var_name, ]
    explained_var = round(var_data[, "Explained variance"] * 100, digits = 1) # paso a pct y redondeo por estética
    p_val = round(var_data[, "Pr(>F)"], digits = 3) # redondeo por estética
    final_label_text = paste0(var_name, " (", explained_var, ", p = ", p_val, ")")
    gg_tbl = gg_tbl %>%
      mutate(label = replace(label, label == var_name, final_label_text))
  }
} 

ait_dbrda_plot = ggplot() + 
  geom_point(data = gg_tbl[point_mask, ],
             aes(x = X, y = Y, color = colour_by, shape = shape_by, size=8, alpha=0.55)) +
  geom_point(data = gg_tbl[point_mask, ],
             aes(x = X, y = Y, color = colour_by, shape = shape_by, size=1)) +
  #geom_segment(data = gg_tbl[biplot_mask, ], 
  #             aes(x = 0, y = 0, xend = X, yend = Y), 
  #             arrow = arrow(length = unit(0.5,"cm"))) +
  #geom_label_repel(data = gg_tbl[biplot_mask, ], # cambia el comando a geom_label_repel()
  #                 aes(x = X, y = Y, label = label),
  #                 fill = "white", # rellena la caja de blanco
  #                 xlim = c(-Inf, Inf), ylim = c(-Inf, Inf)) +
  scale_color_manual(values = chosen4_cols) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1, legend.position = "none") +
  labs(title ="AITCHISON dbRDA",
       x = "dbRDA1", y = "dbRDA2") #+
  #scale_x_continuous(expand = expansion(mult = 0.1))
gg_tbl

# Stitch everything together
three_joined = (ait_nmds_plot | ait_pcoa_plot)/(ait_dbrda_plot | ait_nmds_plot_L)
three_joined
ggsave(filename = "AIT_joined.tif", plot = three_joined, device = "tiff",
       width = 10, height = 8, dpi = 300)

#### ORDINATE UNIFRAC ####
uni_nmds = ordinate(ph, method = "NMDS", distance = dist_uni)
uni_nmds$stress
uni_nmds_plot_L = plot_ordination(ph, uni_nmds, type = "samples", color = "type_f", shape = "sampling_point", 
                title = "UniFrac NMDS") +
  scale_color_manual(values = chosen4_cols) + 
  geom_point(size=6, alpha=0.55) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1)
uni_nmds_plot= plot_ordination(ph, uni_nmds, type = "samples", color = "type_f", shape = "sampling_point", 
                                  title = "UniFrac NMDS") +
  scale_color_manual(values = chosen4_cols) + 
  geom_point(size=6, alpha=0.55) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1, legend.position = "none")


uni_pcoa = ordinate(ph, method = "PCoA", distance = dist_uni)
resis_pcoa_matrix <- as.matrix(dist(uni_pcoa$vectors))
stress_res = sum((resis_pcoa_matrix - as.matrix(dist_uni)^2) / sum(as.matrix(dist_uni))^2)
stress_res
uni_pcoa_plot = plot_ordination(ph, uni_pcoa, type = "samples", color = "type_f", shape = "sampling_point", 
                title = "UniFrac PCoA") +
  scale_color_manual(values = chosen4_cols) + 
  geom_point(size=6, alpha=0.55) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1, legend.position = "none")


metadata = as.data.frame(data.matrix(ph@sam_data)) # Necesario hacer esto de antemano para que vegan reconozca el objeto
uni_dbrda = dbrda(dist_uni ~ type_f + spf + type_f:spf, data = metadata) #calcular el dbrda
scores_dbrda = scores(uni_dbrda) #sacar los valores de los puntos
bip = scores_dbrda$biplot #y de los vectores
sites_dbrda = as.data.frame(scores_dbrda$sites) # sacar los valores de los sitios

gg_tbl = tibble(sites_dbrda$dbRDA1, sites_dbrda$dbRDA2, 
                ph@sam_data$type_f, ph@sam_data$spf, 
                1:nrow(ph@sam_data), row.names(ph@sam_data))
colnames(gg_tbl) = c("X", "Y", "colour_by", "shape_by", "order_by", "rows")
gg_tbl = column_to_rownames(gg_tbl, var = "rows")

gg_tbl = rownames_to_column(gg_tbl, var = "label") %>% # cambiamos los rownames a columna para no perder la info de qué muestra es cuál
  add_column(plot_type = "sites", .before = "label") %>% # añadimos una columna para saber qué datos se tienen que usar en qué plot
  add_row(label = row.names(bip)[1], plot_type = "biplot", #añadimos la fila de la flecha de type_g
          X = bip[row.names(bip)[1], colnames(bip)[1]], Y = bip[row.names(bip)[1], colnames(bip)[2]]) %>%
  add_row(label = row.names(bip)[2], plot_type = "biplot", #añadimos la fila de la flecha de Year
          X = bip[row.names(bip)[2], colnames(bip)[1]], Y = bip[row.names(bip)[2], colnames(bip)[2]]) %>%
  add_row(label = row.names(bip)[3], plot_type = "biplot", #añadimos la fila de la flecha de Year
          X = bip[row.names(bip)[3], colnames(bip)[1]], Y = bip[row.names(bip)[2], colnames(bip)[2]])

point_mask = gg_tbl$plot_type == "sites"
biplot_mask = gg_tbl$plot_type == "biplot"

perma_model = as.data.frame(anova.cca(uni_dbrda, by = NULL)) 
perma_margin = as.data.frame(anova.cca(uni_dbrda, by = "terms", model = "direct")) 
perma_all = rbind(perma_model[1, ], perma_margin)
perma_all[ , "Total variance"] <- perma_all["Model", 2] + perma_all["Residual", 2]
perma_all[ , "Explained variance"] <- perma_all[ , 2] / perma_all[ , "Total variance"]

for (var_name in row.names(perma_all)) {
  if (var_name %in% gg_tbl$label){
    var_data = perma_all[var_name, ]
    explained_var = round(var_data[, "Explained variance"] * 100, digits = 1) # paso a pct y redondeo por estética
    p_val = round(var_data[, "Pr(>F)"], digits = 3) # redondeo por estética
    final_label_text = paste0(var_name, " (", explained_var, ", p = ", p_val, ")")
    gg_tbl = gg_tbl %>%
      mutate(label = replace(label, label == var_name, final_label_text))
  }
} 

uni_dbrda_plot = ggplot() + 
  geom_point(data = gg_tbl[point_mask, ],
             aes(x = X, y = Y, color = colour_by, shape = shape_by, size=8, alpha=0.55)) +
  geom_point(data = gg_tbl[point_mask, ],
             aes(x = X, y = Y, color = colour_by, shape = shape_by, size=1)) +
  scale_color_manual(values = chosen4_cols) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1, legend.position = "none") +
  labs(title ="UniFrac dbRDA",
       x = "dbRDA1", y = "dbRDA2")
gg_tbl

# Stitch everything together
three_joined = (uni_nmds_plot | uni_pcoa_plot)/(uni_dbrda_plot | uni_nmds_plot_L)
three_joined
ggsave(filename = "UNI_joined.tif", plot = three_joined, device = "tiff",
       width = 10, height = 8, dpi = 300)

#### ORDINATE WUNIFRAC ####
wuni_nmds = ordinate(ph, method = "NMDS", distance = dist_wuni)
wuni_nmds$stress
wuni_nmds_plot_L = plot_ordination(ph, wuni_nmds, type = "samples", color = "type_f", shape = "sampling_point", 
                title = "wuniFrac NMDS") +
  scale_color_manual(values = chosen4_cols) + 
  geom_point(size=6, alpha=0.55) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1)
wuni_nmds_plot = plot_ordination(ph, wuni_nmds, type = "samples", color = "type_f", shape = "sampling_point", 
                                   title = "wuniFrac NMDS") +
  scale_color_manual(values = chosen4_cols) + 
  geom_point(size=6, alpha=0.55) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1, legend.position = "none")

wuni_pcoa = ordinate(ph, method = "PCoA", distance = dist_wuni)
resis_pcoa_matrix <- as.matrix(dist(wuni_pcoa$vectors))
stress_res = sum((resis_pcoa_matrix - as.matrix(dist_wuni)^2) / sum(as.matrix(dist_wuni))^2)
stress_res
wuni_pcoa_plot = plot_ordination(ph, wuni_pcoa, type = "samples", color = "type_f", shape = "sampling_point", 
                title = "wuniFrac PCoA") +
  scale_color_manual(values = chosen4_cols) + 
  geom_point(size=6, alpha=0.55) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1, legend.position = "none")

metadata = as.data.frame(data.matrix(ph@sam_data)) # Necesario hacer esto de antemano para que vegan reconozca el objeto
wuni_dbrda = dbrda(dist_wuni ~ type_f + spf + type_f:spf, data = metadata) #calcular el dbrda
scores_dbrda = scores(wuni_dbrda) #sacar los valores de los puntos
bip = scores_dbrda$biplot #y de los vectores
sites_dbrda = as.data.frame(scores_dbrda$sites) # sacar los valores de los sitios

gg_tbl = tibble(sites_dbrda$dbRDA1, sites_dbrda$dbRDA2, 
                ph@sam_data$type_f, ph@sam_data$spf, 
                1:nrow(ph@sam_data), row.names(ph@sam_data))
colnames(gg_tbl) = c("X", "Y", "colour_by", "shape_by", "order_by", "rows")
gg_tbl = column_to_rownames(gg_tbl, var = "rows")

gg_tbl = rownames_to_column(gg_tbl, var = "label") %>% # cambiamos los rownames a columna para no perder la info de qué muestra es cuál
  add_column(plot_type = "sites", .before = "label") %>% # añadimos una columna para saber qué datos se tienen que usar en qué plot
  add_row(label = row.names(bip)[1], plot_type = "biplot", #añadimos la fila de la flecha de type_g
          X = bip[row.names(bip)[1], colnames(bip)[1]], Y = bip[row.names(bip)[1], colnames(bip)[2]]) %>%
  add_row(label = row.names(bip)[2], plot_type = "biplot", #añadimos la fila de la flecha de Year
          X = bip[row.names(bip)[2], colnames(bip)[1]], Y = bip[row.names(bip)[2], colnames(bip)[2]]) %>%
  add_row(label = row.names(bip)[3], plot_type = "biplot", #añadimos la fila de la flecha de Year
          X = bip[row.names(bip)[3], colnames(bip)[1]], Y = bip[row.names(bip)[2], colnames(bip)[2]])

point_mask = gg_tbl$plot_type == "sites"
biplot_mask = gg_tbl$plot_type == "biplot"

perma_model = as.data.frame(anova.cca(wuni_dbrda, by = NULL)) 
perma_margin = as.data.frame(anova.cca(wuni_dbrda, by = "terms", model = "direct")) 
perma_all = rbind(perma_model[1, ], perma_margin)
perma_all[ , "Total variance"] <- perma_all["Model", 2] + perma_all["Residual", 2]
perma_all[ , "Explained variance"] <- perma_all[ , 2] / perma_all[ , "Total variance"]

for (var_name in row.names(perma_all)) {
  if (var_name %in% gg_tbl$label){
    var_data = perma_all[var_name, ]
    explained_var = round(var_data[, "Explained variance"] * 100, digits = 1) # paso a pct y redondeo por estética
    p_val = round(var_data[, "Pr(>F)"], digits = 3) # redondeo por estética
    final_label_text = paste0(var_name, " (", explained_var, ", p = ", p_val, ")")
    gg_tbl = gg_tbl %>%
      mutate(label = replace(label, label == var_name, final_label_text))
  }
} 

wuni_dbrda_plot = ggplot() + 
  geom_point(data = gg_tbl[point_mask, ],
             aes(x = X, y = Y, color = colour_by, shape = shape_by, size=8, alpha=0.55)) +
  geom_point(data = gg_tbl[point_mask, ],
             aes(x = X, y = Y, color = colour_by, shape = shape_by, size=1)) +
  scale_color_manual(values = chosen4_cols) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1, legend.position = "none") +
  labs(title ="wuniFRAC dbRDA",
       x = "dbRDA1", y = "dbRDA2") 
gg_tbl

# Stitch everything together
three_joined = (wuni_nmds_plot | wuni_pcoa_plot)/(wuni_dbrda_plot | wuni_nmds_plot_L)
three_joined
ggsave(filename = "WUNI_joined.tif", plot = three_joined, device = "tiff",
       width = 10, height = 8, dpi = 300)

#### ORDINATE GUNIFRAC ####
guni_nmds = ordinate(ph, method = "NMDS", distance = dist_guni)
guni_nmds$stress
guni_nmds_plot_L = plot_ordination(ph, guni_nmds, type = "samples", color = "type_f", shape = "sampling_point", 
                title = "GUniFrac NMDS") +
  scale_color_manual(values = chosen4_cols) + 
  geom_point(size=6, alpha=0.55) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1)
guni_nmds_plot = plot_ordination(ph, guni_nmds, type = "samples", color = "type_f", shape = "sampling_point", 
                                   title = "GUniFrac NMDS") +
  scale_color_manual(values = chosen4_cols) + 
  geom_point(size=6, alpha=0.55) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1, legend.position = "none")

guni_pcoa = ordinate(ph, method = "PCoA", distance = dist_guni)
resis_pcoa_matrix <- as.matrix(dist(guni_pcoa$vectors))
stress_res = sum((resis_pcoa_matrix - as.matrix(dist_guni)^2) / sum(as.matrix(dist_guni))^2)
stress_res
guni_pcoa_plot = plot_ordination(ph, guni_pcoa, type = "samples", color = "type_f", shape = "sampling_point", 
                title = "GUniFrac PCoA") +
  scale_color_manual(values = chosen4_cols) + 
  geom_point(size=6, alpha=0.55) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1, legend.position = "none")

metadata = as.data.frame(data.matrix(ph@sam_data)) # Necesario hacer esto de antemano para que vegan reconozca el objeto
guni_dbrda = dbrda(dist_guni ~ type_f + spf + type_f:spf, data = metadata) #calcular el dbrda
scores_dbrda = scores(guni_dbrda) #sacar los valores de los puntos
bip = scores_dbrda$biplot #y de los vectores
sites_dbrda = as.data.frame(scores_dbrda$sites) # sacar los valores de los sitios

gg_tbl = tibble(sites_dbrda$dbRDA1, sites_dbrda$dbRDA2, 
                ph@sam_data$type_f, ph@sam_data$spf, 
                1:nrow(ph@sam_data), row.names(ph@sam_data))
colnames(gg_tbl) = c("X", "Y", "colour_by", "shape_by", "order_by", "rows")
gg_tbl = column_to_rownames(gg_tbl, var = "rows")

gg_tbl = rownames_to_column(gg_tbl, var = "label") %>% # cambiamos los rownames a columna para no perder la info de qué muestra es cuál
  add_column(plot_type = "sites", .before = "label") %>% # añadimos una columna para saber qué datos se tienen que usar en qué plot
  add_row(label = row.names(bip)[1], plot_type = "biplot", #añadimos la fila de la flecha de type_g
          X = bip[row.names(bip)[1], colnames(bip)[1]], Y = bip[row.names(bip)[1], colnames(bip)[2]]) %>%
  add_row(label = row.names(bip)[2], plot_type = "biplot", #añadimos la fila de la flecha de Year
          X = bip[row.names(bip)[2], colnames(bip)[1]], Y = bip[row.names(bip)[2], colnames(bip)[2]]) %>%
  add_row(label = row.names(bip)[3], plot_type = "biplot", #añadimos la fila de la flecha de Year
          X = bip[row.names(bip)[3], colnames(bip)[1]], Y = bip[row.names(bip)[2], colnames(bip)[2]])

point_mask = gg_tbl$plot_type == "sites"
biplot_mask = gg_tbl$plot_type == "biplot"

perma_model = as.data.frame(anova.cca(guni_dbrda, by = NULL)) 
perma_margin = as.data.frame(anova.cca(guni_dbrda, by = "terms", model = "direct")) 
perma_all = rbind(perma_model[1, ], perma_margin)
perma_all[ , "Total variance"] <- perma_all["Model", 2] + perma_all["Residual", 2]
perma_all[ , "Explained variance"] <- perma_all[ , 2] / perma_all[ , "Total variance"]

for (var_name in row.names(perma_all)) {
  if (var_name %in% gg_tbl$label){
    var_data = perma_all[var_name, ]
    explained_var = round(var_data[, "Explained variance"] * 100, digits = 1) # paso a pct y redondeo por estética
    p_val = round(var_data[, "Pr(>F)"], digits = 3) # redondeo por estética
    final_label_text = paste0(var_name, " (", explained_var, ", p = ", p_val, ")")
    gg_tbl = gg_tbl %>%
      mutate(label = replace(label, label == var_name, final_label_text))
  }
} 

guni_dbrda_plot = ggplot() + 
  geom_point(data = gg_tbl[point_mask, ],
             aes(x = X, y = Y, color = colour_by, shape = shape_by, size=8, alpha=0.55)) +
  geom_point(data = gg_tbl[point_mask, ],
             aes(x = X, y = Y, color = colour_by, shape = shape_by, size=1)) +
  scale_color_manual(values = chosen4_cols) + 
  theme_bw() + coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1, legend.position = "none") +
  labs(title ="GUniFRAC dbRDA",
       x = "dbRDA1", y = "dbRDA2")

gg_tbl

three_joined = (guni_nmds_plot | guni_pcoa_plot)/(guni_dbrda_plot | guni_nmds_plot_L)
three_joined
ggsave(filename = "GUNI_joined.tif", plot = three_joined, device = "tiff",
       width = 10, height = 8, dpi = 300)