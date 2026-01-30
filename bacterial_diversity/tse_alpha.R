### ALPHA DIVERSITY  ####

#### LIBRERÍAS ####
library("mia")
library("phyloseq")
library("dplyr")
library("MicEco")
library("scater")
library("patchwork")
library("tibble")
library("FSA")
library("ggsignif")
#### DIRECTORIOS ####
script_dir = getwd()
setwd("../data/r_backup")
data_dir = getwd()
setwd("../../results/alpha")
res_dir = getwd()

#### DATA PREPROCESS ####
setwd(data_dir)
ph = readRDS("phystuff_1")

# Once I have the phyloseq objects, I want to fix the naming of various things
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

phtse = convertFromPhyloseq(ph)

# Y ya de paso, defino la paleta aquí al principio
tab10_cols = c("#ff7f0e", "#1f77b4")
chosen4_cols = c("#B22C2C", "#85B22C", "#B26F2C", "#2C85B2")


##### MAKE ALL PLOTS ####
setwd(res_dir)

## Plots
alpha_indices = c("observed", "shannon", "coverage", "relative", "log_modulo_skewness", "faith")
phtse = addAlpha(phtse, # Usa el objeto pstse
                 assay.type = "counts", # Usa los datos guardados en assays -> counts
                 index = alpha_indices, # Calcula este índice de alfa diversidad
                 name = alpha_indices)
View(phtse@colData)

data = as.data.frame(phtse@colData)
View(data)

# Plot without site separation
plots = list()
i = 1
for (alpha_type in alpha_indices) {
  #print(alpha_type)
  plots[[i]] = ggplot(data, aes(x = type_g, y = .data[[alpha_type]], fill = type_g)) +
    geom_jitter(size = 1) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) +
    geom_signif(comparisons = list(c("Plastic", "Surrounding env.")), map_signif_level = TRUE) +
    theme(title = element_text(size = 12),
          axis.title.x = element_blank()) +
    scale_fill_manual(values = tab10_cols) + 
    labs(fill='Sample type')
  print(alpha_type)
  i = i+1
}
plots[[1]]

plots = lapply(plots, "+", 
               theme(axis.title.x = element_blank(),
                     axis.ticks.y = element_blank()))
all_plots = ((plots[[1]] | plots[[2]] | plots[[3]]) / 
               (plots[[4]] | plots[[5]] | plots[[6]])) +
  plot_layout(guides = "collect") + 
  plot_annotation(title = "Different measures of microbiome α diverstity")
all_plots
# Save the plot
ggsave("plas-vs-env-general.png", plot = all_plots)

# Make a plot separating by site
plots = list()
i = 1
for (alpha_type in alpha_indices) {
  #print(alpha_type)
  plots[[i]] = ggplot(data, aes(x = sampling_point, y = .data[[alpha_type]], fill = type_g)) +
    geom_jitter(size = 1) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) +
    theme(title = element_text(size = 12),
          axis.title.x = element_blank()) +
    scale_fill_manual(values = tab10_cols) + 
    labs(fill='Sample type')
  print(alpha_type)
  i = i+1
}
plots[[1]]

plots = lapply(plots, "+", 
                theme(axis.title.x = element_blank(),
                      axis.ticks.y = element_blank()))
all_plots = ((plots[[1]] | plots[[2]] | plots[[3]]) / 
                 (plots[[4]] | plots[[5]] | plots[[6]])) +
  plot_layout(guides = "collect") + 
  plot_annotation(title = "Different measures of microbiome α diverstity")
all_plots
ggsave("All_alpha_tg.png", plot = all_plots)

# Plot without separating by site but separating by sample type
plots = list()
i = 1
for (alpha_type in alpha_indices) {
  #print(alpha_type)
  plots[[i]] = ggplot(data, aes(x = type_f, y = .data[[alpha_type]], fill = type_f)) +
    geom_jitter(size = 1) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) +
    #geom_tukey(where = "whisker") + #this simply desn't work, I belive it has to do with the fact its running on a loop
    theme(title = element_text(size = 12),
          axis.title.x = element_blank()) +
    scale_fill_manual(values = chosen4_cols) + 
    labs(fill='Sample type')
  print(alpha_type)
  i = i+1
}
plots[[1]]

plots = lapply(plots, "+", 
               theme(axis.title.x = element_blank(),
                     axis.ticks.y = element_blank()))
all_plots = ((plots[[1]] | plots[[2]] | plots[[3]]) / 
               (plots[[4]] | plots[[5]] | plots[[6]])) +
  plot_layout(guides = "collect") + 
  plot_annotation(title = "Different measures of microbiome α diverstity")
all_plots
ggsave("plas-vs-env-tf.png", plot = all_plots)

# And now a plot for everything vs everything differentiating by sampling site
plots = list()
i = 1
for (alpha_type in alpha_indices) {
  #print(alpha_type)
  plots[[i]] = ggplot(data, aes(x = sampling_point, y = .data[[alpha_type]], fill = type_f)) +
    geom_jitter(size = 1) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) +
    #geom_tukey(where = "whisker") + #this simply desn't work, I belive it has to do with the fact its running on a loop
    theme(title = element_text(size = 12),
          axis.title.x = element_blank()) +
    scale_fill_manual(values = chosen4_cols) + 
    labs(fill='Sample type')
  print(alpha_type)
  i = i+1
}
plots[[1]]

plots = lapply(plots, "+", 
               theme(axis.title.x = element_blank(),
                     axis.ticks.y = element_blank()))
all_plots = ((plots[[1]] | plots[[2]] | plots[[3]]) / 
               (plots[[4]] | plots[[5]] | plots[[6]])) +
  plot_layout(guides = "collect") + 
  plot_annotation(title = "Different measures of microbiome α diverstity")
all_plots
ggsave("All_alpha-tf.png", plot = all_plots)

# Type_g but with the ard-lakes differentiation
plots = list()
i = 1
for (alpha_type in alpha_indices) {
  #print(alpha_type)
  plots[[i]] = ggplot(data, aes(x = spf, y = .data[[alpha_type]], fill = type_g)) +
    geom_jitter(size = 1) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) +
    theme(title = element_text(size = 12),
          axis.title.x = element_blank()) +
    scale_fill_manual(values = tab10_cols) + 
    labs(fill='Sample type')
  print(alpha_type)
  i = i+1
}
plots[[1]]

plots = lapply(plots, "+", 
               theme(axis.title.x = element_blank(),
                     axis.ticks.y = element_blank()))
all_plots = ((plots[[1]] | plots[[2]] | plots[[3]]) / 
               (plots[[4]] | plots[[5]] | plots[[6]])) +
  plot_layout(guides = "collect") + 
  plot_annotation(title = "Different measures of microbiome α diverstity")
all_plots
ggsave("2sites_alpha_tg.png", plot = all_plots)

# Finally, lake grouping with type_f
plots = list()
i = 1
for (alpha_type in alpha_indices) {
  plots[[i]] = ggplot(data, aes(x = spf, y = .data[[alpha_type]], fill = type_f)) +
    geom_jitter(size = 1) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) +
    theme(title = element_text(size = 12),
          axis.title.x = element_blank()) +
    scale_fill_manual(values = chosen4_cols) + 
    labs(fill='Sample type')
  print(alpha_type)
  i = i+1
}
plots[[1]]

plots = lapply(plots, "+", 
               theme(axis.title.x = element_blank(),
                     axis.ticks.y = element_blank()))
all_plots = ((plots[[1]] | plots[[2]] | plots[[3]]) / 
               (plots[[4]] | plots[[5]] | plots[[6]])) +
  plot_layout(guides = "collect") + 
  plot_annotation(title = "Different measures of microbiome α diverstity")
all_plots
ggsave("lake-tf-spf.png", plot = all_plots)


#### Stats #####
# First, plas vs surrounding env. as is 
sink("tg-no_site.txt")
print("")
sink()
for (alpha_type in alpha_indices) {
  sink("tg-no_site.txt", append = TRUE)
  print(paste0("#################### ", alpha_type, " ####################"))
  wilcoxon_val <- wilcox.test(as.formula(paste0(alpha_type," ~ type_g")), data = colData(phtse))
  print(wilcoxon_val)
  sink()
}

# Then, plas vs surrounding env but differentiating by sampling site
# First I need to make a new variable specific to this question
phtse@colData$tg.sp = paste(phtse@colData$sampling_point, phtse@colData$type_g, sep = " - ")   

# And now we can run the test itself
sink("tg_sp.txt")
print("")
sink()
for (alpha_type in alpha_indices) {
  sink("tg_sp.txt", append = TRUE)
  print(paste0("#################### ", alpha_type, " ####################"))
  kruskal_val = kruskal.test(as.formula(paste0(alpha_type," ~ tg.sp")), data = colData(phtse))
  print(kruskal_val)
  print(" ")
  dunn_val = dunnTest(as.formula(paste0(alpha_type," ~ tg.sp")), data = colData(phtse), method = "bh")
  print(dunn_val)
  sink()
}

# Now differentiating by sample type but not site
sink("tf-no_site.txt")
print("")
sink()
for (alpha_type in alpha_indices) {
  sink("tf-no_site.txt", append = TRUE)
  print(paste0("#################### ", alpha_type, " ####################"))
  kruskal_val = kruskal.test(as.formula(paste0(alpha_type," ~ type_f")), data = colData(phtse))
  print(kruskal_val)
  print(" ")
  dunn_val = dunnTest(as.formula(paste0(alpha_type," ~ type_f")), data = colData(phtse), method = "bh")
  print(dunn_val)
  sink()
}

# And lastly, separating by everything separable
sink("tf_sp.txt")
print("")
sink()
for (alpha_type in alpha_indices) {
  sink("tf_sp.txt", append = TRUE)
  print(paste0("#################### ", alpha_type, " ####################"))
  kruskal_val = kruskal.test(as.formula(paste0(alpha_type," ~ tf.sp")), data = colData(phtse))
  print(kruskal_val)
  print(" ")
  dunn_val = dunnTest(as.formula(paste0(alpha_type," ~ tf.sp")), data = colData(phtse), method = "bh")
  print(dunn_val)
  sink()
}

# Repeating the tests with the lake grouping. Type_g
phtse@colData$tg.spf = paste(phtse@colData$spf, phtse@colData$type_g, sep = " - ")   

sink("spf_tg.txt")
print("")
sink()
for (alpha_type in alpha_indices) {
  sink("spf_tg.txt", append = TRUE)
  print(paste0("#################### ", alpha_type, " ####################"))
  kruskal_val = kruskal.test(as.formula(paste0(alpha_type," ~ tg.spf")), data = colData(phtse))
  print(kruskal_val)
  print(" ")
  dunn_val = dunnTest(as.formula(paste0(alpha_type," ~ tg.spf")), data = colData(phtse), method = "bh")
  print(dunn_val)
  sink()
}

# Repeating with lake grouping, type_f
sink("spf_tf.txt")
print("")
sink()
for (alpha_type in alpha_indices) {
  sink("spf_tf.txt", append = TRUE)
  print(paste0("#################### ", alpha_type, " ####################"))
  kruskal_val = kruskal.test(as.formula(paste0(alpha_type," ~ tf.spf")), data = colData(phtse))
  print(kruskal_val)
  print(" ")
  dunn_val = dunnTest(as.formula(paste0(alpha_type," ~ tf.spf")), data = colData(phtse), method = "bh")
  print(dunn_val)
  sink()
}

##### EACH SITE BY ITSELF ####
# If we look at the graphs and results, there are clear visible differences in them that don't translate into clear differences when 
# comparing by polymers. This could be caused by the (also) clear differences in the scale of Ard vs lakes, so let's separate both
# to make our comparisions

ard_tse <- phtse[ , phtse$spf == "Ardley"]
ard_tse@colData


# Type_g
sink("ARD-tg-no_site.txt")
print("")
sink()
for (alpha_type in alpha_indices) {
  sink("ARD-tg-no_site.txt", append = TRUE)
  print(paste0("#################### ", alpha_type, " ####################"))
  wilcoxon_val <- wilcox.test(as.formula(paste0(alpha_type," ~ type_g")), data = colData(ard_tse))
  print(wilcoxon_val)
  sink()
}

# Type_f
sink("ARD_tf.txt")
print("")
sink()
for (alpha_type in alpha_indices) {
  sink("ARD_tf.txt", append = TRUE)
  print(paste0("#################### ", alpha_type, " ####################"))
  kruskal_val = kruskal.test(as.formula(paste0(alpha_type," ~ type_f")), data = colData(ard_tse))
  print(kruskal_val)
  print(" ")
  dunn_val = dunnTest(as.formula(paste0(alpha_type," ~ type_f")), data = colData(ard_tse), method = "bh")
  print(dunn_val)
  sink()
}

#### LAKES
lake_tse = phtse[, phtse$spf != "Ardley"]

# Type_g 
sink("Lakes-tg.txt")
print("")
sink()
for (alpha_type in alpha_indices) {
  sink("Lakes-tg.txt", append = TRUE)
  print(paste0("#################### ", alpha_type, " ####################"))
  wilcoxon_val <- wilcox.test(as.formula(paste0(alpha_type," ~ type_g")), data = colData(lake_tse))
  print(wilcoxon_val)
  sink()
}

# Type_f
sink("Lakes-tf.txt")
print("")
sink()
for (alpha_type in alpha_indices) {
  sink("Lakes-tf.txt", append = TRUE)
  print(paste0("#################### ", alpha_type, " ####################"))
  kruskal_val = kruskal.test(as.formula(paste0(alpha_type," ~ type_f")), data = colData(lake_tse))
  print(kruskal_val)
  print(" ")
  dunn_val = dunnTest(as.formula(paste0(alpha_type," ~ type_f")), data = colData(lake_tse), method = "bh")
  print(dunn_val)
  sink()
}



