### ALPHA DIVERSITY OF THE RESISTOME ####
# Much like everything resistome before, I want to study MGEs and ARGs by separate
# Otherwise the script is the same I use to calculte the alpha diversity of the bacterial diversity
# with a few minor tweaks

library("mia")
library("phyloseq")
library("dplyr")
library("MicEco")
library("scater")
library("patchwork")
library("tibble")
library("microbiome")
library("FSA")
#library("MetBrewer")




# directorios útiles
script_dir = getwd()
setwd("../data/r_data")
data_dir = getwd()
setwd("../../results/alfa_res")
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

resis_tse = makeTreeSummarizedExperimentFromPhyloseq(resis)
transp_tse = makeTreeSummarizedExperimentFromPhyloseq(transp)

alpha_indices = c("observed", "shannon", "coverage", "relative", "log_modulo_skewness")#, "faith")
for (alpha_type in alpha_indices) {
  #print(alpha_type)
  if (alpha_type == "observed") {
    print(alpha_type)
    resis_tse = estimateRichness(resis_tse, assay.type = "counts",
                                 index = alpha_type, name = alpha_type)
    transp_tse = estimateRichness(transp_tse, assay.type = "counts",
                                  index = alpha_type, name = alpha_type)
  } else if (alpha_type == "relative") {
    resis_tse = estimateDominance(resis_tse, assay.type = "counts", 
                                  index= alpha_type, name = alpha_type)
    transp_tse = estimateDominance(transp_tse, assay.type = "counts", 
                                  index= alpha_type, name = alpha_type)
  } else {
    print(alpha_type)
    resis_tse = estimateDiversity(resis_tse, assay.type = "counts",
                                  index = alpha_type, name = alpha_type)
    transp_tse = estimateDiversity(transp_tse, assay.type = "counts",
                                  index = alpha_type, name = alpha_type)
  }
}

data = as.data.frame(resis_tse@colData)
data_t = as.data.frame(transp_tse@colData)

tab10_cols = c("#ff7f0e", "#1f77b4")
#kandinsky = met.brewer("Wissing", 4)
kandinsky = c("#B22C2C", "#85B22C", "#B26F2C", "#2C85B2")

plots = list()
plots_t = list()

i = 1
for (alpha_type in alpha_indices) {
  #print(alpha_type)
  plots[[i]] = ggplot(data, aes(x = sampling_point, y = .data[[alpha_type]], fill = type_g)) +
    geom_jitter(size = 1) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) + 
    theme(title = element_text(size = 12)) +
    scale_fill_manual(values = tab10_cols) + 
    labs(fill='Sample type')
  
  plots_t[[i]] = ggplot(data_t, aes(x = sampling_point, y = .data[[alpha_type]], fill = type_g)) +
    geom_jitter(size = 1) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) + 
    theme(title = element_text(size = 12)) +
    scale_fill_manual(values = tab10_cols) + 
    labs(fill='Sample type')
  
  print(alpha_type)
  i = i+1
}
#plots[[1]]
plots = lapply(plots, "+", 
                theme(axis.title.x = element_blank(),
                      axis.ticks.y = element_blank()))
plots_t = lapply(plots_t, "+", 
                 theme(axis.title.x = element_blank(),
                       axis.ticks.y = element_blank()))
all_plots = ((plots[[1]] | plots[[2]] | plots[[3]]) / 
                 (plots[[4]] | plots[[5]])) +
  plot_layout(guides = "collect") + 
  plot_annotation(title = "Different measures of ARG α diverstity")

all_plots_t = ((plots_t[[1]] | plots_t[[2]] | plots_t[[3]]) / 
               (plots_t[[4]] | plots_t[[5]])) +
  plot_layout(guides = "collect") + 
  plot_annotation(title = "Different measures of MGE α diverstity")
# Print everything
png(file = "ARG_all_alpha.png", width = 8, height = 6, units = "in", res = 300)
all_plots
invisible(dev.off())

png(file = "MGE_all_alpha.png", width = 8, height = 6, units = "in", res = 300)
all_plots_t
invisible(dev.off())

# print each individual plot. alpha_type to set the name, i to loop trhough the list. 
i = 1
for (alpha_type in alpha_indices) {
  png(file = paste0(alpha_type, "_ARG_index.png"), width = 8, height = 6, units = "in", res = 300)
  print(plots[[i]])
  dev.off()
  png(file = paste0(alpha_type, "_MGE_index.png"), width = 8, height = 6, units = "in", res = 300)
  print(plots_t[[i]])
  dev.off()
  #invisible(dev.off())
  i = i + 1
}
### Save the results of the comparisions! ###
# First, mulch vs soil as is 
sink("ARG-alpha_mulch-vs-soil.txt")
print("")
sink()
for (alpha_type in alpha_indices) {
  sink("ARG-alpha_mulch-vs-soil.txt", append = TRUE)
  print(paste0("#################### ", alpha_type, " ####################"))
  wilcoxon_val <- wilcox.test(as.formula(paste0(alpha_type," ~ type_g")), data = colData(resis_tse))
  print(wilcoxon_val)
  sink()
}

sink("MGE-alpha_mulch-vs-soil.txt")
print("")
sink()
for (alpha_type in alpha_indices) {
  sink("MGE-alpha_mulch-vs-soil.txt", append = TRUE)
  print(paste0("#################### ", alpha_type, " ####################"))
  wilcoxon_val <- wilcox.test(as.formula(paste0(alpha_type," ~ type_g")), data = colData(transp_tse))
  print(wilcoxon_val)
  sink()
}



sink("ARG-alpha_all-vs-all.txt")
print("")
sink()
for (alpha_type in alpha_indices) {
  sink("ARG-alpha_all-vs-all.txt", append = TRUE)
  print(paste0("#################### ", alpha_type, " ####################"))
  kruskal_val = kruskal.test(as.formula(paste0(alpha_type," ~ tf.sp")), data = colData(resis_tse))
  print(kruskal_val)
  print(" ")
  dunn_val = dunnTest(as.formula(paste0(alpha_type," ~ tf.sp")), data = colData(resis_tse), method = "bh")
  print(dunn_val)
  sink()
}

sink("MGE-alpha_all-vs-all.txt")
print("")
sink()
for (alpha_type in alpha_indices) {
  sink("MGE-alpha_all-vs-all.txt", append = TRUE)
  print(paste0("#################### ", alpha_type, " ####################"))
  kruskal_val = kruskal.test(as.formula(paste0(alpha_type," ~ tf.sp")), data = colData(transp_tse))
  print(kruskal_val)
  print(" ")
  dunn_val = dunnTest(as.formula(paste0(alpha_type," ~ tf.sp")), data = colData(transp_tse), method = "bh")
  print(dunn_val)
  sink()
}

#### SEPARATING POLYMER TYPES ####
plots = list()
plots_t = list()

i = 1
for (alpha_type in alpha_indices) {
  #print(alpha_type)
  plots[[i]] = ggplot(data, aes(x = sampling_point, y = .data[[alpha_type]], fill = type_f)) +
    geom_jitter(size = 1) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) + 
    theme(title = element_text(size = 12)) +
    scale_fill_manual(values = kandinsky) + 
    labs(fill='Sample type')
  
  plots_t[[i]] = ggplot(data_t, aes(x = sampling_point, y = .data[[alpha_type]], fill = type_f)) +
    geom_jitter(size = 1) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) + 
    theme(title = element_text(size = 12)) +
    scale_fill_manual(values = kandinsky) +
    labs(fill='Sample type')
  
  print(alpha_type)
  i = i+1
}
#plots[[1]]
plots = lapply(plots, "+", 
               theme(axis.title.x = element_blank(),
                     axis.ticks.y = element_blank()))
plots_t = lapply(plots_t, "+", 
                 theme(axis.title.x = element_blank(),
                       axis.ticks.y = element_blank()))
all_plots = ((plots[[1]] | plots[[2]] | plots[[3]]) / 
               (plots[[4]] | plots[[5]])) +
  plot_layout(guides = "collect") + 
  plot_annotation(title = "Different measures of ARG α diverstity")

all_plots_t = ((plots_t[[1]] | plots_t[[2]] | plots_t[[3]]) / 
                 (plots_t[[4]] | plots_t[[5]])) +
  plot_layout(guides = "collect") + 
  plot_annotation(title = "Different measures of MGE α diverstity")
# Print everything
png(file = "ARG_sep_alpha.png", width = 8, height = 6, units = "in", res = 300)
all_plots
invisible(dev.off())

png(file = "MGE_sep_alpha.png", width = 8, height = 6, units = "in", res = 300)
all_plots_t
invisible(dev.off())

# print each individual plot. alpha_type to set the name, i to loop through the list. 
i = 1
for (alpha_type in alpha_indices) {
  png(file = paste0(alpha_type, "_ARG_index_sep.png"), width = 8, height = 6, units = "in", res = 300)
  print(plots[[i]])
  dev.off()
  png(file = paste0(alpha_type, "_MGE_index_sep.png"), width = 8, height = 6, units = "in", res = 300)
  print(plots_t[[i]])
  dev.off()
  #invisible(dev.off())
  i = i + 1
}
### Save the results of the comparisions! ###
# As we've already done all vs all, we are just missing polymers v polymers 
sink("ARG-alpha_typef.txt")
print("")
sink()
for (alpha_type in alpha_indices) {
  sink("ARG-alpha_typef.txt", append = TRUE)
  print(paste0("#################### ", alpha_type, " ####################"))
  kruskal_val = kruskal.test(as.formula(paste0(alpha_type," ~ type_f")), data = colData(resis_tse))
  print(kruskal_val)
  print(" ")
  dunn_val = dunnTest(as.formula(paste0(alpha_type," ~ type_f")), data = colData(resis_tse), method = "bh")
  print(dunn_val)
  sink()
}

sink("MGE-alpha_typef.txt")
print("")
sink()
for (alpha_type in alpha_indices) {
  sink("MGE-alpha_typef.txt", append = TRUE)
  print(paste0("#################### ", alpha_type, " ####################"))
  kruskal_val = kruskal.test(as.formula(paste0(alpha_type," ~ type_f")), data = colData(transp_tse))
  print(kruskal_val)
  print(" ")
  dunn_val = dunnTest(as.formula(paste0(alpha_type," ~ type_f")), data = colData(transp_tse), method = "bh")
  print(dunn_val)
  sink()
}

##### EACH SITE BY ITSELF ####
# If we look at the graphs and results, there are clear visible differences in them that don't translate into clear differences when 
# comparing by polymers. This could be caused by the (also) clear differences in the scale of Ard vs lakes, so let's separate both
# to make our comparisions

ard_resis = subset_samples(resis, sampling_point == "Ardley")
ard_resis = makeTreeSummarizedExperimentFromPhyloseq(ard_resis)
ard_transp = subset_samples(transp, sampling_point == "Ardley")
ard_transp = makeTreeSummarizedExperimentFromPhyloseq(ard_transp)
# Just for tidyness sake
setwd("Ard")
# Now for graphs
for (alpha_type in alpha_indices) {
  #print(alpha_type)
  if (alpha_type == "observed") {
    print(alpha_type)
    ard_resis = estimateRichness(ard_resis, assay.type = "counts",
                                 index = alpha_type, name = alpha_type)
    ard_transp = estimateRichness(ard_transp, assay.type = "counts",
                                  index = alpha_type, name = alpha_type)
  } else if (alpha_type == "relative") {
    ard_resis = estimateDominance(ard_resis, assay.type = "counts", 
                                  index= alpha_type, name = alpha_type)
    ard_transp = estimateDominance(ard_transp, assay.type = "counts", 
                                   index= alpha_type, name = alpha_type)
  } else {
    print(alpha_type)
    ard_resis = estimateDiversity(ard_resis, assay.type = "counts",
                                  index = alpha_type, name = alpha_type)
    ard_transp = estimateDiversity(ard_transp, assay.type = "counts",
                                   index = alpha_type, name = alpha_type)
  }
}

data = as.data.frame(ard_resis@colData)
data_t = as.data.frame(ard_transp@colData)
plots = list()
plots_t = list()

i = 1
for (alpha_type in alpha_indices) {
  #print(alpha_type)
  plots[[i]] = ggplot(data, aes(x = type_g, y = .data[[alpha_type]], fill = type_g)) +
    geom_jitter(size = 1) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) + 
    theme(title = element_text(size = 12)) +
    scale_fill_manual(values = tab10_cols) + 
    labs(fill='Sample type')
  
  plots_t[[i]] = ggplot(data_t, aes(x = type_g, y = .data[[alpha_type]], fill = type_g)) +
    geom_jitter(size = 1) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) + 
    theme(title = element_text(size = 12)) +
    scale_fill_manual(values = tab10_cols) + 
    labs(fill='Sample type')
  
  print(alpha_type)
  i = i+1
}
#plots[[1]]
plots = lapply(plots, "+", 
               theme(axis.title.x = element_blank(),
                     axis.ticks.y = element_blank()))
plots_t = lapply(plots_t, "+", 
                 theme(axis.title.x = element_blank(),
                       axis.ticks.y = element_blank()))
all_plots = ((plots[[1]] | plots[[2]] | plots[[3]]) / 
               (plots[[4]] | plots[[5]])) +
  plot_layout(guides = "collect") + 
  plot_annotation(title = "Different measures of ARG α diverstity")

all_plots_t = ((plots_t[[1]] | plots_t[[2]] | plots_t[[3]]) / 
                 (plots_t[[4]] | plots_t[[5]])) +
  plot_layout(guides = "collect") + 
  plot_annotation(title = "Different measures of MGE α diverstity")
# Print everything
png(file = "ARG_ardley_all_alpha.png", width = 8, height = 6, units = "in", res = 300)
all_plots
invisible(dev.off())

png(file = "MGE_ardley_all_alpha.png", width = 8, height = 6, units = "in", res = 300)
all_plots_t
invisible(dev.off())

# print each individual plot. alpha_type to set the name, i to loop trhough the list. 
i = 1
for (alpha_type in alpha_indices) {
  png(file = paste0(alpha_type, "_ARG_ardley_index.png"), width = 8, height = 6, units = "in", res = 300)
  print(plots[[i]])
  dev.off()
  png(file = paste0(alpha_type, "_MGE_ardley_index.png"), width = 8, height = 6, units = "in", res = 300)
  print(plots_t[[i]])
  dev.off()
  #invisible(dev.off())
  i = i + 1
}
### Save the results of the comparisions! ###
# First, mulch vs soil as is 
sink("ARG_ardley-alpha_plastic-vs-control.txt")
print("")
sink()
for (alpha_type in alpha_indices) {
  sink("ARG_ardley-alpha_plastic-vs-control.txt", append = TRUE)
  print(paste0("#################### ", alpha_type, " ####################"))
  wilcoxon_val <- wilcox.test(as.formula(paste0(alpha_type," ~ type_g")), data = colData(ard_resis))
  print(wilcoxon_val)
  sink()
}

sink("MGE_ardley-alpha_plastic-vs-control.txt")
print("")
sink()
for (alpha_type in alpha_indices) {
  sink("MGE_ardley-alpha_plastic-vs-control.txt", append = TRUE)
  print(paste0("#################### ", alpha_type, " ####################"))
  wilcoxon_val <- wilcox.test(as.formula(paste0(alpha_type," ~ type_g")), data = colData(ard_transp))
  print(wilcoxon_val)
  sink()
}
# and now separating by polymer type
plots = list()
plots_t = list()

i = 1
for (alpha_type in alpha_indices) {
  #print(alpha_type)
  plots[[i]] = ggplot(data, aes(x = type_g, y = .data[[alpha_type]], fill = type_f)) +
    geom_jitter(size = 1) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) + 
    theme(title = element_text(size = 12)) +
    scale_fill_manual(values = kandinsky) +
    labs(fill='Sample type')
  
  plots_t[[i]] = ggplot(data_t, aes(x = type_g, y = .data[[alpha_type]], fill = type_f)) +
    geom_jitter(size = 1) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) + 
    theme(title = element_text(size = 12)) +
    scale_fill_manual(values = kandinsky) +
    labs(fill='Sample type')
  
  print(alpha_type)
  i = i+1
}
#plots[[1]]
plots = lapply(plots, "+", 
               theme(axis.title.x = element_blank(),
                     axis.ticks.y = element_blank()))
plots_t = lapply(plots_t, "+", 
                 theme(axis.title.x = element_blank(),
                       axis.ticks.y = element_blank()))
all_plots = ((plots[[1]] | plots[[2]] | plots[[3]]) / 
               (plots[[4]] | plots[[5]])) +
  plot_layout(guides = "collect") + 
  plot_annotation(title = "Different measures of ARG α diverstity")

all_plots_t = ((plots_t[[1]] | plots_t[[2]] | plots_t[[3]]) / 
                 (plots_t[[4]] | plots_t[[5]])) +
  plot_layout(guides = "collect") + 
  plot_annotation(title = "Different measures of MGE α diverstity")
# Print everything
png(file = "ARG_ardley__sep_alpha.png", width = 8, height = 6, units = "in", res = 300)
all_plots
invisible(dev.off())

png(file = "MGE_ardley_sep_alpha.png", width = 8, height = 6, units = "in", res = 300)
all_plots_t
invisible(dev.off())

# print each individual plot. alpha_type to set the name, i to loop through the list. 
i = 1
for (alpha_type in alpha_indices) {
  png(file = paste0(alpha_type, "_ARG_ardley_index_sep.png"), width = 8, height = 6, units = "in", res = 300)
  print(plots[[i]])
  dev.off()
  png(file = paste0(alpha_type, "_MGE_ardley_index_sep.png"), width = 8, height = 6, units = "in", res = 300)
  print(plots_t[[i]])
  dev.off()
  #invisible(dev.off())
  i = i + 1
}
### Save the results of the comparisions! ###
sink("ARG_ardley-alpha_typef.txt")
print("")
sink()
for (alpha_type in alpha_indices) {
  sink("ARG_ardley-alpha_typef.txt", append = TRUE)
  print(paste0("#################### ", alpha_type, " ####################"))
  kruskal_val = kruskal.test(as.formula(paste0(alpha_type," ~ type_f")), data = colData(ard_resis))
  print(kruskal_val)
  print(" ")
  dunn_val = dunnTest(as.formula(paste0(alpha_type," ~ type_f")), data = colData(ard_resis), method = "bh")
  print(dunn_val)
  sink()
}

sink("MGE_ardley-alpha_typef.txt")
print("")
sink()
for (alpha_type in alpha_indices) {
  sink("MGE_ardley-alpha_typef.txt", append = TRUE)
  print(paste0("#################### ", alpha_type, " ####################"))
  kruskal_val = kruskal.test(as.formula(paste0(alpha_type," ~ type_f")), data = colData(ard_transp))
  print(kruskal_val)
  print(" ")
  dunn_val = dunnTest(as.formula(paste0(alpha_type," ~ type_f")), data = colData(ard_transp), method = "bh")
  print(dunn_val)
  sink()
}


## AND NOW LAKES 
lake_resis = subset_samples(resis, sampling_point == "Lakes")
lake_resis = makeTreeSummarizedExperimentFromPhyloseq(lake_resis)
lake_transp = subset_samples(transp, sampling_point == "Lakes")
lake_transp = makeTreeSummarizedExperimentFromPhyloseq(lake_transp)
# Just for tidyness sake
setwd("../Lakes")
# Now for graphs
for (alpha_type in alpha_indices) {
  #print(alpha_type)
  if (alpha_type == "observed") {
    print(alpha_type)
    lake_resis = estimateRichness(lake_resis, assay.type = "counts",
                                 index = alpha_type, name = alpha_type)
    lake_transp = estimateRichness(lake_transp, assay.type = "counts",
                                  index = alpha_type, name = alpha_type)
  } else if (alpha_type == "relative") {
    lake_resis = estimateDominance(lake_resis, assay.type = "counts", 
                                  index= alpha_type, name = alpha_type)
    lake_transp = estimateDominance(lake_transp, assay.type = "counts", 
                                   index= alpha_type, name = alpha_type)
  } else {
    print(alpha_type)
    lake_resis = estimateDiversity(lake_resis, assay.type = "counts",
                                  index = alpha_type, name = alpha_type)
    lake_transp = estimateDiversity(lake_transp, assay.type = "counts",
                                   index = alpha_type, name = alpha_type)
  }
}

data = as.data.frame(lake_resis@colData)
data_t = as.data.frame(lake_transp@colData)
plots = list()
plots_t = list()

i = 1
for (alpha_type in alpha_indices) {
  #print(alpha_type)
  plots[[i]] = ggplot(data, aes(x = type_g, y = .data[[alpha_type]], fill = type_g)) +
    geom_jitter(size = 1) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) + 
    theme(title = element_text(size = 12)) +
    scale_fill_manual(values = tab10_cols) + 
    labs(fill='Sample type')
  
  plots_t[[i]] = ggplot(data_t, aes(x = type_g, y = .data[[alpha_type]], fill = type_g)) +
    geom_jitter(size = 1) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) + 
    theme(title = element_text(size = 12)) +
    scale_fill_manual(values = tab10_cols) + 
    labs(fill='Sample type')
  
  print(alpha_type)
  i = i+1
}
#plots[[1]]
plots = lapply(plots, "+", 
               theme(axis.title.x = element_blank(),
                     axis.ticks.y = element_blank()))
plots_t = lapply(plots_t, "+", 
                 theme(axis.title.x = element_blank(),
                       axis.ticks.y = element_blank()))
all_plots = ((plots[[1]] | plots[[2]] | plots[[3]]) / 
               (plots[[4]] | plots[[5]])) +
  plot_layout(guides = "collect") + 
  plot_annotation(title = "Different measures of ARG α diverstity")

all_plots_t = ((plots_t[[1]] | plots_t[[2]] | plots_t[[3]]) / 
                 (plots_t[[4]] | plots_t[[5]])) +
  plot_layout(guides = "collect") + 
  plot_annotation(title = "Different measures of MGE α diverstity")
# Print everything
png(file = "ARG_lakes_all_alpha.png", width = 8, height = 6, units = "in", res = 300)
all_plots
invisible(dev.off())

png(file = "MGE_lakes_all_alpha.png", width = 8, height = 6, units = "in", res = 300)
all_plots_t
invisible(dev.off())

# print each individual plot. alpha_type to set the name, i to loop trhough the list. 
i = 1
for (alpha_type in alpha_indices) {
  png(file = paste0(alpha_type, "_ARG_lakes_index.png"), width = 8, height = 6, units = "in", res = 300)
  print(plots[[i]])
  dev.off()
  png(file = paste0(alpha_type, "_MGE_lakes_index.png"), width = 8, height = 6, units = "in", res = 300)
  print(plots_t[[i]])
  dev.off()
  #invisible(dev.off())
  i = i + 1
}
### Save the results of the comparisions! ###
# First, mulch vs soil as is 
sink("ARG_lakes-alpha_plastic-vs-control.txt")
print("")
sink()
for (alpha_type in alpha_indices) {
  sink("ARG_lakes-alpha_plastic-vs-control.txt", append = TRUE)
  print(paste0("#################### ", alpha_type, " ####################"))
  wilcoxon_val <- wilcox.test(as.formula(paste0(alpha_type," ~ type_g")), data = colData(lake_resis))
  print(wilcoxon_val)
  sink()
}

sink("MGE_lakes-alpha_plastic-vs-control.txt")
print("")
sink()
for (alpha_type in alpha_indices) {
  sink("MGE_lakes-alpha_plastic-vs-control.txt", append = TRUE)
  print(paste0("#################### ", alpha_type, " ####################"))
  wilcoxon_val <- wilcox.test(as.formula(paste0(alpha_type," ~ type_g")), data = colData(lake_transp))
  print(wilcoxon_val)
  sink()
}
# and now separating by polymer type
plots = list()
plots_t = list()

i = 1
for (alpha_type in alpha_indices) {
  #print(alpha_type)
  plots[[i]] = ggplot(data, aes(x = type_g, y = .data[[alpha_type]], fill = type_f)) +
    geom_jitter(size = 1) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) + 
    theme(title = element_text(size = 12)) +
    scale_fill_manual(values = kandinsky) +
    labs(fill='Sample type')
  
  plots_t[[i]] = ggplot(data_t, aes(x = type_g, y = .data[[alpha_type]], fill = type_f)) +
    geom_jitter(size = 1) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) + 
    theme(title = element_text(size = 12)) +
    scale_fill_manual(values = kandinsky) +
    labs(fill='Sample type')
  
  print(alpha_type)
  i = i+1
}
#plots[[1]]
plots = lapply(plots, "+", 
               theme(axis.title.x = element_blank(),
                     axis.ticks.y = element_blank()))
plots_t = lapply(plots_t, "+", 
                 theme(axis.title.x = element_blank(),
                       axis.ticks.y = element_blank()))
all_plots = ((plots[[1]] | plots[[2]] | plots[[3]]) / 
               (plots[[4]] | plots[[5]])) +
  plot_layout(guides = "collect") + 
  plot_annotation(title = "Different measures of ARG α diverstity")

all_plots_t = ((plots_t[[1]] | plots_t[[2]] | plots_t[[3]]) / 
                 (plots_t[[4]] | plots_t[[5]])) +
  plot_layout(guides = "collect") + 
  plot_annotation(title = "Different measures of MGE α diverstity")
# Print everything
png(file = "ARG_lakes_sep_alpha.png", width = 8, height = 6, units = "in", res = 300)
all_plots
invisible(dev.off())

png(file = "MGE_lakes_sep_alpha.png", width = 8, height = 6, units = "in", res = 300)
all_plots_t
invisible(dev.off())

# print each individual plot. alpha_type to set the name, i to loop through the list. 
i = 1
for (alpha_type in alpha_indices) {
  png(file = paste0(alpha_type, "_ARG_lakes_index_sep.png"), width = 8, height = 6, units = "in", res = 300)
  print(plots[[i]])
  dev.off()
  png(file = paste0(alpha_type, "_MGE_lakes_index_sep.png"), width = 8, height = 6, units = "in", res = 300)
  print(plots_t[[i]])
  dev.off()
  #invisible(dev.off())
  i = i + 1
}
### Save the results of the comparisions! ###
sink("ARG_lakes-alpha_typef.txt")
print("")
sink()
for (alpha_type in alpha_indices) {
  sink("ARG_lakes-alpha_typef.txt", append = TRUE)
  print(paste0("#################### ", alpha_type, " ####################"))
  kruskal_val = kruskal.test(as.formula(paste0(alpha_type," ~ type_f")), data = colData(lake_resis))
  print(kruskal_val)
  print(" ")
  dunn_val = dunnTest(as.formula(paste0(alpha_type," ~ type_f")), data = colData(lake_resis), method = "bh")
  print(dunn_val)
  sink()
}

sink("MGE_lakes-alpha_typef.txt")
print("")
sink()
for (alpha_type in alpha_indices) {
  sink("MGE_lakes-alpha_typef.txt", append = TRUE)
  print(paste0("#################### ", alpha_type, " ####################"))
  kruskal_val = kruskal.test(as.formula(paste0(alpha_type," ~ type_f")), data = colData(lake_transp))
  print(kruskal_val)
  print(" ")
  dunn_val = dunnTest(as.formula(paste0(alpha_type," ~ type_f")), data = colData(lake_transp), method = "bh")
  print(dunn_val)
  sink()
}
