##### Getting the core microbiota  #######
library("phyloseq")
library("MicEco")
library("microViz") # ojo
library("dplyr")
library("MiscMetabar")
library("eulerr")
library("patchwork")
library("microbiome")
library("viridis")
library("tibble")
library("ggplot2")
library("reshape2")
library("treemapify")
library("robustbase")
library("FSA")
#library(devtools)
#install_github("Russel88/MicEco")

script_dir = getwd()
setwd("../data/r_data")
data_dir = getwd()
setwd("../../results/core_res")
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
#resis = phyloseq(MGE, META)

setwd(res_dir)
resis@sam_data$tg.sa = as.factor(paste(resis@sam_data$sampling_pint, resis@sam_data$type_g, sep = "-"))
resis@sam_data$rep = row.names(resis@sam_data)

# This is quite nasty, but it ensures data are grouped together nicely when graphing 
resis@sam_data$rep[resis@sam_data$rep == "argl_1"] = "S1s-1"
resis@sam_data$rep[resis@sam_data$rep == "argl_2"] = "S1s-2"
resis@sam_data$rep[resis@sam_data$rep == "argl_3"] = "S1s-3"
resis@sam_data$rep[resis@sam_data$rep == "argl_4"] = "S1w-1"
resis@sam_data$rep[resis@sam_data$rep == "argl_5"] = "S1w-2"
resis@sam_data$rep[resis@sam_data$rep == "argl_6"] = "S1w-3"
resis@sam_data$rep[resis@sam_data$rep == "argl_7"] = "S1e-1"
resis@sam_data$rep[resis@sam_data$rep == "argl_8"] = "S1e-2"
resis@sam_data$rep[resis@sam_data$rep == "argl_9"] = "S1e-3"
resis@sam_data$rep[resis@sam_data$rep == "argl_10"] = "S1p-1"
resis@sam_data$rep[resis@sam_data$rep == "argl_11"] = "S1p-2"
resis@sam_data$rep[resis@sam_data$rep == "argl_12"] = "S1p-3"
resis@sam_data$rep[resis@sam_data$rep == "argl_13"] = "S2s-1"
resis@sam_data$rep[resis@sam_data$rep == "argl_14"] = "S2w-1"
resis@sam_data$rep[resis@sam_data$rep == "argl_15"] = "S2e-1"
resis@sam_data$rep[resis@sam_data$rep == "argl_16"] = "S2e-2"
resis@sam_data$rep[resis@sam_data$rep == "argl_17"] = "S2p-1"
resis@sam_data$rep[resis@sam_data$rep == "argl_18"] = "S2p-2"
resis@sam_data$rep[resis@sam_data$rep == "argl_19"] = "S3s-1"
resis@sam_data$rep[resis@sam_data$rep == "argl_20"] = "S3s-2"
resis@sam_data$rep[resis@sam_data$rep == "argl_21"] = "S3s-3"
resis@sam_data$rep[resis@sam_data$rep == "argl_22"] = "S3e-1"
resis@sam_data$rep[resis@sam_data$rep == "argl_23"] = "S3e-2"
resis@sam_data$rep[resis@sam_data$rep == "argl_24"] = "S3e-3"
resis@sam_data$rep[resis@sam_data$rep == "argl_25"] = "S3p-1"
resis@sam_data$rep[resis@sam_data$rep == "argl_26"] = "S3p-2"
resis@sam_data$rep[resis@sam_data$rep == "argl_27"] = "S3p-3"

row.names(resis@sam_data) = resis@sam_data$rep
sample_names(resis) = row.names(resis@sam_data)
new_order = sort(sample_names(resis))
new_order
resis = ps_reorder(resis, new_order)


#### GETTING THE CORE MICROBIOTA ITSELF #####
# I would define a function that extracts cores but I do not know how to use properly 
# define subset_samples terms within a function soo ugly code

water_core_ps = subset_samples(resis, type_f == "water")
water_core_ps@otu_table = water_core_ps@otu_table[apply(water_core_ps@otu_table, 1, function(row) !0 %in% row), ]
water_core = row.names(water_core_ps@otu_table)
#rm(water_resis)
print(water_core)

EPS_core_ps = subset_samples(resis, type_f == "EPS")
EPS_core_ps@otu_table = EPS_core_ps@otu_table[apply(EPS_core_ps@otu_table, 1, function(row) !0 %in% row), ]
EPS_core = row.names(EPS_core_ps@otu_table) 
rm(EPS_resis)
print(EPS_core)

soil_core_ps = subset_samples(resis, type_f == "soil")
soil_core_ps@otu_table = soil_core_ps@otu_table[apply(soil_core_ps@otu_table, 1, function(row) !0 %in% row), ] # THIS FAILS BECAUSE THERE'S NO CORE!
soil_core = row.names(soil_core_ps@otu_table)
#rm(soil_resis)
print(soil_core)

PUR_core_ps = subset_samples(resis, type_f == "PUR")
PUR_core_ps@otu_table = PUR_core_ps@otu_table[apply(PUR_core_ps@otu_table, 1, function(row) !0 %in% row), ]
PUR_core = row.names(PUR_core_ps@otu_table) 
#rm(PUR_resis)
print(PUR_core)

full_core = c(water_core, EPS_core, PUR_core)
print(full_core)

#core_resis = prune_taxa(full_core, resis)

# Now I want to define a new vector in which the plastic cores and the water core are separated
all_cores = c()
all_cores[["EPS"]] = EPS_core
all_cores[["PUR"]] = PUR_core
all_cores[["water"]] = water_core

##### RESULTS ####
# I'm going to make a very simple loop to generate the results for plastics and soils by separate
for (sam_type in c("EPS", "PUR", "water")) {
  ## GETTING THE STATS OF THE CORE ##
  # Compositional counts are needed for generating the table we want to make, as well as the
  # relative abundance within the global composition data heatmap and barplots and the within core
  # treemap
  area_resis = subset_samples(resis, type_f == sam_type)
  compos_resis = transform(area_resis, transform = "compositional") # get global compos stats
  global_compos = as.data.frame(compos_resis@otu_table)[c(all_cores[[sam_type]]), ] * 100 #core global compos stats, in %
  global_counts = as.data.frame(area_resis@otu_table)[c(all_cores[[sam_type]]), ] #core global counts just in case
  sam_core = as.data.frame(compos_resis@otu_table)[c(all_cores[[sam_type]]), ] #get just the core
  cols_summed = colSums(global_counts) # get the core total for each sample
  
  for (col in colnames(sam_core)){ # transform the counts to within %
    sam_core[[col]] = sam_core[[col]] / cols_summed[[col]]
    }
  #sam_core = sam_core * 100
  # start saving each table as a .csv
  write.csv(global_compos, paste0(sam_type, "--arg-global_pct.csv"))
  write.csv(sam_core, paste0(sam_type, "--arg-core_pct.csv"))
  write.csv(global_counts, paste0(sam_type, "--arg-core_counts.csv"))
  
  ### FIGURES ###
  ## COMPOSITIONAL BARPLOT
  # First get the values themselves
  all_summed = colSums(as.data.frame(area_resis@otu_table)) # total counts per sample
  sam_core_pcts = (cols_summed/all_summed) * 100 # get the total pct of the whole core, per sample
  inverse_core_pcts = ((sam_core_pcts/sam_core_pcts) * 100) - sam_core_pcts #weird hack to circunvent having to recreate the table headers
  # Then put them in a format readily availible for ggplot
  sam_names = row.names(as.data.frame(sam_core_pcts))
  sam_points = c(rep(sam_names, 2))
  condition = c(rep("all", length(sam_names)), rep("core", length(sam_names)))
  vals = c(c(as.numeric(inverse_core_pcts)), c(as.numeric(sam_core_pcts)))
  data = data.frame(sam_points, condition, vals)
  # Plot
  global_barplot = ggplot(data, aes(fill=condition, y=vals, x=sam_points)) + 
    geom_bar(position="stack", stat="identity") +
    labs(title = paste0(sam_type, " full core global pct")) + 
    theme(axis.text.x = element_text(angle=45))
  png(paste0(sam_type, "_arg_global_barplot.png"), width = 8, height = 6, units = "in", res = 300)
  print(global_barplot)
  invisible(dev.off())
  ## Treemap of core composition. 
  # I'm going to agglomerate the data by sampling point for this part, so I will need to redo some datasets I already have
  type_zotu_pcts = rowMeans(sam_core)
  type_zotu_pcts = round(type_zotu_pcts, digits = 3)
  type_zotu_pcts = as.data.frame(type_zotu_pcts)
  type_zotu_pcts$zotu = row.names(type_zotu_pcts)
  #print(type_zotu_pcts)
  type_treemap = ggplot(type_zotu_pcts, aes(area = type_zotu_pcts, fill = zotu, 
                                            label = paste(zotu, type_zotu_pcts, sep = "\n"))) +
    geom_treemap() + 
    geom_treemap_text(colour = "white", place = "centre", size = 20, reflow = TRUE)
  png(paste0(sam_type,"-arg-global_treemap.png"), width = 8, height = 6, units = "in", res = 300)
  print(type_treemap)
  invisible(dev.off())
  
}

# Now I want to get some statistical results of whether the core is more prevalent in mucl or in soil
#resis = transform(resis, transform = "rclr")
sink("ARG-core_pct-kruskal.txt") # Pre-clean it in case to prevent append from doing weird stuff
print("")
sink()

# obtain the total of the pcts of the core (EPS in this case), per sample
EPS_core_summed = EPS_core_ps@otu_table %>%
                  as.data.frame() %>%
                  colSums()
# obtain all (EPS in this case) pcts, per sample
all_EPS_summed = subset_samples(resis, type_f == "EPS")@otu_table %>%
                 as.data.frame() %>%
                 colSums()

EPS_core_pcts = (EPS_core_summed/all_EPS_summed) * 100

# Now the same for water and PUR
PUR_core_summed = PUR_core_ps@otu_table %>%
  as.data.frame() %>%
  colSums()
all_PUR_summed = subset_samples(resis, type_f == "PUR")@otu_table %>%
  as.data.frame() %>%
  colSums()
PUR_core_pcts = (PUR_core_summed/all_PUR_summed) * 100 

water_core_summed = water_core_ps@otu_table %>%
  as.data.frame() %>%
  colSums()
all_water_summed = subset_samples(resis, type_f == "water")@otu_table %>%
  as.data.frame() %>%
  colSums()
water_core_pcts = (water_core_summed/all_water_summed) * 100 

data = c(EPS_core_pcts, PUR_core_pcts, water_core_pcts)
groups = factor(rep(1:3, c(length(EPS_core_pcts), length(PUR_core_pcts), length(water_core_pcts))),
                labels = c("EPS", "PUR", "Water"))
kruskal.test(data, groups) #result = no difference. Which is weird, so I'm going to test by putting all plastic pcts together
plas_core_pcts = c(EPS_core_pcts, PUR_core_pcts)
wilcox.test(plas_core_pcts, water_core_pcts)

sink("ARG-core_pct-kruskal.txt", append = TRUE)
print("#################### Wilcox test to test whether the core microbiota is significantly more present on a specific sample type ####################")
wilcoxon_val = wilcox.test(plas_core_pcts, water_core_pcts)
print(wilcoxon_val)
print("#################### Kruskal test to test whether the core microbiota is significantly more present on a specific sample type ####################")
kruskal_val = kruskal.test(data, groups)
print(kruskal_val)
sink()

## Now let's make some boxplots to represent this difference
# let's begin with the df of everything by separate
pcts_df = data.frame(data, groups)
tab10_cols = c("#ff7f0e", "#1f77b4")
chosen3_cols = c("#B22C2C", "#85B22C", "#2C85B2")


all = ggplot(data = pcts_df, mapping = aes(x = groups, y = data, fill = groups)) +
               geom_boxplot(alpha = 0.8) +
               stat_summary(fun = "mean", geom = "point", shape = 4,
                            size = 4, color = "black") +
               scale_fill_manual(values = chosen3_cols)

all
# Now the df of plastics all together
groups = factor(rep(1:2, c(length(plas_core_pcts), length(water_core_pcts))),
                           labels = c("Plastic", "Water"))
data = c(plas_core_pcts, water_core_pcts)
binary_df = data.frame(data, groups)
two = ggplot(data = binary_df, mapping = aes(x = groups, y = data, fill = groups)) +
  geom_boxplot(alpha = 0.8) +
  stat_summary(fun = "mean", geom = "point", shape = 4,
               size = 4, color = "black") + 
  scale_fill_manual(values = tab10_cols)

two
final_pct_graph = two + all + plot_layout(guides = "collect", axes = "collect")

png(file = "RESIS_pcts_compared.png", width = 8, height = 6, units = "in", res = 300)
final_pct_graph
invisible(dev.off())


### MOBILOME(?) ANALYSIS
## Let's go back to the MGEs and do the same
setwd(res_dir)
transp@sam_data$tg.sa = as.factor(paste(transp@sam_data$sampling_pint, transp@sam_data$type_g, sep = "-"))
transp@sam_data$rep = row.names(transp@sam_data)

# This is quite nasty, but it ensures data are grouped together nicely when graphing 
transp@sam_data$rep[transp@sam_data$rep == "argl_1"] = "S1s-1"
transp@sam_data$rep[transp@sam_data$rep == "argl_2"] = "S1s-2"
transp@sam_data$rep[transp@sam_data$rep == "argl_3"] = "S1s-3"
transp@sam_data$rep[transp@sam_data$rep == "argl_4"] = "S1w-1"
transp@sam_data$rep[transp@sam_data$rep == "argl_5"] = "S1w-2"
transp@sam_data$rep[transp@sam_data$rep == "argl_6"] = "S1w-3"
transp@sam_data$rep[transp@sam_data$rep == "argl_7"] = "S1e-1"
transp@sam_data$rep[transp@sam_data$rep == "argl_8"] = "S1e-2"
transp@sam_data$rep[transp@sam_data$rep == "argl_9"] = "S1e-3"
transp@sam_data$rep[transp@sam_data$rep == "argl_10"] = "S1p-1"
transp@sam_data$rep[transp@sam_data$rep == "argl_11"] = "S1p-2"
transp@sam_data$rep[transp@sam_data$rep == "argl_12"] = "S1p-3"
transp@sam_data$rep[transp@sam_data$rep == "argl_13"] = "S2s-1"
transp@sam_data$rep[transp@sam_data$rep == "argl_14"] = "S2w-1"
transp@sam_data$rep[transp@sam_data$rep == "argl_15"] = "S2e-1"
transp@sam_data$rep[transp@sam_data$rep == "argl_16"] = "S2e-2"
transp@sam_data$rep[transp@sam_data$rep == "argl_17"] = "S2p-1"
transp@sam_data$rep[transp@sam_data$rep == "argl_18"] = "S2p-2"
transp@sam_data$rep[transp@sam_data$rep == "argl_19"] = "S3s-1"
transp@sam_data$rep[transp@sam_data$rep == "argl_20"] = "S3s-2"
transp@sam_data$rep[transp@sam_data$rep == "argl_21"] = "S3s-3"
transp@sam_data$rep[transp@sam_data$rep == "argl_22"] = "S3e-1"
transp@sam_data$rep[transp@sam_data$rep == "argl_23"] = "S3e-2"
transp@sam_data$rep[transp@sam_data$rep == "argl_24"] = "S3e-3"
transp@sam_data$rep[transp@sam_data$rep == "argl_25"] = "S3p-1"
transp@sam_data$rep[transp@sam_data$rep == "argl_26"] = "S3p-2"
transp@sam_data$rep[transp@sam_data$rep == "argl_27"] = "S3p-3"

row.names(transp@sam_data) = transp@sam_data$rep
sample_names(transp) = row.names(transp@sam_data)
new_order = sort(sample_names(transp))
new_order
transp = ps_reorder(transp, new_order)

water_core_ps = subset_samples(transp, type_f == "water")
water_core_ps@otu_table = water_core_ps@otu_table[apply(water_core_ps@otu_table, 1, function(row) !0 %in% row), ]
water_core = row.names(water_core_ps@otu_table)
#rm(water_transp)
print(water_core)

EPS_core_ps = subset_samples(transp, type_f == "EPS")
EPS_core_ps@otu_table = EPS_core_ps@otu_table[apply(EPS_core_ps@otu_table, 1, function(row) !0 %in% row), ]
EPS_core = row.names(EPS_core_ps@otu_table) 
#rm(EPS_transp)
print(EPS_core)

soil_core_ps = subset_samples(transp, type_f == "soil")
soil_core_ps@otu_table = soil_core_ps@otu_table[apply(soil_core_ps@otu_table, 1, function(row) !0 %in% row), ] # THIS FAILS BECAUSE THERE'S NO CORE!
soil_core = row.names(soil_core_ps@otu_table)
#rm(soil_transp)
print(soil_core)

PUR_core_ps = subset_samples(transp, type_f == "PUR")
PUR_core_ps@otu_table = PUR_core_ps@otu_table[apply(PUR_core_ps@otu_table, 1, function(row) !0 %in% row), ]
PUR_core = row.names(PUR_core_ps@otu_table) 
#rm(PUR_transp)
print(PUR_core)

# So, there's only one MGE in the core of both EPS, PUR (intl3 for both) and water (trb-C). No need to make any graph, then. 
# But, just in case, let's make all the calculations related to the pct

# EPS
EPS_core_summed = EPS_core_ps@otu_table %>%
  as.data.frame() %>%
  colSums()
all_EPS_summed = subset_samples(resis, type_f == "EPS")@otu_table %>%
  as.data.frame() %>%
  colSums()
EPS_core_pcts = (EPS_core_summed/all_EPS_summed) * 100

# PUR
PUR_core_summed = PUR_core_ps@otu_table %>%
  as.data.frame() %>%
  colSums()
all_PUR_summed = subset_samples(resis, type_f == "PUR")@otu_table %>%
  as.data.frame() %>%
  colSums()
PUR_core_pcts = (PUR_core_summed/all_PUR_summed) * 100 
# Water
water_core_summed = water_core_ps@otu_table %>%
  as.data.frame() %>%
  colSums()
all_water_summed = subset_samples(resis, type_f == "water")@otu_table %>%
  as.data.frame() %>%
  colSums()
water_core_pcts = (water_core_summed/all_water_summed) * 100 

# Now, the statistics
data = c(EPS_core_pcts, PUR_core_pcts, water_core_pcts)
groups = factor(rep(1:3, c(length(EPS_core_pcts), length(PUR_core_pcts), length(water_core_pcts))),
                labels = c("EPS", "PUR", "Water"))
kruskal.test(data, groups) #result = no difference
plas_core_pcts = c(EPS_core_pcts, PUR_core_pcts)
wilcox.test(plas_core_pcts, water_core_pcts) #result = no difference

sink("MGE-core_pct-kruskal.txt") # Pre-clean it in case to prevent append from doing weird stuff
print("")
sink()
sink("MGE-core_pct-kruskal.txt", append = TRUE)
print("#################### Wilcox test to test whether the core microbiota is significantly more present on a specific sample type ####################")
wilcoxon_val = wilcox.test(plas_core_pcts, water_core_pcts)
print(wilcoxon_val)
print("#################### Kruskal test to test whether the core microbiota is significantly more present on a specific sample type ####################")
kruskal_val = kruskal.test(data, groups)
print(kruskal_val)
sink()

pcts_df = data.frame(data, groups)

all = ggplot(data = pcts_df, mapping = aes(x = groups, y = data, fill = groups)) +
  geom_boxplot(alpha = 0.8) +
  stat_summary(fun = "mean", geom = "point", shape = 4,
               size = 4, color = "black") +
  scale_fill_manual(values = chosen3_cols)

all
# Now the df of plastics all together
groups = factor(rep(1:2, c(length(plas_core_pcts), length(water_core_pcts))),
                labels = c("Plastic", "Water"))
data = c(plas_core_pcts, water_core_pcts)
binary_df = data.frame(data, groups)
two = ggplot(data = binary_df, mapping = aes(x = groups, y = data, fill = groups)) +
  geom_boxplot(alpha = 0.8) +
  stat_summary(fun = "mean", geom = "point", shape = 4,
               size = 4, color = "black") + 
  scale_fill_manual(values = tab10_cols)

two
final_pct_graph = two + all + plot_layout(guides = "collect", axes = "collect")

png(file = "TRANSP_pcts_compared.png", width = 8, height = 6, units = "in", res = 300)
final_pct_graph
invisible(dev.off())


## There's one last thing to be done: generate the .csv files with the compositional counts for the second supplementary
## So let's reuse the but of code of the ARG part that does just that

all_cores = c()
all_cores[["EPS"]] = EPS_core
all_cores[["PUR"]] = PUR_core
all_cores[["water"]] = water_core
for (sam_type in c("EPS", "PUR", "water")) {
  ## GETTING THE STATS OF THE CORE ##
  # Compositional counts are needed for generating the table we want to make, as well as the
  # relative abundance within the global composition data heatmap and barplots and the within core
  # treemap
  area_transp = subset_samples(transp, type_f == sam_type)
  compos_transp = transform(area_transp, transform = "compositional") # get global compos stats
  global_compos = as.data.frame(compos_transp@otu_table)[c(all_cores[[sam_type]]), ] * 100 #core global compos stats, in %
  global_counts = as.data.frame(area_transp@otu_table)[c(all_cores[[sam_type]]), ] #core global counts just in case
  sam_core = as.data.frame(compos_transp@otu_table)[c(all_cores[[sam_type]]), ] #get just the core
  cols_summed = colSums(global_counts) # get the core total for each sample
  
  for (col in colnames(sam_core)){ # transform the counts to within %
    sam_core[[col]] = sam_core[[col]] / cols_summed[[col]]
  }
  #sam_core = sam_core * 100
  # start saving each table as a .csv
  write.csv(global_compos, paste0(sam_type, "--mge-global_pct.csv"))
  write.csv(sam_core, paste0(sam_type, "--mge-core_pct.csv"))
  write.csv(global_counts, paste0(sam_type, "--mge-core_counts.csv"))
}
















