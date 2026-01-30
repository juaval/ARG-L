####### SET THINGS
# Objective of the script: get a .json of ARGs unique to specific sets and generate venn graphs of said sets
 
##### JSON MAKING ####
library("phyloseq")
library("MicEco")
library("dplyr")
library("MiscMetabar")
library("ComplexUpset")
library("patchwork")
library("tibble")
library("jsonlite")
library("eulerr")

# directorios útiles
script_dir = getwd()
setwd("../data/clean_data")
data_dir = getwd()
setwd("../../results/Venn_res/all")
all_res = getwd()
setwd("../ard")
ard_res = getwd()
setwd("../ion")
ion_res = getwd()
setwd("../uru")
uru_res = getwd()

setwd(data_dir)
meta_tab = read.csv("../metadata/metadata_resistome.csv", sep = ",", row.names = 1)
meta_tab[order(row.names(meta_tab)), ] #por consistencia con luego las abundancias
META = sample_data(meta_tab)
arg_freqs_tab = read.csv("count_arg.csv", sep = ",", header = TRUE, row.names = 1)
mge_freqs_tab = read.csv("count_mge.csv", sep = ",", header = TRUE, row.names = 1)
ARG = otu_table(arg_freqs_tab, taxa_are_rows = TRUE)
MGE = otu_table(mge_freqs_tab, taxa_are_rows = TRUE)
resis = phyloseq(ARG, META)
transp = phyloseq(MGE, META)

set.seed(1312)
##### GENERAL RESULTS #####
setwd(all_res)
arg_typeg = ps_venn(resis, group = "type_g", plot = FALSE)
arg_typeg = toJSON(arg_typeg, pretty = TRUE)
write(arg_typeg, "arg_typeg_all.json")
mge_typeg = ps_venn(transp, group = "type_g", plot = FALSE)
mge_typeg = toJSON(mge_typeg, pretty = TRUE)
write(mge_typeg, "mge_typeg_all.json")
# repeat for type_f, used to make tables
arg_typef = ps_venn(resis, group = "type_f", plot = FALSE)
arg_typef = toJSON(arg_typef, pretty = TRUE)
write(arg_typef, "arg_typef_all.json")
mge_typef = ps_venn(transp, group = "type_f", plot = FALSE)
mge_typef = toJSON(mge_typef, pretty = TRUE)
write(mge_typef, "mge_typef_all.json")

rm(arg_typef)
rm(arg_typeg)
rm(mge_typeg)
rm(mge_typef)
# Make plots
arg_sets = read_json("arg_typeg_all.json", simplifyVector = TRUE)
mge_sets = read_json("mge_typeg_all.json", simplifyVector = TRUE)

arg_plas = as.numeric(length(arg_sets$plastic))
arg_sh = as.numeric(length(arg_sets$control__plastic))
arg_c = as.numeric(length(arg_sets$control))
# Regular plot
lala = euler(c("Plastic" = arg_plas, "Surrounding environment" = arg_c, "Plastic&Surrounding environment" = arg_sh))
png("Global resistome.png", width = 10, height = 6, units = "in", res = 300)
plot(lala, counts = TRUE, font=2, cex=1, alpha=0.8,
     quantities = list(type = c("percent", "counts"),
                       col = c("black", "black","white"),
                       cex = 1.5
     ),
     labels = list(col = "black",
                   cex = 1.5
     ),
     edges = list(col = c("gray5","gray5","white"), 
                  lex = 3
     ),
     fills = c("#ff7f0e","#1f77b4","gray20"))
invisible(dev.off())
# Empty plot
png("Global resistome-empty.png", width = 10, height = 6, units = "in", res = 300) #these empy ones make adjusting the graphs in other software much, much easier
plot(lala, counts = TRUE, font=2, cex=1, alpha=0.8,
     quantities = list(type = c("percent", "counts"),
                       col = c("black", "black","white"),
                       cex = 0
     ),
     labels = list(col = "black",
                   cex = 0
     ),
     edges = list(col = c("gray5","gray5","white"), 
                  lex = 3
     ),
     fills = c("#ff7f0e","#1f77b4","gray20"))
invisible(dev.off())

mge_plas = as.numeric(length(mge_sets$plastic))
mge_sh = as.numeric(length(mge_sets$control__plastic))
mge_c = 0 # there are no exclusivities in this set

lala = euler(c("Plastic" = mge_plas, "Plastic&Surrounding environment" = mge_sh, "Surrounding environment" = 0))
png("Global mobilome.png", width = 10, height = 6, units = "in", res = 300)
plot(lala, counts = TRUE, font=2, cex=1, alpha=0.8,
     quantities = list(type = c("percent", "counts"),
                       col = c("black", "white"),
                       cex = 1.5
     ),
     labels = list(col = "black",
                   cex = c(1.5, 0)
     ),
     edges = list(col = c("gray5","gray5"), 
                  lex = 3
     ),
     fills = c("#ff7f0e", "gray20"))
invisible(dev.off())
png("Global mobilome-empty.png", width = 10, height = 6, units = "in", res = 300)
plot(lala, counts = TRUE, font=2, cex=1, alpha=0.8,
     quantities = list(type = c("percent", "counts"),
                       col = c("black", "white"),
                       cex = 0
     ),
     labels = list(col = "black",
                   cex = 0
     ),
     edges = list(col = c("gray5","gray5"), 
                  lex = 3
     ),
     fills = c("#ff7f0e", "gray20"))
invisible(dev.off())


##### PLACE SPECIFIC RESULTS #####
corr_df = data.frame(varname = unique(resis@sam_data$sampling_point), dirname = c(ard_res, ion_res, uru_res))
corr_df = corr_df %>% remove_rownames %>% column_to_rownames(var="varname")
rownames(corr_df)

resis@sam_data

for (area in rownames(corr_df)) {
  # Create the subset
  area_resis = NA
  area_trans = NA
  area_resis = subset_samples(resis, sampling_point == area)
  area_transp = subset_samples(transp, sampling_point == area)
  print(paste0("Working on ", unique(area_resis@sam_data$sampling_point)))
  #Go to the subset's specific directory to store its results
  setwd(corr_df[area, ])
  
  # Repeat everything
  arg_typeg = ps_venn(area_resis, group = "type_g", plot = FALSE)
  arg_typeg = toJSON(arg_typeg, pretty = TRUE)
  write(arg_typeg, paste0("arg_typeg_", area, ".json"))
  mge_typeg = ps_venn(area_transp, group = "type_g", plot = FALSE)
  mge_typeg = toJSON(mge_typeg, pretty = TRUE)
  write(mge_typeg, paste0("mge_typeg_", area, ".json"))
  
  # repeat for type_f, used to make tables
  arg_typef = ps_venn(area_resis, group = "type_f", plot = FALSE)
  arg_typef = toJSON(arg_typef, pretty = TRUE)
  write(arg_typef, paste0("arg_typef_", area, ".json"))
  mge_typef = ps_venn(area_transp, group = "type_f", plot = FALSE)
  mge_typef = toJSON(mge_typef, pretty = TRUE)
  write(mge_typef, paste0("mge_typef_", area, ".json"))
  
  rm(arg_typef)
  rm(arg_typeg)
  rm(mge_typeg)
  rm(mge_typef)
  
  arg_sets = read_json(paste0("arg_typeg_", area, ".json"), simplifyVector = TRUE)
  mge_sets = read_json(paste0("mge_typeg_", area, ".json"), simplifyVector = TRUE)
  
  arg_plas = as.numeric(length(arg_sets$plastic))
  arg_sh = as.numeric(length(arg_sets$control__plastic))
  arg_c = as.numeric(length(arg_sets$control))
  
  lala = euler(c("Plastic" = arg_plas, "Soil" = arg_c, "Plastic&Soil" = arg_sh))
  lolo = plot(lala, counts = TRUE, font=2, cex=1, alpha=0.8,
              quantities = list(type = c("percent", "counts"),
                                col = c("black", "black","white"),
                                cex = 1.5),
              labels = list(col = "black",
                            cex = 1.5),
              edges = list(col = c("gray5","gray5","white"),
                           lex = 3),
              fills = c("#ff7f0e","#1f77b4","gray20"))
  png(paste0(area, "-resistome.png"), width = 10, height = 6, units = "in", res = 300)
  print(lolo)
  invisible(dev.off())
  
  lolo = plot(lala, counts = TRUE, font=2, cex=1, alpha=0.8,
              quantities = list(type = c("percent", "counts"),
                                col = c("black", "black","white"),
                                cex = 0),
              labels = list(col = "black",
                            cex = 0),
              edges = list(col = c("gray5","gray5","white"),
                           lex = 3),
              fills = c("#ff7f0e","#1f77b4","gray20"))
  png(paste0(area, "-resistome-empty.png"), width = 10, height = 6, units = "in", res = 300)
  print(lolo)
  invisible(dev.off())
  
  mge_plas = as.numeric(length(mge_sets$plastic))
  mge_sh = as.numeric(length(mge_sets$control__plastic))
  mge_c = as.numeric(length(mge_sets$control))
  
  lala = euler(c("Plastic" = mge_plas, "Soil" = mge_c, "Plastic&Soil" = mge_sh))
  lolo = plot(lala, counts = TRUE, font=2, cex=1, alpha=0.8,
              quantities = list(type = c("percent", "counts"),
                                col = c("black", "black","white"),
                                cex = 1.5),
              labels = list(col = "black",
                            cex = 1.5),
              edges = list(col = c("gray5","gray5","white"),
                           lex = 3),
              fills = c("#ff7f0e","#1f77b4","gray20"))
  png(paste0(area, "-mobilome.png"), width = 10, height = 6, units = "in", res = 300)
  print(lolo)
  invisible(dev.off())
  lolo = plot(lala, counts = TRUE, font=2, cex=1, alpha=0.8,
              quantities = list(type = c("percent", "counts"),
                                col = c("black", "black","white"),
                                cex = 0),
              labels = list(col = "black",
                            cex = 0),
              edges = list(col = c("gray5","gray5","white"),
                           lex = 3),
              fills = c("#ff7f0e","#1f77b4","gray20"))
  png(paste0(area, "-mobilome-empty.png"), width = 10, height = 6, units = "in", res = 300)
  print(lolo)
  invisible(dev.off())
}

