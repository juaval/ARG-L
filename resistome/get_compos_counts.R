#### MAKE COMPOSITIONAL COUNTS #### 
# The idea with this is just to make a quick script which spits out the compositional version of the count table for paper purposes
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


script_dir = getwd()
setwd("../data/r_data")
data_dir = getwd()


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







resis = transform(resis@otu_table, transform = "compositional")
transp = transform(transp@otu_table, transform = "compositional")

resis = resis * 100
resis
write.csv(resis, "compos_ARG.csv")
write.csv(transp, "compos_MGE.csv")