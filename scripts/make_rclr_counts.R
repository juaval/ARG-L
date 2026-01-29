##### Getting the core microbiota  #######
library("phyloseq")
library("microbiome")


script_dir = getwd()
setwd("../data/r_data")
data_dir = getwd()
setwd("../../../metabarcoding/data/r_backup")
rback_dir = getwd()

setwd(rback_dir)
ph = readRDS("phystuff_1")
ph = transform(ph, transform = "rclr")
setwd(data_dir)
write.csv(ph@otu_table, "rclr_zotu_counts.csv")
