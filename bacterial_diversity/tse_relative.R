library("mia")
library("phyloseq")
library("dplyr")
library("MicEco")
library("scater")
library("patchwork")
library("tibble")
library("microbiome")
library("vegan")
library("miaViz")

script_dir = getwd()
setwd("../data/r_backup")
rback_dir = getwd()
setwd("../../results/mia_results/relative")
res_dir = getwd()

setwd(rback_dir)
ph = readRDS("phystuff_1")
ph = ps_tax_clean(ph)
ph@phy_tree = phangorn::midpoint(ph@phy_tree)

ph@sam_data$type_f[ph@sam_data$type_f == "Suelo"] = "Soil"
ph@sam_data$type_f[ph@sam_data$type_f == "Agua"] = "Water"

ph@sam_data$tf.sp = as.factor(paste(ph@sam_data$sampling_point, ph@sam_data$type_f))

#ph@sam_data$sam_area[ph@sam_data$sam_area == "inv"] = "S1"
#ph@sam_data$sam_area[ph@sam_data$sam_area == "sapo"] = "S2"
#ph@sam_data$sam_area[ph@sam_data$sam_area == "alm"] = "S3"
#ph@sam_data$tg.sa = as.factor(paste(ph@sam_data$sam_area, ph@sam_data$type_g, sep = "-"))



phtse = convertFromPhyloseq(ph)

##### GRÁFICO DE BARRAS #####
# No hace falta mucho más para generar un gráfico de barras en condiciones
setwd(res_dir)
lala = plotAbundance(phtse, rank = "phylum", assay.type = "counts", 
              features = "tf.sp", order_rank_by = "abund", order_sample_by = "tf.sp")
wrap_plots(lala, ncol = 1)

lala[[1]] <- plotAbundance(phtse, assay.type = "counts", rank = "phylum",
                           features = "tf.sp", order_sample_by = "tf.sp",
                           add_legend = TRUE) [[1]]
# \donttest{
png(file = "Relative_phylum.png", width = 8, height = 6, units = "in", res = 300)
wrap_plots(lala, ncol = 1, heights = c(1,0.05))
invisible(dev.off())
#### SACAR LAS ABUNDANCIAS ####
phylums = tax_glom(ph, taxrank = "phylum")
phylums@tax_table
lolo = abundances(phylums, transform = "compositional")
taxons_tab =  as.data.frame(phylums@tax_table)
row.names(lolo) = as.data.frame(phylums@tax_table)$phylum
colnames(lolo) = phylums@sam_data$tf.sp
write.csv(lolo, "relative_phylum.csv")

##### A NIVEL DE GÉNERO ####
# Getting top taxa on a Phylum level
phtse = transformAssay(phtse, method="relabundance")
phtse_genus <- agglomerateByRank(phtse, rank ="genus", onRankOnly=TRUE)
top_taxa <- getTop(phtse_genus,top = 50, assay.type = "relabundance")

# Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
genus_renamed <- lapply(rowData(phtse)$genus,
                         function(x){if (x %in% top_taxa) {x} else {"Other"}})
rowData(phtse)$genus <- as.character(genus_renamed)

# Compositional barplot
top_gen_plot = plotAbundance(phtse, assay.type="counts", rank = "genus",
                             order_rank_by="abund", order_sample_by = "tf.sp", 
                             add_legend = FALSE, features = "tf.sp")
top_gen_plot[[1]] = plotAbundance(phtse, assay.type = "counts", rank = "genus",
                                  features = "tf.sp", order_sample_by = "tf.sp",
                                  add_legend = FALSE)[[1]]
 
png(file = "Relative_genus.png", width = 12, height = 9, units = "in", res = 300)
wrap_plots(top_gen_plot, ncol = 1, heights = c(1,0.05))
invisible(dev.off())
top_gen_plot = plotAbundance(phtse, assay.type="counts", rank = "genus",
                             order_rank_by="abund", order_sample_by = "tf.sp", 
                             add_legend = TRUE)
png(file = "Legend_genus.png", width = 36, height = 27, units = "in", res = 300)
top_gen_plot
invisible(dev.off())


genusi = tax_glom(ph, taxrank = "genus")
genusi@tax_table
lolo = abundances(genusi, transform = "compositional")
taxons_tab =  as.data.frame(genusi@tax_table)
row.names(lolo) = as.data.frame(genusi@tax_table)$genus
colnames(lolo) = phylums@sam_data$tf.sp
write.csv(lolo, "relative_genus.csv")


#### GENERATE THE TABLES
# The only thing missing is to generate the tables for the phyla and the genera to add to the supplementary
# material. This can be done in a very roundabout way without having to write a parser making use of a couple 
# of libraries that have this functionality buried in the documentation
library("MiscMetabar")
library("ComplexUpset")
library("jsonlite")

setwd(res_dir)

# Let's first do the phyla tables
ph_phyla = tax_glom(ph, taxrank = "phylum", NArm = FALSE)
ph_phyla = ps_tax_clean(ph_phyla)
ph_genus = tax_glom(ph, taxrank = "genus", NArm = FALSE)
ph_genus = ps_tax_clean(ph_genus)
# Fix the phylum names
lala = ph_phyla@tax_table
lala = as.data.frame(lala)
taxa_names(ph_phyla) = lala$genus
# and fix the genus names
lala = ph_genus@tax_table
lala = as.data.frame(lala)
taxa_names(ph_genus) = lala$genus

# Generate the jsons
phyla_typeg = ps_venn(ph_phyla, group = "type_g", plot = FALSE)
phyla_typeg = toJSON(phyla_typeg, pretty = TRUE)
write(phyla_typeg, "phyla_typeg_all.json")
genus_typeg = ps_venn(ph_genus, group = "type_g", plot = FALSE)
genus_typeg = toJSON(genus_typeg, pretty = TRUE)
write(genus_typeg, "genus_typeg_all.json")
# repeat for type_f, used to make tables
phyla_typef = ps_venn(ph_phyla, group = "type_f", plot = FALSE)
phyla_typef = toJSON(phyla_typef, pretty = TRUE)
write(phyla_typef, "phyla_typef_all.json")
genus_typef = ps_venn(ph_genus, group = "type_f", plot = FALSE)
genus_typef = toJSON(genus_typef, pretty = TRUE)
write(genus_typef, "genus_typef_all.json")

rm(zotu_typef)
rm(zotu_typeg)
rm(genus_typeg)
rm(genus_typef)