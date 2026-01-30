#### ABUNDANCIA DIFERENCIAL ####
# Abundancia diferencial hace referencia a una pregunta muy simplona:
# ¿hay algún zOTU que sea más abundante en un tipo de muestra que en otra?
# Eso puede llevar a pensar que para qué nos hacemos esta pregunta, si ya sabemos
# que hay diferencias entre un tipo de muestra, lo hemos comprobado a nivel de alfa y beta diversidad.
# La clave reside en el hecho de que este tipo de test va zOTU a zOTU (o género a género, o ud. taxonómica a ud. taxonómica)
# asi que nos dice qué zOTUs son los responsables de las diferencias observadas.
# O, en el caso de que no lleguemos a observar diferencias, que diferencias existen a nivel de especie pero sin que lleguen
# a trascender al nivel de la comunidad 

# Ahora bien, al igual que sucede con la beta diversidad, no hay un método único para comprobar que ud. taxonómica es diferencialmente abundante
# hay muchos, muchísimos. Pero, a diferencia de lo que sucede con la beta diversidad, no todos están pensados para funcionar
# con datos de microbioma (hasta donde yo sepa, solo 4 lo hacen, los que uso aquí), y la comparación de resultados es bastante
# inmediata

#### LIBRERÍAS ####
library("mia")
library("phyloseq")
library("ALDEx2")
library("maaslin3")
library("MicrobiomeStat")
library("GUniFrac")
library("MicEco")
library("dplyr")
library("tibble")
library("stringr")
library(patchwork)
library(knitr)
#### DIRECTORIOS ####
script_dir = getwd()
setwd("../data/r_backup")
data_dir = getwd()
setwd("../../results/DA_markers")
result_dir = getwd()


#### PREPARACIÓN PREVIA DE LOS DATOS ####
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

# Drop water samples
ph = subset_samples(ph, type_f != "Water")


# Quick comparision of how will grouping by genera affect the dimensionality of the dataset
dim(ph@otu_table)[1]

# Quickly check how many zOTUs need their genera fixed
nas = sum(is.na(as.data.frame(ph@tax_table)$genus)) 
print(paste0("There's ", nas, " zOTUs with unassigned genus, ", 
             round((nas/dim(ph@otu_table)[1]) *100, digits = 3) , "% of the total"))
ph = ps_tax_clean(ph)
nas = sum(is.na(as.data.frame(ph@tax_table)$genus))
print(paste0("There's ", nas, " zOTUs with unassigned genus, ", 
             round((nas/dim(ph@otu_table)[1]) *100, digits = 3) , "% of the total"))

# From here on, I'm switching to mia 
phtse = convertFromPhyloseq(ph) 

# Rarefaction to deal with sequencing depth differences 
rarefaction_depth = assay(phtse, "counts") %>% 
  colSums() %>% 
  min() 
phtse = rarefyAssay(phtse, # So now, assay(phtse, 2) will be the rarefied counts
                     assay.type = "counts",
                     sample = rarefaction_depth,
                     replace = FALSE,
                     name = "rarefied")

# And agglomerate into genus
phgen = agglomerateByRank(phtse, rank = "genus")
phgen

#### ALDEx2 ####
setwd(result_dir)
set.seed(1312) 
aldex_data = aldex.clr(reads = assay(phtse, 2), 
                       conds = colData(phtse)$type_g, 
                       mc.samples = 128, 
                       denom = "all", 
                       verbose = TRUE
                       )

aldex_tt <- aldex.ttest(aldex_data,  
                        paired.test = FALSE, 
                        verbose = TRUE)

aldex_effect <- aldex.effect(aldex_data, CI = TRUE, verbose = TRUE)
aldex_out <- data.frame(aldex_tt, aldex_effect)

# Lets visualize aldex results
png("aldex2_results_zOTUs.png", width = 10, height = 6, units = "in", res = 300)
par(mfrow = c(1, 2))
aldex.plot(
  aldex_out, 
  type = "MA", 
  test = "welch", 
  xlab = "Log-ratio abundance",
  ylab = "Difference",
  cutoff = 0.05
)
aldex.plot(
  aldex_out, 
  type = "MW", 
  test = "welch",
  xlab = "Dispersion",
  ylab = "Difference",
  cutoff = 0.05
)
invisible(dev.off())

# Let's check which zOTUs appeared as significant
aldex_out_filt = aldex_out |> 
  rownames_to_column(var = "zOTU") |>
  filter(we.eBH <= 0.05)  |> 
  dplyr::select(zOTU, we.eBH, wi.eBH, effect, overlap) |>
  kable() 
aldex_out_filt


#### MAASLIN3 ####
maaslin3_out <- maaslin3(
  input_data = phtse,
  output = "maaslin_res_zOTU",
  formula = "~ type_g", 
  normalization = "TSS", 
  transform = "LOG",     
  verbosity = "ERROR"
)
maaslin3_out <- maaslin3_out[["fit_data_abundance"]][["results"]] 
maaslin3_out |>
  filter(qval_joint <= 0.05) |>
  kable() # Mind that this table has way more columns, so it won't seem quite right from the get go

#### LINDA ####
zotu_tab <- as.data.frame(assay(phtse, 2)) 
meta <- as.data.frame(colData(phtse)) %>% select(type_g) 

linda_res <- linda(feature.dat = zotu_tab, 
                   meta.dat = meta, 
                   formula = '~type_g', 
                   alpha = 0.05, 
                   prev.filter = 0, # Ojo a esto. No estamos aplicando ningún filtro más allá de la rarefacción, lo que afecta a los resultados 
                   verbose = TRUE)
linda.plot(linda_res, 
           variables.plot = "type_gplastic",
           legend = TRUE, 
           width = 16, height = 16, alpha = 0.05,
           titles = "Soil vs Plastic", 
           directory = "linda_res_zOTU")

#### ZICOSEQ ####
zotu_tab <- as.matrix(assay(phtse,2)) # ZicoSeq won't work with the rarefied data, as it gives an error relative to the meta which I'm wholly unable to decipher 
meta = as.data.frame(colData(phtse))
zicoseq_out <- ZicoSeq(feature.dat = as.matrix(assay(phtse)),
                       meta.dat = as.data.frame(colData(phtse)),
                       grp.name = "type_g",
                       feature.dat.type = "count",
                       return.feature.dat = TRUE,
                       prev.filter = 0,
                       mean.abund.filter = 0,
                       max.abund.filter = 0,
                       perm.no = 999,
                       outlier.pct = 0.01)
zicoseq_res <- cbind.data.frame(p.raw = zicoseq_out$p.raw,
                                p.adj.fdr = zicoseq_out$p.adj.fdr)
png("zicoseq_results_all.png", width = 12, height = 12, units = "in", res = 300)
ZicoSeq.plot(ZicoSeq.obj = zicoseq_out,
             pvalue.type = 'p.adj.fdr', cutoff = 0.05, text.size = 2)
invisible(dev.off())


#### JUNTAR TODOS LOS RESULTADOS ####
sum_table = rownames_to_column(aldex_out, "genus") %>%
  select(genus, aldex2 = wi.eBH)
sum_table = full_join(sum_table, select(maaslin3_out, genus = feature, maaslin3 = qval_joint) %>%
                        mutate(genus = str_remove(genus, "X")), by = "genus")
sum_table = full_join(sum_table, select(rownames_to_column(linda_res$output$`type_gSurrounding env.`, "genus"),
                                        genus, linda = padj), by = "genus")
sum_table = full_join(sum_table, select(rename(as.data.frame(rownames_to_column(zicoseq_res)), genus = rowname), 
                                        genus, zicoseq = p.adj.fdr), by = "genus")
filt_table = mutate(sum_table, across(c(aldex2, maaslin3, linda, zicoseq), ~ .x <= 0.05),
                    across(-genus, function(x) ifelse(is.na(x), FALSE, x)),
                    score = rowSums(across(c(aldex2, maaslin3, linda, zicoseq))))
common = filter(filt_table, score == 4)
common_3 = filter(filt_table, score == 3)

write.csv(sum_table, "fdr_pvals_DA_zOTU.csv")
write.csv(filt_table, "TF_DA_zOTU.csv")
write.csv(common, "only_true_zOTU.csv")
write.csv(common_3, "3true_zOTU.csv")



#### GENUS LEVEL ####
#### ALDEx2 ####
aldex_data = aldex.clr(reads = assay(phgen, 2), 
                       conds = colData(phgen)$type_g, 
                       mc.samples = 128, 
                       denom = "all", 
                       verbose = TRUE
)

aldex_tt <- aldex.ttest(aldex_data,  
                        paired.test = FALSE, 
                        verbose = TRUE)

aldex_effect <- aldex.effect(aldex_data, CI = TRUE, verbose = TRUE)
aldex_out <- data.frame(aldex_tt, aldex_effect)

# Lets visualize aldex results
png("aldex2_results_genus.png", width = 10, height = 6, units = "in", res = 300)
par(mfrow = c(1, 2))
aldex.plot(
  aldex_out, 
  type = "MA", 
  test = "welch", 
  xlab = "Log-ratio abundance",
  ylab = "Difference",
  cutoff = 0.05
)
aldex.plot(
  aldex_out, 
  type = "MW", 
  test = "welch",
  xlab = "Dispersion",
  ylab = "Difference",
  cutoff = 0.05
)
invisible(dev.off())

# Let's check which zOTUs appeared as significant
aldex_out_filt = aldex_out |> 
  rownames_to_column(var = "zOTU") |>
  filter(we.eBH <= 0.05)  |> 
  dplyr::select(zOTU, we.eBH, wi.eBH, effect, overlap) |>
  kable() 
aldex_out_filt


#### MAASLIN3 ####
maaslin3_out <- maaslin3(
  input_data = phgen,
  output = "maaslin_res_genus",
  formula = "~ type_g", 
  normalization = "TSS", 
  transform = "LOG",     
  verbosity = "ERROR"
)
maaslin3_out <- maaslin3_out[["fit_data_abundance"]][["results"]] 
maaslin3_out |>
  filter(qval_joint <= 0.05) |>
  kable() # Mind that this table has way more columns, so it won't seem quite right from the get go

#### LINDA ####
zotu_tab <- as.data.frame(assay(phgen, 2)) 
meta <- as.data.frame(colData(phgen)) %>% select(type_g) 

linda_res <- linda(feature.dat = zotu_tab, 
                   meta.dat = meta, 
                   formula = '~type_g', 
                   alpha = 0.05, 
                   prev.filter = 0, # Ojo a esto. No estamos aplicando ningún filtro más allá de la rarefacción, lo que afecta a los resultados 
                   verbose = TRUE)
linda.plot(linda_res, 
           variables.plot = "type_gplastic",
           legend = TRUE, 
           width = 16, height = 16, alpha = 0.05,
           titles = "Soil vs Plastic", 
           directory = "linda_res_genus")

#### ZICOSEQ ####
zotu_tab <- as.matrix(assay(phgen,2)) # ZicoSeq won't work with the rarefied data, as it gives an error relative to the meta which I'm wholly unable to decipher 
meta = as.data.frame(colData(phgen))
zicoseq_out <- ZicoSeq(feature.dat = as.matrix(assay(phgen)),
                       meta.dat = as.data.frame(colData(phgen)),
                       grp.name = "type_g",
                       feature.dat.type = "count",
                       return.feature.dat = TRUE,
                       prev.filter = 0,
                       mean.abund.filter = 0,
                       max.abund.filter = 0,
                       perm.no = 999,
                       outlier.pct = 0.01)
zicoseq_res <- cbind.data.frame(p.raw = zicoseq_out$p.raw,
                                p.adj.fdr = zicoseq_out$p.adj.fdr)
png("zicoseq_results_genus.png", width = 12, height = 12, units = "in", res = 300)
ZicoSeq.plot(ZicoSeq.obj = zicoseq_out,
             pvalue.type = 'p.adj.fdr', cutoff = 0.05, text.size = 2)
invisible(dev.off())


#### JUNTAR TODOS LOS RESULTADOS ####
sum_table = rownames_to_column(aldex_out, "genus") %>%
  select(genus, aldex2 = wi.eBH)
sum_table = full_join(sum_table, select(maaslin3_out, genus = feature, maaslin3 = qval_joint) %>%
                        mutate(genus = str_remove(genus, "X")), by = "genus")
sum_table = full_join(sum_table, select(rownames_to_column(linda_res$output$`type_gSurrounding env.`, "genus"),
                                        genus, linda = padj), by = "genus")
sum_table = full_join(sum_table, select(rename(as.data.frame(rownames_to_column(zicoseq_res)), genus = rowname), 
                                        genus, zicoseq = p.adj.fdr), by = "genus")
filt_table = mutate(sum_table, across(c(aldex2, maaslin3, linda, zicoseq), ~ .x <= 0.05),
                    across(-genus, function(x) ifelse(is.na(x), FALSE, x)),
                    score = rowSums(across(c(aldex2, maaslin3, linda, zicoseq))))
common = filter(filt_table, score == 4)
common_3 = filter(filt_table, score == 3)

write.csv(sum_table, "fdr_pvals_DA_genus.csv")
write.csv(filt_table, "TF_DA_genus.csv")
write.csv(common, "only_true_genus.csv")
write.csv(common_3, "3true_genus.csv")