## OBJETIVO DE ESTE SCRIPT: llegados a este punto ya debiera tener un objeto phyloseq con las frecuencias de los ASV, la clasificación 
## taxonómica y los metadatos dentro. Esto me permite llevar a cabo los dos primeros filtrados:
## 1. El filtrado de secuencias contaminantes
## 2. El filtrado de solterones
## Una vez tenga estos dos filtrados podré generar una tabla actualizada de frecuencias y un archivo de secuencias FASTA definitivo.
## Con ellos podré calcular las distancias eucliadianas entre muestras, y con las distancias euclidianas, un arbol taxonómico: 
## nuestro objeto phyloseq será definitivo al añadirle el árbol. Ahora bien, calcular distancias euclidianas a partir de datos de microbioma
## requiere tirar de DESeq2, lo cual no es tan sencillo como pudiera parecer.
library(data.table)
library(phyloseq)
#library(DESeq2)
library(ggplot2)
library(ggbeeswarm)
library(ggrepel)
library(vegan)
library(tidyverse)
library(seqateurs)
library(Biostrings)
library(phangorn)
library(DECIPHER)
library(decontam)

# Como siempre, lo primero: establecer los directorios por comodidad
script_dir = getwd()
setwd("../data/r_backup")
rback_dir = getwd()
setwd("../tree_building")
tree_dir = getwd() #en este directorio iré guardando los archivos relacionados con crear un árbol
setwd("../../results/phy_results")
phyr_dir = getwd() #aquí iran los resultados nuevos


# Cargar nuestros archivos de interés, empezando por phystuff
setwd(rback_dir)
phystuff_raw = readRDS("phystuff_raw.rds")
#phystuff_decont = readRDS("phystuff_decont.rds")

# Una visualización rápida: nº de reads por muestra, colores por tipología
df = as.data.frame(sample_data(phystuff_raw)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize = sample_sums(phystuff_raw)
df = df[order(df$LibrarySize),]
df$Index = seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=tf.sp)) + geom_point()

# Identificar contaminantes por su frecuencia
contamdf.freq <- isContaminant(phystuff_raw, method="frequency", conc="ADN") #si son contaminantes, saldrán como TRUE
table(contamdf.freq$contaminant)
contams_id = which(contamdf.freq$contaminant)
#plot_frequency(phystuff_raw, taxa_names(phystuff_raw)[c(1049,1050)], conc="ADN") + 
#  xlab("DNA Concentration (QBit intensity)")
# Y ya quitarlos
phystuff_decont = prune_taxa(!contamdf.freq$contaminant, phystuff_raw)
saveRDS(phystuff_decont, "phystuff_decont.rds")

# Aquí ya filtro para quitarme solterones. El criterio de entrada es: 
# - SI ESTÁ IDENTIFICADO A NIVEL DE CLASE:
#   - Mínimo 3 copias en total
# - SI NO ESTÁ IDENTIFICADO A NIVEL DE CLASE
#   - Mínimo 20 copias
# Nuevamente: muestras antárticas, conservar lo máximo posible. Eso no quita que los cloroplastos y las mitocondrias
# las quite y punto

# Primero, quito aquellas secuencias no identificadas a nivel de clase.
ps0 = subset_taxa(phystuff_decont, 
                  domain == "Bacteria" 
                  & !is.na(phylum) & !is.na(class)
                  & class != "Chloroplast" & family != "Mitochondria") #ps0 = phystuff pero sin NAs
tdt = data.table(setDT(as.data.frame(tax_table(ps0))), 
                 TotalCounts = taxa_sums(ps0), SV = taxa_names(ps0))
ggplot(tdt, aes(TotalCounts)) + geom_histogram(bins = 50) + theme_bw() + 
  ggtitle("Histogram of Total Counts")

tdt[(TotalCounts == 1), .N] #nº de solterones
tdt[(TotalCounts == 2), .N] #nº de dobletones

# taxa cumulative sum
taxcumsum = tdt[, .N, by = TotalCounts]
setkey(taxcumsum, TotalCounts)
taxcumsum[, CumSum := cumsum(N)]

# Graficar la suma del total de ASVs vs la abundancia de cada ASV
pCumSum = ggplot(taxcumsum, aes(TotalCounts, CumSum)) + geom_point() + theme_bw() + 
  xlab("Filtering Threshold") + ylab("zOTU Filtered")

png("../../results//phy_results/test_abundance.png", width = 8, height = 8, units = "in", res = 300)
gridExtra::grid.arrange(pCumSum, pCumSum + xlim(0, 500), 
                        pCumSum + xlim(0, 100), pCumSum + xlim(0, 50), nrow = 2, 
                        top = "zOTUs that would be filtered vs. minimum taxa counts threshold")
invisible(dev.off())

# Ahora una tabla de prevalencia. La hago un poco por arrastre, pq, nuevamente, lo que aquí llegue en muestras antárticas me da igual, no 
# voy a filtrar por prevalencia
mdt = fast_melt(ps0) #esto convierte el objeto phyloseq en una data.table de formato largo (cada fila es una observación)
#mdt = fast_melt(phystuff_decont)
# nueva tabla, con la prevalencia, el nº total de apariciones en todas las muestras y el máximo entre muestras
prevdt = mdt[, list(Prevalence = sum(count > 0), TotalCounts = sum(count), 
                    MaxCounts = max(count)), by = taxaID]
#un histograma con la misma idea
ggplot(prevdt, aes(Prevalence)) + geom_histogram(bins = 50) + theme_bw() +
  ggtitle("Histogram of Taxa Prevalence")
#singletons y doubletons
prevdt[(TotalCounts == 1), .N]
prevdt[(TotalCounts == 2), .N]
#histograma con el nº máximo de apariciones por muestra
ggplot(prevdt, aes(MaxCounts)) + geom_histogram(bins = 100) + xlim(0, 500) + theme_bw() + 
  ggtitle("Histogram of Maximum TotalCounts")
#y repito
prevdt[(MaxCounts == 1), .N]
prevdt[(MaxCounts == 2), .N]

# Cumsum de taxones
prevcumsum = prevdt[, .N, by = Prevalence]
setkey(prevcumsum, Prevalence)
prevcumsum[, CumSum := cumsum(N)]

# Cumsum de los ASV vs la prevalencia
pPrevCumSum = ggplot(prevcumsum, aes(Prevalence, CumSum)) + geom_point(size = 2, alpha = 0.5) + 
  theme_bw() + xlab("Filtering Threshold") + ylab("zOTUs Filtered") + 
  ggtitle("zOTUs that would be filtered vs. minimum sample count threshold")

png("../../results/phy_results/Filter-taxa_cumsum.png", width = 8, height = 8, units = "in", res = 300)
pPrevCumSum
invisible(dev.off())

# Y ahora lo mismo, pero a lo bestia: en vez de cumsum, scatterplot
png("../../results//phy_results/Filter-taxa_total.png", width = 8, height = 8, units = "in", res = 300)
ggplot(prevdt, aes(Prevalence, TotalCounts)) + geom_point(size = 2, alpha = 0.5) + 
  scale_y_log10() + theme_bw() + xlab("Prevalence [No. Samples]") + ylab("TotalCounts [Taxa]")
invisible(dev.off())


# Por último, puedo hacer esto mismo pero con un nivel más de complejidad: diferenciado por filo
addPhylum = unique(copy(mdt[, list(taxaID, phylum)])) #tablita corta que tiene el nº de ASV (índice) y el filo

# Junto vía .join para poder diferenciar las prevalencias por filo
setkey(prevdt, taxaID)
setkey(addPhylum, taxaID)
prevdt = addPhylum[prevdt]
setkey(prevdt, phylum)

png("../../results/phy_results/Filter_taxa_total_phylum.png", width = 8, height = 8, 
    units = "in", res = 300)
ggplot(prevdt, aes(Prevalence, TotalCounts, color = phylum)) + 
  geom_point(size = 1, alpha = 0.5) + scale_y_log10() + theme_bw() + 
  facet_wrap(~phylum, nrow = 4) + theme(legend.position="none") + 
  xlab("Prevalence [No. Samples]") + ylab("Total Abundance")
invisible(dev.off())

# Aunque para la muestra, por ser de naturaleza antártica, venía dado, ya he podido sacar mis puntos de filtrado a base de observar 
# las gráficas
preval_thres = 0 # above this prevalence in all samples ( = at least 1 sample presents it)
total_thres = 3 # above this number of copies per sample (= no sample has less than 3 copies of it)
max_thres = 2 # above this maximum on at least one sample (= at least one sample has two copies)

# Y filtrar en base a abundancia
keepTaxa = prevdt[(Prevalence > preval_thres & 
                     TotalCounts > total_thres & 
                     MaxCounts > max_thres), taxaID] #taxones a conservar

ps1 = prune_taxa(keepTaxa, ps0) #ps1 = phystuff, sin contaminantes, NAs o singletons/doubletons

# Hasta ahora he filtrado en base a los criterios para muestras identificadas al nivel de clase. Sin embargo, 
# he ignorado por completo las muestras no identificadas a nivel de clase (las cuales debían ser > 20)
ps0_alt = subset_taxa(phystuff_decont, is.na(class) & domain == "Bacteria") #quédate todo lo que no esté identificado a nivel de clase
#OJO QUE TAMBIÉN ME QUITO LOS NO IDENTIFICADOS A NIVEL DE DOMINIO
ps0_alt = prune_taxa(taxa_sums(ps0_alt)>= 20, ps0_alt) #mantén sólo lo superior a 20

ps1 = merge_phyloseq(ps1, ps0_alt) #junta ambos


write.table(as.data.table(otu_table(ps1), keep.rownames = T), file = "../../results/phy_results/zOTU_counts_filtered.csv", 
            sep = ";", quote = F, row.names = F, col.names = T)
write.table(as.data.table(tax_table(ps1), keep.rownames = T), file = "../../results/phy_results/zOTU_taxonomy_filtered.csv", 
            sep = ";", quote = F, row.names = F, col.names = T)

saveRDS(ps1, "phystuff_1.rds")

## ¿Qué me queda? Tras la asignación taxonómica definitiva, me quedará sólamente una serie de secuencias
## Como en estas muestras prima el conservar el máximo nº de ASVs, aquí no vuelvo a filtrar, pero tampoco afectará
## (siempre puedo re-filtrar aun teniendo árbol, sólo pierdo más tiempo de cómputo calculando el árbol)
final_seqs = refseq(ps1)
alignment = AlignSeqs(final_seqs, 
                      anchor = NA,
                      processors = NULL,
                      iterations = 3)
# Y exporto el alineamiento
setwd(tree_dir)
phang_align = phyDat(as(alignment, "matrix"), type = "DNA")
write.phyDat(phang_align, file = "alignment.fasta", format = "fasta") #Ojo aquí a dónde guardo el archivo

rm(phang_align)
rm(alignment)
rm(final_seqs)

# Ahora hago el árbol
fasttree_dir = "/usr/bin/fasttreeMP" #hardcodeado a lo bestia, ojo
system2(fasttree_dir, args = c("-nt", "-gtr", "-log", "-gamma",  "alignment.fasta"), stdout = "treefasted.tre")
#raxml_ng_dir = "/home/pak/Desktop/raxml-ng/raxml-ng" #hardcodeado a lo bestia, ojo
#system2(raxml_ng_dir, args = c("--evaluate", "--force", "--seed 1234", "--log progress",
#                               "--msa /home/pak/Desktop/p_osp_lakes/data/tree_building/alignment.fasta", "--msa-format FASTA",
#                               "--model GTR+G", "--brlen scaled", 
#                               "--tree /home/pak/Desktop/p_osp_lakes/data/tree_building/treefasted.tre", "--prefix eval_tree"))

# Y ya por fin, termino el objeto phyloseq
TREE = read_tree("eval_tree.raxml.bestTree")
TREE = read_tree("treefasted.tre")
phystuff_1 = merge_phyloseq(ps1, TREE)
setwd(rback_dir)
saveRDS(phystuff_1, "phystuff_1")




