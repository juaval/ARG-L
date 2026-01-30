## OBJETIVO DE ESTE SCRIPT: generar, partiendo de los resultados de los archivos de dada, un objeto phyloseq con el que continuar trabajando
## Esto va a conllevar dos cosas a la vez: no solo generar un objeto phyloseq, sino 
library(phyloseq)
library(Biostrings)

# Antes de nada, como siempre, el juego del directorio
setwd("../results/dada_results")
dadar_dir = getwd()
setwd("../../data/metadata")
meta_dir = getwd()

# Cargar los archivos con la info. Debiera empezar en meta_dir
meta_tab = read.csv("metadatos.csv", sep = ",", row.names = 1)
meta_tab[order(row.names(meta_tab)), ] #por consistencia con luego las abundancias
#meta_tab$id = as.character(meta_tab$id) #un pequeño apaño para asegurar que sean iguales
setwd(dadar_dir)
zotu_freqs_tab = read.csv("zOTU_count.csv", sep = ";", header = TRUE, row.names = 1)
zotu_freqs_tab[order(colnames(zotu_freqs_tab)), ] #ídem
colnames(zotu_freqs_tab) = as.character(colnames(zotu_freqs_tab)) #3/4 de lo mismo
zotu_taxa_tab = read.csv("zOTUs_taxonomy.csv", header = TRUE, row.names = 1)
zotu_seqs = readDNAStringSet("zOTUs_fasta.fa")  #Esto debiera venir formateado de entrada

# Empiezo a meterlo en phyloseq paso a paso
ZOTU = otu_table(zotu_freqs_tab, taxa_are_rows = TRUE)
TAX = tax_table(zotu_taxa_tab)
META = sample_data(meta_tab)
SEQS = refseq(zotu_seqs)
#esto se cepilla el nombre de las columnas Y DE LAS FILAS de ASV_taxa_tab, asi que lo parcheo
colnames(TAX) = colnames(zotu_taxa_tab)
rownames(TAX) = rownames(zotu_taxa_tab)
print(unique(taxa_names(ZOTU) == taxa_names(TAX)))
rownames(ZOTU) = rownames(TAX)

## ¡¡¡¡PRIMER MERGE!!!!
phystuff = phyloseq(ZOTU, TAX, META, SEQS)
saveRDS(phystuff, "../../data/r_backup/phystuff_raw.rds")
# para chequear que funciona, tontísimamente: un barplot
plot_bar(phystuff, fill = "domain") # si se baja de dominio a rangos inferiores, el cangrejo es exterminado









