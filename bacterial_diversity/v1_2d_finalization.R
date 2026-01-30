## Este script busca coger los objetos .RDS generados a lo largo del script 2 y convertirlos a un formato
## que pueda cargar de rebote a phyloseq directo, y ahí ya funcionar

library(dada2)
library(DECIPHER)
#library(phyloseq)

# Antes de nada, arreglar el path
script_dir = getwd()
setwd("../data/metadata")
meta_dir = getwd()
setwd("../r_backup")
backup_dir = getwd()

# Primero, cargar todos los objetos que voy a querer convertir
#length_tab = readRDS("lenght_tab.rds")
asv_tab = readRDS("chimless_tab_merge.rds")
tax_info = readRDS("tax_info_d")
#track_filt = read.csv("../metadata/track_filt.csv", row.names = 1)
#track_dada = read.csv("../metadata/dada_track.csv", row.names = 1)
#info_tab = read.csv("../metadata/metadata.csv", sep = ";")

# Esto aquí no debiera estar, pero hasta moverlo de sitio: vamos a terminar la tabla de tracking
row.names(track_filt) = row.names(track_dada) #esto sólo es posible pq tienen el mismo orden, lo suyo sería generar las columnas con un nombre en condiciones de entrada
full_track = cbind(track_filt, track_dada)
full_track["pct_retained"] = (full_track$nonchim/full_track$reads.in) * 100
write.csv(full_track, "../metadata/full_track.csv")

# Estos debieran ser todos. Ahora, con estos datos hemos de generar 3 nuevos objetos en formato csv:
# 1. Archivo de texto con las secuencias de los ASV en formato fasta
# 2. Un .csv con las frecuencias relativas de cada ASV por muestra
# 3. Una tabla con las clasificaciones taxonómicas de cada ASV

# 1. FASTA de las secuencias de los ASV
asv_seqs = colnames(asv_tab) #saca la secuencia en si
asv_headers = paste(">ASV", 1:ncol(asv_tab), sep = "_") #ASV_1, ASV_2, ASV_3, etc
asv_fasta = c(rbind(asv_headers, asv_seqs))
write(asv_fasta, file = file.path("../../results", "ASVs_fasta.fa"))

# 2. csv con las frecuencias relativas. Véase, la tabla asv_tab, pero transpuesta
t_asv = t(asv_tab)
row.names(t_asv) = sub(">", "", asv_headers)
write.csv(t_asv, "../../results/ASV_count.csv", sep = ";", quote = TRUE, col.names = NA)

# 3. Tabla con las clasificaciones taxonómicas de cada ASV
# Esta parte es más graciosa que las anteriores, puesto que el objeto que tengo en su formato actual
# no me sirve: primero he de pasarlo a un formato que dada entiende

# Estos son los rangos que generaría dada
ranks = c("domain", "phylum", "class", "order", "family", "genus", "species") #la especie no debiera tener ninguna
asv_tax = t(sapply(as.table(tax_info), function(x) { #transpón al resultado de aplicar a tax_info
  m = match(ranks, x$rank) # te quedas con lo que coincida en mis rangos y los del objeto taxa original
  taxa = x$taxon[m] # filtras para que los nuevos taxones sean los que hayan coincidido
  taxa[startsWith(taxa, "unclassified_")] <- NA # si hay alguno sin clasificar, lo paso a NA (lo que pondría dada)
  taxa #me devuelves el nuevo objeto taxa generado
}))

colnames(tax_info) = ranks #nombres de columnas = ranks
rownames(tax_info) = gsub(pattern=">", replacement="", x=asv_headers) # cambia > por nada en los nombres de cada ASV, está el > ahí por tema formato de .fa
write.csv(tax_info, "../../results/ASVs_taxonomy.csv", sep = ";", quote=TRUE, col.names=NA)


