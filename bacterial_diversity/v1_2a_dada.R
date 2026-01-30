library(dada2)
library(Biostrings)


## El script para poder funcionar debe correr desde un directorio en el cual se encuentren las secuencias fw y rv ya 
## filtradas por calidad y recortadas. En este caso, dicho directorio es ../data/cutadapt/filtered --> CAMBIARLO
## Esto conlleva una consideración importante: las secuencias que entran en este punto del pipeline pudiendo haber perdido
## el solapamiento de 12 bases necesario para ser combinadas, por lo que siempre es necesario previsualizar el solapamiento
## antes de empezar a correr, no sea que mueran en el paso de merge tras horas corriendo código

setwd("../data/filtered_data")
filt_path = getwd()
setwd("../filtered_qi")
qi_path = getwd()
setwd("../metadata")
meta_path = getwd()
setwd("../../results")
result_path = getwd()

# Antes de nada, una manera sencilla de comprobar en cada momento cómo está discurriendo el proceso es consultar el % de secuencias
# que voy reteniendo: debiera conservar un nº similar al indicado por figaro (si esque lo corro) tras filtrar, y que este número cambie
# lo MÍNIMO en cada paso. Para ello, tiro de tabla:

setwd(meta_path)
track_tab_mp = read.csv("filt_step_mp.csv", row.names = 1) # abro la tabla de seguimiento del filtrado
track_tab_qi = read.csv("filt_step_qi.csv", row.names = 1)
pct_tab_mp = track_tab_mp #copio la tabla, me interesa solo la estructura
pct_tab_qi = track_tab_qi
pct_tab_mp$reads.in = 100 #al principio = 100% de lecturas
pct_tab_qi$reads.in = 100
pct_tab_mp$reads.out = (track_tab_mp$reads.out * pct_tab_mp$reads.in) / track_tab_mp$reads.in #calculo el siguiente pct
pct_tab_qi$reads.out = (track_tab_qi$reads.out * pct_tab_qi$reads.in) / track_tab_qi$reads.in #calculo el siguiente pct
print(mean(pct_tab_mp$reads.out)) #y esto para tener una primera idea de de dónde partimos
print(mean(pct_tab_qi$reads.out))

# Ya podemos empezar a machacar código. Primero, aprender los modelos de error de cada grupo de secuencias y representarlos
# para poder inspeccionar si están bien o mal. Esto requiere de, previamente, diferenciar archivos fw de rv
fw_seqs = sort(list.files(filt_path, pattern = "*_1", full.names = TRUE))
rv_seqs = sort(list.files(filt_path, pattern = "*_2", full.names = TRUE))
sample_names = sapply(strsplit(basename(fw_seqs), "_1.fastq"), `[`, 1)
sample_names = paste0("s_", as.character(sample_names))

# segunda parte de este bloque
fw_seqs_qi = sort(list.files(qi_path, pattern = "*_1", full.names = TRUE))
rv_seqs_qi = sort(list.files(qi_path, pattern = "*_2", full.names = TRUE))
qi_names = sapply(strsplit(basename(fw_seqs_qi), "_1.fastq"), `[`, 1)
qi_names = paste0("s_", as.character(qi_names))


# Ya errores
error_F = learnErrors(fw_seqs, multithread=TRUE, MAX_CONSIST = 30, verbose = TRUE)
saveRDS(error_F, "../r_backup/error_F.rds") # Voy guardando pq es código lento y así puedo retomar
png("../../results/fw_error.png", width = 1280, height = 960, unit = "px") #guarda el gráfico de ajuste del modelo de error
plotErrors(error_F, nominalQ = TRUE)
dev.off() #cierra el controlador gráfico
error_R = learnErrors(rv_seqs, multithread=TRUE, MAX_CONSIST = 30, verbose = TRUE)
saveRDS(error_R, "../r_backup/error_R.rds")
png("../../results/rw_error.png", width = 1280, height = 960, unit = "px")
plotErrors(error_R, nominalQ = TRUE)
dev.off() 

# Errores, segunda parte
qi_error_F = learnErrors(fw_seqs_qi, multithread=TRUE, MAX_CONSIST = 30, verbose = TRUE)
saveRDS(qi_error_F, "../r_backup/qi_error_F.rds") # Voy guardando pq es código lento y así puedo retomar
png("../../results/qi_fw_error.png", width = 1280, height = 960, unit = "px") #guarda el gráfico de ajuste del modelo de error
plotErrors(qi_error_F, nominalQ = TRUE)
dev.off() #cierra el controlador gráfico
qi_error_R = learnErrors(rv_seqs_qi, multithread=TRUE, MAX_CONSIST = 30, verbose = TRUE)
saveRDS(qi_error_R, "../r_backup/qi_error_R.rds")
png("../../results/qi_rw_error_1.png", width = 1280, height = 960, unit = "px")
plotErrors(qi_error_R, nominalQ = TRUE)
dev.off() 

# Segundo, aplicar el algoritmo dada perse para terminar de limpiar las secuencias.
dada_fw = dada(fw_seqs, err = error_F, multithread = TRUE, verbose = TRUE, pool = "pseudo")
saveRDS(dada_fw, "../r_backup/dada_fw.rds")
dada_rv = dada(rv_seqs, err = error_R, multithread = TRUE, verbose = TRUE, pool = "pseudo")
saveRDS(dada_rv, "../r_backup/dada_rv.rds")

# Idem para qiagen
dada_fw_qi = dada(fw_seqs_qi, err = qi_error_F, multithread = TRUE, verbose = TRUE, pool = "pseudo")
saveRDS(dada_fw_qi, "../r_backup/dada_fw_qi.rds")
dada_rv_qi = dada(rv_seqs_qi, err = qi_error_R, multithread = TRUE, verbose = TRUE, pool = "pseudo")
saveRDS(dada_rv_qi, "../r_backup/dada_rv_qi.rds")

#Nuevamente, nos interesa saber cómo está yendo el proceso, asi que queremos saber cómo vamos a nivel de nº de secuencias
# Para ello, necesitamos poder saber cuántas secuencias quedan
dada_fw = readRDS("../r_backup/dada_fw.rds")
dada_fw_qi = readRDS("../r_backup/dada_fw_qi.rds")
dada_rv = readRDS("../r_backup/dada_rv.rds")
dada_rv_qi = readRDS("../r_backup/dada_rv_qi.rds")
get_N = function(x) sum(getUniques(x)) #para contar
#la original
track_tab_mp = cbind(track_tab_mp, sapply(dada_fw, get_N), sapply(dada_rv, get_N))
colnames(track_tab_mp) = c("original", "filtered", "dada_fw", "dada_rv") #lo rehago para que no quede feo
pct_tab_mp["dada_fw"] = (track_tab_mp$dada_fw * pct_tab_mp$reads.in) / track_tab_mp$original
print(mean(pct_tab_mp$dada_fw))
pct_tab_mp["dada_rv"] = (track_tab_mp$dada_rv * pct_tab_mp$reads.in) / track_tab_mp$original 
print(mean(pct_tab_mp$dada_rv))
#las separadas
track_tab_qi = cbind(track_tab_qi, sapply(dada_fw_qi, get_N), sapply(dada_rv_qi, get_N))
colnames(track_tab_qi) = c("original", "filtered", "dada_fw", "dada_rv") #lo rehago para que no quede feo
pct_tab_qi["dada_fw"] = (track_tab_qi$dada_fw * pct_tab_qi$reads.in) / track_tab_qi$original
print(mean(pct_tab_qi$dada_fw))
pct_tab_qi["dada_rv"] = (track_tab_qi$dada_rv * pct_tab_qi$reads.in) / track_tab_qi$original 
print(mean(pct_tab_qi$dada_rv))



# Tercero, juntar las secuencias fw y rv para formar contigs/ASVs. Es aquí donde se va a ver
# si el corte y filtrado de secuencias ha sido adecuado: si no hay solapamiento de 12 nts,
# se va bastante al carajo. ESTO SE PUEDE CAMBIAR A ESTE NIVEL SI NO HAY ALTERNATIVA
merged_seqs = mergePairs(dada_fw, fw_seqs, dada_rv, rv_seqs, verbose = TRUE)
saveRDS(merged_seqs, "../r_backup/merged_seqs.rds")

merged_seqs_qi = mergePairs(dada_fw_qi, fw_seqs_qi, dada_rv_qi, rv_seqs_qi, verbose = TRUE)
saveRDS(merged_seqs_qi, "../r_backup/merged_seqs_qi.rds")


#yyy repetir el seguimiento
merged_seqs = readRDS("../r_backup/merged_seqs.rds")
merged_seqs_qi = readRDS("../r_backup/merged_seqs_qi.rds")

track_tab_mp = cbind(track_tab_mp, sapply(merged_seqs, get_N))
colnames(track_tab_mp) = c("original", "filtered", "dada_fw", "dada_rv", "merged") #lo rehago para que no quede feo
pct_tab_mp["merged"] = (track_tab_mp$merged * pct_tab_mp$reads.in) / track_tab_mp$original
print(mean(pct_tab_mp$merged))

track_tab_qi = cbind(track_tab_qi, sapply(merged_seqs_qi, get_N))
colnames(track_tab_qi) = c("original", "filtered", "dada_fw", "dada_rv", "merged") #lo rehago para que no quede feo
pct_tab_qi["merged"] = (track_tab_qi$merged * pct_tab_qi$reads.in) / track_tab_qi$original
print(mean(pct_tab_qi$merged))

# Cuarto: construir una tabla con las ASVs. 
# Esto nos permitirá un primer filtrado de solterones
asv_tab = makeSequenceTable(merged_seqs)
length_tab = table(nchar(getSequences(asv_tab))) #esto sirve para ver qué ASVs son únicos en base a una baja frecuencia
png("../../results/length_dist.png", width = 1280, height = 960, unit = "px")
plot(length_tab)
dev.off() 

# De nuevo, para qiagen
asv_tab_qi= makeSequenceTable(merged_seqs_qi)
length_tab_qi = table(nchar(getSequences(asv_tab_qi))) #esto sirve para ver qué ASVs son únicos en base a una baja frecuencia
png("../../results/length_dist_qi.png", width = 1280, height = 960, unit = "px")
plot(length_tab_qi)
dev.off()

# Quinto: purgar los "solterones". La verdad es que si uno se guiase sólamente por la frecuencia de las secuencias inferiría
# que el tamaño medio es 375, pero la realidad no es así: el amplicón mide 569 pbs, por lo que, para ser exactos, ninguna secuencia alcanza su tamaño
# (esto puede que no sea cierto, pues no termino de comprender cómo encaja el primer reverse aquí). Es por esta razón que, tras el pico en 375, hay un segundo
# pico a partir del 536: el pico de nuestro amplicón de verdad. La pena es que por la baja calidad de la parte final de las secuencias es muy, muy reducido
# por lo que perderé muchas secuencias en este punto. Pero me las creo todas
asv_tab = asv_tab[,nchar(colnames(asv_tab)) %in% 401:429] #ahora esto lo estoy estimando en base a correr en Rstudio bloque a bloque, lo que ha de cambiar necesariamente
asv_tab_qi = asv_tab_qi[,nchar(colnames(asv_tab_qi)) %in% 401:429]


# Y AQUÍ HAGO UN MERGE
dim(asv_tab)
rownames(asv_tab)
dim(asv_tab_qi)
rownames(asv_tab_qi)

merge_tab = mergeSequenceTables(asv_tab, asv_tab_qi) # unnamed arguments assumed to be sequence tables

dim(merge_tab)

# Sexto: quitar quimeras. Utilizo los parámetros básicos de entrada, luego me lo pienso más
merge_tab.nochim = removeBimeraDenovo(merge_tab, method = "consensus", multithread = TRUE, verbose = TRUE)
saveRDS(merge_tab.nochim, "../r_backup/chimless_tab_merge.rds")
(sum(merge_tab.nochim))/sum(merge_tab) #esto nos dice cuántas secuencias han sobrevivido el merge y cuántas no

# Séptimo: para comprobar qué ha ido ocurriendo en cada paso del proceso vía cuántas secuencias han 
# ido permaneciendo. Esta info está en distintos objetos que se han ido creando a lo largo del proceso es necesario 
# chequear la cantidad de "elementos" que los componen. Esto es lo que hace este bloque

# Ahora, antes de nada: juntar ambas tablas de seguimiento, tipo a tipo

# Esto es un chapuzote enorme, pero al mergear pierde el nombre de las filas ??
track_tab_mp["sam"] = row.names(track_tab_mp)
track_tab_qi["sam"] = row.names(track_tab_qi)
pct_tab_mp["sam"] = row.names(pct_tab_mp)
pct_tab_qi["sam"] = row.names(pct_tab_qi)
track_tab = merge(track_tab_mp, track_tab_qi, all = TRUE)
pct_tab = merge(pct_tab_mp, pct_tab_qi, all = TRUE)
track_tab = track_tab[order(track_tab$sam), ]
pct_tab = pct_tab[order(pct_tab$sam), ]
rownames(pct_tab) = pct_tab$sam
rownames(track_tab) = track_tab$sam
track_tab = subset(track_tab, select = -c(sam))
pct_tab = subset(pct_tab, select = -c(sam))

merge_tab.nochim = readRDS("../r_backup/chimless_tab_merge.rds")
merge_tab.nochim = merge_tab.nochim[order(row.names(merge_tab.nochim)), ] #esto hay que hacerlo también aquí porque se cambia el orden al mergear asv_tab's

track_tab = cbind(track_tab, rowSums(merge_tab.nochim))
colnames(track_tab) = c("original", "filtered", "dada_fw", "dada_rv", "merged", "chimless") #lo rehago para que no quede feo
pct_tab["chimless"] = (track_tab$chimless * pct_tab$reads.in) / track_tab$original
print(mean(pct_tab$chimless)) 
colnames(pct_tab) = c("original", "filtered","denoised_fw", "denoised_rv", "merged", "nonchim")
#rownames(track_tab) = sample_names
setwd(meta_path)
write.csv(track_tab, "full_track.csv", sep = ";")
write.csv(pct_tab, "pct_track.csv", sep = ";")
png("../../results/seqs_retained.png", width = 1280, height = 960, unit = "px")
plot(pct_tab)
dev.off() 

# Por último, antes de salir de dada y empezar a usar otras librerías, lo más sencillo va a ser guardar el output de dos formas:
# 1. Un archivo .fa con los asvs en formato fasta
# 2. Un archivo en formato csv con el nº de cada asv por muestra
#AQUÍ NECESITO JUNTAR SAMPLE_NAMES Y QI_NAMES
all_names = append(sample_names, qi_names) 
row.names(merge_tab.nochim) = all_names #por quitarme futuros dolores de cabeza
asv_seqs = colnames(merge_tab.nochim)
asv_headers = paste(">ASV", 1:ncol(merge_tab.nochim), sep = "_") #ASV_1, ASV_2, ASV_3, etc
# Entrelazarlos
asv_fasta = c(rbind(asv_headers, asv_seqs))
# Guardarlo
setwd(result_path)
write(asv_fasta, file = file.path("ASVs_fasta.fa"))
# Misma historia para un archivo con las cuentas
t_asv = t(merge_tab.nochim)
row.names(t_asv) = sub(">", "", asv_headers)
write.table(t_asv, "ASV_count.csv", sep = ";", quote = F, col.names = NA)

# Y ya debiera estar