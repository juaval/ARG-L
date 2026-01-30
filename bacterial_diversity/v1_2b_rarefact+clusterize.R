## OBJETIVO DE ESTE SCRIPT ##
## Llegados a este punto tenemos una serie X de ASVs, cada cual con n repeticiones por muestra. Sin embargo,
## la naturaleza del algoritmo de dada impide cazar falsos positivos en dichas muestras que provengan de errores
## de PCR en cualquiera de las múltiples PCRs que el método conlleva. Una solución simple a este problema es 
## clusterizar las ASVs (¿pasando a zOTUs?) y juntar las n repeticiones de todas las ASVs que acaben juntas en 
## su zOTU correspondiente. Esto, así como un análisis de correlación para determinar si realmente un falso 
## positivo puede considerarse como tal, nos lo da hecho lulu

## Antes de ello, un aparte: no tiene sentido ponernos a clusterizar si sabemos de antemano que el resultado
## no es representativo. Por ello, antes (¡y después!) de clusterizar se comprobará la representatividad de la
## secuenciación con una curva de rarefacción


library(vegan)
library(lulu)
library(stringr)

script_dir = getwd()
setwd("../data/r_backup")
rback_dir = getwd()
setwd("../../results/dada_results")
results_dir = getwd()

setwd(results_dir)
count_tab = read.table("ASV_count.csv", header=T, row.names=1,
                        check.names=F, sep=",")
ASV_matchlist = read.table("../../data/zOTU_building/ASV_matchlist.txt", header=F,
                           as.is=TRUE, stringsAsFactors=FALSE)

# Recordatorio: la curva de rarefacción nos dice cuántas ASVs/OTUs observamos cuando consideramos nuevas lecturas
# de la misma muestra. Si tenemos suficiente profundidad de secuenciación, veremos una meseta: incluso si metiésemos
# nuevas lecturas se corresponderían con ASV/OTUs ya presentes, la muestra es representativa en cuanto a diversidad

png("rarefact_ASVs.png", width = 1280, height = 960, unit = "px")
rarecurve(t(count_tab), step=1000, lwd=2, ylab="ASVs", label=T)
# Representar una línea en la profundidad de secuencia de la muestra con la menor cantidad de secuencias
abline(v=(min(rowSums(t(count_tab))))) 
dev.off()

## LULU DE AQUÍ A ABAJO ##
# Previamente, he creado la matchlist que pide el paquete vía usar blastn
lulu_res = lulu(count_tab, ASV_matchlist, minimum_ratio_type = "min", minimum_ratio = 1,
                minimum_match = 84, minimum_relative_cooccurence = 0.95) #tarda un pico en correr
setwd(rback_dir)
saveRDS(lulu_res, "lulu_res")
#lulu_res = readRDS("lulu_res")
setwd(results_dir)
png("rarefact_zOTUs.png", width = 1280, height = 960, unit = "px")
rarecurve(t(lulu_res$curated_table), step=1000, lwd=2, ylab="zOTUs", label=T)
# Representar una línea en la profundidad de secuencia de la muestra con la menor cantidad de secuencias
abline(v=(min(rowSums(t(lulu_res$curated_table))))) 
dev.off()




# El resultado final devuelve una lista de ASVs definitivas, una tabla de cuentas en orden alfabético (que no en orden de abundancia)
# y con los nombres del ASV "parental". Así pues, me parece que lo correcto sería rehacer el fasta de las ASVs definitivas para tener
# uno en el que sólo permanezcan las "ASVs" definitivas y un archivo de cuentas en el que las ASVs estén ordenadas por abundancia.
# Y, para terminar de rizar el rizo: que dejen de llamarse ASVs y pasen a llamarse "zOTUs"

# Si quiero que estén ordenados los archivos .fasta final y los zOTUs vean reflejado en su nombre el orden, lo primero es ordenar
zOTU_tab = lulu_res$curated_table
zOTU_tab["total"] = rowSums(zOTU_tab) #suma todas las columnas
zOTU_tab = zOTU_tab[order(-zOTU_tab$total), ] #ordena según el total de la fila
zOTU_tab = subset(zOTU_tab, select = -c(total)) #borra el total de la fila
zOTU_tab
rm(lulu_res)
rm(ASV_matchlist)

# Ahora ya puedo generar la equivalencia ASV_x a zOTU_y
zotu_headers = paste("zOTU", 1:nrow(zOTU_tab), sep = "_") #primero, los nombres definitivos
parent_asv_headers = rownames(zOTU_tab)
trans_df = rbind(zotu_headers, parent_asv_headers)
setwd(results_dir)
original_fasta = readLines("ASVs_fasta.fa")

# Con esto obtendré un vector con las secuencias de los ASVs ordenadas en su nuevo orden de zOTUs por frecuencia
final_seqs = c()
#count = 0 #this is here for debugging purposes
for (header in parent_asv_headers) {
  asv_num = as.integer(str_split(header, "_")[[1]][2])
  print(header)
  asv_seq = original_fasta[(asv_num*2)]
  final_seqs = append(final_seqs, asv_seq)
  #count = count + 1
  #if (count == 5){break}
}
# Y ya hacer el archivo en sí
zotu_fasta = c(rbind(zotu_headers, final_seqs))

write(zotu_fasta, file = file.path("zOTUs_fasta.fa"))

# Y ya terminar de apañar las filas del archivo de cuentas
row.names(zOTU_tab) = zotu_headers
write.table(zOTU_tab, "zOTU_count.csv", sep = ";", quote = F, col.names = NA)
write.table(trans_df, "zOTU_to_ASV_equivalence.csv", sep = ";", quote = F, col.names = NA)

# Me facilita bastante la vida, de cara a hacer la asignación taxonómica con DECIPHER, generar un último objeto con el mismo
# formato que tiene dada2 en su tutorial
zOTU_tab = t(zOTU_tab)  
colnames(zOTU_tab) = final_seqs
setwd(rback_dir)
saveRDS(zOTU_tab, "zotu_to_taxo")
