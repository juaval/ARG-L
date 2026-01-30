#script rápido y sucio para hacer la asignación taxonómica
library(dada2)
library(DECIPHER)

# poner el wd en orden. Debiera estar en scripts ahora mismo
script_dir = getwd()
setwd("../data/r_backup")
rback_dir = getwd()

# Ahora, abrir la tabla de asv sin quimeras
zotu_tab = readRDS("zotu_to_taxo")

# Cargar el clasificador
load("../databases/SILVA_SSU_r138_2019.RData")

# cambiar el tipo de objeto de las secuencias de los ASV a DNAStringSet
dna_seqs = DNAStringSet(getSequences(zotu_tab))

# Y clasificar <- This kills the crab
tax_info = IdTaxa(test = dna_seqs, trainingSet = trainingSet, 
                  strand = "both", processors = NULL, verbose = TRUE)


# Me guardo el objeto tipo taxa porsiaca
saveRDS(tax_info, "tax_info_zotu")
#tax_info = readRDS("tax_info_zotu")

# Y lo rearmo para que sea como los que genera dada
ranks = c("domain", "phylum", "class", "order", "family", "genus", "species") #la especie no debiera tener ninguna
zotu_tax = t(sapply(tax_info, function(x) { #transpón al resultado de aplicar a tax_info
  m = match(ranks, x$rank) # te quedas con lo que coincida en mis rangos y los del objeto taxa original
  taxa = x$taxon[m] # filtras para que los nuevos taxones sean los que hayan coincidido
  taxa[startsWith(taxa, "unclassified_")] <- NA # si hay alguno sin clasificar, lo paso a NA (lo que pondría dada)
  taxa #me devuelves el nuevo objeto taxa generado
}))
saveRDS(zotu_tax, "tax_info_zotu")
#zotu_tax = readRDS("tax_info_zotu")
colnames(zotu_tax) = ranks #nombres de columnas = ranks
zotu_headers = paste("zOTU", 1:nrow(zotu_tax), sep = "_")
rownames(zotu_tax) = zotu_headers # cambia > por nada en los nombres de cada ASV, está el > ahí por tema formato de .fa
write.csv(zotu_tax, "../../results/dada_results/zOTUs_taxonomy.csv", sep = ";", quote=TRUE, col.names=NA)

# Ahora voy a probar a hacer la asignación de especies con el comodísimo add-on de addSpecies. Si funciona, chachi, si no, ni
# me molesto.
zotu_tax.plus = zotu_tax
rownames(zotu_tax.plus) = getSequences(zotu_tab)
zotu_tax.plus = addSpecies(taxtab = zotu_tax.plus, refFasta = "../databases/silva_species_assignment_v138.1.fa.gz", verbose = TRUE, tryRC = TRUE)
