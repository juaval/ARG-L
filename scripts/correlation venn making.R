####### SET THINGS
# Objective of the script: parse the correlation results to retain only the most relevant and then 
# generate venn diagrams of said results
 
##### JSON MAKING ####
library("phyloseq")
library("MicEco")
library("dplyr")
library("MiscMetabar")
library("ComplexUpset")
library("patchwork")
library("tibble")
library("jsonlite")
library("eulerr")

# directorios útiles
script_dir = getwd()
setwd("../results/Corr_res")
data_dir = getwd()
setwd("venn_corrs")
res_dir = getwd()

setwd(data_dir)
zap = read.csv("ZAP_corrs.csv") # Zotu - Arg in Plastics
zac = read.csv("ZAC_corrs.csv") # Zotu - Arg in Control 
zmp = read.csv("ZMP_corrs.csv") # Zotu - Mge in Plastics
zmc = read.csv("ZMC_corrs.csv") # Zotu - Mge in Control
amp = read.csv("AMP_corrs.csv") # Arg - Mge in Plastics
amc = read.csv("AMC_corrs.csv") # Arg - Mge in Control

# Let's filter the datasets to only keep the actual interesting results. A for loop doesn't work as is,
# as they turn into atomic vectors when inside the for loop
#print(dim(zap))
zap = zap[zap$coef > 0,]
zap = zap[zap$p.val <= 0.01, ]
#print(dim(zap))
zac = zac[zac$coef > 0,]
zac = zac[zac$p.val <= 0.01, ]
zmp = zmp[zmp$coef > 0,]
zmp = zmp[zmp$p.val <= 0.01, ]
zmc = zmc[zmc$coef > 0,]
zmc = zmc[zmc$p.val <= 0.01, ]
amp = amp[amp$coef > 0,]
amp = amp[amp$p.val <= 0.01, ]
amc = amc[amc$coef > 0,]
amc = amc[amc$p.val <= 0.01, ]

# Now I need a way to check whether each particular correlation is present in the opposing dataset (plastics vs controls)
# I can do this in a quite hacky way
zap$combo = paste0(zap$var1, " - ", zap$var2)
# df_list = list(zap, zac, zmp, zmc, amp, amc) I wish I had the time to make this work
# lapply(df_list, function(x) {x$combo = paste0(x$var1, " - ", x$var2)})
zmp$combo = paste0(zmp$var1, " - ", zmp$var2)
amp$combo = paste0(amp$var1, " - ", amp$var2)
zac$combo = paste0(zac$var1, " - ", zac$var2)
zmc$combo = paste0(zmc$var1, " - ", zmc$var2)
amc$combo = paste0(amc$var1, " - ", amc$var2)

# And now its just a matter of doing set operations. Thankfully dplyr exists
zap_shared = as.numeric(length(intersect(zap$combo, zac$combo)))
zap_plas = as.numeric(length(setdiff(zap$combo, zac$combo)))
zap_cont = as.numeric(length(setdiff(zac$combo, zap$combo)))

zmp_shared = as.numeric(length(intersect(zmp$combo, zmc$combo)))
zmp_plas = as.numeric(length(setdiff(zmp$combo, zmc$combo)))
zmp_cont = as.numeric(length(setdiff(zmc$combo, zmp$combo)))

amp_shared = as.numeric(length(intersect(amp$combo, amc$combo)))
amp_plas = as.numeric(length(setdiff(amp$combo, amc$combo)))
amp_cont = as.numeric(length(setdiff(amc$combo, amp$combo)))

# Let's start making the diagrams themselves
setwd(res_dir)
lala = euler(c("Plastic" = zap_plas, "Surrounding environment" = zap_cont, "Plastic&Surrounding environment" = zap_shared))
png("ZAP.png", width = 10, height = 6, units = "in", res = 300)
plot(lala, counts = TRUE, font=2, cex=1, alpha=0.8,
     quantities = list(type = c("percent", "counts"),
                       col = c("black", "black","white"),
                       cex = 1.5
     ),
     labels = list(col = "black",
                   cex = 1.5
     ),
     edges = list(col = c("gray5","gray5","white"), 
                  lex = 3
     ),
     fills = c("#ff7f0e","#1f77b4","gray20"))
invisible(dev.off())
# Empty variant
png("ZAP-empty.png", width = 10, height = 6, units = "in", res = 300) #these empty ones make adjusting the graphs in other software much, much easier
plot(lala, counts = TRUE, font=2, cex=1, alpha=0.8,
     quantities = list(type = c("percent", "counts"),
                       col = c("black", "black","white"),
                       cex = 0
     ),
     labels = list(col = "black",
                   cex = 0
     ),
     edges = list(col = c("gray5","gray5","white"), 
                  lex = 3
     ),
     fills = c("#ff7f0e","#1f77b4","gray20"))
invisible(dev.off())
#ZMP
lala = euler(c("Plastic" = zmp_plas, "Surrounding environment" = zmp_cont, "Plastic&Surrounding environment" = zmp_shared))
png("ZMP.png", width = 10, height = 6, units = "in", res = 300)
plot(lala, counts = TRUE, font=2, cex=1, alpha=0.8,
     quantities = list(type = c("percent", "counts"),
                       col = c("black", "black","white"),
                       cex = 1.5
     ),
     labels = list(col = "black",
                   cex = 1.5
     ),
     edges = list(col = c("gray5","gray5","white"), 
                  lex = 3
     ),
     fills = c("#ff7f0e","#1f77b4","gray20"))
invisible(dev.off())
# Empty variant
png("ZMP-empty.png", width = 10, height = 6, units = "in", res = 300) #these empty ones make adjusting the graphs in other software much, much easier
plot(lala, counts = TRUE, font=2, cex=1, alpha=0.8,
     quantities = list(type = c("percent", "counts"),
                       col = c("black", "black","white"),
                       cex = 0
     ),
     labels = list(col = "black",
                   cex = 0
     ),
     edges = list(col = c("gray5","gray5","white"), 
                  lex = 3
     ),
     fills = c("#ff7f0e","#1f77b4","gray20"))
invisible(dev.off())
#AMP
lala = euler(c("Plastic" = amp_plas, "Surrounding environment" = amp_cont, "Plastic&Surrounding environment" = amp_shared))
png("AMP.png", width = 10, height = 6, units = "in", res = 300)
plot(lala, counts = TRUE, font=2, cex=1, alpha=0.8,
     quantities = list(type = c("percent", "counts"),
                       col = c("black", "black","white"),
                       cex = 1.5
     ),
     labels = list(col = "black",
                   cex = 1.5
     ),
     edges = list(col = c("gray5","gray5","white"), 
                  lex = 3
     ),
     fills = c("#ff7f0e","#1f77b4","gray20"))
invisible(dev.off())
# Empty variant
png("AMP-empty.png", width = 10, height = 6, units = "in", res = 300) #these empty ones make adjusting the graphs in other software much, much easier
plot(lala, counts = TRUE, font=2, cex=1, alpha=0.8,
     quantities = list(type = c("percent", "counts"),
                       col = c("black", "black","white"),
                       cex = 0
     ),
     labels = list(col = "black",
                   cex = 0
     ),
     edges = list(col = c("gray5","gray5","white"), 
                  lex = 3
     ),
     fills = c("#ff7f0e","#1f77b4","gray20"))
invisible(dev.off())