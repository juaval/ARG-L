# Trying to make a venn diagram with proportional set sizes using ggvenn

##### TRACING GENERA ACROSS SAMPLES  #######
library("dplyr")
library("eulerr")
library("patchwork")
library("tibble")
library("jsonlite")
library("tidyr")
library("ggplot2")

script_dir = getwd()
setwd("../results/phy_results/venn/all")
all_dir = getwd()
setwd("../inv")
inv_dir = getwd()
setwd("../sapo")
sapo_dir = getwd()
setwd("../beach")
beach_dir = getwd()

##### ALL ########
setwd(all_dir)
all_zotu = read_json("typeg_zotus.json", simplifyVector = TRUE)
all_plas = as.numeric(length(all_zotu$plastic))
all_sh = as.numeric(length(all_zotu$plastic__soil))
all_c = as.numeric(length(all_zotu$soil))

lala = euler(c("Plastic" = all_plas, "Soil" = all_c, "Plastic&Soil" = all_sh))
png("Regular all.png", width = 10, height = 6, units = "in", res = 300)
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
png("ALL-empty.png", width = 10, height = 6, units = "in", res = 300)
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




setwd(inv_dir)
inv_zotu = read_json("typeg_zotus_inv.json", simplifyVector = TRUE)
setwd(sapo_dir)
sapo_zotu = read_json("typeg_zotus_sapo.json", simplifyVector = TRUE)
setwd(beach_dir)
beach_zotu = read_json("typeg_zotus_alm.json", simplifyVector = TRUE)

###### INV #####
setwd(inv_dir)
inv_plas = as.numeric(length(inv_zotu$plastic))
inv_sh = as.numeric(length(inv_zotu$plastic__soil))
inv_c = as.numeric(length(inv_zotu$soil))

lala = euler(c("Plastic" = inv_plas, "Soil" = inv_c, "Plastic&Soil" = inv_sh))
png("Tracing inv.png", width = 10, height = 6, units = "in", res = 300)
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
png("Tracing inv-empty.png", width = 10, height = 6, units = "in", res = 300)
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


##### SAPO #####
setwd(sapo_dir)
sapo_plas = as.numeric(length(sapo_zotu$plastic))
sapo_sh = as.numeric(length(sapo_zotu$plastic__soil))
sapo_c = as.numeric(length(sapo_zotu$soil))
route1 = as.numeric(length(intersect(sapo_zotu$plastic, inv_zotu$plastic__soil)))

lala = euler(c("Plastic" = sapo_plas, "Soil" = sapo_c, "Greenhouse Shared" = 0, 
               "Plastic&Soil" = sapo_sh, "Plastic&Greenhouse Shared" = route1),
             shape = "ellipse")
png("Tracing sapo.png", width = 10, height = 6, units = "in", res = 300)
plot(lala,counts = TRUE, font=2, cex=1, alpha=0.8,
     quantities = list(type = c("percent", "counts"), 
                       col = c("gray5", "gray5","white", "gray5"),
                       cex = 1.5
                       ),
     labels = list(col = "gray5",
                   cex = 1.5
                   ),
     edges = list(col = c("gray5", "gray5", "black", "gray5"),
                  lex = c(3, 3, 1.5,3)
                  ),
     fills = c("#ff7f0e","#1f77b4","firebrick", "gray20"))
dev.off()
png("Tracing sapo-EMPTY.png", width = 10, height = 6, units = "in", res = 300)
plot(lala,counts = TRUE, font=2, cex=1, alpha=0.8,
     quantities = list(type = c("percent", "counts"), 
                       col = c("gray5", "gray5","white", "gray5"),
                       cex = 0
     ),
     labels = list(col = "gray5",
                   cex = 0
     ),
     edges = list(col = c("gray5", "gray5", "black", "gray5"),
                  lex = c(3, 3, 1.5,3)
     ),
     fills = c("#ff7f0e","#1f77b4","firebrick", "gray20"))
dev.off()

##### BEACH ####
setwd(beach_dir)
beach_plas = as.numeric(length(beach_zotu$plastic))
beach_sh = as.numeric(length(beach_zotu$plastic__soil))
beach_c = as.numeric(length(beach_zotu$soil))
route1 = intersect(sapo_zotu$plastic, inv_zotu$plastic__soil)
route1 = as.numeric(length(intersect(beach_zotu$plastic, route1)))
route2 = as.numeric(length(intersect(beach_zotu$plastic, inv_zotu$plastic__soil)))
route3 = as.numeric(length(intersect(beach_zotu$plastic, sapo_zotu$plastic)))

lala = euler(c("Plastic" = beach_plas, "Soil" = beach_c, 
               "Main route" = 0, "Greenhouse shared" = 0, "Sapo plastic" = 0,
               "Plastic&Soil" = beach_sh, "Plastic&Main route" = route1, 
               "Plastic&Greenhouse shared" = route2, "Plastic&Sapo plastic" = route3),
             shape = "ellipse")
png("Tracing beach.png", width = 10, height = 6, units = "in", res = 300)
plot(lala, counts = TRUE, font=2, cex=1, alpha=0.8,
     quantities = list(type = c("percent", "counts"), 
                       col = c("black", "black","white", "black", "black", "black"),
                       cex = 1.5
                       ),
     labels = list(col = "black",
                   cex = 1.5
                   ),
     edges = list(col = c("gray5", "gray5", "black", "black", "black", "gray5"), 
                  lex = c(3, 3, 1.5, 1.5, 1.5, 3)
                  ),
     fills = c("#ff7f0e","#1f77b4", "firebrick", "firebrick", "firebrick","gray20")
     )
invisible(dev.off())
png("Tracing beach-EMPTY.png", width = 10, height = 6, units = "in", res = 300)
plot(lala, counts = TRUE, font=2, cex=1, alpha=0.8,
     quantities = list(type = c("percent", "counts"), 
                       col = c("black", "black","white", "black", "black", "black"),
                       cex = 0
     ),
     labels = list(col = "black",
                   cex = 0
     ),
     edges = list(col = c("gray5", "gray5", "black", "black", "black", "gray5"), 
                  lex = c(3, 3, 1.5, 1.5, 1.5, 3)
     ),
     fills = c("#ff7f0e","#1f77b4", "firebrick", "firebrick", "firebrick","gray20")
)
invisible(dev.off())