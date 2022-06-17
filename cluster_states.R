# packages
library(sf)
library(dplyr)
library(vegan)
library(fastcluster)
library(distances)

# Donnees

vegpot <- read.csv("data/0_Vegpot.csv")

sp_ba <- readRDS("data/sp_mat_ba_oct2020.RDS")

# pep_pe <- st_read("data/PEP_GDB/PEP.gdb", layer = "STATION_PE")
# 
# # VEGPOT = prendre la valeur de la derniere mesure
# pep_pe <- pep_pe %>% 
#   group_by(ID_PE) %>%
#   arrange(NO_MES, .by_group = TRUE) %>% 
#   mutate(VEG_POT = last(na.omit(VEG_POT))) %>%
#   ungroup()
# 
# # enlever la PEP qui n'ont pas de VEGPOT et garder seulement les 26 vegpot qui sont dans artemis
# 
# pep_pe <- pep_pe %>% 
#   filter(!is.na(VEG_POT)) %>%
#   filter(VEG_POT %in% vegpot$VegPotName)
# 
# # enlever qq colonnes
# pep_pe <- pep_pe %>%
#   select(ID_PE:PERTURB, VEG_POT, GR_ESS,TYPE_ECO, COUV_ARBO, GRESP_INDI, TYPE_FOR)
# 
# saveRDS(pep_pe, "data/pep_pe.RDS")

pep_pe <- readRDS("data/pep_pe.RDS")
pep_vegpot <- pep_pe %>%
  select(ID_PE, ID_PE_MES, VEG_POT)


# garder seulement les PEP qui sont dans pep_pe
MySpecies <- colnames(sp_ba[,-c(1:3,59)])

sp_ba <- sp_ba %>%
  filter(ID_PE %in% pep_pe$ID_PE) %>%
  left_join(pep_vegpot, by = c("ID_PE", "ID_PE_MES"))


sort(colSums(sp_ba[,MySpecies]))



# Distance ####

# Hellinger distance matrice
dist_hell <- dist(decostand(sp_ba[,MySpecies], "hel"))

# log-chord distances

dist_chord <- (decostand(log1p(sp_ba[,MySpecies]), "normalize"))
dist_chord <- dist(dist_chord)


# hclust ####

#### #### #### #### #### #### #### ####
# Ward clustering on log chord ####
#### #### #### #### #### #### #### ####

grWard_chord <- hclust(dist_chord, method = "ward.D2")


ss <- names(which(colSums(sp_ba[,MySpecies])>2000))
seqy = seq(0, 1, len = length(ss))

# Dendrogramme ####
xgr = grWard_chord
gr = cutree(xgr, k = 18)
gr_ord <- xgr$order

png("res/dendro_ward_chord.png", width = 10, height = 6, units = "in", res = 300)
layout(matrix(1:2), heights = c(.6,.4))
par(mar = c(0,3,1,1))
plot(xgr, hang = -1,xaxs = "i", labels = FALSE, main = "Ward clustering on log-Chord distance - 18 groups", 
     xlab="", ylab = "", sub="", cex.axis = .6,
     las = 1)
rect.hclust(xgr, k = 18, border = 2:20) 

par(mar = c(1,3,0,1))
image(as.matrix(sp_ba[gr_ord,rev(ss)]), axes = F, 
      col = hcl.colors(12, "Viridis"))
text(-.02, rev(seqy), ss, xpd = NA, 
     cex = .5, adj = 1, font = 3)
dev.off()

# Boxplot de la composition par groupe ####

png('res/ward_chord_sp_boxplot.png', width = 13, height = 8.5, 
    res = 300, units = "in")
par(mfrow = c(4,5), mar = c(1.5,3.4,1,.5))

for(i in 1:18){
  boxplot(sp_ba[gr==i, ss], 
          horizontal = T, cex.axis = .65, las = 1, xaxt = "n", 
          outwex = 1, pch = 20, cex = .1, outcol = "grey", col = "red3") 
  axis(1, cex.axis = .5, tick = F, line = -1)
  mtext(paste('Groupe', i), 3, line = 0, cex = .6)
}

dev.off()

# Distribution des groupes par vegpot ####
sp_ba_gr <- cbind(sp_ba, gr = gr)

gr_vegpot = table(sp_ba_gr$VEG_POT, sp_ba_gr$gr)

write.csv(gr_vegpot, "vegpot_gr_ward_chord.csv")

#### #### #### #### #### #### #### ####
# Ward clustering on hellinger ####
#### #### #### #### #### #### #### ####

grWard_hell <- hclust(dist_hell, method = "ward.D2")

ss <- names(which(colSums(sp_ba[,MySpecies])>2000))
seqy = seq(0, 1, len = length(ss))

# Dendrogramme ####
xgr = grWard_hell
gr = cutree(xgr, k = 16)
gr_ord <- xgr$order

png("res/dendro_ward_hell.png", width = 10, height = 6, units = "in", res = 300)
layout(matrix(1:2), heights = c(.6,.4))
par(mar = c(0,3,1,1))
plot(xgr, hang = -1,xaxs = "i", labels = FALSE, main="Ward clustering on Hellinger distance - 16 groups", 
     xlab="", ylab = "", sub="", cex.axis = .6,
     las = 1)
rect.hclust(xgr, k = 16, border = 2:18) 

par(mar = c(1,3,0,1))
image(as.matrix(sp_ba[gr_ord,rev(ss)]), axes = F, 
      col = hcl.colors(12, "Viridis"))
text(-.02, rev(seqy), ss, xpd = NA, 
     cex = .5, adj = 1, font = 3)
dev.off()

# Boxplot de la composition par groupe ####

png('res/ward_hell_sp_boxplot.png', width = 13, height = 8.5, 
    res = 300, units = "in")
par(mfrow = c(4,5), mar = c(1.5,3.4,1,.5))

for(i in 1:18){
  boxplot(sp_ba[gr==i, ss], 
          horizontal = T, cex.axis = .65, las = 1, xaxt = "n", 
          outwex = 1, pch = 20, cex = .1, outcol = "grey", col = "red3") 
  axis(1, cex.axis = .5, tick = F, line = -1)
  mtext(paste('Groupe', i), 3, line = 0, cex = .6)
}

dev.off()

# Distribution des groupes par vegpot ####
sp_ba_gr <- cbind(sp_ba, gr = gr)

gr_vegpot = table(sp_ba_gr$VEG_POT, sp_ba_gr$gr)

write.csv(gr_vegpot, "vegpot_gr_ward_hell.csv")



