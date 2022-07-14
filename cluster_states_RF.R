#///////////////////////////////////////////////////////////
#   Modification du code de MH Brice cluster_states.R
#  Loop to produce clustering (2 methods) and plots for a range of clusters `[minClus, MaxClus]`
# RandomForest model produced for each number of clusters, ntrees= nArb
# No parallel so it can be very slow (cluster distances + random forest loop)
#///////////////////////////////////////////////////////////

rm(list = ls())
gc()

set.seed(111)  #  repeat

# packages
library(sf)
library(dplyr)
library(vegan)
library(fastcluster)
library(distances)
require(randomForest)

# Paramètres pour boucle et RF
nArb <- 750
minClus <- 13
MaxClus <- 21

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


ss <- names(which(colSums(sp_ba[,MySpecies]) > 2000))
seqy = seq(0, 1, len = length(ss))



#### #### #### #### #### #### #### ####
# Loop Ward clustering on log chord ####
#### #### #### #### #### #### #### ####

for (i in MaxClus:minClus){         #    i <- 16
  
xgr = grWard_chord
gr = cutree(xgr, k = i)
gr_ord <- xgr$order

# Titles and filenames depending on loop index
name1 <- paste0("resLoop/dendro_ward_chord", i,".png")
name2 <- paste0('resLoop/ward_chord_sp_boxplot', i,".png")
name3 <- paste0('resLoop/vegpot_gr_ward_chord', i,".csv")
name4 <- paste0("resLoop/rf/rfTestWardChord", i,".RData")
main1 <- paste0("resLoop/Ward clustering on log-Chord distance - ", i," groups")


# Dendrogramme ####
png(name1, width = 10, height = 6, units = "in", res = 300)
layout(matrix(1:2), heights = c(.6,.4))
par(mar = c(0,3,1,1))
plot(xgr, hang = -1,xaxs = "i", labels = FALSE, main = main1, 
     xlab="", ylab = "", sub="", cex.axis = .6,
     las = 1)
rect.hclust(xgr, k = i, border = 2:20) 

par(mar = c(1,3,0,1))
image(as.matrix(sp_ba[gr_ord,rev(ss)]), axes = F, 
      col = hcl.colors(12, "Viridis"))
text(-.02, rev(seqy), ss, xpd = NA, 
     cex = .5, adj = 1, font = 3)
dev.off()

# Boxplot de la composition par groupe ####

png(name2, width = 13, height = 8.5, 
    res = 300, units = "in")
par(mfrow = c(4,5), mar = c(1.5,3.4,1,.5))
if(i==21){
  par(mfrow = c(5,5), mar = c(1.5,3.4,1,.5))
}

for(j in 1:i){
  boxplot(sp_ba[gr==j, ss], 
          horizontal = T, cex.axis = .65, las = 1, xaxt = "n", 
          outwex = 1, pch = 20, cex = .1, outcol = "grey", col = "red3") 
  axis(1, cex.axis = .5, tick = F, line = -1)
  mtext(paste('Groupe', j), 3, line = 0, cex = .6)
}

dev.off()

# Distribution des groupes par vegpot ####
sp_ba_gr <- cbind(sp_ba, gr = gr)

gr_vegpot = table(sp_ba_gr$VEG_POT, sp_ba_gr$gr)

write.csv(gr_vegpot, name3)

# Prepare data for RandomForest


# RF files

 arbres <- nrow(sp_ba_gr)
# ntr <- 1500
# sz <- 5000
#

# sample variable
# train <- sample(1:nrow(sp_ba_gr), arbres)
# trainY <- as.factor(sp_ba_gr$gr[train])
# trainX <-  sp_ba_gr[train,ss]
# test <- sp_ba_gr[-train, c("gr", ss)]
# rm(train)

trainY <- as.factor(sp_ba_gr$gr)
trainX <-  sp_ba_gr[, ss]



# RandoForest 
rfTestPar <- randomForest(trainX, trainY, ntree = nArb, do.trace = 50) #,sampsize=sz)

# rfTestPar
# varImpPlot(rfTestPar)
# plot(rfTestPar)

#Keep results par cluster number
dfI <- data.frame(rfTestPar$importance)
dfI$spe <- row.names(dfI)
dfIm <- dfI[order(-dfI$MeanDecreaseGini),]
dfIm$value <- paste0(dfIm$spe, " (", round(dfIm$MeanDecreaseGini, digits =0), ")")
dfIm$spe <- NULL
dfIm$MeanDecreaseGini <- NULL
dfIm <- as.data.frame(t(dfIm))
names(dfIm) <- paste0("varImportance_", c(1:22))
row.names(dfIm) <- paste0("Imp_MeanGini_",i)



dfCat <- data.frame(table(trainY))
dfCat$err <- rfTestPar$confusion[,ncol(rfTestPar$confusion)]
dfCat$value2 <- paste0(dfCat$Freq, " (", round(dfCat$err * 100, digits = 2), ")")
dfCat$Freq <- NULL
dfCat$err <- NULL
names(dfCat)[2] <- paste0("Freq_Err", i)

 sort(colSums(sp_ba[,MySpecies]))

save(rfTestPar, file = name4)

conf <- rfTestPar$confusion[,-ncol(rfTestPar$confusion)]
oob <- round((1 - (sum(diag(conf))/sum(conf))) * 100,digit =2)
maxErr <- round(100 * max((rfTestPar$confusion[,ncol(rfTestPar$confusion)])), digit = 2)
# dfErr <- data.frame(rfTestPar$confusion[,ncol(rfTestPar$confusion)])
# names(dfErr) <- paste0("ErrClass_", i)
# dfErr$clus <- c(1:i)
# dfErr <- dfErr[, c(2,1)]


# Assemble info for each loop value
if(i == MaxClus){
  err <- oob
  MAXerr<- maxErr
  dfImp <- dfIm
  dfCate <- dfCat
  # dfErro <- dfErr
} else{
  err <- c(err, oob)
  MAXerr <- c(MAXerr,maxErr)
  dfImp <- rbind(dfImp, dfIm)
  dfCate <- left_join(dfCate, dfCat)
  # dfErro <- left_join(dfErro,dfErr)
}



} # End loop ward log chord ####

#Prepare Final DF par clustering method
clust <- c(MaxClus:minClus)

df <- data.frame(clust, err, MAXerr)

dfCate$trainY <- NULL
dfCateT <-  as.data.frame(t(dfCate))
names(dfCateT) <- paste0("Freq&ERR(%)_", c(1:MaxClus))
df1 <- cbind(df, dfCateT)
df3 <- cbind(df1, dfImp)
row.names(df3) <- c(MaxClus:minClus)
name5 <- paste0("resLoop/rf/errRrfTestWardChord",".csv")
write.csv(df3, name5, row.names = FALSE )


rm(df,df1 , dfCateT, dfCate, dfCat,dfI, dfImp)

rm(rfTestPar, xgr, grWard_chord, trainX, trainY, dist_chord)


#### #### #### #### #### #### #### ####
# Ward clustering on hellinger ####
#### #### #### #### #### #### #### ####

grWard_hell <- hclust(dist_hell, method = "ward.D2")

ss <- names(which(colSums(sp_ba[,MySpecies])>2000))
seqy = seq(0, 1, len = length(ss))



for (i in MaxClus:minClus){     #  i <- 11

xgr = grWard_hell
gr = cutree(xgr, k = i)
gr_ord <- xgr$order

name1 <- paste0("resLoop/dendro_Hell", i,".png")
name2 <- paste0('resLoop/Hell_sp_boxplot', i,".png")
name3 <- paste0('resLoop/vegpot_gr_Hell', i,".csv")
name4 <- paste0("resLoop/rf/rfTestHell", i,".RData")
main1 <- paste0("resLoop/Ward clustering on Hell distance - ", i," groups")


# Dendrogramme ####
png(name1, width = 10, height = 6, units = "in", res = 300)
layout(matrix(1:2), heights = c(.6,.4))
par(mar = c(0,3,1,1))
plot(xgr, hang = -1,xaxs = "i", labels = FALSE, main= main1, 
     xlab="", ylab = "", sub="", cex.axis = .6,
     las = 1)
rect.hclust(xgr, k = i, border = 2:18) 

par(mar = c(1,3,0,1))
image(as.matrix(sp_ba[gr_ord,rev(ss)]), axes = F, 
      col = hcl.colors(12, "Viridis"))
text(-.02, rev(seqy), ss, xpd = NA, 
     cex = .5, adj = 1, font = 3)
dev.off()

# Boxplot de la composition par groupe ####

png(name2, width = 13, height = 8.5, 
    res = 300, units = "in")
par(mfrow = c(4,5), mar = c(1.5,3.4,1,.5))
if(i==21){
  par(mfrow = c(5,5), mar = c(1.5,3.4,1,.5))
}

for(j in 1:i){
  boxplot(sp_ba[gr== j, ss], 
          horizontal = T, cex.axis = .65, las = 1, xaxt = "n", 
          outwex = 1, pch = 20, cex = .1, outcol = "grey", col = "red3") 
  axis(1, cex.axis = .5, tick = F, line = -1)
  mtext(paste('Groupe', j), 3, line = 0, cex = .6)
}

dev.off()

# Distribution des groupes par vegpot ####
sp_ba_gr <- cbind(sp_ba, gr = gr)

gr_vegpot = table(sp_ba_gr$VEG_POT, sp_ba_gr$gr)

write.csv(gr_vegpot, name3)

# Prepare RandmForest
arbres <- nrow(sp_ba_gr)
# ntr <- 1500
# sz <- 5000
#

# sample variable
# train <- sample(1:nrow(sp_ba_gr), arbres)
# trainY <- as.factor(sp_ba_gr$gr[train])
# trainX <-  sp_ba_gr[train,ss]
# test <- sp_ba_gr[-train, c("gr", ss)]
# rm(train)

trainY <- as.factor(sp_ba_gr$gr)
trainX <-  sp_ba_gr[, ss]

#RandomForests
rfTestPar <- randomForest(trainX, trainY, ntree = nArb, do.trace =50) #,sampsize=sz)

# rfTestPar
# varImpPlot(rfTestPar)
# plot(rfTestPar)

#Keep info per cluster number
dfI <- data.frame(rfTestPar$importance)
dfI$spe <- row.names(dfI)
dfIm <- dfI[order(-dfI$MeanDecreaseGini),]
dfIm$value <- paste0(dfIm$spe, " (", round(dfIm$MeanDecreaseGini, digits =0), ")")
dfIm$spe <- NULL
dfIm$MeanDecreaseGini <- NULL
dfIm <- as.data.frame(t(dfIm))
names(dfIm) <- paste0("varImportance_", c(1:22))
row.names(dfIm) <- paste0("Imp_MeanGini_",i)



dfCat <- data.frame(table(trainY))
dfCat$err <- rfTestPar$confusion[,ncol(rfTestPar$confusion)]
dfCat$value2 <- paste0(dfCat$Freq, " (", round(dfCat$err * 100, digits = 2), ")")
dfCat$Freq <- NULL
dfCat$err <- NULL
names(dfCat)[2] <- paste0("Freq_Err", i)
                       
sort(colSums(sp_ba[,MySpecies]))

save(rfTestPar, file = name4)

conf <- rfTestPar$confusion[,-ncol(rfTestPar$confusion)]
oob <- round((1 - (sum(diag(conf))/sum(conf))) * 100,digit =2)
maxErr <- round(100 * max((rfTestPar$confusion[,ncol(rfTestPar$confusion)])), digit = 2)
# dfErr <- data.frame(rfTestPar$confusion[,ncol(rfTestPar$confusion)])
# names(dfErr) <- paste0("ErrClass_", i)
# dfErr$clus <- c(1:i)
# dfErr <- dfErr[, c(2,1)]


# Combine info for each loop value
if(i == MaxClus){
  err <- oob
  MAXerr<- maxErr
  dfImp <- dfIm
  dfCate <- dfCat
  # dfErro <- dfErr
} else{
  err <- c(err, oob)
  MAXerr <- c(MAXerr,maxErr)
  dfImp <- rbind(dfImp, dfIm)
  dfCate <- left_join(dfCate, dfCat)
  # dfErro <- left_join(dfErro,dfErr)
}



}    # End loop Hellinger ####

#Prepare Final DF par clustering method
clust <- c(MaxClus:minClus)

df <- data.frame(clust, err, MAXerr)

dfCate$trainY <- NULL
dfCateT <-  as.data.frame(t(dfCate))
names(dfCateT) <- paste0("Freq&ERR(%)_", c(1:MaxClus))
df1 <- cbind(df, dfCateT)
df2 <- cbind(df1, dfImp)
row.names(df2) <- c(MaxClus:minClus)
name5 <- paste0("resLoop/rf/errRrfTestWardHell",".csv")
write.csv(df2, name5, row.names = FALSE )

rm(df,df1 , dfCateT, dfCate, dfCat,dfI, dfImp)

rm(rfTestPar, xgr, grWard_chord, trainX, trainY, dist_chord)

