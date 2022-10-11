#****************************
##  Produce Groupement essences
#   Selon schemes Norme de stratification écoforestière
# https://mffp.gouv.qc.ca/documents/forets/inventaire/norme-stratification.pdf
# Pages 11-24
# Inputs are stacks from landis with  %biomass per species [0-1].
# On en a besoin d'info de la carte ecoforestière - classe hydrique)
# Code written for 17 species simulated in Landis (see objet speLa)
# Modifying species will demand modify code as it considers only our species
# colones number should be switched by colones names
# Results are csv files with grouping integer codes par pixel (integer code), scenario and year
#***********************************
#

rm(list=ls())
gc()

#***********************************

require(raster)
require(utils)
require(dplyr)
require(tidyr)



#-----------------------
# RScriptDir <- "C:\\Users\\jpascual\\Documents\\R\\R-4.1.1\\bin\\x64\\Rscript.exe"
# 
# rasterOptions(tmpdir = "U:/YBoulanger/Naturalité/Maps/anci/stacks/TEMP/")
# write("R_USER = U:/YBoulanger/Naturalité/Maps/anci/stacks/TEMP/", file=file.path(Sys.getenv('R_USER'), '.Renviron'))



outDir <- paste0("U:/YBoulanger/Naturalité/Maps/anci/stacks/TEMP/",Sys.Date())
dir.create(outDir)
setwd(outDir)

# Dir for intermediate and final results
dir.create("eta1")
dir.create("eta2")
dir.create("Final")





#************************************************
#------------------------------------------------------------
#     Read raster veg_pot and csv integer codes
#------------------------------------------------------------
#************************************************

Dir <- paste0("U:\\YBoulanger\\Naturalité\\Maps\\anci")



# Class Hydrique from carte ecoforestière. Needed for groupements

rasClaHid <- raster("P:/A_2022/landisToGrESS/ras/ClaHydric.tif")
#  0  xerique et mesique
#  1 sub et hydrique

dfClaHid <- as.data.frame(values(rasClaHid))
names(dfClaHid) <- "claHyd"
rm (rasClaHid)



#************************************************
#------------------------------------------------------------
#     Read raster ID pixels
#------------------------------------------------------------
#************************************************

codRas <- raster("P:/A_2022/landisToGrESS/ras/rasId.tif") #Raster codInteger pixel
dfId <- as.data.frame(values(codRas))
names(dfId) <- "codRas"
rm(codRas)


#||||||||||||||||||||||||||||||||||||

# SCENARIOS LANDIS avec noms
#HARV
scenHar1 <- c("BudwormBaselineFire", "BudwormBaselineFireBaselineHarvest")
scenHar2 <- c("GrowthBudwormProjectedFire", "GrowthBudwormProjectedFireBaselineHarvest")

scenHarv <- c("NoHar", "YeHar")

#CLIMA
scenCli <- c("baseline", "RCP45", "RCP85")
scenClim <- c("basel", "RCP45", "RCP85")



speLa <- c("ABIE.BAL", "ACER.RUB", "ACER.SAH", "BETU.ALL", "BETU.PAP", "FAGU.GRA",
           "LARI.LAR", "PICE.GLA", "PICE.MAR", "PICE.RUB", "PINU.BAN", "PINU.RES",
           "PINU.STR", "POPU.TRE", "QUER.RUB", "THUJ.SPP.ALL", "TSUG.CAN")

speQc <- c("SAB", "ERR", "ERS", "BOJ", "BOP", "HEG",
           "MEL", "EPB", "EPN", "EPR", "PIG", "PIR",
           "PIB", "PET", "CHR", "THO", "PRU")

# To name speces groups  2 character MFFP codes
speQc2 <- c("Sb", "Eo", "Es", "Bj", "Bp", "Hg",
           "Ml", "Eb", "En", "Eu", "Pg", "Pr",
           "Pb", "Pt", "Cr", "To", "Pu")

years <- seq(0,130,10)



# Code integer for all groupements; faster and less space needed to save integer than 5-8 characters code 
# 
# vecCods <- c("Sb","Eo", "Es", "Bj", "Bp", "Hg", "Ml", "Eb", "En", "Eu", "Pg", "Pr",  "Pb",  "Pt", "Cr", "To", "Pu",      
#              "Er", "Ep",  "Pi","Se", "Ft", "Fi","Fh", "Rx", "Fx")

# all softwood codes in our case
vecCodRe <- c("En","Eb", "Sb","Pg",  "Ml", "Eu", "Pr", "Pu", "Pb","To",       
              "Ep",  "Pi","Se",  "Rx")
# all hardwood codes in our case
vecCodFe <- c("Pt","Bp","Eo","Bj", "Es", "Hg",    "Cr",      
             "Er",  "Ft", "Fi","Fh", "Fx")

# All groups dataframe with groupemetn integer code
Rcods <- expand.grid(uno = vecCodRe, dos = vecCodRe)
Rcods$G <- paste0(Rcods$uno, " ", Rcods$dos)
Fcods <- expand.grid(uno = vecCodFe, dos = vecCodFe)
Fcods$G <- paste0(Fcods$uno, " ", Fcods$dos)

RCod <- Rcods$G      # softwoods
FCod <- Fcods$G     # Hardwoods

MRCods <- expand.grid(uno = RCod, dos = vecCodFe)
MRCods$G <- paste0(MRCods$uno, " ", MRCods$dos)
MRCod <- MRCods$G # mixte mainly softwood

MFCods <- expand.grid(uno = FCod, dos = vecCodRe)
MFCods$G <- paste0(MFCods$uno, " ", MFCods$dos)
MFCod <- MFCods$G  # mixte mainly hardwoods

Groupem <- c(RCod, MRCod, MFCod, FCod)
codGroup <- seq(1:length(Groupem))
linkGroup <- data.frame(Groupem, codGroup)    # 4708 posiible groupings in our case (17 species)
write.csv(linkGroup, "P:/A_2022/landisToGrESS/data/intGroup_Groupement.csv", row.names =  FALSE )

rm(vecCodRe, vecCodFe,Rcods, Fcods, RCod, FCod,MRCods,MRCod ,MFCods, MFCod, Groupem)


codYear <- letters[seq(from = 1, to = 14)] # code for years 0-130
#-----------------------------------------------
#-------------------- Lecture stack percentage + maxAge
#-----------------------------------------------




#/-/-/-/-/-/-/-/-/-/-
#Prepare parallel

require(doSNOW)
require(parallel)
require(foreach)


clusterN <- 18

# cl = makeCluster(clusterN, rscript=RScriptDir, type='SOCK')
cl = makeCluster(clusterN,  type='SOCK')
registerDoSNOW(cl)
#/-/-/-/-/-/-/-/-/-/-




#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------    PARALLEL    --------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------






 foreach(i = seq_along(scenCli)) %:%
 
  
  foreach(j = seq_along(scenHar1))%:%
 
  
  foreach(k = seq_along(years)) %:%
   
  
  foreach(rep = 1:5) %dopar% {     #     i <- 1 ;   j <- 1; k <- 3; rep <- 1
    # gc()
    # for (rep in 1:5){
    
    require(raster)
    require(reshape2)
    require(data.table)
    require(plyr)
    require(dplyr)
    
    
    
    if(scenCli[i] == "baseline"){scenHar <- scenHar1} else{scenHar <- scenHar2}
    # read merged stack
    staPerc <- stack(paste0(Dir,"/stacks/merged","\\","stack_percBiom_","_",rep, scenCli[i],"_",scenHar[j],"_", years[k],".tif"))
   
    names(staPerc)[1:17] <- speQc2  # code mffp 2 character
    # Delete bioTot et maxAge
     staPerc <- staPerc[[1:17]]
     # summary(values(staPerc)[[1:2]]) 
     
     
    df21 <- data.frame(values(staPerc))     # raster values into df
    
    rm(staPerc)
    df21$IDras <- dfId$codRas        # Create ID par pixel
    df21 <- df21[, c(18, 1:17)]
    
    #Addition de code Hydrique (0 et 1 (sub et hydryque) pour bien classer FT ou FH-ERR et BOJ)
    df21$claHyd <- dfClaHid$claHyd
    df21[is.na(df21$claHyd), "claHyd"] <- 0   #cheap mais on ne perd pas information ! Si pas d'information on considere non humide (0)
    # summary(df21)
    
    # creer colonnes pour Bj, Eo,  en terrain hydriques et sub et en xerique-mesiques
    
    
    df21$BjXM <- 0
    df21$BjSH <- 0
    
    df21[df21$claHyd == 0, "BjXM"] <- df21[df21$claHyd == 0, "Bj"]
    df21[df21$claHyd == 1, "BjSH"] <- df21[df21$claHyd == 1, "Bj"]
    
    df21$EoXM <- 0
    df21$EoSH <- 0
    
    df21[df21$claHyd == 0, "EoXM"] <- df21[df21$claHyd == 0, "Eo"]
    df21[df21$claHyd == 1, "EoSH"] <- df21[df21$claHyd == 1, "Eo"]
    
    
    
    
    # Étape 0
    
    #//////////////////////////////////////////////////
    # R' F ou M
   #///////////////////////////////////////////////////
    # Somme RE
    df21$RE <- df21[,"Sb"] + df21[,"Ml"]+ df21[,"Eb"] + df21[,"En"] + df21[,"Eu"] + df21[,"Pg"] + 
      df21[,"Pr"] + df21[,"Pb"] + df21[,"To"] + df21[,"Pu"]
    # Somme FE
    df21$FE <- df21[,"Eo"] + df21[,"Es"]+ df21[,"Bj"] + df21[,"Bp"] + df21[,"Hg"] + df21[,"Pt"] + 
      df21[,"Cr"]
    
    # if data from landis RE !is.na
    df21 <- df21[!is.na(df21$RE) | !is.na(df21$FE), ]
    
    # df21$tipCouv <- NULL
    # 
    # df21$tipCouv <- "M"
    df21[which(df21$RE > 0.75), "tipCouv"] <- "R"
    df21[which(df21$FE > 0.749), "tipCouv"] <- "F"
    df21[which(df21$RE >= 0.50 &  df21$RE <= 0.75), "tipCouv"] <- "MR"
    df21[which(df21$RE > 0.25 &  df21$RE < 0.50), "tipCouv"] <- "MF"
    # table(df21$tipCouv)
    
    #Registre de code RASter et typCouverture
    IdTipCouv <- df21[,c(1,26)]
    
    
  
    
    
    # creer colonnes avec % par combinaison essences 
    
    df21$Er <- df21$Es + df21$Eo
    
    df21$Ep <- df21$En + df21$Eu
    
    df21$Pi <- df21$Pg + df21$Pr + df21$Pb
    
    # creer colonnes avec % par association  d'essences 
    
    
    df21$Se <- df21$Eb + df21$Sb
    
    df21$Ft <- df21$Cr + df21$Es + df21$EoXM + df21$BjXM + df21$Hg
    
    df21$Fi <- df21$Pt + df21$Bp
    
    df21$Fh <-  df21$BjSH + df21$EoSH    #  à reflechir pour caracteriser hydriques et separer PEB de PET
    
    # On ventile le DF en RE (schéma 4), FE (schéma 5) er ME (Schéma 6)
    
    #    R   Schéma 4  (R) et partie du 6 (MR)
    
    df21R <- df21[df21$tipCouv == "R" | df21$tipCouv == "MR",]  #   we use it now  (schéma 4 and 5)
    df21 <- df21[!(df21$tipCouv == "R"),]  # we keep for schéma 5 and part of 6) (risqué??)
    
    
    
    # Étape 1 when sp > 0.75
    df21R$Etape1<- apply(df21R[, c(2, 8:14,17:18)] > (0.75 * df21R$RE),
                         1, function(x) ifelse(any(x), paste(colnames(df21R)[c(2, 8:14,17:18)][x], collapse=', '), 'None'))
    
    df21R_E1 <- df21R[!(df21R$Etape1 =="None"), c(1,34)]
    df21R_E1$Groupem <- paste0(df21R_E1$Etape1, " ", df21R_E1$Etape1)
    
    #---------------------
    df21R_E1 <-  df21R_E1[,c(1,3)]
    df21R_E1 <- left_join(df21R_E1, linkGroup)
    df21R_E1 <- df21R_E1[, c(1,3)]
    #--------------------
    
     fileNamee <- paste0 ("eta1/R_", scenClim[i], "_", scenHarv[j], "_", codYear[k], "_",  rep, "_E1.csv")
     fwrite(df21R_E1, file = fileNamee)
     rm(df21R_E1)
    
    df21R <- df21R[df21R$Etape1 =="None",]
    
    
    # Étape 2 when sp > 0.5 and sp > 0.25
    df21R$Etape2_1 <- apply(df21R[, c(2, 8:14,17:18)] > (0.50 * df21R$RE) & df21R[, c(2, 8:14,17:18)] <= (0.75 * df21R$RE),
                            1, function(x) ifelse(any(x), paste(colnames(df21R)[c(2, 8:14,17:18)][x], collapse=', '), 'None'))
    df21R$Etape2_2 <- apply(df21R[, c(2, 8:14,17:18)] > (0.25 * df21R$RE) & df21R[, c(2, 8:14,17:18)] <= (0.50 * df21R$RE),
                            1, function(x) ifelse(any(x), paste(colnames(df21R)[c(2, 8:14,17:18)][x], collapse=', '), 'None'))
    
    df21R_E2_E <- df21R[!(df21R$Etape2_1 == "None") & !(df21R$Etape2_2 =="None"), c(1,35:36)]
    df21R_E2_E$Groupem <- paste0(df21R_E2_E$Etape2_1, " ", df21R_E2_E$Etape2_2)
    
    #---------------------
    df21R_E2_E <-  df21R_E2_E[,c(1,4)]
    df21R_E2_E <- left_join(df21R_E2_E, linkGroup)
    df21R_E2_E <- df21R_E2_E[, c(1,3)]
    #--------------------
    
    
    df21R  <- df21R[(df21R$Etape2_1 == "None" | df21R$Etape2_2 =="None"),]
    
    # Étape 2 (suite) when sp > 0.5 and groupe essences > 0.25
    # df21R$Etape2_2G <- apply(df21R[, c(28:30)] > (0.25 * df21R$RE) & df21R[, c(28:30)] <= (0.50 * df21R$RE), 1, function(x) ifelse(any(x), paste(colnames(df21R)[c(28:30)][x], collapse=', '), 'None'))
    
   
    df21R$Etape2_2Ge <-  colnames(df21R[c(28:30)])[apply(df21R[c(28:30)],1,which.max)]
    df21R$Etape2_2GeV <- apply(X=df21R[c(28:30)], MARGIN=1, FUN=max)
    df21R$Etape2_2G <- "None"
    df21R[(df21R$Etape2_2GeV  > (0.25 * df21R$RE)) & (df21R$Etape2_2GeV  <= (0.50 * df21R$RE)), "Etape2_2G"  ] <- 
      df21R[(df21R$Etape2_2GeV  > (0.25 * df21R$RE)) & (df21R$Etape2_2GeV  <= (0.50 * df21R$RE)), "Etape2_2Ge"  ]
    
    df21R <- df21R[, c(1:36,39)]
    
    df21R[!(df21R$Etape2_1 == "None") & (df21R$Etape2_2G == "None"), "Etape2_2G"] <- "Rx"
    
    
    df21R_E2_G <- df21R[!(df21R$Etape2_1 == "None") & !(df21R$Etape2_2G =="None"),c(1,35,37)]
    #
    df21R_E2_G$Groupem <- paste0(df21R_E2_G$Etape2_1, " ", df21R_E2_G$Etape2_2G)
    #
    #---------------------
    df21R_E2_G <-  df21R_E2_G[,c(1,4)]
    df21R_E2_G <- left_join(df21R_E2_G, linkGroup)
    df21R_E2_G <- df21R_E2_G[, c(1,3)]
    #--------------------
    
    df21R_E2 <- rbind(df21R_E2_E, df21R_E2_G)
    
    fileNamee <- paste0 ("eta1/R_", scenClim[i], "_", scenHarv[j], "_", codYear[k], "_",  rep, "_E2.csv")
    fwrite(df21R_E2, file = fileNamee)
    rm(df21R_E2, df21R_E2_E, df21R_E2_G)
    
   
    
    
    df21R <- df21R[(df21R$Etape2_1 == "None" | df21R$Etape2_2G =="None"),]
    
    # Étape 3 groupe d'essences > 0.75  and one essence > 25 and > other in the group
    
    df21R$Etape3_1<- apply(df21R[, c(28:30)] > (0.75 * df21R$RE), 1, function(x) ifelse(any(x), paste(colnames(df21R)[c(28:30)][x], collapse=', '), 'None'))
    
    #//////////// Ep
    d21R_Ep <- df21R[df21R$Etape3_1 == "Ep",]
    d21R_Ep$Etape3_2<- apply(d21R_Ep[, c(10:11)] > (0.5 * d21R_Ep$Ep), 1, function(x) ifelse(any(x), paste(colnames(d21R_Ep)[c(10:11)][x], collapse=', '), 'None'))
    d21R_Ep <- d21R_Ep[,c(1, 38:39)]
    d21R_Ep[d21R_Ep$Etape3_2 == "None", "Etape3_2"] <- "Ep"
    d21R_Ep$Groupem <- paste0(d21R_Ep$Etape3_1, " ", d21R_Ep$Etape3_2)
    
    #//////////// Se
    d21R_Se <- df21R[df21R$Etape3_1 == "Se",]
    d21R_Se$Etape3_2<- apply(d21R_Se[, c(2,9)] > (0.5 * d21R_Se$Se), 1, function(x) ifelse(any(x), paste(colnames(d21R_Se)[c(2,9)][x], collapse=', '), 'None'))
    d21R_Se <- d21R_Se[,c(1, 38:39)]
    d21R_Se[d21R_Se$Etape3_2 == "None", "Etape3_2"] <- "Se"
    d21R_Se$Groupem <- paste0(d21R_Se$Etape3_1, " ", d21R_Se$Etape3_2)
    
    #//////////// Pi
    d21R_Pi <- df21R[df21R$Etape3_1 == "Pi",]
    d21R_Pi$Etape3_2e <- colnames(d21R_Pi[c(12:14)])[apply(d21R_Pi[c(12:14)],1,which.max)]
    d21R_Pi$Etape3_2eV <- apply(X = d21R_Pi[c(12:14)], MARGIN=1, FUN=max)
    d21R_Pi$Etape3_2 <- "None"
    d21R_Pi[d21R_Pi$Etape3_2eV > (0.25 * d21R_Pi$RE), "Etape3_2"] <-  d21R_Pi[d21R_Pi$Etape3_2eV > (0.25 * d21R_Pi$RE), "Etape3_2e"]
  
    
    d21R_Pi <- d21R_Pi[,c(1, 38,41)]
    d21R_Pi[d21R_Pi$Etape3_2 == "None", "Etape3_2"] <- "Pi"
    d21R_Pi$Groupem <- paste0(d21R_Pi$Etape3_1, " ", d21R_Pi$Etape3_2a) 
    
    #-------------------------
    df21R_E3<- rbind(d21R_Ep[, c(1,4)], d21R_Se[, c(1,4)], d21R_Pi[, c(1,4)])
    df21R_E3 <- left_join(df21R_E3, linkGroup)
    df21R_E3 <- df21R_E3[, c(1,3)]
    #--------------------------
    fileNamee <- paste0 ("eta1/R_", scenClim[i], "_", scenHarv[j], "_", codYear[k], "_",  rep,  "_E3.csv")
    fwrite(df21R_E3, file = fileNamee)
    rm(df21R_E3)
    
     rm(d21R_Pi, d21R_Se, d21R_Ep)                     
  
     df21R <- df21R[!(df21R$Etape3_1 %in% c("Ep","Pi", "Se")),]
     
     
     
     # Étape 4 RE groupe d'essences > 0.5  and one essence autre groupe > 25 and > other in the group
    
     df21R$Etape4_1<- apply(df21R[, c(28:30)] > (0.5 * df21R$RE), 1, function(x) ifelse(any(x), paste(colnames(df21R)[c(28:30)][x], collapse=', '), 'None'))
     
     #//////////// EP
     d21R_Ep <- df21R[df21R$Etape4_1 == "Ep",]
     d21R_Ep$Etape4_2<- apply(d21R_Ep[, c(2,8:9,12:14,17:18)] > (0.25 * d21R_Ep$RE), 1, function(x) ifelse(any(x), paste(colnames(d21R_Ep)[c(2,8:9, 12:14,17:18)][x], collapse=', '), 'None'))
     #  c(2, 8:14,17:18)
     d21R_Ep$Etape4_3<- apply(d21R_Ep[, c(29,30)] > (0.25 * d21R_Ep$RE), 1, function(x) ifelse(any(x), paste(colnames(d21R_Ep)[c(29:30)][x], collapse=', '), 'None'))
     d21R_Ep <- d21R_Ep[,c(1, 39:41)]
     
     d21R_Ep[!(d21R_Ep$Etape4_2 == "None"),"Groupem"] <- paste0(d21R_Ep[!(d21R_Ep$Etape4_2 == "None"),"Etape4_1"],
                                                                                                  " ", d21R_Ep[ !(d21R_Ep$Etape4_2 == "None"),"Etape4_2"] )
     d21R_Ep[(d21R_Ep$Etape4_2 == "None") & !(d21R_Ep$Etape4_3 == "None") ,"Groupem"] <- paste0(d21R_Ep[(d21R_Ep$Etape4_2 == "None") & !(d21R_Ep$Etape4_3 == "None") ,"Etape4_1"],
                                                                  " ", d21R_Ep[(d21R_Ep$Etape4_2 == "None") & !(d21R_Ep$Etape4_3 == "None") ,"Etape4_3"] )
     d21R_Ep[is.na(d21R_Ep$Groupem), "Groupem"] <-  paste0(d21R_Ep[is.na(d21R_Ep$Groupem),"Etape4_1"],     " ", "Rx")
    
     
     #//////////// Se
     d21R_Se <- df21R[df21R$Etape4_1 == "Se",]
     d21R_Se$Etape4_2<- apply(d21R_Se[, c(8, 10:14,17:18)] > (0.25 * d21R_Se$RE), 1, function(x) ifelse(any(x), paste(colnames(d21R_Se)[c(8, 10:14,17:18)][x], collapse=', '), 'None'))
     #  c(2, 8:14,17:18)
     d21R_Se$Etape4_3<- apply(d21R_Se[, c(28,29)] > (0.25 * d21R_Se$RE), 1, function(x) ifelse(any(x), paste(colnames(d21R_Se)[c(28,29)][x], collapse=', '), 'None'))
     d21R_Se <- d21R_Se[,c(1, 39:41)]
     
     d21R_Se[!(d21R_Se$Etape4_2 == "None"),"Groupem"] <- paste0(d21R_Se[!(d21R_Se$Etape4_2 == "None"),"Etape4_1"],
                                                                  " ", d21R_Se[ !(d21R_Se$Etape4_2 == "None"),"Etape4_2"] )
     d21R_Se[(d21R_Se$Etape4_2 == "None") & !(d21R_Se$Etape4_3 == "None") ,"Groupem"] <- paste0(d21R_Se[(d21R_Se$Etape4_2 == "None") & !(d21R_Se$Etape4_3 == "None") ,"Etape4_1"],
                                                                                                  " ", d21R_Se[(d21R_Se$Etape4_2 == "None") & !(d21R_Se$Etape4_3 == "None") ,"Etape4_3"] )
     d21R_Se[is.na(d21R_Se$Groupem), "Groupem"] <-  paste0(d21R_Se[is.na(d21R_Se$Groupem),"Etape4_1"],     " ", "Rx")
     
     #//////////// Pi
     d21R_Pi <- df21R[df21R$Etape4_1 == "Pi",]
     d21R_Pi$Etape4_2<- apply(d21R_Pi[, c(2,8:11,17:18)] > (0.25 * d21R_Pi$RE), 1, function(x) ifelse(any(x), paste(colnames(d21R_Pi)[ c(2,8:11,17:18)][x], collapse=', '), 'None'))
     #  c(2, 8:14,17:18)
     d21R_Pi$Etape4_3<- apply(d21R_Pi[, c(28,30)] > (0.25 * d21R_Pi$RE), 1, function(x) ifelse(any(x), paste(colnames(d21R_Pi)[c(28,30)][x], collapse=', '), 'None'))
     d21R_Pi <- d21R_Pi[,c(1, 39:41)]
     
     d21R_Pi[!(d21R_Pi$Etape4_2 == "None"),"Groupem"] <- paste0(d21R_Pi[!(d21R_Pi$Etape4_2 == "None"),"Etape4_1"],
                                                                  " ", d21R_Pi[ !(d21R_Pi$Etape4_2 == "None"),"Etape4_2"] )
     d21R_Pi[(d21R_Pi$Etape4_2 == "None") & !(d21R_Pi$Etape4_3 == "None") ,"Groupem"] <- paste0(d21R_Pi[(d21R_Pi$Etape4_2 == "None") & !(d21R_Pi$Etape4_3 == "None") ,"Etape4_1"],
                                                                                                  " ", d21R_Pi[(d21R_Pi$Etape4_2 == "None") & !(d21R_Pi$Etape4_3 == "None") ,"Etape4_3"] )
     d21R_Pi[is.na(d21R_Pi$Groupem), "Groupem"] <-  paste0(d21R_Pi[is.na(d21R_Pi$Groupem),"Etape4_1"],     " ", "Rx")
     
     #----------------------------
     df21R_E4 <- rbind(d21R_Ep[, c(1,5)], d21R_Se[, c(1,5)], d21R_Pi[, c(1,5)])
     df21R_E4 <- left_join(df21R_E4, linkGroup)
     df21R_E4 <- df21R_E4[, c(1,3)]
     #----------------------------
     fileNamee <- paste0 ("eta1/R_", scenClim[i], "_", scenHarv[j], "_", codYear[k], "_",  rep,  "_E4.csv")
     fwrite(df21R_E4, file = fileNamee)
     rm(df21R_E4)
     
     
     rm(d21R_Pi, d21R_Se, d21R_Ep)                     
     
     df21R <- df21R[!(df21R$Etape4_1 %in% c("Ep","Pi", "Se")),]
     
    
     
     
     
     # Étape 5 groupe d'essences ou essence  > 0.25  
     
     
     df21R$Etape5_e <-  colnames(df21R[c(2, 8:14,17:18)])[apply(df21R[c(2, 8:14,17:18)],1,which.max)]
     df21R$Etape5_eV <- apply(X=df21R[c(2, 8:14,17:18)], MARGIN=1, FUN=max)
     df21R$Etape5_g <-  colnames(df21R[c(28:30)])[apply(df21R[c(28:30)],1,which.max)]
     df21R$Etape5_gV <- apply(X=df21R[c(28:30)], MARGIN=1, FUN=max)
     
     df21R <- df21R[,c(1,24,40:43)]
     
     
   
     df21R[ df21R$Etape5_eV <=  (0.25 * df21R$RE), "Etape5_e"] <- "None"
     df21R[ df21R$Etape5_gV <=  (0.25 * df21R$RE), "Etape5_g"] <- "None"
     
     df21R$Groupem <- "Rx Rx"
     df21R[!(df21R$Etape5_g == "None") , "Groupem"] <- paste0("Rx", " ", df21R[!(df21R$Etape5_g == "None"),  "Etape5_g"]) 
     df21R[!(df21R$Etape5_e == "None") , "Groupem"] <- paste0("Rx", " ", df21R[!(df21R$Etape5_e == "None"),  "Etape5_e"]) 
     
     
     
     #----------------------------
     df21R_E5 <- df21R[, c(1,7)]
     df21R_E5 <- left_join(df21R_E5, linkGroup)
     df21R_E5 <- df21R_E5[, c(1,3)]
     #----------------------------
     fileNamee <- paste0 ("eta1/R_", scenClim[i], "_", scenHarv[j], "_", codYear[k], "_",  rep,  "_E5.csv")
     fwrite(df21R_E5, file = fileNamee)
     rm(df21R_E5)
     
     
     #  Lecture files étapes RE , combine and save
     
     pat <- paste0("R_", scenClim[i], "_", scenHarv[j], "_", codYear[k], "_",  rep )
     fil <- list.files("eta1", full.names = TRUE)
     files <- fil[grep(pat, fil)]
     df <- rbind.fill(lapply(files, fread, header=TRUE))
     file.remove(files)
     
     
     df1 <- left_join(IdTipCouv, df)
     
     dfR <- df1[df1$tipCouv == "R",]
     
     #SAVE RE in eta2   (check NAs)
     
     fileNamee <- paste0 ("eta2/FILE_", scenClim[i], "_", scenHarv[j], "_", codYear[k], "_",  rep,  "_RE.csv")
     fwrite(dfR, file = fileNamee)
     rm(dfR, df)
     
     
     ##  ADD Hardwood code(3rd) to MR 
     dfMR <- df1[df1$tipCouv == "MR",]
     dfMR1 <- left_join(dfMR, linkGroup)
     dfMR2 <- left_join(dfMR1, df21)
     
     #---------------------
     #Schema 6
     #----------------------   
     # Étape 2 when sp RE > 0.5 RE
     dfMR2$Etape2_1<- apply(dfMR2[, c(6:10, 18:19)] > (0.50 * dfMR2$FE), 1, function(x) ifelse(any(x), paste(colnames(dfMR2)[c(6:10, 18:19)][x], collapse=', '), 'None'))
     # dfMR2$Etape2_2<- apply(dfMR2[, c(29, 33:35)] > (0.50 * dfMR2$FE), 1, function(x) ifelse(any(x), paste(colnames(dfMR2)[c(29, 33:35)][x], collapse=', '), 'None'))
     
     dfMR2$Etape2_2e <- colnames(dfMR2[c(29, 33:35)])[apply(dfMR2[c(29, 33:35)],1,which.max)]
     dfMR2$Etape2_2eV <- apply(X = dfMR2[c(29, 33:35)], MARGIN=1, FUN=max)
     dfMR2$Etape2_2 <- "None"
     dfMR2[dfMR2$Etape2_2eV >  (0.50 * dfMR2$FE), "Etape2_2" ] <- 
       dfMR2[dfMR2$Etape2_2eV >  (0.50 * dfMR2$FE), "Etape2_2e" ]
     dfMR2 <- dfMR2[,c(1:36,39)]
     
     
     names(dfMR2)[4] <- "Groupem1"
     
     dfMR2$Groupem <- paste0(dfMR2$Groupem1, " ", "Fx")
     dfMR2[!(dfMR2$Etape2_2 == "None"), "Groupem"] <- paste0( dfMR2[!(dfMR2$Etape2_2 == "None"), "Groupem1"], " ",
                                                              dfMR2[!(dfMR2$Etape2_2 == "None"), "Etape2_2"])
     dfMR2[!(dfMR2$Etape2_1 == "None"), "Groupem"] <- paste0( dfMR2[!(dfMR2$Etape2_1 == "None"), "Groupem1"], " ",
                                                              dfMR2[!(dfMR2$Etape2_1 == "None"), "Etape2_1"])
     dfMR3 <- dfMR2[, c(1,2,38)] 
     dfMR4 <- left_join(dfMR3, linkGroup)
     #idRAS, tipCouv, codGroup
     dfMR4 <-  dfMR4[, c(1,2,4)]
     
     fileNamee <- paste0 ("eta2/FILE_", scenClim[i], "_", scenHarv[j], "_", codYear[k], "_",  rep,  "_MR.csv")
     fwrite(dfMR4, file = fileNamee)
     rm(dfMR4, dfMR3, dfMR2, dfMR1, dfMR)
     rm(df1, df21R)

     
     
     
     #---------------------
     #Schema 5 (et partie du 6 FE)
     #----------------------
     df21 <- df21[!(df21$tipCouv =="MR"),]
     
     
     
     
     # Étape 1 FE when sp FE > 0.75
     df21$Etape1 <- apply(df21[, c(3:7, 15:16)] > (0.75 * df21$FE), 1, function(x) ifelse(any(x), paste(colnames(df21)[c(3:7, 15:16)][x], collapse=', '), 'None'))
     
     df21F_E1 <- df21[!(df21$Etape1 =="None"), c(1,34)]
     df21F_E1$Groupem <- paste0(df21F_E1$Etape1, " ", df21F_E1$Etape1)
     
     #---------------------
     df21F_E1 <-  df21F_E1[,c(1,3)]
     df21F_E1 <- left_join(df21F_E1, linkGroup)
     df21F_E1 <- df21F_E1[, c(1,3)]
     #--------------------
     
     fileNamee <- paste0 ("eta1/F_", scenClim[i], "_", scenHarv[j], "_", codYear[k], "_",  rep,  "_E1.csv")
     fwrite(df21F_E1, file = fileNamee)
     rm(df21F_E1)
     
     df21F <- df21[df21$Etape1 =="None",]
     
     
     # Étape 2 FE when sp > 0.5 and sp < 0.75
     df21F$Etape2_1 <- apply(df21F[, c(3:7, 15:16)] > (0.50 * df21F$FE) & df21F[, c(3:7, 15:16)] <= (0.75 * df21F$FE), 1, function(x) ifelse(any(x), paste(colnames(df21F)[c(3:7, 15:16)][x], collapse=', '), 'None'))
     df21F$Etape2_2 <- apply(df21F[, c(3:7, 15:16)] > (0.25 * df21F$FE) & df21F[, c(3:7, 15:16)] <= (0.50 * df21F$FE), 1, function(x) ifelse(any(x), paste(colnames(df21F)[c(3:7, 15:16)][x], collapse=', '), 'None'))
     
    
     df21F_E2_E <- df21F[!(df21F$Etape2_1 == "None") & !(df21F$Etape2_2 =="None"), c(1,35:36)]
     
     #erreur si nrow = 0
     if(nrow(df21F_E2_E > 0)){
     df21F_E2_E$Groupem <- paste0(df21F_E2_E$Etape2_1, " ", df21F_E2_E$Etape2_2)
     
     #---------------------
     df21F_E2_E <-  df21F_E2_E[,c(1,4)]
     df21F_E2_E <- left_join(df21F_E2_E, linkGroup)
     df21F_E2_E <- df21F_E2_E[, c(1,3)]
     #--------------------
     }
     
     df21F  <- df21F[(df21F$Etape2_1 == "None" | df21F$Etape2_2 =="None"),]
     
     # Étape 2 FE(suite) when sp > 0.5 and groupe > 0.25
     # df21F$Etape2_2G <- apply(df21F[, c(27, 31:33)] > (0.25 * df21F$FE) & df21F[, c(27, 31:33)] <= (0.50 * df21F$FE), 1, function(x) ifelse(any(x), paste(colnames(df21F)[c(27, 31:33)][x], collapse=', '), 'None'))
     
     df21F$Etape2_2e <- colnames(df21F[c(27, 31:33)])[apply(df21F[c(27, 31:33)],1,which.max)]
     df21F$Etape2_2eV <- apply(X = df21F[c(27, 31:33)], MARGIN=1, FUN=max)
     df21F$Etape2_2G <- "None"
     df21F[(df21F$Etape2_2eV > (0.25 * df21F$FE) & df21F$Etape2_2eV <= (0.50 * df21F$FE)), "Etape2_2G" ] <- 
       df21F[(df21F$Etape2_2eV > (0.25 * df21F$FE) & df21F$Etape2_2eV <= (0.50 * df21F$FE)), "Etape2_2e" ]
     df21F <- df21F[,c(1:36,39)]
     
     
     
     
     df21F[!( df21F$Etape2_1 == "None") & ( df21F$Etape2_2G == "None"), "Etape2_2G"] <- "Fx"
     
      df21F_E2_G <-  df21F[!( df21F$Etape2_1 == "None") & !( df21F$Etape2_2G =="None"),c(1,35,37)]
     #
      df21F_E2_G$Groupem <- paste0( df21F_E2_G$Etape2_1, " ",  df21F_E2_G$Etape2_2G)
     #
     #---------------------
      df21F_E2_G <-   df21F_E2_G[,c(1,4)]
      df21F_E2_G <- left_join( df21F_E2_G, linkGroup)
      df21F_E2_G <-  df21F_E2_G[, c(1,3)]
     #--------------------
     
      df21F_E2 <- rbind( df21F_E2_E,  df21F_E2_G)
     
     fileNamee <- paste0 ("eta1/F_", scenClim[i], "_", scenHarv[j], "_", codYear[k], "_",  rep, "_E2.csv")
     fwrite( df21F_E2, file = fileNamee)
     rm( df21F_E2,  df21F_E2_E,  df21F_E2_G)
     
     
     
      df21F <-  df21F[( df21F$Etape2_1 == "None" |  df21F$Etape2_2G =="None"),]
      
     
      # Étape 3 FE ombinaisson d'essences > 0.75  and one essence > 25 and > other in the group
      df21F$Etape3_1 <- "None"
       df21F[df21F$Er > (0.75 * df21F$FE), "Etape3_1"] <- "Er"
      # df21F$Etape3_1<- apply( df21F[, c(27, 18)] > (0.75 *  df21F$FE), 1, function(x) ifelse(any(x), paste(colnames( df21F)[c(27,18)][x], collapse=', '), 'None'))
      
      #//////////// Er
      d21F_Er <-  df21F[ df21F$Etape3_1 == "Er",]
      
      # d21F_Er$Etape3_2 <- with(d21F_Er, {
      #   names(d21F_Er)[(Eo > Es) * 3 + (Es > Eo) * 4]
      # })
      # 
      d21F_Er$Etape3_2 <- apply(d21F_Er[, c(3:4)] > (0.50* d21F_Er$Er), 1, function(x) ifelse(any(x), paste(colnames(d21F_Er)[c(3:4)][x], collapse=', '), 'None'))
      d21F_Er <- d21F_Er[,c(1, 38:39)]
      d21F_Er[d21F_Er$Etape3_2 == "None", "Etape3_2"] <- "Er"
      d21F_Er$Groupem <- paste0(d21F_Er$Etape3_1, " ", d21F_Er$Etape3_2)
      
     
      #-------------------------
       df21F_E3<- d21F_Er[, c(1,4)]
       df21F_E3 <- left_join( df21F_E3, linkGroup)
       df21F_E3 <-  df21F_E3[, c(1,3)]
      #--------------------------
      fileNamee <- paste0 ("eta1/F_", scenClim[i], "_", scenHarv[j], "_", codYear[k], "_",  rep, "_E3.csv")
      fwrite( df21F_E3, file = fileNamee)
      rm( df21F_E3)
      
      rm(d21F_Er)                     
      
       df21F <-  df21F[!( df21F$Etape3_1 == "Er"),]
      
      
       # Étape 4 FE combinaison d'essences > 0.5  and one essence autre groupe > 25 and > other in the group
       
       
       df21F$Etape4_1 <- "None"
       df21F[df21F$Er > (0.5 * df21F$FE), "Etape4_1"] <- "Er"
        # df21F$Etape4_1<- apply( df21F[, c(28:30)] > (0.5 *  df21F$RE), 1, function(x) ifelse(any(x), paste(colnames( df21F)[c(28:30)][x], collapse=', '), 'None'))
       
       #//////////// Er
        d21F_Er <-  df21F[ df21F$Etape4_1 == "Er",]
        d21F_Er$Etape4_2 <- apply( d21F_Er[, c(5:7, 15:16)] > (0.25 *  d21F_Er$FE), 1, function(x) ifelse(any(x), paste(colnames( d21F_Er)[c(5:7, 15:16)][x], collapse=', '), 'None'))
       #  c(2, 8:14,17:18)
        d21F_Er$Etape4_3 <- "None"
        
        d21F_Er[d21F_Er$Fi > (0.25 * d21F_Er$FE), "Etape4_3"] <- "Fi"
        
        # d21F_Ep$Etape4_3<- apply( d21F_Ep[, c(29,30)] > (0.25 *  d21F_Ep$RE), 1, function(x) ifelse(any(x), paste(colnames( d21F_Ep)[c(29:30)][x], collapse=', '), 'None'))
        d21F_Er <-  d21F_Er[,c(1, 39:41)]
        
        
        d21F_Er[!( d21F_Er$Etape4_3 == "None"),"Groupem"] <- paste0( d21F_Er[!( d21F_Er$Etape4_3 == "None"),"Etape4_1"],
                                                                     " ",  d21F_Er[ !( d21F_Er$Etape4_3 == "None"),"Etape4_3"] )
        d21F_Er[!( d21F_Er$Etape4_2 == "None"),"Groupem"] <- paste0( d21F_Er[!( d21F_Er$Etape4_2 == "None"),"Etape4_1"],
                                                                  " ",  d21F_Er[ !( d21F_Er$Etape4_2 == "None"),"Etape4_2"] )
        # d21F_Er[( d21F_Er$Etape4_2 == "None") & !( d21F_Er$Etape4_3 == "None") ,"Groupem"] <- paste0( dC_Er[( d21F_Er$Etape4_2 == "None") & !( d21F_Er$Etape4_3 == "None") ,"Etape4_3"] )
        d21F_Er[is.na( d21F_Er$Groupem), "Groupem"] <-  paste0( d21F_Er[is.na( d21F_Er$Groupem),"Etape4_1"],     " ", "Fx")
       
       
       #----------------------------
        df21F_E4 <-  d21F_Er[, c(1,5)]
        df21F_E4 <- left_join( df21F_E4, linkGroup)
        df21F_E4 <-  df21F_E4[, c(1,3)]
       #----------------------------
       fileNamee <- paste0 ("eta1/F_", scenClim[i], "_", scenHarv[j], "_", codYear[k], "_",  rep, "_E4.csv")
       fwrite( df21F_E4, file = fileNamee)
       rm( df21F_E4)
       
       
       rm(  d21F_Er)                     
       
       df21F <-  df21F[!( df21F$Etape4_1 == "Er"),]
       
       
       # Étape 5 FE association d'essences > 0.75  and one essence > 25 and > other in the group
       
       df21F$Etape5_1 <- apply(df21F[, c(31:33)] > (0.75 * df21F$FE), 1, function(x) ifelse(any(x), paste(colnames(df21F)[c(31:33)][x], collapse=', '), 'None'))
       
       #//////////// Fh
       d21F_Fh <- df21F[df21F$Etape5_1 == "Fh",]
       d21F_Fh$Etape5_2e <-  colnames(d21F_Fh[c(21,23)])[apply(d21F_Fh[c(21,23)],1,which.max)]
       d21F_Fh$Etape5_2eV <- apply(X=d21F_Fh[c(21,23)], MARGIN=1, FUN=max)
       d21F_Fh$Etape5_2 <- "None"
       d21F_Fh[d21F_Fh$Etape5_2eV > (0.25 *  d21F_Fh$FE), "Etape5_2"] <- substr(d21F_Fh[d21F_Fh$Etape5_2eV >(0.25 *  d21F_Fh$FE), "Etape5_2e"], 1,2)
       d21F_Fh <-  d21F_Fh[,c(1,40,43)]
       
       d21F_Fh[d21F_Fh$Etape5_2 == "None", "Etape5_2"] <- "Fh"
       d21F_Fh$Groupem <- paste0(d21F_Fh$Etape5_1, " ", d21F_Fh$Etape5_2)
       
       #//////////// Fi
       d21F_Fi <- df21F[df21F$Etape5_1 == "Fi",]
       d21F_Fi$Etape5_2e <-  colnames(d21F_Fi[c(6,15)])[apply(d21F_Fi[c(6,15)],1,which.max)]
       d21F_Fi$Etape5_2eV <- apply(X=d21F_Fi[c(6,15)], MARGIN=1, FUN=max)
       d21F_Fi$Etape5_2 <- "None"
       d21F_Fi[d21F_Fi$Etape5_2eV > (0.25 *  d21F_Fi$FE), "Etape5_2"] <- substr(d21F_Fi[d21F_Fi$Etape5_2eV >(0.25 *  d21F_Fi$FE), "Etape5_2e"], 1,2)
       d21F_Fi <-  d21F_Fi[,c(1,40,43)]
       
       d21F_Fi[d21F_Fi$Etape5_2 == "None", "Etape5_2"] <- "Fi"
       d21F_Fi$Groupem <- paste0(d21F_Fi$Etape5_1, " ", d21F_Fi$Etape5_2)
       
       
       #//////////// Ft
       d21F_Ft <- df21F[df21F$Etape5_1 == "Ft",]
       d21F_Ft$Etape5_2e <-  colnames(d21F_Ft[c( 4,  7, 16, 20,22)])[apply(d21F_Ft[c( 4,  7, 16, 20,22)],1,which.max)]
       d21F_Ft$Etape5_2eV <- apply(X=d21F_Ft[c( 4,  7, 16, 20,22)], MARGIN=1, FUN=max)
       d21F_Ft$Etape5_2 <- "None" 
       d21F_Ft[d21F_Ft$Etape5_2eV > (0.25 *  d21F_Ft$FE), "Etape5_2"] <- substr(d21F_Ft[d21F_Ft$Etape5_2eV >(0.25 *  d21F_Ft$FE), "Etape5_2e"], 1,2)
       d21F_Ft <-  d21F_Ft[,c(1,40,43)]
       
       d21F_Ft[d21F_Ft$Etape5_2 == "None", "Etape5_2"] <- "Ft"
       d21F_Ft$Groupem <- paste0(d21F_Ft$Etape5_1, " ", d21F_Ft$Etape5_2)
       
       
       
       
       #-------------------------
       df21F_E5<- rbind(d21F_Fh[, c(1,4)], d21F_Fi[, c(1,4)], d21F_Ft[, c(1,4)])
       df21F_E5 <- left_join(df21F_E5, linkGroup)
       df21F_E5 <- df21F_E5[, c(1,3)]
       #--------------------------
       fileNamee <- paste0 ("eta1/F_", scenClim[i], "_", scenHarv[j], "_", codYear[k], "_",  rep, "_E5.csv")
       fwrite(df21F_E5, file = fileNamee)
       rm(df21F_E5)
       
       rm(d21F_Fi, d21F_Fh, d21F_Ft)                     
       
       df21F <- df21F[!(df21F$Etape5_1 %in% c("Fi","Fh", "Ft")),]
       
       
       
       # Étape 6 FE association d'essences > 0.50  and one essence autre groupe 25 and > other 
       
       df21F$Etape6_1<- apply(df21F[, c(31:33)] > (0.50 * df21F$FE), 1, function(x) ifelse(any(x), paste(colnames(df21F)[c(31:33)][x], collapse=', '), 'None'))
       
       #//////////// Fh 21, 23
       d21F_Fh <- df21F[df21F$Etape6_1 == "Fh",]
       d21F_Fh$Etape6_2e <-  colnames(d21F_Fh[c(4, 6, 7, 15, 16, 20,22)])[apply(d21F_Fh[c(4, 6, 7, 15, 16, 20,22)],1,which.max)]
       d21F_Fh$Etape6_2eV <- apply(X=d21F_Fh[c(4, 6, 7, 15, 16, 20,22)], MARGIN=1, FUN=max)
       d21F_Fh$Etape6_2 <- "None"
       d21F_Fh[d21F_Fh$Etape6_2eV > (0.25 *  d21F_Fh$FE), "Etape6_2"] <- substr(d21F_Fh[d21F_Fh$Etape6_2eV >(0.25 *  d21F_Fh$FE), "Etape6_2e"], 1,2)
       d21F_Fh <-  d21F_Fh[,c(1:41,44)]
       
       
       d21F_Fh$Etape6_3 <- apply(d21F_Fh[, c(31:32)] > (0.25* d21F_Fh$FE), 1, function(x) ifelse(any(x), paste(colnames(d21F_Fh)[c(31:32)][x], collapse=', '), 'None'))
       d21F_Fh[d21F_Fh$Etape6_3 == "None", "Etape6_3" ] <- "Fx"
       d21F_Fh$Groupem <- "None"
       d21F_Fh[!(d21F_Fh$Etape6_2 == "None"), "Groupem"] <- paste0( d21F_Fh[!(d21F_Fh$Etape6_2 == "None"), "Etape6_1"], " ", 
                                                                     d21F_Fh[!(d21F_Fh$Etape6_2 == "None"), "Etape6_2"])
       d21F_Fh[d21F_Fh$Etape6_2 == "None", "Groupem"] <- paste0( d21F_Fh[d21F_Fh$Etape6_2 == "None", "Etape6_1"], " ", 
                                                                    d21F_Fh[d21F_Fh$Etape6_2 == "None", "Etape6_3"])
       
       d21F_Fh <-  d21F_Fh[,c(1,44)]
       
       #//////////// Fi 6,15
       d21F_Fi <- df21F[df21F$Etape6_1 == "Fi",]
       d21F_Fi$Etape6_2e <-  colnames(d21F_Fi[c(4, 7, 16, 20:23)])[apply(d21F_Fi[c(4, 7, 16, 20:23)],1,which.max)]
       d21F_Fi$Etape6_2eV <- apply(X=d21F_Fi[c(4, 7, 16, 20:23)], MARGIN=1, FUN=max)
       d21F_Fi$Etape6_2 <- "None"
       d21F_Fi[d21F_Fi$Etape6_2eV >(0.25 *  d21F_Fi$FE), "Etape6_2"] <- substr(d21F_Fi[d21F_Fi$Etape6_2eV >(0.25 *  d21F_Fi$FE), "Etape6_2e"], 1,2)
       d21F_Fi <-  d21F_Fi[,c(1:41,44)]
       
       
       d21F_Fi$Etape6_3 <- apply(d21F_Fi[, c(31,33)] > (0.25* d21F_Fi$FE), 1, function(x) ifelse(any(x), paste(colnames(d21F_Fi)[c(31:32)][x], collapse=', '), 'None'))
       d21F_Fi[d21F_Fi$Etape6_3 == "None", "Etape6_3" ] <- "Fx"
       d21F_Fi$Groupem <- "None"
       d21F_Fi[!(d21F_Fi$Etape6_2 == "None"), "Groupem"] <- paste0( d21F_Fi[!(d21F_Fi$Etape6_2 == "None"), "Etape6_1"], " ", 
                                                                    d21F_Fi[!(d21F_Fi$Etape6_2 == "None"), "Etape6_2"])
       d21F_Fi[d21F_Fi$Etape6_2 == "None", "Groupem"] <- paste0( d21F_Fi[d21F_Fi$Etape6_2 == "None", "Etape6_1"], " ", 
                                                                 d21F_Fi[d21F_Fi$Etape6_2 == "None", "Etape6_3"])
       
       d21F_Fi <-  d21F_Fi[,c(1,44)]
       
       
       #//////////// Ft     4,  7, 16, 20,22
       d21F_Ft <- df21F[df21F$Etape6_1 == "Ft",]
       d21F_Ft$Etape6_2e <-  colnames(d21F_Ft[c(6, 15, 21, 23)])[apply(d21F_Ft[c(6, 15, 21, 23)],1,which.max)]
       d21F_Ft$Etape6_2eV <- apply(X=d21F_Ft[c(6, 15, 21, 23)], MARGIN=1, FUN=max)
       d21F_Ft$Etape6_2 <- "None"
       d21F_Ft[d21F_Ft$Etape6_2eV > (0.25 *  d21F_Ft$FE), "Etape6_2"] <- substr(d21F_Ft[d21F_Ft$Etape6_2eV > (0.25 *  d21F_Ft$FE), "Etape6_2e"], 1,2)
       d21F_Ft <-  d21F_Ft[,c(1:41,44)]
       
       
       d21F_Ft$Etape6_3 <- apply(d21F_Ft[, c(32,33)] > (0.25* d21F_Ft$FE), 1, function(x) ifelse(any(x), paste(colnames(d21F_Ft)[c(31:32)][x], collapse=', '), 'None'))
       d21F_Ft[d21F_Ft$Etape6_3 == "None", "Etape6_3" ] <- "Fx"
       d21F_Ft$Groupem <- "None"
       d21F_Ft[!(d21F_Ft$Etape6_2 == "None"), "Groupem"] <- paste0( d21F_Ft[!(d21F_Ft$Etape6_2 == "None"), "Etape6_1"], " ", 
                                                                    d21F_Ft[!(d21F_Ft$Etape6_2 == "None"), "Etape6_2"])
       d21F_Ft[d21F_Ft$Etape6_2 == "None", "Groupem"] <- paste0( d21F_Ft[d21F_Ft$Etape6_2 == "None", "Etape6_1"], " ", 
                                                                 d21F_Ft[d21F_Ft$Etape6_2 == "None", "Etape6_3"])
       
       d21F_Ft <-  d21F_Ft[,c(1,44)]
       
       
       
       #-------------------------
       df21F_E6<- rbind(d21F_Fh[, c(1,2)], d21F_Fi[, c(1,2)], d21F_Ft[, c(1,2)])
       df21F_E6 <- left_join(df21F_E6, linkGroup)
       df21F_E6 <- df21F_E6[, c(1,3)]
       #--------------------------
       fileNamee <- paste0 ("eta1/F_", scenClim[i], "_", scenHarv[j], "_", codYear[k], "_",  rep,"_E6.csv")
       fwrite(df21F_E6, file = fileNamee)
       rm(df21F_E6)
       
       rm(d21F_Fi, d21F_Fh, d21F_Ft)                     
       
       df21F <- df21F[!(df21F$Etape6_1 %in% c("Fi","Fh", "Ft")),]
       
       
       # Étape 7 FE espèce ou groupe > 25 
       
       
       df21F$Etape7_1 <- "Fx"
       df21F$Etape7_2e <-  colnames(df21F[c(3, 4,5,6,7,15,16)])[apply(df21F[c(3, 4,5,6,7,15,16)],1,which.max)]
       df21F$Etape7_2eV <- apply(X= df21F[c(3, 4,5,6,7,15,16)], MARGIN=1, FUN=max)
       
       df21F$Etape7_2 <- "None"
       
       df21F[df21F$Etape7_2eV > (0.25 * df21F$FE), "Etape7_2"] <- substr(df21F[df21F$Etape7_2eV > (0.25 * df21F$FE), "Etape7_2e"], 1,2)
       df21F <-  df21F[,c(1:42,45)]
       
       
       df21F$Etape7_3e <-  colnames(df21F[c(27, 31:33)])[apply(df21F[c(27, 31:33)],1,which.max)]
       df21F$Etape7_3eV <- apply(X= df21F[c(27, 31:33)], MARGIN=1, FUN=max)
       df21F$Etape7_3 <- "None"
       
       df21F[df21F$Etape7_3eV > (0.25 * df21F$FE), "Etape7_3"] <- substr(df21F[df21F$Etape7_3eV > (0.25 * df21F$FE), "Etape7_3e"], 1,2)
       df21F <-  df21F[,c(1:43,46)]
       
       
       df21F[df21F$Etape7_3 == "None", "Etape7_3" ] <- "Fx"
       df21F$Groupem <- "None"
       df21F[!(df21F$Etape7_2 == "None"), "Groupem"] <- paste0( df21F[!(df21F$Etape7_2 == "None"), "Etape7_1"], " ", 
                                                                df21F[!(df21F$Etape7_2 == "None"), "Etape7_2"])
       df21F[df21F$Etape7_2 == "None", "Groupem"] <- paste0( df21F[df21F$Etape7_2 == "None", "Etape7_1"], " ", 
                                                                 df21F[df21F$Etape7_2 == "None", "Etape7_3"])
       
       df21F <-  df21F[,c(1,45)]
       
       #-------------------------
       df21F_E7<-  df21F
       df21F_E7 <- left_join(df21F_E7, linkGroup)
       df21F_E7 <- df21F_E7[, c(1,3)]
       #--------------------------
       fileNamee <- paste0 ("eta1/F_", scenClim[i], "_", scenHarv[j], "_", codYear[k], "_",  rep, "_E7.csv")
       fwrite(df21F_E7, file = fileNamee)
       rm(df21F_E7)
       rm(df21F)
                         
       #   FIN schéma 5
       
       
       
       #  Lecture files étapes FE , combine and save
       
       pat <- paste0("F_", scenClim[i], "_", scenHarv[j], "_", codYear[k], "_",  rep)
       fil <- list.files("eta1", full.names = TRUE)
       files <- fil[grep(pat, fil)]
       df <- rbind.fill(lapply(files, fread, header=TRUE))     ## ERREUR étape 2 Schéma 5
       # dff <- df[is.na(df$codGroup),]
       file.remove(files)
       
       
       df1 <- left_join(IdTipCouv, df)
       
       dfF <- df1[df1$tipCouv == "F",]
       
       #SAVE RE in eta2   (check NAs)
       
       fileNamee <- paste0 ("eta2/FILE_", scenClim[i], "_", scenHarv[j], "_", codYear[k], "_",  rep, "_FE.csv")
       fwrite(dfF, file = fileNamee)
       rm(dfF, df)
       
       
       ##  ADD Hardwood code(3rd) to MR 
       dfMF <- df1[df1$tipCouv == "MF",]
       dfMF1 <- left_join(dfMF, linkGroup)
       dfMF2 <- left_join(dfMF1, df21)
       
       #---------------------
       #Schema 6
       #----------------------   
       # Étape 2 when sp RE > 0.5 RE
       dfMF2$Etape2_1<- apply(dfMF2[, c(5,11:17,20,21)] > (0.50 * dfMF2$RE), 1, function(x) ifelse(any(x), paste(colnames(dfMF2)[c(5,11:17,20,21)][x], collapse=', '), 'None'))
       dfMF2$Etape2_2<- apply(dfMF2[, c(30:32)] > (0.50 * dfMF2$RE), 1, function(x) ifelse(any(x), paste(colnames(dfMF2)[c(30:32)][x], collapse=', '), 'None'))
       
       names(dfMF2)[4] <- "Groupem1"
       
       dfMF2$Groupem <- paste0(dfMF2$Groupem1, " ", "Rx")
       dfMF2[!(dfMF2$Etape2_2 == "None"), "Groupem"] <- paste0( dfMF2[!(dfMF2$Etape2_2 == "None"), "Groupem1"], " ",
                                                                dfMF2[!(dfMF2$Etape2_2 == "None"), "Etape2_2"])
       dfMF2[!(dfMF2$Etape2_1 == "None"), "Groupem"] <- paste0( dfMF2[!(dfMF2$Etape2_1 == "None"), "Groupem1"], " ",
                                                                dfMF2[!(dfMF2$Etape2_1 == "None"), "Etape2_1"])
       dfMF3 <- dfMF2[, c(1,2,39)] 
       dfMF4 <- left_join(dfMF3, linkGroup)
       #idRAS, tipCouv, codGroup
       dfMF4 <-  dfMF4[, c(1,2,4)]
       
       fileNamee <- paste0 ("eta2/FILE_", scenClim[i], "_", scenHarv[j], "_", codYear[k], "_",  rep, "_MF.csv")
       fwrite(dfMF4, file = fileNamee)
       rm(dfMF4, dfMF3, dfMF2, dfMF1, dfMF)
       rm(df1, df21)
       
       
       # Combine RE, MR, FE et MF et save as .csv in /

       
       pat <- paste0("FILE_", scenClim[i], "_", scenHarv[j], "_", codYear[k], "_",  rep)
       fil <- list.files("eta2", full.names = TRUE)
       files <- fil[grep(pat, fil)]
       df <- rbind.fill(lapply(files, fread, header=TRUE))     ## ERREUR étape 2 Schéma 5
       # dff <- df[is.na(df$codGroup),]
       file.remove(files)
       
       fileNamee <- paste0 ("Final/regr_", scenClim[i], "_", scenHarv[j], "_", codYear[k], "_",  rep, ".csv")
       fwrite(df, file = fileNamee)
       
       
       
  }



#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------  END  PARALLEL    -----------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------  
#---------------------
# Close clusters

stopCluster(cl)
showConnections()
showConnections(all = TRUE)
closeAllConnections()



       
     