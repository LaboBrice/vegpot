#****************************
##  Produce Groupement essences
#   Selon cluster clasification log-ward (fastcluster package)
# Inputs are stacks from landis with  %biomass per species [0-1].
# 
# Code written for 17 species simulated in Landis (see objet speLa)
#  Clusterin was done considering 17 species &
# 
# Results are .tifs files (stackes) with cells probabilities of belonging to a cluster....
#***********************************
#

rm(list=ls())
gc()

#***********************************

require(raster)
require(utils)
require(dplyr)
require(tidyr)

#***********************************
#*
#*
#*
#*
setwd("P:/A_2022/DevenirVP/rapports/Rapport_Oct_2022/codeR")

# Load RF objet
clus20WardRF <- get(load("rfTestWardChord20.RData"))
rm(rfTestPar)



#************************************************
#------------------------------------------------------------
#     Read raster ID pixels
#------------------------------------------------------------
#************************************************

codRas <- raster("tifs/rasId.tif") #Raster codInteger pixel
dfId <- as.data.frame(values(codRas))
names(dfId) <- "codRas"
# rm(codRas)
spp17 <- read.csv("data/QC_RIA_SppList17.csv")$mycode

#-----------------------
# RScriptDir <- "C:\\Users\\jpascual\\Documents\\R\\R-4.1.1\\bin\\x64\\Rscript.exe"
# 
# rasterOptions(tmpdir = "U:/YBoulanger/Naturalité/Maps/anci/stacks/TEMP/")
# write("R_USER = U:/YBoulanger/Naturalité/Maps/anci/stacks/TEMP/", file=file.path(Sys.getenv('R_USER'), '.Renviron'))



# Outputs folders


outDir <- paste0("U:/YBoulanger/Naturalité/Maps/anci/stacks/TEMP/",Sys.Date())
dir.create(outDir)
setwd(outDir)

# Dir for intermediate and final results
# dir.create("eta1")
# dir.create("eta2")
dir.create("Final")

Dir <- paste0("U:\\YBoulanger\\Naturalité\\Maps\\anci")  # dir reading stacks



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



speMH <- c("ABIBAL", "ACERUB", "ACESAC", "BETALL", "BETPAP", "FAGGRA",
           "LARLAR", "PICGLA", "PICMAR", "PICRUB", "PINBAN", "PINRES",
           "PINSTR", "POPTRE", "QUERUB", "THUOCC","TSUCAN"  )
     
    


nRepli <- 5
# nRepli <- 2
years <- seq(0,130,10)
# years <- c(0,100,110,120,130)


codYear <- letters[seq(from = 1, to = 14)] # code for years 0-130
# codYear <- c("a", "k", "l", "m", "n")
#-----------------------------------------------
#-------------------- Lecture stack percentage 
#-----------------------------------------------




#/-/-/-/-/-/-/-/-/-/-
#Prepare parallel

require(doSNOW)
require(parallel)
require(foreach)


clusterN <- 6

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
   
  
  foreach(rep = 1:nRepli) %dopar% {     #     i <- 1 ;   j <- 1; k <- 1; rep <- 1
    gc()
    
    
    require(raster)
   
    require(dplyr)
    require(randomForest)
    require(rgdal)
   
    
    
    
    if(scenCli[i] == "baseline"){scenHar <- scenHar1} else{scenHar <- scenHar2}
    # read merged stack
    staPerc <- stack(paste0(Dir,"\\stacks\\merged","\\","stack_percBiom_","_",rep, scenCli[i],"_",scenHar[j],"_", years[k],".tif"))
   # Delete bioTot et maxAge
     staPerc <- staPerc[[1:17]]
    names(staPerc)[1:17] <- speMH  # code mffp 2 character
    
     # summary(values(staPerc)[[1:2]]) 
     
     # predStack <- predict(staPerc,clus20WardRF, type ='response', progress = 'text')
     
     
     
    df <- data.frame(values(staPerc))# raster values into df
    rm(staPerc)
    df0 <- cbind(dfId, df)
    dfTest <- df0[complete.cases(df0),]
    dfPre <- predict(clus20WardRF,dfTest, type ='prob')
    df1 <- cbind(dfTest[,1], as.data.frame(dfPre))
    names(df1)[1] <- "codRas"
    dfRasVal <- left_join(dfId, df1)
    rm(df,df0, df1, dfTest, dfPre)
    
    rasStack <- stack(codRas, codRas,codRas,codRas,codRas,
                     codRas,codRas,codRas,codRas,codRas,
                     codRas,codRas,codRas,codRas,codRas,
                     codRas,codRas,codRas,codRas,codRas)
    # rasStack1 <- subs(codRas,dfRasVal)
    
   
    
    for(ras in 1:20){  # ras <- 1
      values(rasStack[[ras]]) <- as.vector(dfRasVal[,ras+1])
    }
    
    rm(dfRasVal)
    
    # df21$IDras <- dfId$codRas        # Create ID par pixel
    # df21 <- df21[, c(18, 2:21)]
    
    nameOutStack <-  paste0 ("Final/regr_", scenClim[i], "_", scenHarv[j], "_", codYear[k], "_",  rep, ".tif")
    
    writeRaster(rasStack, filename = nameOutStack, format= 'GTiff', overwrite =TRUE )
     rm( rasStack)
    
    
    
    
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



       
     