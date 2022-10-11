#-------------------------------------------
# 
# Converti outputs biomass LANDIS-II into en stacks with biomass percentage
# stack per scenario CLI-Scenario HAR & year (0-130)
#  layers[[1:17]]  % biomass per species ,[[18]] biom total & [[19]] âge Max
#  Merge stack- Code prepared for  simulations study area divided into 5 areas and 5 replications 
#  Manuy steps are for specific problems of having our study area divided into 5 areas with buffers
#-------------------------------------------


#||||||||||||||||||||||||||||||||||||
rm(list=ls())
gc()
#////////////////////////////////

#||||||||||||||||||||||||||||||||||||
require(raster)
library(sp)   
library(rgdal)
library(dplyr)
require(sf)
library(reshape2)
require(fasterize)




#////////////////////////////////
#||||||||||||||||||||||||||||||||||||
inDir <- "U:\\YBoulanger\\LANDIS\\RIAmodv2"
outDir <- paste0("U:\\YBoulanger\\Naturalité\\Maps\\anci\\stacks\\",Sys.Date())

# RScriptDir <- "C:\\Users\\jpascual\\Documents\\R\\R-4.1.1\\bin\\x64\\Rscript.exe"
# tempDir = "U:/YBoulanger/Naturalité/Maps/anci/stacks/TEMP/"
# rasterOptions(tmpdir = tempDir)
# write("R_USER = U:/YBoulanger/Naturalité/Maps/anci/stacks/TEMP/", file=file.path(Sys.getenv('R_USER'), '.Renviron'))



# To reproject if needed
rasterProjection <- "+proj=lcc +lat_1=46 +lat_2=60 +lat_0=44 +lon_0=-68.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

# Read shapefile Sous-Domaines bioclimatiques (used to merge all the 5 zones not needed in other cases)
sDomBio<- st_read(dsn="P:\\RIA\\GIS\\shp\\Lamb\\suDomBioRIAqcLamb.shp", layer ="suDomBioRIAqcLamb")


# 

# raster pattern (landTypes)
ltRA <- raster("U:/Ecofor/_FcRIA/Jesus/QcRIA/QcRIA_5zones_June/QcRiaData/tif/landTypesQcRIA.tif")
ltRAS <- projectRaster(ltRA, crs = rasterProjection)
rm(ltRA)

#||||||||||||||||||||||||||||||||||||

# SCENARIOS

#Harvest    different names for baseline and RCP scenarios
scenHar1 <- c("BudwormBaselineFire", "BudwormBaselineFireBaselineHarvest")
scenHar2 <- c("GrowthBudwormProjectedFire", "GrowthBudwormProjectedFireBaselineHarvest")

#Clima


scenCli <- c("baseline", "RCP45", "RCP85")

# For files names
joker1 <- "scenarios"
joker2 <- "output/biomass"
joker3 <- "output/cohort-stats"

speLa <- c("ABIE.BAL", "ACER.RUB", "ACER.SAH", "BETU.ALL", "BETU.PAP", "FAGU.GRA",
         "LARI.LAR", "PICE.GLA", "PICE.MAR", "PICE.RUB", "PINU.BAN", "PINU.RES",
         "PINU.STR", "POPU.TRE", "QUER.RUB", "THUJ.SPP.ALL", "TSUG.CAN")

speQc <- c("SAB", "ERR", "ERS", "BOJ", "BOP", "HEG",
           "MEL", "EPB", "EPN", "EPR", "PIG", "PIR",
           "PIB", "PET", "CHR", "THO", "PRU")

 years <- seq(0,130,10)
 
 
 nRepl <- 5  #nombre de repliques

 #sousDomaines par zone
 zo123E <- c("1", "2E","3E")
 zo345O <- c("3O", "4O","5O")
 zo45E <- c( "4E","5E")
 zo6E <- c( "6E")
 zo6O <- c( "6O")
 
 # to rasterize SdomBios according to zones names. Avoid to calculate means in bi=ufer areas
 
 zon <- list(zo123E ,zo345O, zo45E, zo6E, zo6O ) # to rasterize SdomBios according to zones names

 zones <- c("123est" ,"345ouest", "45est", "6est", "6ouest" )


#/-/-/-/-/-/-/-/-/-/-
#Prepare parallel

require(doSNOW)
require(parallel)
require(foreach)


clusterN <- 16

# cl = makeCluster(clusterN, rscript=RScriptDir, type='SOCK')
cl = makeCluster(clusterN,  type='SOCK')


registerDoSNOW(cl)
#/-/-/-/-/-/-/-/-/-/-



# Area study divided into 5 zones for LANDIS
# Only one area makes zz loop inutil

 for(zz in seq_along(zones)) {   #   zz <- 2
  
  mask <- raster(paste0("U:/Ecofor/_FcRIA/Jesus/QcRIA/QcRIA_5zones_2020/",zones[zz],"/initial-communities_",zones[zz],".tif"))
  
  #Rasterize sousdomaines par zone (easy to merge; nomean to be obtained in buffer zones)
  
  sDomZon <- sDomBio[sDomBio$SDOM_BIO %in% zon[[zz]], ]
  sDomZon$idd <- 1

  sDomRAS<- fasterize( sDomZon,mask, field="idd")
  
  
   
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------    PARALLEL    --------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------



iniTime <- Sys.time()

  
  foreach(i = seq_along(scenCli)) %:%
 
   foreach(j = seq_along(scenHar1))%:%
   
      foreach(k = seq_along(years))%:%
    
      foreach(rep = 1:nRepl) %dopar%   {  #      k <- 1;    i <- 3 ;   j <- 2; rep <- 1    5 replications
    gc()
    require(dplyr)
    require(raster)
   
     
        
  
  if(scenCli[i] == "baseline"){scenHar <- scenHar1} else{scenHar <- scenHar2}   # differnt names if baseline or RCPs scens 
       
#----- get list of raster to be read
    list1 <- list.files(paste0(inDir,"/",zones[zz],"/",joker1, "/",
                              scenCli[i],"/",  scenHar[j], "/", rep,"/",joker2),
                               pattern =".tif",full.names =TRUE)
      list2 <- list1[grep(paste0("_",years[k],".tif"), list1)]
      list3 <- list2[-grep("TotalBiomass", list2)]   # plus vite de la calculer que de la lire
     
      
#-----       
     
#read rasters in list (loop), stack      
      for(nl in 1:length(list3)){   #  nl <- 1
        
        rasNL <- raster(list3[nl])
        values(mask)  <- values(rasNL)
        
        if(nl == 1) {sta <- stack(mask) }
        if(nl > 1) {sta[[nl]] <- mask }
      }
      
 # get total biomass   
      sta[[length(list3)+1]] <- sum(sta[[1:length(list3)]], na.rm=TRUE)
      names(sta[[length(list3)+1]]) <- "totBiomass"
      
# get % biomass per species  (0-1)     
      for(nl in 1:length(list3)){   #  nl <- 1}
        sta[[nl]] <-  sta[[nl]]/ sta[[length(list3)+1]]
      }
# names layer with species names      
      names(sta)[1:length(list3)] <-  speLa
      
# read rasters with max Age      
      list1 <- list.files(paste0(inDir,"/",zones[zz],"/",joker1, "/",
                                 scenCli[i],"/",  scenHar[j], "/", rep,"/",joker3),
                          pattern =".tif",full.names =TRUE)
      list2 <- list1[grep(paste0("-",years[k],".tif"), list1)]
      # list3 <- list2[grep(paste0(scenCli[i],"_",scenHar[j],"_"), list2)]
      
     maxAge <- raster(list2)
     values(mask)  <- values(maxAge)
# addd to stack
      sta[[nlayers(sta) +1]] <- mask
      
# crop  +  mask       
       staCrop <- crop(sta,sDomRAS)
       staMask <- mask(staCrop,sDomRAS)
       sta17 <- staMask[[length(list3)+1]]
# Write raster       
        tifName <- paste0(outDir, "\\","stack_percBiom_",zones[zz],"_",rep, scenCli[i],"_",scenHar[j],"_", years[k],".tif")
       # tifName <- paste0("C:\\Res\\","stack_percBiom_",zones[zz],"_",rep, scenCli[i],"_",scenHar[j],"_", years[k],".tif")
       writeRaster(sta17, filename= tifName, format="GTiff", overwrite=TRUE)
      }
    


  
fin <- Sys.time() - iniTime
print(fin)
print(zz)

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



rm(sDomBio)


# Code fait pour study area made in 5 areas
#Code merges the n (with n > 1) zones

scenCli <- c( "RCP85")


if(length(zones) > 1){
  
  # outDir <- paste0("U:\\YBoulanger\\NaturalitÃ©\\Maps\\anci\\stacks\\merged")
  # inDir <- paste0("U:\\YBoulanger\\NaturalitÃ©\\Maps\\anci\\stacks\\",Sys.Date()) 
  outDir <- paste0("U:\\YBoulanger\\Naturalité\\Maps\\anci\\stacks\\merged")
  inDir <- paste0("U:\\YBoulanger\\Naturalité\\Maps\\anci\\stacks\\",Sys.Date())
  
  #/-/-/-/-/-/-/-/-/-/-
  #Prepare parallel
  
  require(doSNOW)
  require(parallel)
  require(foreach)
  
  
  clusterN <- 14
  
  # cl = makeCluster(clusterN, rscript=RScriptDir, type='SOCK')
  
  cl = makeCluster(clusterN,  type='SOCK')
  
  registerDoSNOW(cl)
  #/-/-/-/-/-/-/-/-/-/-
  
  
  
  # mask_123E <- extend(mask_123E, ltRAS, value=NA)
  
  
  
  
  #---------------------------------------------------------------
  #---------------------------------------------------------------
  #---------------------------------------------------------------
  #---------------------    PARALLEL    --------------------------
  #---------------------------------------------------------------
  #---------------------------------------------------------------
  #---------------------------------------------------------------
  
  
  
  iniTime <- Sys.time()
  
  
  foreach(i = seq_along(scenCli) ) %:%
    
    foreach(j = seq_along(scenHar1))%:%

    foreach(k = seq_along(years)) %:%

    foreach(rep = 1:nRepl) %dopar%{     #     i <- 1 ;   j <- 1; k <- 1; rep <- 1
      # gc()
      require(dplyr)
      require(raster)
      
    if(scenCli[i] == "baseline"){scenHar <- scenHar1} else{scenHar <- scenHar2}
       # scenHar <- scenHar1
     
      # if(scenCli[i] == "baseline") {scenHar <- scenHar1} else {scenHar <- scenHar2}
      # if(i < 2) {scenHar =scenHar1} #  else {scenHar <- scenHar2}
      
      
      #List files to be read (stacks)
      list1 <- list.files(inDir,
                          pattern =".tif", full.names =TRUE)
      list2 <- list1[grep(paste0("_",years[k],".tif"), list1)]
      list3 <- list2[grep(paste0(rep,scenCli[i],"_",scenHar[j],"_"), list2)]
     
      
      
      #Open inside loop and rename ras1, ras2...
      for(nl in 1:length(list3)){   #  nl <- 1
        rasN <- stack(list3[nl])[[1:17]]
        
        nom <- paste0("ras",nl)
        assign(nom, rasN)
      }
      
      #merge with mosaic, function for cells with values in 2 layers (as masked with Sdom, only  possible in SDOM contours)
      # allRas <- merge(ras1,ras2,ras3,ras4, ras5)
      allRas <- mosaic(ras1,ras2,ras3,ras4, ras5, fun=mean)
      allRas <- stack(allRas)  # not sure
      rm(ras1,ras2,ras3,ras4,ras5, rasN)
      
      # rename stack layers
      names(allRas)[1:17] <-  speLa
     
     
      
      
      #Write merged stacks
      
      tifName <- paste0(outDir, "\\","stack_percBiom_","_",rep, scenCli[i],"_",scenHar[j],"_", years[k],".tif")
      writeRaster(allRas, filename= tifName, format="GTiff", overwrite=TRUE)
      
    }
  
  
  
  
  
  fin <- Sys.time() - iniTime
  print(fin)
  
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
  
   
  
  
  
}
