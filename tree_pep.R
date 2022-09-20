### Formatting PEP tree data from Quebec ####

# The geodatabase containing the PEP tree data (placette-échantillon permanente) is available online at https://www.donneesquebec.ca/recherche/fr/dataset/placettes-echantillons-permanentes-1970-a-aujourd-hui

# download.file("ftp://transfert.mffp.gouv.qc.ca/Public/Diffusion/DonneeGratuite/Foret/DONNEES_FOR_ECO_SUD/Placettes_permanentes/PEP_GDB.zip", destfile = "raw_data/PEP.zip")



# unzip("raw_data/PEP.zip", exdir = "raw_data")
# database downloaded on May the 24th 2018

# Steps for data cleaning
# remove trees which dhp < 90mm and gaules
# remove some states (25 (intru), 44, 45, 46 (dead recruit), 34, 35, 36 (dead forgotten))
# renumber resurected trees (some recruits were given the IDs from dead trees)
# remove second death measures 
# Correct tree ids
# Manage renumbered trees : ETAT==29 and 50,52,54,55,56
# Simplify state to alive, dead, unknown
# Check for realistic growth rate

### PACKAGES ####
library(sf)
library(dplyr)
library(plyr)


### READ GDB LAYERS ####
st_layers("data/PEP_GDB/PEP.gdb")

# Plot coordinates
plot_xy <- st_read("data/PEP_GDB/PEP.gdb", layer = "PLACETTE")


# Plot date of measurements
plot_mes <- st_read("data/PEP_GDB/PEP.gdb", layer = "PLACETTE_MES")

# DHP + cote DENDRO_ARBRES
tree_mes <- st_read("data/PEP_GDB/PEP.gdb", layer = "DENDRO_ARBRES")

# Select premier inventaire

plot_mes <- plot_mes %>%
  filter(VERSION == "1er inv. 1970 à 1974")

# Species code

sps_code <- read.csv2("data/ref_spCode.csv")

### DATA CLEANING ####

# Keep only necessary columns 
# Remove abandonned plots (STATUT_MES %in% c(AB, RE, RL, DE, NT, SR))

table(plot_mes$STATUT_MES)
ID_PE_aband <- plot_mes %>% filter(!is.na(STATUT_MES))

plot_mes <- plot_mes %>%
  filter(!(ID_PE %in% ID_PE_aband$ID_PE)) %>%
  mutate(year_measured = as.integer(format(DATE_SOND, format="%Y"))) %>%
  select(ID_PE, ID_PE_MES, year_measured) 

plot_xy <- plot_xy %>%
  filter(ID_PE %in% plot_mes$ID_PE) %>%
  select(ID_PE, SHAPE) 

tree_mes <- tree_mes %>%
  filter(ID_PE %in% plot_mes$ID_PE) %>%
  select(-c(NO_MES, IN_ESS_NC, IN_1410, DEFOL_MIN, 
                   CL_QUAL:PRIO_RECOL, STADE_DEGR:VMB_HA))

### Change species code

tree_mes$sp_code <- sps_code$spCode[match(tree_mes$ESSENCE, sps_code$qc_code)]


length(which(!(plot_mes$ID_PE_MES %in% unique(tree_mes$ID_PE_MES))))
# 532 PE_MES are missing from tree_data...

# likely disturbed hence no tree

# Merge with plot info

tree_data0 <- plot_mes %>% 
  left_join(tree_mes, by = c("ID_PE", "ID_PE_MES"))

added_ID <- tree_data0 %>% filter(!(ID_PE_MES %in% tree_mes$ID_PE_MES))

levels(tree_data0$ETAT) <- c(levels(tree_data0$ETAT),"AllDead")
tree_data0$ETAT[tree_data0$ID_PE_MES %in% added_ID$ID_PE_MES] <- "AllDead"


### REMOVE ALL GAULES & DHP < 90mm ####
#but not dhp = 0 or NA because can be dead trees

tree_data <- tree_data0 %>% 
  filter(!ETAT %in% c("GV","GM","GA")) %>%
  filter(DHP > 90 | DHP == 0 | is.na(DHP))
  

#### RECLASSIFY STATES ####

#### ETAT == 26 -> state == "harvested"
#tree_data$state[which(tree_data$ETAT == 26)] <- "harvested"
#none

#### ETAT == 10,12,30,32,40,42 -> state == "alive"

tree_data$state[which(tree_data$ETAT %in% c(10,12,30,32,40,42))] <- "alive"

#### ETAT == 14,15,16,23,24,54,55,56 -> state == "dead"
tree_data$state[which(tree_data$ETAT %in% c(14,15,16,23,24,54,55,56) & is.na(tree_data$state))] <- "dead"

tree_data$state[which(tree_data$ETAT == "AllDead")] <- "dead" 

#### ETAT == 23,24 -> state == "mechanical death"
#tree_data$state[which(tree_data$ETAT %in% c(23,24))] <- "mechanical"

#### ETAT == 14,16,54,56 -> state == "stress and biotic death death"
#tree_data$state[which(tree_data$ETAT %in% c(14,15,16,54,55,56))] <- "stress"


#### Join xy ####

tree_data <- left_join(tree_data, plot_xy, by = "ID_PE")

#### Reshape - site x species matrix ####

sp_mat <- tree_data %>% 
  filter(state == "alive") %>%
  tidyr::pivot_wider(id_cols = c(ID_PE, ID_PE_MES, year_measured),
                     names_from = sp_code,
                     names_sort = TRUE,
                     values_from = DHP,
                     values_fn = length,
                     values_fill = 0)

# remettre les PE vides qui ont été enlevés
add_rm <- tree_data %>% 
  ungroup() %>% subset(!(ID_PE_MES %in% sp_mat$ID_PE_MES)) %>% 
  distinct(ID_PE, ID_PE_MES, year_measured)

sp_mat <- sp_mat %>% bind_rows(add_rm) %>%
  arrange(ID_PE, ID_PE_MES, year_measured)
sp_mat[is.na(sp_mat)] <- 0

# remove NA when computing basal area (to keep all stems of a species that where measured)

sp_BA <- tree_data %>% 
  filter(state == "alive") %>%
  tidyr::pivot_wider(id_cols = c(ID_PE, ID_PE_MES, year_measured),
                     names_from = sp_code,
                     names_sort = TRUE,
                     values_from = DHP,
                     values_fn = function(x) sum(pi*(x/(2 * 1000))^2, na.rm = T)*(10000/399.7312),
                     values_fill = 0)


add_rm <- tree_data %>% 
  ungroup() %>% 
  subset(!(ID_PE_MES %in% sp_BA$ID_PE_MES)) %>% 
  distinct(ID_PE, ID_PE_MES, year_measured)

sp_BA <- sp_BA %>% bind_rows(add_rm) %>%
  arrange(ID_PE, ID_PE_MES, year_measured)
sp_BA[is.na(sp_BA)] <- 0


#### SAVE ####

saveRDS(sp_mat, "data/sp_abun_1er_inv_sept22.RDS")
saveRDS(sp_BA, "data/sp_ba_1er_inv_sept22.RDS")


