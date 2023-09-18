########################################################
## Spécifier le directoire et charger les librairies ##
setwd("C:/Users/Valerie/Documents/PESSAC/Analyses/Scripts finaux")

library(xlsx) # Pour sauver les fichiers dans excel
library(XLConnect) 
library(reshape2)
library(lme4)
require(FD)
require(entropart)
library(effects)
library(afex) # for Kenward-Roger model correction
#install.packages("glmmADMB", repos="http://r-forge.r-project.org", type="source")
library(glmmADMB)
library(multcomp) # post-hoc Tukey
library(MASS)
library(Rmisc) # summarySE and multiplot
library(predictmeans) # residplots
library(dunn.test)
library(ade4)
library(adegraphics)
library(vegan)
library(cheddar)
library(igraph)
library(devtools)
library(betapart)
library(ggplot2)
library(gridExtra)
library(gtools)
require(gplots)
library(Hmisc)
library(MuMIn)
library(FSA) # fonction "Subset": enlève complètement le niveau des facteurs omis

########################################################

## Charger les matrices de flux créées av le script "code_réseau trophique" ##
  # Utilisation du package XLConnect; Une liste de matrice est chargée pour chaque année
Fluxmatlist.T0 <- readWorksheetFromFile("fluxmatT0.xlsx",sheet = 1:20,header = TRUE,rownames=1)
Fluxmatlist.T2 <- readWorksheetFromFile("fluxmatT2.xlsx",sheet = 1:20,header = TRUE,rownames=1)
Fluxmatlist.T4 <- readWorksheetFromFile("fluxmatT4.xlsx",sheet = 1:20,header = TRUE,rownames=1)
#########################################

## Transformer les matrices pour qu'elles correspondent au format du package Cheddar ##
  # Cheddar travaille avec trois sortes de fichiers:
   # Un fichier cvs "properties" qui contient les attributs généraux des réseaux (voir l'onglet "Données pour Cheddar_properties")
   # Un fichier cvs "nodes" qui contient l'identité des noeuds, les biomasses, abondances,...
   # Un fichier cvs "trophic.links" qui contient les interactions.
   # A noter: le noms des colonnes est assez rigide dans cheddar et ne correspond pas forcément exactement à se que la colonne contient (par ex "diet.fraction")
  # Plusieurs vignettes explicatives se trouvent sur la page du package dans CRAN qu'il est utile de lire

# Préparer en amont des dossiers vides "obligatoires": un dossier nommé "communities" dans lequel on place un dossier pour chaque parcelle (le nom est libre)

# 1) Création du fichier "properties"

Dataprop <- read.table("data_WebProperties_cheddar.txt",header=TRUE, sep='\t', dec=",")
# Split par parcelle
Dataprop.T0 <-  split(Dataprop[Dataprop$date=="Year0",], Dataprop[Dataprop$date=="Year0",]$plot) 
Dataprop.T2 <-  split(Dataprop[Dataprop$date=="Year2",], Dataprop[Dataprop$date=="Year2",]$plot) 
Dataprop.T4 <-  split(Dataprop[Dataprop$date=="Year4",], Dataprop[Dataprop$date=="Year4",]$plot) 

# Ecriture des données dans les dossiers excel; adapter le path à l'emplacement du fichier
lapply(1:length(Dataprop.T0), function(i) write.csv(Dataprop.T0[[i]], 
                                                         file = paste0("F:/5. RH/4. Post-doctorant/2014 SOFIA V Coudrain/dossier juin 2016/Scripts finaux/TrophcommT0_Cheddar/communities/",Dataprop.T0[[i]]$title,"/properties", ".csv"),
                                                         row.names = FALSE)) 
lapply(1:length(Dataprop.T2), function(i) write.csv(Dataprop.T2[[i]], 
                                                         file = paste0("F:/5. RH/4. Post-doctorant/2014 SOFIA V Coudrain/dossier juin 2016/Scripts finaux/TrophcommT2_Cheddar/communities/",Dataprop.T2[[i]]$title,"/properties", ".csv"),
                                                         row.names = FALSE))
lapply(1:length(Dataprop.T4), function(i) write.csv(Dataprop.T4[[i]], 
                                                         file = paste0("F:/5. RH/4. Post-doctorant/2014 SOFIA V Coudrain/dossier juin 2016/Scripts finaux/TrophcommT4_Cheddar/communities/",Dataprop.T4[[i]]$title,"/properties", ".csv"),
                                                         row.names = FALSE))

# 2) Création du fichier "nodes"

DataWebs <- read.table("data_web.txt",header=TRUE, sep='\t', dec=",") 
IDnames.T0 <- paste("T0", unique(DataWebs$Plot),sep="") # vecteur de noms qui correspond au nom des dossiers excel
IDnames.T2 <- paste("T2", unique(DataWebs$Plot),sep="")
IDnames.T4 <- paste("T4", unique(DataWebs$Plot),sep="")

# Split par parcelle
DataWebs.T0 <- split(DataWebs[DataWebs$Date=="Year0",], DataWebs[DataWebs$Date=="Year0",]$Plot) 
DataWebs.T2 <- split(DataWebs[DataWebs$Date=="Year2",], DataWebs[DataWebs$Date=="Year2",]$Plot) 
DataWebs.T4 <- split(DataWebs[DataWebs$Date=="Year4",], DataWebs[DataWebs$Date=="Year4",]$Plot) 

# Ecriture des données dans les dossiers excel; adapter le path à l'emplacement du fichier
lapply(1:length(DataWebs.T0), function(i) write.csv(DataWebs.T0[[i]], 
                                                    file = paste0("F:/5. RH/4. Post-doctorant/2014 SOFIA V Coudrain/dossier juin 2016/Scripts finaux/TrophcommT0_Cheddar/communities/",IDnames.T0[i],"/nodes", ".csv"),
                                                    row.names = FALSE))
lapply(1:length(DataWebs.T2), function(i) write.csv(DataWebs.T2[[i]], 
                                                    file = paste0("F:/5. RH/4. Post-doctorant/2014 SOFIA V Coudrain/dossier juin 2016/Scripts finaux/TrophcommT2_Cheddar/communities/",IDnames.T2[i],"/nodes", ".csv"),
                                                    row.names = FALSE))
lapply(1:length(DataWebs.T4), function(i) write.csv(DataWebs.T4[[i]], 
                                                    file = paste0("F:/5. RH/4. Post-doctorant/2014 SOFIA V Coudrain/dossier juin 2016/Scripts finaux/TrophcommT4_Cheddar/communities/",IDnames.T4[i],"/nodes", ".csv"),
                                                    row.names = FALSE))



# 3) Création du fichier "trophic.links"
  # On part des matrices de flux

flux_list.T0 <- lapply(Fluxmatlist.T0, function (x) melt(as.matrix(x))) # transformation "wide to long"
flux_list.T0.2 <- list()
for (i in 1:length(flux_list.T0)) {
  flux_list.T0.2[[i]] <- flux_list.T0[[i]][flux_list.T0[[i]]$value != 0,] # enlève les interactions nulles
  rownames(flux_list.T0.2[[i]]) <- seq(length=nrow(flux_list.T0.2[[i]])) # juste pour enlever l'index initial des lignes
  colnames(flux_list.T0.2[[i]]) <- c("resource","consumer","diet.fraction") # les noms standards utilisés par Cheddar
}
flux_list.T0.2

flux_list.T2 <- lapply(Fluxmatlist.T2, function (x) melt(as.matrix(x)))
flux_list.T2.2 <- list()
for (i in 1:length(flux_list.T2)) {
  flux_list.T2.2[[i]] <- flux_list.T2[[i]][flux_list.T2[[i]]$value != 0,] 
  rownames(flux_list.T2.2[[i]]) <- seq(length=nrow(flux_list.T2.2[[i]]))
  colnames(flux_list.T2.2[[i]]) <- c("resource","consumer","diet.fraction")
}
flux_list.T2.2

flux_list.T4 <- lapply(Fluxmatlist.T4, function (x) melt(as.matrix(x)))
flux_list.T4.2 <- list()
for (i in 1:length(flux_list.T4)) {
  flux_list.T4.2[[i]] <- flux_list.T4[[i]][flux_list.T4[[i]]$value != 0,] 
  rownames(flux_list.T4.2[[i]]) <- seq(length=nrow(flux_list.T4.2[[i]]))
  colnames(flux_list.T4.2[[i]]) <- c("resource","consumer","diet.fraction")
}
flux_list.T4.2

# Ecriture des données dans les dossiers excel; adapter le path à l'emplacement du fichier

lapply(1:20, function(i) write.csv(flux_list.T0.2[[i]], 
                                   file = paste0("F:/5. RH/4. Post-doctorant/2014 SOFIA V Coudrain/dossier juin 2016/Scripts finaux/TrophcommT0_Cheddar/communities/",IDnames.T0[i],"/trophic.links", ".csv"),
                                   row.names = FALSE))

lapply(1:20, function(i) write.csv(flux_list.T2.2[[i]], 
                                   file = paste0("F:/5. RH/4. Post-doctorant/2014 SOFIA V Coudrain/dossier juin 2016/Scripts finaux/TrophcommT2_Cheddar/communities/",IDnames.T2[i],"/trophic.links", ".csv"),
                                   row.names = FALSE))

lapply(1:20, function(i) write.csv(flux_list.T4.2[[i]], 
                                   file = paste0("F:/5. RH/4. Post-doctorant/2014 SOFIA V Coudrain/dossier juin 2016/Scripts finaux/TrophcommT4_Cheddar/communities/",IDnames.T4[i],"/trophic.links", ".csv"),
                                   row.names = FALSE))


##########################################################################
#########################################################################

## Maintenant que les "communautés" sont créées, elles peuvent être utilisées
## avec les différentes fonctions du package Cheddar ##

# Charger les communautés #
  # Remarque: comme notre dossier "communities" contient plusieurs communautés, 
  # on a ce que Cheddar appelle une "Collection" de communautés et on utilise les fonctions correspondantes (voir la vignette du package)
CollectionT0 <- LoadCollection("F:/5. RH/4. Post-doctorant/2014 SOFIA V Coudrain/dossier juin 2016/Scripts finaux/TrophcommT0_Cheddar") # path de l'emplacement du dossier "communities"
CollectionT0 # > A collection of 20 communities
CollectionT2 <- LoadCollection("F:/5. RH/4. Post-doctorant/2014 SOFIA V Coudrain/dossier juin 2016/Scripts finaux/TrophcommT2_Cheddar") 
CollectionT2
CollectionT4 <- LoadCollection("F:/5. RH/4. Post-doctorant/2014 SOFIA V Coudrain/dossier juin 2016/Scripts finaux/TrophcommT4_Cheddar") 
CollectionT4


###########################################################################

## Calculs des indices de Bersier et Cohen ##
  
  # se référer aux articles de Bersier et al.(2002, Ecology) et Cohen et al. (2009, Proc. Natl Acad. Sci)
####################################################################

# 1) Indices de Bersier
BersierIndex_T0 <- lapply(CollectionT0, function(x) QuantitativeDescriptors(x, "diet.fraction" )) # "diet.fraction" correspond aux flux 
BersierIndex_T2 <- lapply(CollectionT2, function(x) QuantitativeDescriptors(x, "diet.fraction" ))
BersierIndex_T4 <- lapply(CollectionT4, function(x) QuantitativeDescriptors(x, "diet.fraction" ))

# Création d'un data frame contenat années, parcelles, blocs, traitements
BersierIndex_T0.df.list <- lapply(seq(length(BersierIndex_T0)), function(i){
  data.frame(community = rep(factor(CollectionCPS(CollectionT0)$title[i]),length(BersierIndex_T0[[i]][,1])), 
             Date=rep(factor(CollectionCPS(CollectionT0)$ date[i]),length(BersierIndex_T0[[i]][,1])),Plot=rep(factor(CollectionCPS(CollectionT0)$plot[i]),length(BersierIndex_T0[[i]][,1])),Mod=rep(factor(CollectionCPS(CollectionT0)$mod[i]),length(BersierIndex_T0[[i]][,1])),IDindex = factor(rownames(BersierIndex_T0[[i]])), BersierIndex_T0[[i]])
})

BersierIndex_T2.df.list <- lapply(seq(length(BersierIndex_T2)), function(i){
  data.frame(community = rep(factor(CollectionCPS(CollectionT2)$title[i]),length(BersierIndex_T2[[i]][,1])), 
             Date=rep(factor(CollectionCPS(CollectionT2)$ date[i]),length(BersierIndex_T2[[i]][,1])),Plot=rep(factor(CollectionCPS(CollectionT2)$plot[i]),length(BersierIndex_T2[[i]][,1])),Mod=rep(factor(CollectionCPS(CollectionT2)$mod[i]),length(BersierIndex_T2[[i]][,1])),IDindex = factor(rownames(BersierIndex_T2[[i]])), BersierIndex_T2[[i]])
})

BersierIndex_T4.df.list <- lapply(seq(length(BersierIndex_T4)), function(i){
  data.frame(community = rep(factor(CollectionCPS(CollectionT4)$title[i]),length(BersierIndex_T4[[i]][,1])), 
             Date=rep(factor(CollectionCPS(CollectionT4)$ date[i]),length(BersierIndex_T4[[i]][,1])),Plot=rep(factor(CollectionCPS(CollectionT4)$plot[i]),length(BersierIndex_T4[[i]][,1])),Mod=rep(factor(CollectionCPS(CollectionT4)$mod[i]),length(BersierIndex_T4[[i]][,1])),IDindex = factor(rownames(BersierIndex_T4[[i]])), BersierIndex_T4[[i]])
})

# Tout combiner dans un seul dataframe et le sauvegarder
BersierIndex_T0.df  <- plyr::rbind.fill(BersierIndex_T0.df.list)
BersierIndex_T2.df  <- plyr::rbind.fill(BersierIndex_T2.df.list)
BersierIndex_T4.df  <- plyr::rbind.fill(BersierIndex_T4.df.list)
BersierIndex_all.df <- rbind(BersierIndex_T0.df,BersierIndex_T2.df,BersierIndex_T4.df)
BersierIndex_all.df <- arrange(BersierIndex_all.df, IDindex, Date, Plot) # Mettre un peu d'ordre
BersierIndex_all.df$Block <- rep(NA,length(BersierIndex_all.df[,1])) # ajouter les blocs en utilisant "grep" (voir détails dans le script "Code_réseau trophique")
for (i in c("A01","A02","A03","A05","A06")){ BersierIndex_all.df$Block[grep(i, as.character(BersierIndex_all.df$community), ignore.case=TRUE)] <-  "B1"}
for (i in c("A07","A08","A09","A10","A12")){ BersierIndex_all.df$Block[grep(i, as.character(BersierIndex_all.df$community), ignore.case=TRUE)] <-  "B2"}
for (i in c("A13","A14","A16","A17","A18")){ BersierIndex_all.df$Block[grep(i, as.character(BersierIndex_all.df$community), ignore.case=TRUE)] <-  "B3"}
for (i in c("A20","A21","A22","A23","A24")){ BersierIndex_all.df$Block[grep(i, as.character(BersierIndex_all.df$community), ignore.case=TRUE)] <-  "B4"}
BersierIndex_all.df$Block<- factor(BersierIndex_all.df$Block)

write.xlsx(BersierIndex_all.df,"c:/Users/Valerie/BersierIndex_all.df.xlsx")

# 2) Indices de Cohen

CohenIndT0 <- NvMTriTrophicTable(CollectionT0)
CohenIndT2 <- NvMTriTrophicTable(CollectionT2)
CohenIndT4 <- NvMTriTrophicTable(CollectionT4)

# Les indices 18 à 25 ont des noms en double qui correspondent à la prise en compte ou non des interactions cannibales
  # Changer les noms pour différencier les deux
Indcoh <- attributes(CohenIndT0)$dimnames[[1]] # vecteur d'indices
Indcoh <- c(Indcoh[1:21], paste0(Indcoh[22:25],"_cann")) # les indices 22 à 25 prennent en compte les interactions cannibales 
CohenIndT0.df <- data.frame(CohenInd = factor(Indcoh),CohenIndT0,row.names=1:25)
CohenIndT2.df <- data.frame(CohenInd = factor(Indcoh),CohenIndT2,row.names=1:25)
CohenIndT4.df <- data.frame(CohenInd = factor(Indcoh),CohenIndT4,row.names=1:25)

# Création d'un dataframe des indices de Cohen
longCIT0 <- melt(CohenIndT0.df, id=1)
longCIT2 <- melt(CohenIndT2.df, id=1)
longCIT4 <- melt(CohenIndT4.df, id=1)
longCIT_all <- rbind(longCIT0,longCIT2,longCIT4)
names(longCIT_all)[2] <- "community"

longCIT_all$Date <- longCIT_all$Plot <- longCIT_all$Block <- longCIT_all$Mod <- rep(NA,length(longCIT_all[,1]))

for (i in c("T0","T2","T4")){ 
  if (i=="T0") { longCIT_all$Date[grep(i, as.character(longCIT_all$community), ignore.case=TRUE)] <-  "Year0"}
  else if (i=="T2"){ longCIT_all$Date[grep(i, as.character(longCIT_all$community), ignore.case=TRUE)] <-  "Year2"}
  else if (i=="T4"){ longCIT_all$Date[grep(i, as.character(longCIT_all$community), ignore.case=TRUE)] <-  "Year4"}
}
longCIT_all$Date <- factor(longCIT_all$Date)

for(i in c("A03","A12","A13","A22")){
  longCIT_all$Mod[grep(i, as.character(longCIT_all$community), ignore.case=TRUE)] <-  "CONV"
  for(i in c("A02","A10","A18","A20"))
    longCIT_all$Mod[grep(i, as.character(longCIT_all$community), ignore.case=TRUE)] <-  "RT"
  for(i in c("A06","A07","A17","A23"))
    longCIT_all$Mod[grep(i, as.character(longCIT_all$community), ignore.case=TRUE)] <-  "RT-RR"
  for(i in c("A01","A09","A14","A24"))
    longCIT_all$Mod[grep(i, as.character(longCIT_all$community), ignore.case=TRUE)] <-  "RN"
  for(i in c("A05","A08","A16","A21"))  
    longCIT_all$Mod[grep(i, as.character(longCIT_all$community), ignore.case=TRUE)] <-  "RR-PER"
}
longCIT_all$Mod <- factor(longCIT_all$Mod, levels = c("CONV","RN","RT","RT-RR","RR-PER"))

for (i in c("A01","A02","A03","A05","A06")){ longCIT_all$Block[grep(i, as.character(longCIT_all$community), ignore.case=TRUE)] <-  "B1"}
for (i in c("A07","A08","A09","A10","A12")){ longCIT_all$Block[grep(i, as.character(longCIT_all$community), ignore.case=TRUE)] <-  "B2"}
for (i in c("A13","A14","A16","A17","A18")){ longCIT_all$Block[grep(i, as.character(longCIT_all$community), ignore.case=TRUE)] <-  "B3"}
for (i in c("A20","A21","A22","A23","A24")){ longCIT_all$Block[grep(i, as.character(longCIT_all$community), ignore.case=TRUE)] <-  "B4"}
longCIT_all$Block<- factor(longCIT_all$Block)

longCIT_all$Plot <- substr(as.character(longCIT_all$community),3,6) # Découpage du nom des communautés pour obtenir le nom des parcelles
longCIT_all$Plot<- factor(longCIT_all$Plot)

write.xlsx(longCIT_all,"c:/Users/Valerie/longCIT_all.xlsx")




















