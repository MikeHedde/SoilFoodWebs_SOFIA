# -----------------------------------------------------------------------------
# PROJECT:
#    Linking soil foodweb and soil functions dynamics under
#    changing agricultural practices managment
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# STEP:
#    0.1   Computing the fluxes within the foodweb
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# LIBRARIES:
library(rJava) # n?cessaire pour certains packages
library(reshape2)
library(lme4)
require(FD)
require(entropart)
library(effects)
library(afex) # for Kenward-Roger model correction
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
library(FSA) # fonction "Subset": enl?ve compl?tement le niveau des facteurs omis
install_version("XLConnectJars", version = "0.2-12", repos = "http://cran.us.r-project.org") 
install_version("XLConnect", version = "0.2-12", repos = "http://cran.us.r-project.org")
library(xlsx) # Pour sauver les fichiers dans excel
library(XLConnect) 

librarian::shelf(tidyr, dplyr)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Loading trophic preference matrix 
# -----------------------------------------------------------------------------
entryweb <- read.table("data/raw-data/adjacency_matrix.txt", h = T, dec = ",")
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Mod?le utilisant l'estimation des biomasses des bact?ries et champignons
# -----------------------------------------------------------------------------
# Parameters definition
  nomsp   <- names(entryweb)[-1]    # List of nodes
  nbgp    <- length(nomsp)          # Total number of nodes
  spID    <- 1:29                   # Nodes' ID 
  web     <- tibble(entryweb) %>%
                  pivot_longer(cols = !ID, names_to = "variable") %>%
                  dplyr::filter(value > 0)
  nbinter <- length(web$value)      # Total number of edges
  res     <- data.frame(pred = match(web$ID, nomsp), prey = match(web$variable, nomsp), pred_pref = web$value) # remplacement des noms par les nr d'ID
  res$numinter <- 1:nbinter         # edges' ID

# Loading bacteria and fungi dataset
  Ref     <- 0.65 # valeur moyenne ratio biomasse bact sur champ dans SoilService, ? modifier pour regarder sensibilit? des r?sultats
  micro   <- tibble(read.table("data/raw-data/biom_micro.txt", h = T, dec = ",")) %>%
    mutate(valconv = Ref/mean(Ratio16S18S)) %>% # c'est le ratio beta/phi, voir fichier doc d'explication
    mutate(BACT = valconv * Ratio16S18S * WetW_MiOr/(1 + valconv * Ratio16S18S)) %>%
    mutate(FUNG = WetW_MiOr - BACT,
           a = 1,
           metaind = 0) %>%
    pivot_longer(cols = c(BACT, FUNG), names_to = "node", values_to = "biomasstot") %>%
    pivot_longer(cols = c(BactADN_m2, FungiADN_m2), values_to = "abundance") %>%
    select(Date, Plot, Block, Mod, node, biomasstot, abundance, a, metaind)

# Loading invertebrate biomass dataset
  datatot <- read.table("data/raw-data/biom_other.txt", h = T, dec = ",") %>%
                  filter(!is.na(indmass)) %>%
                  mutate(a = assimcoeff/100) %>%
                  select(Date, Plot, Block, Mod, node, biomasstot, abundance, a, metaind)
  
# Merge the two datasets
  biom <- bind_rows(datatot, micro) 

# Créer des listes de données pour chaque année et parcelle
  dataplot_list <- biom %>%
    filter(Date %in% c("Year0", "Year2", "Year4")) %>%
    group_by(Plot, Date) %>%
    nest()
  
# Load trophic flux model code
  source("analyses/1. trophic flux model.R")
  
# Compute fluxes for each plot and each year
  Fluxfct(dataplot_list[[3]][1], node, biomasstot, abundance, a, metaind)  #continuer
  


# -----------------------------------------------------------------------------
# Caract?ristiques des r?seaux
# -----------------------------------------------------------------------------
  # Extraction des noms des noeuds de chaque parcelle
  names.T0 <- names.T2 <- names.T4 <- list() 
  Fun_names <- function(dl, db) {
    for (i in 1:length(dataB.T0)) {
      dl[[i]] <-db[[i]]$name[which(db[[i]]$biomass>0 | is.na(db[[i]]$biomass))]  
    }
    return(dl)
  }
  
      names.T0 <- Fun_names(names.T0, dataB.T0)
      names.T2 <- Fun_names(names.T2, dataB.T2)
      names.T4 <- Fun_names(names.T4, dataB.T4)

  # int?gration des flux dans une matrice d'interactions
  weblist.T0 <- fun_web(F_list_T0)
  weblist.T2 <- fun_web(F_list_T2)
  weblist.T4 <- fun_web(F_list_T4)
      # beaucoup de flux n?gatifs (pourquoi?) mais la plupart n?gligeables, sauf pour le site T4A16:
      # lapply(weblist.T4, function(x) which(x< -10^-3)==T)

  # Remplacement des num?ros de lignes et colonnes par les noeuds correspondant
  # les matrices "weblist" contiennent tous les noeuds irrespectivement de leur pr?sence dans la parcelle
  # Il faut donc ?liminer les noeuds absents (biomasse de 0 dans la parcelle)
  fun_index <- function(db, wb, nm) {
    index <- replicate(20, list()) 
    for (i in 1:20) {
      index[[i]] <-which(db[[i]]$biomass>0|is.na(db[[i]]$biomass)) # s?lection de l'index des noeuds pr?sents dans le r?seaux, sans oublier les plantes et d?tritus (|is.na(db[[i]]$biomass)) 
      wb[[i]] <-wb[[i]][index[[i]],index[[i]]]                     # on ne garde que les noeuds dont la biomasse n'est pas nulle
      rownames(wb[[i]]) <- colnames(wb[[i]]) <-  nm[[i]]           # 
    }
    wb
  }

wmatT0_list <- fun_index(db=dataB.T0, wb=weblist.T0, nm=names.T0)
wmatT2_list <- fun_index(db=dataB.T2, wb=weblist.T2, nm=names.T2)
wmatT4_list <- fun_index(db=dataB.T4, wb=weblist.T4, nm=names.T4)
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Sauvegarde des matrices dans excel 
# -----------------------------------------------------------------------------  
      # La fonction "writeWorksheetToFile" du package "XLConnect" est tr?s utile pour ?a
      # On sp?cifie le nom du fichier qui va ?tre cr??, puis les donn?es, puis le nom de l'onglet, puis le nom de la premi?re colonne (ne pas oublier de sp?cifier une case vide d'abord "")
      # ignorer les warnings (?)

plotID <- as.character(unique(datatot[,"Plot"])) # vecteur du nom des parcelles

for(i in 1:20){    
  writeWorksheetToFile("fluxmatT0.xlsx", wmatT0_list[[i]],
                       sheet = paste("T0",plotID[[i]],sep=""),rownames=c("",names.T0[[i]]))
}

for(i in 1:20){    
  writeWorksheetToFile("fluxmatT2.xlsx",wmatT2_list[[i]],
                       sheet = paste("T2",plotID[[i]],sep=""),rownames=c("",names.T2[[i]]))
}

for(i in 1:20){    
  writeWorksheetToFile("fluxmatT4.xlsx",wmatT4_list[[i]],
                       sheet = paste("T4",plotID[[i]],sep=""),rownames=c("",names.T4[[i]]))
}


  ## Les matrices peuvent ensuite ?tre charg?es pour ne pas avoir ? refaire tous les calculs de flux ##
  # de nouveau en utilisant une fonction du package "XLConnect"
Fluxmatlist.T0 <- readWorksheetFromFile("fluxmatT0.xlsx", sheet = 1:20, h = T, rownames = 1)
Fluxmatlist.T2 <- readWorksheetFromFile("fluxmatT2.xlsx", sheet = 1:20, h = T, rownames = 1)
Fluxmatlist.T4 <- readWorksheetFromFile("fluxmatT4.xlsx", sheet = 1:20, h = T, rownames = 1)
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Calcul des niveaux trophiques 
# -----------------------------------------------------------------------------
TL_T0 <- list(); TL_T0 <- lapply(Fluxmatlist.T0, GetTL2)
TL_T2 <- list(); TL_T2 <- lapply(Fluxmatlist.T2, GetTL2)
TL_T4 <- list(); TL_T4 <- lapply(Fluxmatlist.T4, GetTL2)

## Cr?ation d'un dataframe contenant les niveaux trophiques ##
  ## c'est un peu du bricolage...

plotID <- as.character(unique(datatot[,"Plot"])) # vecteur du nom des parcelles

plotID.T0 <- paste0("T0",plotID)
plotID.T2 <- paste0("T2",plotID)
plotID.T4 <- paste0("T4",plotID)

dfTLT0 <- lapply(TL_T0, function(x) data.frame(TrL = x)) 
dfTLT02 <- list()
for(i in 1:length(plotID.T0)){
  dfTLT02[[i]] <- data.frame(commID=rep(plotID.T0[i],length(dfTLT0[[i]][,1])),Date=rep("Year0",length(dfTLT0[[i]][,1])), Plot=rep(plotID[i],length(dfTLT0[[i]][,1])), dfTLT0[[i]],GroupID=row.names(dfTLT0[[i]]))
}

dfTLT2 <- lapply(TL_T2, function(x) data.frame(TrL = x))
dfTLT22 <- list()
for(i in 1:length(plotID.T2)){
  dfTLT22[[i]] <- data.frame(commID=rep(plotID.T2[i],length(dfTLT2[[i]][,1])),Date=rep("Year2",length(dfTLT2[[i]][,1])), Plot=rep(plotID[i],length(dfTLT2[[i]][,1])), dfTLT2[[i]],GroupID=row.names(dfTLT2[[i]]))
}

dfTLT4 <- lapply(TL_T4, function(x) data.frame(TrL = x))
dfTLT42 <- list()
for(i in 1:length(plotID.T4)){
  dfTLT42[[i]] <- data.frame(commID=rep(plotID.T4[i],length(dfTLT4[[i]][,1])),Date=rep("Year4",length(dfTLT4[[i]][,1])), Plot=rep(plotID[i],length(dfTLT4[[i]][,1])), dfTLT4[[i]],GroupID=row.names(dfTLT4[[i]]))
}

# Coercion de toutes les parcelles
  # fonction "rbind.fill" du package "plyr" permet de merger des listes
dfTLall <- rbind(rbind.fill(dfTLT02), rbind.fill(dfTLT22), rbind.fill(dfTLT42))
  # Utilisation de la fonction "grep" pour ajouter la colonne des blocs et des pratiques 
    # Cette fonction rep?re les ?l?ments (i) dans la colonne sp?cifi?e (dfTLall$commID) et rempli la colonne (dfTLall$Block) avec la valeur sp?cifi?e ("B1")
dfTLall$Block <-  dfTLall$Mod <- rep(NA,length(dfTLall$commID)) # d'abord ajouter des colonnes vides
for (i in c("A01","A02","A03","A05","A06")){ dfTLall$Block[grep(i, as.character(dfTLall$commID), ignore.case=TRUE)] <-  "B1"}
for (i in c("A07","A08","A09","A10","A12")){ dfTLall$Block[grep(i, as.character(dfTLall$commID), ignore.case=TRUE)] <-  "B2"}
for (i in c("A13","A14","A16","A17","A18")){ dfTLall$Block[grep(i, as.character(dfTLall$commID), ignore.case=TRUE)] <-  "B3"}
for (i in c("A20","A21","A22","A23","A24")){ dfTLall$Block[grep(i, as.character(dfTLall$commID), ignore.case=TRUE)] <-  "B4"}
dfTLall$Block<- factor(dfTLall$Block) # transforme les caract?res en facteurs

for (i in c("A03","A12","A13","A22")){ dfTLall$Mod[grep(i, as.character(dfTLall$commID), ignore.case=TRUE)] <-  "CONV"}
for (i in c("A01","A09","A14","A24")){ dfTLall$Mod[grep(i, as.character(dfTLall$commID), ignore.case=TRUE)] <-  "RN"}
for (i in c("A02","A10","A18","A20")){ dfTLall$Mod[grep(i, as.character(dfTLall$commID), ignore.case=TRUE)] <-  "RT"}
for (i in c("A06","A07","A17","A23")){ dfTLall$Mod[grep(i, as.character(dfTLall$commID), ignore.case=TRUE)] <-  "RT-RR"}
for (i in c("A05","A08","A16","A21")){ dfTLall$Mod[grep(i, as.character(dfTLall$commID), ignore.case=TRUE)] <-  "RR-PER"}
dfTLall$Mod<- factor(dfTLall$Mod, levels=c("CONV","RN","RT","RT-RR","RR-PER"))

# Sauvegarde dans excel
write.xlsx(dfTLall,"F:/5. RH/4. Post-doctorant/2014 SOFIA V Coudrain/dossier juin 2016/Scripts finaux/dfTLall.xlsx")
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# Repr?sentation des r?seaux #
# -----------------------------------------------------------------------------
    # On repart des dataframes de base (utilis?s pour le calcul des flux)
    # Soit refaits, soit sauvegard?s et recharg?s

# Index num?raire des noeuds
index.T0<-lapply(dataB.T0,function(x) which(x$biomass>0|is.na(x$biomass)))
index.T2<-lapply(dataB.T2,function(x) which(x$biomass>0|is.na(x$biomass)))
index.T4<-lapply(dataB.T4,function(x) which(x$biomass>0|is.na(x$biomass)))

# Fonction pour repr?senter les r?seaux
  # les arguments de la fonction reprennent plusieurs des ?l?ments calcul?s pr?c?demment
PlotWeb <- function(dataB, index, TL, webTL, mylab="Mod", mycol, coltext){    
  biom <- dataB$biomass[index]
  nomweb <- nomsp[index]
  biom[which(is.na(biom))] <- 10^4    # estimation de la biomasse des plantes et d?tritus: ? adapter si possible d'?tre plus r?aliste
  Sweb <- length(index)
  webTL <- webTL[nomweb, nomweb] # CONTROLER QUE L'ORDRE DES NOMS DE WEBTL ET INDEX CORRESPOND!
  TL <- TL[nomweb]
  g = matrix(0, nrow = Sweb, ncol=3)
  g[,3]<-0.01*log10(biom*10) # grandeur des cercle: peut ?tre adapt?e
  g[,2]<-TL/4.5 #sum(TL)
  TLround<-round(TL)
  for (i in 1:max(round(TL))) {
    a <-TLround==i
    b <-1:Sweb
    b <-b[a]
    xaxis <-(1:length(b))/length(b)
    g[a,1] <-xaxis+ 0.5 - sum(xaxis)/length(b)
  }
  symbols(g[,1], g[,2], circles=g[,3], inches=FALSE, fg="grey70", bg = mycol, xlab="",ylab="",
          ylim=c(0.1,1),bty="n",xaxt="n",yaxt="n") 
  for (i in 1:Sweb){
    for (j in 1:Sweb){
      if (webTL[i,j]>0){ 
        arrows(g[i,1],g[i,2],g[j,1],g[j,2], lwd=(log10(webTL[i,j]*10^2)), 
               col="grey50",length=0.1)}
    }
  } 
  
  text(g[,1],g[,2],labels=nomweb,col= coltext,cex=0.7,font=2)
  mtext(mylab,3,cex=0.7)
}



## Quelques petites fonctions qui ?vitent les r?p?titions ##
  # tlab est ? adapter au titre d?sir?

plotwCONV <- function(dataB,index,TL,webTL,tlab="DATE CONV", mycol, coltext){
  for(i in c( 3, 10, 11, 18)){
    PlotWeb(dataB[[i]], index[[i]], TL[[i]], webTL[[i]],tlab, mycol, coltext)
  }
}


plotwRN <- function(dataB,index,TL,webTL,tlab="DATE RN", mycol, coltext){
  for(i in c( 1, 8, 12, 20)){
    PlotWeb(dataB[[i]], index[[i]], TL[[i]], webTL[[i]],tlab, mycol, coltext)
  }
}

plotwRT <- function(dataB,index,TL,webTL,tlab="DATE RT", mycol, coltext){
  for(i in c(2, 9, 15, 16)){
    PlotWeb(dataB[[i]], index[[i]], TL[[i]], webTL[[i]],tlab, mycol, coltext)
  }
}

plotwRTRR <- function(dataB,index,TL,webTL,tlab="DATE RT-RR", mycol, coltext){
  for(i in c(5,  6, 14, 19)){
    PlotWeb(dataB[[i]], index[[i]], TL[[i]], webTL[[i]],tlab, mycol, coltext)
  }
}

plotwRRPER <- function(dataB,index,TL,webTL,tlab="DATE RR-PER", mycol, coltext){
  for(i in c(4,  7, 13, 17)){
    PlotWeb(dataB[[i]], index[[i]], TL[[i]], webTL[[i]],tlab, mycol, coltext)
  }
}


#x11() 
par(mfrow=c(3,4))
plotwCONV(dataB = dataB.T0, index = index.T0, TL = TL_T0, webTL = Fluxmatlist.T0, tlab = "", mycol = "#CCCCCC", coltext = "black")
plotwCONV(dataB.T2,index.T2,TL_T2,Fluxmatlist.T2,"T2 CONV")
plotwCONV(dataB.T4,index.T4,TL_T4,Fluxmatlist.T4,"", mycol = "skyblue4", coltext = "white")
par(mfrow=c(1,1))

x11() 
par(mfrow=c(3,4))
plotwRN(dataB.T0,index.T0,TL_T0,Fluxmatlist.T0,"T0 RN")
plotwRN(dataB.T2,index.T2,TL_T2,Fluxmatlist.T2,"T2 RN")
plotwRN(dataB.T4,index.T4,TL_T4,Fluxmatlist.T4,"", mycol = "slateblue4", coltext = "white")
par(mfrow=c(1,1))

x11() 
par(mfrow=c(3,4))
plotwRT(dataB.T0,index.T0,TL_T0,Fluxmatlist.T0,"T0 RT")
plotwRT(dataB.T2,index.T2,TL_T2,Fluxmatlist.T2,"T2 RT")
plotwRT(dataB.T4,index.T4,TL_T4,Fluxmatlist.T4,"T", mycol = "forestgreen", coltext = "white")
par(mfrow=c(1,1))

x11() 
par(mfrow=c(3,4))
plotwRTRR(dataB.T0,index.T0,TL_T0,Fluxmatlist.T0,"T0 RT-RR")
plotwRTRR(dataB.T2,index.T2,TL_T2,Fluxmatlist.T2,"T2 RT-RR")
plotwRTRR(dataB.T4,index.T4,TL_T4,Fluxmatlist.T4,"", mycol = "darkgreen", coltext = "white")
par(mfrow=c(1,1))

x11() 
par(mfrow=c(3,4))
plotwRRPER(dataB.T0,index.T0,TL_T0,Fluxmatlist.T0,"T0 RR-PER")
plotwRRPER(dataB.T2,index.T2,TL_T2,Fluxmatlist.T2,"T2 RR-PER")
plotwRRPER(dataB.T4,index.T4,TL_T4,Fluxmatlist.T4,"T4 RR-PER")
par(mfrow=c(1,1))


## Comparaison entre modalit?s, un graph par parcelle ##

x11() 
par(mfrow = c(5,4),oma = c(5,4,0,0) + 0.1, mar = c(0,0,1,1) + 0.1) # r?duire les espaces entre graphs
plotwCONV(dataB.T0,index.T0,TL_T0,Fluxmatlist.T0,"T0 CONV", mycol = "#003300")
plotwRN(dataB.T0,index.T0,TL_T0,Fluxmatlist.T0,"T0 RN")
plotwRT(dataB.T0,index.T0,TL_T0,Fluxmatlist.T0,"T0 RT")
plotwRTRR(dataB.T0,index.T0,TL_T0,Fluxmatlist.T0,"T0 RT-RR")
plotwRRPER(dataB.T0,index.T0,TL_T0,Fluxmatlist.T0,"T0 RR-PER")

x11() 
par(mfrow = c(5,4),oma = c(5,4,0,0) + 0.1, mar = c(0,0,1,1) + 0.1)
plotwCONV(dataB.T2,index.T2,TL_T2,Fluxmatlist.T2,"T2 CONV")
plotwRN(dataB.T2,index.T2,TL_T2,Fluxmatlist.T2,"T2 RN")
plotwRT(dataB.T2,index.T2,TL_T2,Fluxmatlist.T2,"T2 RT")
plotwRTRR(dataB.T2,index.T2,TL_T2,Fluxmatlist.T2,"T2 RT-RR")
plotwRRPER(dataB.T2,index.T2,TL_T2,Fluxmatlist.T2,"T2 RR-PER")


x11() 
par(mfrow = c(5,4),oma = c(5,4,0,0) + 0.1, mar = c(0,0,1,1) + 0.1)
plotwCONV(dataB.T4,index.T4,TL_T4,Fluxmatlist.T4,"T4 CONV")
plotwRN(dataB.T4,index.T4,TL_T4,Fluxmatlist.T4,"T4 RN")
plotwRT(dataB.T4,index.T4,TL_T4,Fluxmatlist.T4,"T4 RT")
plotwRTRR(dataB.T4,index.T4,TL_T4,Fluxmatlist.T4,"T4 RT-RR")
plotwRRPER(dataB.T4,index.T4,TL_T4,Fluxmatlist.T4,"T4 RR-PER")
# -----------------------------------------------------------------------------