# -----------------------------------------------------------------------------
# PROJECT:
#    Linking soil foodweb and soil functions dynamics under
#    changing agricultural practices managment
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# STEP:
#    0.2   Computing foodweb indices
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# NOTES:
# Script fourni par Elisa le 16.08.15, calcul de quelques indices de réseau
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# LIBRARIES:
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
# -----------------------------------------------------------------------------

# directory
setwd("F:/5. RH/4. Post-doctorant/2014 SOFIA V Coudrain/dossier juin 2016/Scripts finaux")

# -----------------------------------------------------------------------------
# FUNCTIONS
# -----------------------------------------------------------------------------
    
    # 1. FONCTIONS NIVEAUX TROPHIQUES
    GetIND <- function(web){  # fonction pour calculer le nombre, la moyenne et le max des niveaux trophiques
      
      tweb <- t(web)
      
      ## La somme des lignes doit être 1
      rs <- rowSums(tweb)
      for(i in 1:length(tweb[,1]))
        tweb[i,tweb[i,]>0] = tweb[i,tweb[i,]>0]/rs[i]
      
      nb.TL <- try(solve(diag(length(tweb[,1])) - tweb), T)
      
      if(class(nb.TL)=="try-error")
        nbTL <- rep(NA, length(tweb[,1]))
      
      if(class(nb.TL)!="try-error")
        nbTL <- rowSums(nb.TL)
      
      return(res <- list(nbTL, mean(nbTL),max(nbTL)))
      
    }
    
    
    # 2. CHEMIN FONG:BACT
    Path <- function(web, num){    #calcule pour chaque groupe la proportion de sa biomasse qui est
      #originaire du groupe source "num" (num = position du gpe dans la matrice web)
      tweb <- t(web)
      Vectpath <- rep(0, length(web[,1]))  #definition du vecteur pour la source des effets indirects
      Vectpath[num] <- 1
      
      ## Somme des lignes doit être 1
      rs <- rowSums(tweb)
      for(i in 1:length(tweb[,1]))
        tweb[i,tweb[i,]>0] = tweb[i,tweb[i,]>0]/rs[i]
      nb.TL <- try(solve(diag(length(tweb[,1])) - tweb), T)
      if(class(nb.TL)=="try-error")
        nbTL <- rep(NA, length(tweb[,1]))
      if(class(nb.TL)!="try-error")
        nbTL <- nb.TL%*%Vectpath
      nbTL
    }
    
    PropPath <- function(web, numP1, numP2, biom){   # calcule ratio des chemins bacterien sur fongique
      #biom = vecteur des biomasses des groupes présents dans web
      path1<-Path(web, numP1)
      path2<-Path(web, numP2)
      proppath<-(sum(biom * path1)) / (sum(biom * path2))
      proppath
    }
    
   
    # 3. DIVERSITE DE LIENS
    Nblinks<- function(web) {    #nombre de liens dans le réseau
      flux<-as.vector(web)
      flux<-flux[flux!=0]
      nblink<-length(flux)
      
      nblink
    }
    
    Divlinks<- function(web) {    #diversité de Shannon des liens
      flux<-as.vector(web)
      flux<-flux[flux!=0]
      Ftot<-sum(flux)
      div<-0
      for (i in 1:length(flux))   div<-div - flux[i]/Ftot*log2(flux[i]/Ftot)
      
      div
    }
    
    Evenlinks<- function(web) {    #evenness des liens
      even<-Divlinks(web)/Nblinks(web)
      
      even
    }
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Charger les matrices 
# -----------------------------------------------------------------------------
    # de flux
    Fluxmatlist.T0 <- readWorksheetFromFile("fluxmatT0.xlsx", sheet = 1:20, h = T, rownames = 1)
    Fluxmatlist.T2 <- readWorksheetFromFile("fluxmatT2.xlsx", sheet = 1:20, h = T, rownames = 1)
    Fluxmatlist.T4 <- readWorksheetFromFile("fluxmatT4.xlsx", sheet = 1:20, h = T, rownames = 1)

    # la matrice générale pour récolter les données des sites 
    basedat <- read.table("data_basefile.txt", h = T, sep = '\t', dec = ",")
    basedat <- Subset(basedat, Date!="Year3")
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Calcul d'indices
# -----------------------------------------------------------------------------

  # Niveaux trophiques
  IndTL_T0 <- list(); IndTL_T0 <- lapply(Fluxmatlist.T0, GetIND)
  IndTL_T2 <- list(); IndTL_T2 <- lapply(Fluxmatlist.T2, GetIND)
  IndTL_T4 <- list(); IndTL_T4 <- lapply(Fluxmatlist.T4, GetIND)

        # placer les mmean et max TL dans un df
        TLmeandatT0 <- TLmaxdatT0 <- c() 
        for(i in 1:length(IndTL_T0)){
          TLmeandatT0 <- c(TLmeandatT0,IndTL_T0[[i]][[2]])
          TLmaxdatT0 <- c(TLmaxdatT0,IndTL_T0[[i]][[3]])
        }
        TLmeandatT2 <- TLmaxdatT2 <- c() 
        for(i in 1:length(IndTL_T2)){
          TLmeandatT2 <- c(TLmeandatT2,IndTL_T2[[i]][[2]])
          TLmaxdatT2 <- c(TLmaxdatT2,IndTL_T2[[i]][[3]])
        }
        TLmeandatT4 <- TLmaxdatT4 <- c() 
        for(i in 1:length(IndTL_T4)){
          TLmeandatT4 <- c(TLmeandatT4,IndTL_T4[[i]][[2]])
          TLmaxdatT4 <- c(TLmaxdatT4,IndTL_T4[[i]][[3]])
        }
        
        dfTLmaxmaxT0 <- data.frame(TLmean = TLmeandatT0, TLmax = TLmaxdatT0)
        dfTLmaxmaxT2 <- data.frame(TLmean = TLmeandatT2, TLmax = TLmaxdatT2)
        dfTLmaxmaxT4 <- data.frame(TLmean = TLmeandatT4, TLmax = TLmaxdatT4)
        dfTLmaxmax <- data.frame(basedat[,1:5], rbind(dfTLmaxmaxT0, dfTLmaxmaxT2, dfTLmaxmaxT4))
        
        write.xlsx(dfTLmaxmax,"F:/5. RH/4. Post-doctorant/2014 SOFIA V Coudrain/dossier juin 2016/Scripts finaux/dfTLmaxmax.xlsx")
        
  # Rapport entre voie fongique et voie bactérienne

          # D'abord obtenir le vecteur de biomasses ##
              # Attention à ce que l'ordre des groupes corresponde avec la matrice de flux
          
          # ORDRER les groupes des matrices de flux dans l'ordre alphabétique!
          test[,order(names(test))]
          
          fluxmat_ordT0 <- lapply(Fluxmatlist.T0, function(x) x[order(names(x)), order(names(x))])
          fluxmat_ordT2 <- lapply(Fluxmatlist.T2, function(x) x[order(names(x)), order(names(x))])
          fluxmat_ordT4 <- lapply(Fluxmatlist.T4, function(x) x[order(names(x)), order(names(x))])
          
          dataWEB <- read.table("data_web.txt", h = TRUE, sep = '\t', dec=",")
          dataWEB_list <- lapply(split(dataWEB, dataWEB$Date), function (x) split(x, x$Plot))
          dataWEB_listT0 <- dataWEB_list[[1]]
          dataWEB_listT2 <- dataWEB_list[[2]]
          dataWEB_listT4 <- dataWEB_list[[3]]
          
          # Petit test pour être sûr que les noms concordent
          NameVecT4_list <- lapply(dataWEB_listT4, function(x) factor(x$node))
          NamefluxT4_list <- lapply(fluxmat_ordT4, function(x) factor(names(x)))
          equtest <- list()
          for(i in 1:20){
            equtest[[i]] <- all(NameVecT4_list[[i]] == NamefluxT4_list[[i]])
            }
          unlist(equtest) # ok
          
          BiomVecT0_list <- lapply(dataWEB_listT0, function(x) x$biomasstot)
          BiomVecT2_list <- lapply(dataWEB_listT2, function(x) x$biomasstot)
          BiomVecT4_list <- lapply(dataWEB_listT4, function(x) x$biomasstot)
          
              # Calcul des chemins
              PropPathBF_T0 <- c()
              for(i in 1:length(fluxmat_ordT0)){
                web <- fluxmat_ordT0[[i]]
                numP1 <- which(names(fluxmat_ordT0[[i]])=="BACT")
                numP2 <- which(names(fluxmat_ordT0[[i]])=="FUNG")
                biom <- BiomVecT0_list[[i]]
                biom[is.na(biom)] <- 0
                PropPathBF_T0 <- c(PropPathBF_T0, PropPathBF(web, numP1, numP2, biom))                   
                dfPropPathBF_T0 <- data.frame(BFpath = PropPathBF_T0)
              }
              
              PropPathBF_T2 <- c()
              for(i in 1:length(fluxmat_ordT2)){
                web   <- fluxmat_ordT2[[i]]
                numP1 <- which(names(fluxmat_ordT2[[i]])=="BACT")
                numP2 <- which(names(fluxmat_ordT2[[i]])=="FUNG")
                biom  <- BiomVecT2_list[[i]]
                biom[is.na(biom)] <- 0
                PropPathBF_T2 <- c(PropPathBF_T2, PropPathBF(web, numP1, numP2, biom))                   
                dfPropPathBF_T2 <- data.frame(BFpath=PropPathBF_T2)
              }
              
              PropPathBF_T4 <- c()
              for(i in 1:length(fluxmat_ordT4)){
                web <- fluxmat_ordT4[[i]]
                numP1 <- which(names(fluxmat_ordT4[[i]])=="BACT")
                numP2 <- which(names(fluxmat_ordT4[[i]])=="FUNG")
                biom <- BiomVecT4_list[[i]]
                biom[is.na(biom)] <- 0
                PropPathBF_T4 <- c(PropPathBF_T4,PropPathBF(web, numP1, numP2, biom))                  
                dfPropPathBF_T4 <- data.frame(BFpath=PropPathBF_T4)
              }
              
              PropPathBF_dftot <- data.frame(basedat[,1:5],rbind(dfPropPathBF_T0,dfPropPathBF_T2,dfPropPathBF_T4))


              # Rapport entre voie verte et voie brune

              PropPathBG_T0 <- c()
              for(i in 1:length(fluxmat_ordT0)){
                web <- fluxmat_ordT0[[i]]
                numP2 <- which(names(fluxmat_ordT0[[i]])=="Plantparts")
                numP1 <- which(names(fluxmat_ordT0[[i]])=="Detritus")
                biom <- BiomVecT0_list[[i]]
                biom[is.na(biom)] <- 0
                PropPathBG_T0 <- c(PropPathBG_T0, PropPathBF(web, numP1, numP2, biom))                   
                dfPropPathBG_T0 <- data.frame(BGpath = PropPathBG_T0)
              }
              
              PropPathBG_T2 <- c()
              for(i in 1:length(fluxmat_ordT2)){
                web <- fluxmat_ordT2[[i]]
                numP2 <- which(names(fluxmat_ordT2[[i]])=="Plantparts")
                numP1 <- which(names(fluxmat_ordT2[[i]])=="Detritus")
                biom <- BiomVecT2_list[[i]]
                biom[is.na(biom)] <- 0
                PropPathBG_T2 <- c(PropPathBG_T2, PropPathBF(web, numP1, numP2, biom))                   
                dfPropPathBG_T2 <- data.frame(BGpath = PropPathBG_T2)
              }
              
              PropPathBG_T4 <- c()
              for(i in 1:length(fluxmat_ordT4)){
                web <- fluxmat_ordT4[[i]]
                numP2 <- which(names(fluxmat_ordT4[[i]])=="Plantparts")
                numP1 <- which(names(fluxmat_ordT4[[i]])=="Detritus")
                biom <- BiomVecT4_list[[i]]
                biom[is.na(biom)] <- 0
                PropPathBG_T4 <- c(PropPathBG_T4, PropPathBF(web, numP1, numP2, biom))                   
                dfPropPathBG_T4 <- data.frame(BGpath = PropPathBG_T4)
              }
              
              
              PropPath_dftot <- data.frame(basedat[, 1:5],rbind(dfPropPathBG_T0, dfPropPathBG_T2, dfPropPathBG_T4), rbind(dfPropPathBF_T0,dfPropPathBF_T2,dfPropPathBF_T4))


                      # Figures ratio voies de flux d'énergie
                      
                      mydat <- PropPath_dftot[PropPath_dftot$Mod != "RR-PER",]
                      g <- ggplot(data = mydat, aes(x = Date, y = log(BGpath), group=interaction(Date, Mod)))+
                        geom_boxplot(aes(fill = factor(mydat$Mod))) +
                        ylab("Brown to green pathways ratio (log)")
                      
                      g
                      
                      g <- ggplot(data = mydat, aes(x = Date, y = log(BFpath), group=interaction(Date, Mod)))+
                        geom_boxplot(aes(fill = factor(mydat$Mod))) +
                        ylab("Bacterial to fungal pathways ratio (log)")
                      
                      g



  # DIVERSITE DE LIENS
    NbLT0 <- data.frame(NbL=unlist(lapply(fluxmat_ordT0,function(x) Nblinks(x))))
    NbLT2 <- data.frame(NbL=unlist(lapply(fluxmat_ordT2,function(x) Nblinks(x))))
    NbLT4 <- data.frame(NbL=unlist(lapply(fluxmat_ordT4,function(x) Nblinks(x))))
    
    DivLT0 <- data.frame(DivL=unlist(lapply(fluxmat_ordT0,function(x) Divlinks(x))))
    DivLT2 <- data.frame(DivL=unlist(lapply(fluxmat_ordT2,function(x) Divlinks(x))))
    DivLT4 <- data.frame(DivL=unlist(lapply(fluxmat_ordT4,function(x) Divlinks(x))))
    
    EvenLT0 <- data.frame(EvenL=unlist(lapply(fluxmat_ordT0,function(x) Evenlinks(x))))
    EvenLT2 <- data.frame(EvenL=unlist(lapply(fluxmat_ordT2,function(x) Evenlinks(x))))
    EvenLT4 <- data.frame(EvenL=unlist(lapply(fluxmat_ordT4,function(x) Evenlinks(x))))
    
    df_totdiv <- data.frame(PropPath_dftot, rbind(NbLT0, NbLT2, NbLT4), TLmean = dfTLmaxmax$TLmean,
                            rbind(DivLT0, DivLT2, DivLT4), rbind(EvenLT0, EvenLT2, EvenLT4))
    
    write.xlsx(df_totdiv,"F:/5. RH/4. Post-doctorant/2014 SOFIA V Coudrain/dossier juin 2016/Scripts finaux/df_totdiv.xlsx")

    
    
###################################
# Quelques relations #
###################################
# Lire les données sauvées précédemment

IndicesET <- read.table("data_indicesET.txt",sep="\t",header=T, dec=",")
str(IndicesET)
IndicesET$Mod <- factor(IndicesET$Mod, levels=c("CONV","RN","RT","RT-RR","RR-PER"))

# Quelques transformations utiles pour les graphs
IndicesET_wide <- melt(IndicesET, id.vars = c(1:5), measure.vars = c(6:11))
IndicesET_basedat <- cbind(basedat[,c(1:5,11,17,23,28,30,32,35,38,41,44,47,50,53)],IndicesET[,-c(1:5)])
Indbasedat_wide <- melt(IndicesET_basedat, id.vars = c(1:18), measure.vars = c(19:24))

# Quelques graphs

ggplot(data=IndicesET_wide, aes(x=value,group=variable)) +
  geom_histogram(fill="#880011") +  
  facet_grid(Date~variable, scales = "free") +
  theme_bw()+
  theme(panel.grid = element_blank())

ggplot(data=IndicesET_wide, aes(x=Mod, y=value,group=Mod, color=Mod)) +
  geom_boxplot(alpha=.4) +
  geom_point(alpha=.4, size=3) +
  facet_grid(variable~Date, scales = "free") +
  theme_bw()+
  theme(panel.grid = element_blank())

names(Indbasedat_wide)

# Relations avec les variables continues
plotfct1 <- function(df, v="C.N_0.5"){
  ggplot(data=df, aes_string(x=v, y="value")) +
    geom_point(size=3, aes(col=Mod)) +
    geom_smooth(alpha=0) +
    facet_grid(variable~Date, scales = "free") +
    theme_bw()+
    theme(panel.grid = element_blank())
}

# Inverser x et y pour les variables de fonctions
plotfct2 <- function(df, v="Nitrif_0.5"){
  ggplot(data=df, aes_string(y=v, x="value")) +
    geom_point(size=3, aes(col=Mod)) +
    geom_smooth(alpha=0) +
    facet_grid(Date~variable, scales = "free") +
    theme_bw()+
    theme(panel.grid = element_blank())
}


plotfct(Indbasedat_wide, v="C.N_0.5")
plotfct2(Indbasedat_wide)


