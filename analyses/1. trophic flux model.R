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
# NOTES:
# Modèle de flux trophiques créé par Elisa Thébault, code adapté par Valérie Coudrain
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
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# FUNCTIONS
# -----------------------------------------------------------------------------
# 1. Calcul des flux dans le r?seau
Fluxfct <- function(DBtest) {
  G<-vector(mode="numeric", length = nbinter)    
  B<-vector(mode="numeric", length = nbinter)
  A<-matrix(0, nrow = nbinter, ncol = nbinter)
  diag(A)<-1
  for (i in 1:nbgp) {
    predi<-(subset(res, res$pred==i))
    iprey<-predi$prey
    if(length(iprey)>0){
      if (length(iprey)>1){
        nobiom<-which(iprey==1|iprey==2)
        Sum_nB<-0
        if (length(nobiom)>0) {
          Sum_nB<-sum(predi$pred_pref[nobiom])
          G[predi$numinter[nobiom]]<-predi$pred_pref[nobiom]
        }
        withbiom<-which(iprey!=1&iprey!=2)
        if (length(withbiom)>0){
          Sum_wB<-sum(predi$pred_pref[withbiom]*DBtest$biomass[iprey[withbiom]])
          G[predi$numinter[withbiom]]<-(1-Sum_nB)*predi$pred_pref[withbiom]*DBtest$biomass[iprey[withbiom]]/Sum_wB  
          if (Sum_wB==0) { G[predi$numinter[withbiom]]<-predi$pred_pref[withbiom]*DBtest$biomass[iprey[withbiom]]  
          G[predi$numinter[nobiom]]<-predi$pred_pref[nobiom]/Sum_nb } 
        }
      }
      else  G[predi$numinter]<-1
      B[predi$numinter] <- G[predi$numinter]*DBtest$abundance[i]*DBtest$metaind[i]/(DBtest$a[i])
      
      preyi<-(subset(res,res$prey==i)) 
      ipred<-preyi$pred
      A[predi$numinter,preyi$numinter]<- A[predi$numinter,preyi$numinter] - G[predi$numinter]/(DBtest$a[i])
    }
  }
  F<-solve(A,B)         #F contient l'estimation des flux
  return(F)
}

# 2. integration des flux dans une matrice d'interactions
fun_web <- function(fl) {
  webmat <- replicate(20, matrix(0, nrow=nbgp, ncol=nbgp), simplify = F) # liste de matrices vides pour entrer les flux; le nombre de matrice correspond au nombre de parcelles (par ann?e)
  for (j in 1:20) {
    for (i in 1:(nbinter-1)) {
      webmat[[j]][res$prey[i], res$pred[i]] <- fl[[j]][i]    
      if (fl[[j]][i]<0) {print(fl[[j]][i])                # Les flux n?gatifs sont imprim?s 
        webmat[[j]][res$prey[i], res$pred[i]]<- 0       # ? sp?cifier si on veut ?liminer les flux n?gatifs, mais le probl?me est que ?a "cache" les probl?mes et donne une valeur artificielle de 0
      }
    }
  }
  webmat
}

# 3. calcul des niveaux trophiques
GetTL2 <- function(web){  
  tweb <- t(web) # transposer le r?seau 
  rs <- rowSums(tweb) # La somme des lignes doit ?tre 1
  for(i in 1:length(tweb[,1]))
    tweb[i,tweb[i,]>0] = tweb[i,tweb[i,]>0]/rs[i]
  
  nb.TL <- try(solve(diag(length(tweb[,1])) - tweb), T) 
  
  if(class(nb.TL)=="try-error")
    nbTL <- rep(NA, length(tweb[,1]))
  
  if(class(nb.TL)!="try-error")
    nbTL <- rowSums(nb.TL)
  
  nbTL
}
# -----------------------------------------------------------------------------
