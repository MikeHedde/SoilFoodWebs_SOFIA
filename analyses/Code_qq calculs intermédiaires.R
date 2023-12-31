########################################################
## Sp�cifier le directoire et charger les librairies ##
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
require(gplots) # attention: conflits av predictmeans
library(Hmisc)
library(MuMIn)
library(FSA) # fonction "Subset": enl�ve compl�tement le niveau des facteurs omis

########################################################

## Calcul des valeurs des noeuds � partir des donn�es par taxon ##
  
  # Important: Au sein de chaque noeud se trouve des taxons de masses et m�tabolismes diff�rents
  # qui sont plus ou moins repr�sent�s au sein du noeud
  # la biomasse et le m�tabolisme par individu pour un noeud doit donc �tre calcul�
  # comme une moyenne pond�r�e (par l'abondance des taxons)!  

DataWebs <- read.table("data_par taxon.txt",header=TRUE, sep='\t', dec=",")
DataWebs$community <- factor(paste0(DataWebs$Date,DataWebs$Plot))
DataWebs_plotlist <- split(DataWebs,DataWebs$community)

library(data.table) # permet d'arranger les donn�es et calculer la moyenne pond�r�e

names(DataWebs_plotlist)
dt1 <- dt2 <- dt3 <- dt.all <- list()  
for (i in 1:length(DataWebs_plotlist)){
  dt1[[i]] <- as.data.table(DataWebs_plotlist[[i]])
  dt2[[i]] <- as.data.frame(dt1[[i]][, lapply(as.list(.SD)[c("M","MRInd","AssimilEff")], weighted.mean, w = N), 
                                          by = list(Date,Plot,Block,Mod,TrophID)])
  dt3[[i]] <- as.data.frame(dt1[[i]][, lapply(as.list(.SD)["N"], sum),       # Somme des abondance par noeud
                                          by = list(Date,Plot,Block,Mod,TrophID)])
  dt.all[[i]] <- data.frame(dt2[[i]][,-5],node=dt2[[i]][,5],N=dt3[[i]][,"N"]) # Un peu d'ordre
  dt.all[[i]] <- data.frame (dt.all[[i]][,c(1:4,8,9,5:7)])
}
dt.all

dt.all.df <- plyr::rbind.fill(dt.all)
write.xlsx(dt.all.df,"c:/Users/Valerie/dt.all.df2.xlsx")

















