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

# Pour cette analyse on utilise les poids secs en accord avec les études existantes ##

  # Article utile : Reuman et al. 2009 ("Allometry of Body Size and Abundance in 166 Food Webs", ADVANCES IN ECOLOGICAL RESEARCH VOL. 41)

###########################################################

## Charger les données; le fichier ne contient pas de données pour les microorganismes (non disponibles) ##

Data_allom <- read.table("data_allometry.txt",header=TRUE, sep='\t', dec=",")
Data_allom$Mod <- factor(Data_allom$Mod, levels=c("CONV","RN","RT","RT-RR","RR-PER"))
# On ajoute les log des abondances et poids qui sont utilisés pour les relations allométriques
Data_allom$Log10N <- log10(Data_allom$N)
Data_allom$Log10M <- log10(Data_allom$M)

# Juste un test de l'effet bloc
with(Data_allom[Data_allom$Date=="Year0",],kruskal.test(Log10N,Block))
with(Data_allom[Data_allom$Date=="Year0",],kruskal.test(Log10M,Block))
with(Data_allom[Data_allom$Date=="Year2",],kruskal.test(Log10N,Block))
with(Data_allom[Data_allom$Date=="Year2",],kruskal.test(Log10M,Block))
with(Data_allom[Data_allom$Date=="Year4",],kruskal.test(Log10N,Block))
with(Data_allom[Data_allom$Date=="Year4",],kruskal.test(Log10M,Block))

# Test de la relation linéaire entre LogM et logN
  # Reuman et al. suggère différents tests dont jarquebera, shapiro et lilliefors

# Représentation --> très similaire entre traitements
p1 <- ggplot(Data_allom,aes(x=Log10M,y=Log10N,col=Mod)) + geom_point(size=0.5) + facet_grid(~Date)+
  geom_smooth(alpha=0)
p2 <- ggplot(Data_allom,aes(x=Log10M,y=Log10N,col=Mod)) + geom_point(size=0.5) + facet_grid(~Date)+
  stat_smooth(method = "lm",alpha=0)
p3 <- ggplot(Data_allom,aes(x=Log10M,y=Log10N,col=Mod)) + geom_point(size=0.5) + facet_grid(~Date)+
  stat_smooth(method = "lm", formula = y ~ poly(x, n=2, raw=TRUE),alpha=0)
x11();multiplot(p1,p2,p3)

## Histograms de la distribution des poids ##
ggplot(Data_allom, aes(Log10M, colour = Mod,fill=Mod)) +
  geom_histogram(binwidth =0.8,aes(y = ..density..))+
  facet_wrap(Date~Mod,ncol=5)
range(Data_allom$Log10M)
ggplot(Data_allom, aes(Log10M, fill = Mod)) +
  geom_density(alpha = 0.2,adjust=0.5) + xlim(-7, 2) +
  facet_wrap(~ Date)+
  theme_bw()

# Fonction pour tester la normalité des résidus
require(fBasics)
normfuns <- list(jarqueberaTest=jarqueberaTest, shapiroTest=shapiroTest, lillieTest=lillieTest)

fnormtest <- function(d) {
  norm_test <- res_reg <- list()  
  for (i in c("CONV","RN","RT","RT-RR","RR-PER")){
    res_reg[[i]] <- residuals(lm(Log10N~Log10M, data=d[d$Mod==i,]))
    norm_test[[i]] <- lapply(normfuns, function(f) f(res_reg[[i]]))                             
  }
  return(norm_test) 
}

res_listT0 <- fnormtest(Data_allom[Data_allom$Date=="Year0",])
res_listT2 <- fnormtest(Data_allom[Data_allom$Date=="Year2",])
res_listT4 <- fnormtest(Data_allom[Data_allom$Date=="Year4",])

getpval <- function(myp) {
  m <- matrix(NA, nrow=length(myp), ncol=3)
  for (i in (1:length(myp))){
    m[i,] <- sapply(1:3,function(x) myp[[i]][[x]]@test$p.value)
  }
  dfpval <- data.frame(names(myp), m)
  colnames(dfpval) <- c("Mod", names(normfuns))
  return(dfpval)
}

getpval(res_listT0)
getpval(res_listT2)
getpval(res_listT4)

# La linéarité des résidus est plutôt bien supportée pour les années 2012 et 2014 mais pas 2010
  # Problème: la pente ne se mesure pas au niveau du traitement mais au niveau de la parcelle
    
commlmtest <- function(dfv, plot=TRUE) {
  lmmod <- shapt <- list()
  x11();par(mfrow=c(2,4))
  for (i in unique(dfv$Plot)){
    lmmod[[i]] <- lm(Log10N~Log10M, data=dfv[dfv$Plot==i,])
    shapt[[i]] <- shapiro.test(residuals(lmmod[[i]]))
    if(plot==TRUE){plot(lmmod[[i]], which=c(1:2),main=paste(i,unique(dfv$Date)))}    
  }
  # return(lapply(lmmod,function(x) anova(x)))
  shapt
}

commlmtest(dfv=Data_allom[Data_allom$Date=="Year0" & Data_allom$Mod=="CONV",], plot=T)
commlmtest(dfv=Data_allom[Data_allom$Date=="Year0" & Data_allom$Mod=="RN",], plot=T)
commlmtest(dfv=Data_allom[Data_allom$Date=="Year0" & Data_allom$Mod=="RT",], plot=T)
commlmtest(dfv=Data_allom[Data_allom$Date=="Year0" & Data_allom$Mod=="RT-RR",], plot=T)
commlmtest(dfv=Data_allom[Data_allom$Date=="Year0" & Data_allom$Mod=="RR-PER",], plot=T)

commlmtest(dfv=Data_allom[Data_allom$Date=="Year2" & Data_allom$Mod=="CONV",], plot=T)
commlmtest(dfv=Data_allom[Data_allom$Date=="Year2" & Data_allom$Mod=="RN",], plot=T)
commlmtest(dfv=Data_allom[Data_allom$Date=="Year2" & Data_allom$Mod=="RT",], plot=T)
commlmtest(dfv=Data_allom[Data_allom$Date=="Year2" & Data_allom$Mod=="RT-RR",], plot=T)
commlmtest(dfv=Data_allom[Data_allom$Date=="Year2" & Data_allom$Mod=="RR-PER",], plot=T)

commlmtest(dfv=Data_allom[Data_allom$Date=="Year4" & Data_allom$Mod=="CONV",], plot=T)
commlmtest(dfv=Data_allom[Data_allom$Date=="Year4" & Data_allom$Mod=="RN",], plot=T)
commlmtest(dfv=Data_allom[Data_allom$Date=="Year4" & Data_allom$Mod=="RT",], plot=T)
commlmtest(dfv=Data_allom[Data_allom$Date=="Year4" & Data_allom$Mod=="RT-RR",], plot=T)
commlmtest(dfv=Data_allom[Data_allom$Date=="Year4" & Data_allom$Mod=="RR-PER",], plot=T)


# Beaucoup moins linéaire du coup...
      # ... rélexion nécessaire sur la pertinence de la comparaison des pentes 
# tests avec des modèles gls, mais ça n'améliore rien


#####################################################
# Si on décide tout de même de comparer les pentes:

## Extraction des pentes et des intercepts N vs M ##

SlopeFun <- function(dsl) {
  dfslope <- ddply(dsl, .(Date, Plot, Block, Mod), summarise,  slope = NA, intercept = NA) # création d'un dataframe vide
  lmmod <- list() # liste vide qui contiendra les modèles
  slparam <- interparam <- numeric() # vecteur qui contiendra les valeurs de pentes et d'intercepts
  for (i in unique(dsl$Plot)){
    lmmod[[i]] <- lm(Log10N~Log10M, data=dsl[dsl$Plot==i,])
    slparam[i] <- coefficients(lmmod[[i]])[[2]]
    interparam[i] <- coefficients(lmmod[[i]])[[1]]
  }
  dfslope$slope <- slparam
  dfslope$intercept <- interparam
  return (dfslope)
}

dfslopeT0 <- SlopeFun(Data_allom[Data_allom$Date=="Year0",])
dfslopeT2 <- SlopeFun(Data_allom[Data_allom$Date=="Year2",])
dfslopeT4 <- SlopeFun(Data_allom[Data_allom$Date=="Year4",])

dfslope_all <- rbind(dfslopeT0,dfslopeT2,dfslopeT4)
range(dfslope_all$slope)
# [1] -1.486945 -1.029553  Pente plus forte qu'attendu (-0.75 à -1), pourquoi?
  # Peut-être parce que les groupes sont trop "grossiers", la pente diminue si on se place au niveau de l'espèce (Sechi et al. 2014)
    # L'abondance des protozoaires est peut-être sur-estimée par la méthode utilisée

# Sauvegarde dans excel
write.xlsx(dfslope_allY,"c:/Users/Valerie/dfslope_allY.xlsx")

hist(dfslope_allY$slope) # distribution des valeurs
# Quelques modèles
mslTall <- lm(slope~Mod*Date,data=dfslope_allY)
predictmeans::residplot(mslTall) # Le package "predictmeans" permet de visualiser les résidus avec "residplot" (alternative à plot(mod)); "predictmeans::" à cause du conflit avec le package gplots, on peut aussi détacher gplots: detach(package:gplots)
mslTall2 <- lm(slope~Mod+Date,data=dfslope_allY)
residplot(mslTall2)
mslT0 <- lm(slope~Mod,data=dfslopeT0)
residplot(mslT0)
mslT2 <- lm(slope~Mod,data=dfslopeT2)
residplot(mslT2)
mslT4 <- lm(slope~Mod,data=dfslopeT4)
residplot(mslT4)
lapply(list(mslTall,mslTall2,mslT0,mslT2,mslT4),function(x) anova(x))

# Rien de significatif

###################################################"

## Exploration des relations entre pente et les diverses variables mesurées ##

# Chargement des données
soil.dat <- read.table("data_basefile.txt",header=TRUE, sep='\t', dec=",")
soil.dat  <- Subset(soil.dat,Date!="Year3",drop=T)
soil.dat$Mod <- factor(soil.dat$Mod, levels=c("CONV","RN","RT","RT-RR","RR-PER"))
# dfslope_allY <- read.table("dfslope_allY.txt",header=TRUE, sep='\t', dec=",")
# dfslope_allY$Mod <- factor(dfslope_allY$Mod, levels=c("CONV","RN","RT","RT-RR","RR-PER"))

# Tout mettre dans un seul df
SlSodf <- cbind(dfslope_allY,soil.dat[,-c(2:5)])

# Après c'est un peu du bricolage pour définir une fonction qui permet d'automatiser un peu les choses

slvarFun <- function(d,v=13){ # v=13 est aléatoire, c'est juste pour indiquer une variable par défaut (et éviter que le code bug...)
  d <-  melt(d[,c(1:6,v)], id=1:6, na.rm=TRUE) # Astuce pour que le nom de la colonne s'appelle "value" indépendemment de la variable sélectionnée; vérifier que ID correspondent bien aux variables "descriptives" (années, parcelles, blocs, pratiques)
  mslTall <- lm(slope~value*Date,data=d)
  mslTall2 <- lm(slope~ Date + value,data=d)
  mslTall3 <- lm(slope~value,data=d)
  mslT0 <- lm(slope~value,data=d[d$Date=="Year0",])
  mslT2 <- lm(slope~value,data=d[d$Date=="Year2",])
  mslT4 <- lm(slope~value,data=d[d$Date=="Year4",])
  return(flmslres <- list(mslTall,mslTall2,mslTall3,mslT0,mslT2,mslT4))
}

names(SlSodf) # indice de toutes les variables

lapply(res_C.N_0.5 <- slvarFun(SlSodf,v=13), function(x) anova(x))
lapply(res_C.N_0.5, function(x) {x11(); residplot(x)}) 

# ... avec toutes les variables à tester ...

# La relation avec la nitrification semble significative (interactions année*nitrif)
lapply(res_Nitrif_0.5 <- slvarFun(SlSodf,v=30), function(x) anova(x)) # Nitrate_0.5 -->  interaction sign!
residplot(res_Nitrif_0.5[[1]])
ggplot(SlSodf,aes(x=Nitrif_0.5, y=slope,col=Date)) + geom_point() +stat_smooth(method = "lm", formula = y ~x,alpha=0)
  # Apparemment en 2014 la pente est d'autant plus négative que la nitrification est élevée
    # difficile à intérpéter sans avoir les données microorganismes... 
       # Moins de noeuds de niveau trophique élevé (= pente plus raide) --> plus forte nitrification???












