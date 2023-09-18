#########################################################

# Il s'agit d'une ébauche, les indices à considérer doivent être
  # réfléchis et les modèles adaptés en fonction de la distribution des variables
# Penser aussi que les indices peuvent être influencés par la taille du réseau (nb de noeuds)

###############################################################

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

## Charger les matrices de flux créées av le script "code_Cheddar" ##
dfBersierInd <- read.table("data_BersierInd.txt",header=TRUE, sep='\t', dec=",") # 
dfBersierInd$Mod <- factor(dfBersierInd$Mod,levels=c("CONV", "RN", "RT", "RT-RR", "RR-PER"))
dfCohenInd <- read.table("data_CohenInd.txt",header=TRUE, sep='\t', dec=",") # 
dfCohenInd$Mod <- factor(dfCohenInd$Mod,levels=c("CONV", "RN", "RT", "RT-RR", "RR-PER"))


## Créer un format "wide" des données qui peut être utile dans ceraines analyses ##
  # Pour les indices de Bersier, on choisit les valeurs à considérer: ici "weighted"
    # mais "unweighted" ou "qualitative" pourrait aussi être considérés
      # Voir l'article de Bersier et al. 2002 pour une explication
dfBersierInd_wide <- dcast(dfBersierInd, community + Date+ Plot+ Block+Mod ~ BersierInd, value.var="Weighted")
dfCohenInd_wide <- dcast(dfCohenInd, community + Date+ Block+ Plot+Mod ~ CohenInd, value.var="value")

## Pour voir quelques corrélations ##
library(PerformanceAnalytics)
names(dfBersierInd)
chart.Correlation(dfBersierInd_wide[,c(6,10,15,16,17,18,22,23,24)], histogram=TRUE, pch=19)
# Certaines corrélation sont très fortes

# Fontcions pour automatiser un peu le calcul des modèles
  # relation indices - traitements
lmIndFun <- function(df,ix="Con") {
  mod1 <- lm(df[,ix] ~ Block + Mod*Date, data=df)
  mod2 <- lm(df[,ix]~ Block + Mod + Date, data=df)
  mod3 <- lm(df[df$Date=="Year0",ix]~ Block + Mod, data=df[df$Date=="Year0",])
  mod4 <- lm(df[df$Date=="Year2",ix]~ Block + Mod, data=df[df$Date=="Year2",])
  mod5 <- lm(df[df$Date=="Year4",ix]~ Block + Mod, data=df[df$Date=="Year4",])
   res  <- list(mod1,mod2,mod3,mod4,mod5)
}

x11(); par(mfrow=c(3,4));lapply(lmIndFun(dfBersierInd_wide),function(x) plot(x, which = 1:2))
lapply(lmIndFun(dfBersierInd_wide),function(x) anova(x))

x11(); par(mfrow=c(3,4));lapply(lmIndFun(dfCohenInd_wide, ix="S"),function(x) plot(x, which = 1:2))
lapply(lmIndFun(dfCohenInd_wide, ix="S"),function(x) anova(x))
#ggplot(dfCohenInd_wide, aes(x=Mod, y=S, group=Mod)) + geom_boxplot() + facet_wrap(~Date) 
ggplot(dfCohenInd_wide,aes(x=Date, y=S,col=Mod))+  geom_boxplot() + geom_jitter() # +stat_smooth(alpha=0,aes(group = Mod))


# ... tester avec les différents indices ... #

#################################################################

## Constatation : La normalité de la distribution des résidus n'est pas souvent respectée ##
## Une solution est d'utiliser un test non-paramétrique ##

## Le test de Friedman permet de prendre en compte le design en bloc (randomized complete block design)
## Post-hoc: posthoc.friedman.nemenyi.test(y), package PMCMR
library(PMCMR) 
#########################################

# Tester chaque année séparément

dfBersierInd_wide_T0 <- split(dfBersierInd_wide,dfBersierInd_wide$Date)[[1]]
dfBersierInd_wide_T2 <- split(dfBersierInd_wide,dfBersierInd_wide$Date)[[2]]
dfBersierInd_wide_T4 <- split(dfBersierInd_wide,dfBersierInd_wide$Date)[[3]]
dfCohenInd_wide_T0 <- split(dfCohenInd_wide,dfCohenInd_wide$Date)[[1]]
dfCohenInd_wide_T2 <- split(dfCohenInd_wide,dfCohenInd_wide$Date)[[2]]
dfCohenInd_wide_T4 <- split(dfCohenInd_wide,dfCohenInd_wide$Date)[[3]]

# Une fonction pour automatiser le test, pas forcément idéale
op <- options()
bvarnames <- names(dfBersierInd_wide)[-c(1:5)]
FMtest_B <- function(df){
  FMtl <- phFMtl <- plotsig <- list()
  for (bvar in bvarnames) {
    print(paste(bvar,round(friedman.test(df[,bvar] ~ Mod | Block, data=df)$p.value,3)))
    if(friedman.test(df[,bvar] ~ Mod | Block, data=df)$p.value <0.055){
      FMt <- friedman.test(df[,bvar] ~ Mod | Block, data=df)
      phFMt <- posthoc.friedman.nemenyi.test(df[,bvar] ~ Mod | Block, data=df)
      plotsig[[bvar]] <- ggplot(melt(df[,c(bvar,"Mod")]), aes(y= value, x= Mod)) + geom_boxplot() + scale_y_continuous(bvar)
      FMtl[[bvar]] <- FMt
      phFMtl[[bvar]] <- phFMt
    }
  }
  
  print(list(FMtl,phFMtl))
  if(length(plotsig)!=0){ multiplot(plotlist = plotsig,cols = 2 )}
}

FMtestBT0 <- FMtest_B(dfBersierInd_wide_T0) 
FMtestBT2 <- FMtest_B(dfBersierInd_wide_T2)
FMtestBT4 <- FMtest_B(dfBersierInd_wide_T4) 

bvarnames <- names(dfCohenInd_wide)[-c(1:5)] # choisir également uniquement les indices pertinents...
FMtestCT0 <- FMtest_B(dfCohenInd_wide_T0) 
FMtestCT2 <- FMtest_B(dfCohenInd_wide_T2)
FMtestCT4 <- FMtest_B(dfCohenInd_wide_T4) 


##########################################################

## Relation to state and function variables ##

##########################################################

soil.dat <- read.table("data_basefile.txt",header=TRUE, sep='\t', dec=",")
soil.dat <- Subset(soil.dat, Date!="Year3"& PlotID!="T4A16")
length(soil.dat[,1]);length(dfBersierInd_wide[,1]);length(dfCohenInd_wide[,1])
soil.dat_T0 <- split(soil.dat,soil.dat$Date)[[1]]
soil.dat_T2 <- split(soil.dat,soil.dat$Date)[[2]]
soil.dat_T4 <- split(soil.dat,soil.dat$Date)[[3]]

# Grouper les indices
Indexdat <- cbind(dfBersierInd_wide,dfCohenInd_wide[,-c(1:5)])
Indexdat_T0 <- split(Indexdat,Indexdat$Date)[[1]]
Indexdat_T2 <- split(Indexdat,Indexdat$Date)[[2]]
Indexdat_T4 <- split(Indexdat,Indexdat$Date)[[3]]

## Fonction pour automatiser un peu la procédure
  ## Mais rester conscient que le modèle linéaire n'est pas forcément le plus adapté
    ## éventuellement transformer certains indices (sqrt, log) avant

Indvarnames <- names(Indexdat)[-c(1:5)] # Choisir les indices 

Indtest <- function(df1,df2){     # df1: dataframe contenant les valeurs des indices; df2: dataframe contenant les variables x
  df2 <- melt(df2)                # nécessaire pour que le titre de la colonne soit "value". Il y a sans doute un moyen de sélectionner la colonne directement dans les arguments de la fct mais j'ai galéré av les erreurs d'environnement
  Indml  <- plotsig <- pval_l <- list()      # listes qui contiendront les modèles significatifs, les graphs pour les relations significatives et les p-values
  for (indvar in Indvarnames) {
    Indm <-lm(df1[,indvar]~df2$value)     # Le bloc n'est pas pris en compte, je ne vois pas l'utilité d'inclure un effet bloc pour les relations entre indices et variables quantitatives; l'année pourrait être ignorée aussi
    pval_l[[indvar]] <- round(anova(Indm)$Pr[1],3) #  
    Indml[[indvar]] <- Indm               # liste de modèles
    if(anova(Indm)$Pr[1] <0.055){         # pour les relations significatives...
      Indmsig <- Indm                       
      df3 <- melt(data.frame(df1[,c(1:5)],indvar = df1[,indvar])) # reformatage nécessaire pour ggplot2
      plotsig[[indvar]] <- ggplot(data.frame(df2,value2=df3[,length(df3)]), aes(y= value2, x= value)) + geom_point() + scale_y_continuous(indvar) + stat_smooth(method="lm")
      #Indml[[indvar]] <- Indmsig # si on veut uniquement retourner les modèles significatifs --> enlever "Indml[[indvar]] <- Indm"
    }  
  }
  if(length(plotsig)!=0){ multiplot(plotlist = plotsig,cols = 2 )} # représentation des relations sign si il y en a
  print(unlist(pval_l)) # impression des p-values
  return(Indml) # 
}

# Même chose avec transformation sqrt

SQRTIndtest <- function(df1,df2){      
  df2 <- melt(df2)                
  Indml  <- plotsig <- pval_l <- list()       
  for (indvar in Indvarnames) {
    Indm <-lm(sqrt(df1[,indvar])~df2$value)     
    pval_l[[indvar]] <- round(anova(Indm)$Pr[1],3)   
    Indml[[indvar]] <- Indm
    if(anova(Indm)$Pr[1] <0.055){         
      Indmsig <- Indm                       
      df3 <- melt(data.frame(df1[,c(1:5)],indvar = df1[,indvar])) 
      df3[,length(df3)] <- sqrt(df3[,length(df3)]) # rajouter la racine carrée
      plotsig[[indvar]] <- ggplot(data.frame(df2,value2=df3[,length(df3)]), aes(y= value2, x= value)) + geom_point() + scale_y_continuous(indvar) + stat_smooth(method="lm")
      #Indml[[indvar]] <- Indmsig # si on veut uniquement retourner les modèles significatifs --> enlever "Indml[[indvar]] <- Indm"
    }  
  }
  if(length(plotsig)!=0){ multiplot(plotlist = plotsig,cols = 2 )} 
  print(unlist(pval_l))
  return(Indml) 
}


# Même chose avec transformation log10

LOGIndtest <- function(df1,df2){      
  df2 <- melt(df2)                
  Indml  <- plotsig <- pval_l <- list()       
  for (indvar in Indvarnames) {
    Indm <-lm(log10(df1[,indvar])~df2$value)     
    pval_l[[indvar]] <- round(anova(Indm)$Pr[1],3)   
    Indml[[indvar]] <- Indm
    if(anova(Indm)$Pr[1] <0.055){         
      Indmsig <- Indm                       
      df3 <- melt(data.frame(df1[,c(1:5)],indvar = df1[,indvar])) 
      df3[,length(df3)] <- log10(df3[,length(df3)]) # rajouter la racine carrée
      plotsig[[indvar]] <- ggplot(data.frame(df2,value2=df3[,length(df3)]), aes(y= value2, x= value)) + geom_point() + scale_y_continuous(indvar) + stat_smooth(method="lm")
      #Indml[[indvar]] <- Indmsig # si on veut uniquement retourner les modèles significatifs --> enlever "Indml[[indvar]] <- Indm"
    }  
  }
  if(length(plotsig)!=0){ multiplot(plotlist = plotsig,cols = 2 )} 
  print(unlist(pval_l))
  return(Indml) 
}

# Eventuellement considérer le problème des tests multiples. 

names(soil.dat)

Indtest_C.N_0.5 <- Indtest(Indexdat,soil.dat[,11])
SQRTIndtest_C.N_0.5 <- SQRTIndtest(Indexdat,soil.dat[,11])
LOGIndtest_C.N_0.5 <- LOGIndtest(Indexdat,soil.dat[11])
lapply(Indtest_C.N_0.5, function(x) anova(x))


Indtest_Nitrif_0.5 <- Indtest(Indexdat,soil.dat[,28])
SQRTIndtest_Nitrif_0.5 <- SQRTIndtest(Indexdat,soil.dat[,28])
LOGIndtest_Nitrif_0.5 <- LOGIndtest(Indexdat,soil.dat[,28])
lapply(Indtest_Nitrif_0.5, function(x) anova(x))
x11();par(mfrow=c(3,4));lapply(Indtest_Nitrif_0.5, function(x) plot(x,which=c(1)));par(mfrow=c(1,1))
x11();par(mfrow=c(3,4));lapply(Indtest_Nitrif_0.5, function(x) plot(x,which=c(2)));par(mfrow=c(1,1))
x11();par(mfrow=c(3,4));lapply(SQRTIndtest_Nitrif_0.5, function(x) plot(x,which=c(1)));par(mfrow=c(1,1))
x11();par(mfrow=c(3,4));lapply(SQRTIndtest_Nitrif_0.5, function(x) plot(x,which=c(2)));par(mfrow=c(1,1))
x11();par(mfrow=c(3,4));lapply(LOGIndtest_Nitrif_0.5, function(x) plot(x,which=c(1)));par(mfrow=c(1,1))
x11();par(mfrow=c(3,4));lapply(LOGIndtest_Nitrif_0.5, function(x) plot(x,which=c(2)));par(mfrow=c(1,1))

  # Peut-être que transformer également les variables sol aiderait...

## Digression ##

# La variable causale et la variable réponse change en fonction de ce qu'on considère 
  # Les variables d'état sont plutôt causales, alors que les variables d'activité sont plutôt réponse...
  # Ca ne change pas forcément grand chose dans le cas d'une régression linéaire simple, mais peut-être
  # qu'il serait au final plus simple de faire simplement un test de corrélation de spearman entre les variables ?!

# Test de corrélation. Les variables x et y peuvent être adaptées dans le graph en fonction des variables cause/réponse
  # Les graphs assument des relations linéaires; on peut enlever "method="lm"" pour obtenir un courbe smooth ou inclure un autre modèle
CORIndtest <- function(df1,df2){    
  df2 <- melt(df2)                
  Indml  <- plotsig <- pval_l <- list()  
  for (indvar in Indvarnames) {
    Indm <-cor.test(df2$value,df1[,indvar], method = "spearman")     
    pval_l[[indvar]] <- round(Indm$p.value,3) #  
    Indml[[indvar]] <- Indm
    if(Indm$p.value < 0.055){        
      Indmsig <- Indm                       
      df3 <- melt(data.frame(df1[,c(1:5)],indvar = df1[,indvar])) # reformatage nécessaire pour ggplot2
      plotsig[[indvar]] <- ggplot(data.frame(df2,value2=df3[,length(df3)]), aes(y= value, x= value2)) + geom_point() + scale_x_continuous(indvar) + stat_smooth(method="lm")
      #Indml[[indvar]] <- Indmsig # si on veut uniquement retourner les modèles significatifs --> enlever "Indml[[indvar]] <- Indm"
    }  
  }
  if(length(plotsig)!=0){ multiplot(plotlist = plotsig,cols = 2 )} # représentation des relations sign si il y en a
  print(unlist(pval_l)) # impression des p-values
  return(Indml) # 
}

CORIndtest_Nitrif_0.5 <- CORIndtest(Indexdat,soil.dat[,28])
# les résultats restent similaires
CORIndtest_Mineraliz_0.5 <- CORIndtest(Indexdat,soil.dat[,30]) # pas sign
CORIndtest_Organisation_0.5 <- CORIndtest(Indexdat,soil.dat[,32])
  # nombreuses relations significatives, mais pas franchement évidentes sur les graphs...
CORIndtest_Organisation_0.5_T0 <- CORIndtest(Indexdat[Indexdat$Date=="Year0",],soil.dat[soil.dat$Date=="Year0",32])
CORIndtest_Organisation_0.5_T2 <- CORIndtest(Indexdat[Indexdat$Date=="Year2",],soil.dat[soil.dat$Date=="Year2",32])
CORIndtest_Organisation_0.5_T4 <- CORIndtest(Indexdat[Indexdat$Date=="Year4",],soil.dat[soil.dat$Date=="Year4",32])
# Intéressant: pas de relations en 2010 et des relations en 2012 et 2014
# Ne pas oublier que les indices sont corrélés entre eux!
summary(modOrg <-lm(soil.dat[,32]~Indexdat[,"Gen"]*Indexdat[,"Date"]))    
anova(modOrg)
par(mfrow=c(2,2)); plot(modOrg); par(mfrow=c(1,1)) # bof...
ggplot(modOrg$model,aes(x=modOrg$model[,2],y=modOrg$model[,1], col=modOrg$model[,3]))+geom_point()+geom_smooth(method="lm")
  # Apparemment surtout en 2012, avec des valeurs d'organisations bien supérieures (?!)

CORIndtest_Respsol_0.5 <- CORIndtest(Indexdat,soil.dat[,35]) # de nouveau très significatif globalement
CORIndtest_Respsol_0.5_T0 <- CORIndtest(Indexdat[Indexdat$Date=="Year0",],soil.dat[soil.dat$Date=="Year0",35])
CORIndtest_Respsol_0.5_T2 <- CORIndtest(Indexdat[Indexdat$Date=="Year2",],soil.dat[soil.dat$Date=="Year2",35])
CORIndtest_Respsol_0.5_T4 <- CORIndtest(Indexdat[Indexdat$Date=="Year4",],soil.dat[soil.dat$Date=="Year4",35])
summary(modResp <-lm(soil.dat[,35]~Indexdat[,"Con"]*Indexdat[,"Date"]))    
ggplot(modResp$model,aes(x=modResp$model[,2],y=modResp$model[,1], col=modResp$model[,3]))+geom_point()+geom_smooth(method="lm")

CORIndtest_Phosphatase_0.5 <- CORIndtest(Indexdat,soil.dat[,50]) # 
CORIndtest_Phosphatase_0.5_T0 <- CORIndtest(Indexdat[Indexdat$Date=="Year0",],soil.dat[soil.dat$Date=="Year0",50])
CORIndtest_Phosphatase_0.5_T2 <- CORIndtest(Indexdat[Indexdat$Date=="Year2",],soil.dat[soil.dat$Date=="Year2",50])
CORIndtest_Phosphatase_0.5_T4 <- CORIndtest(Indexdat[Indexdat$Date=="Year4",],soil.dat[soil.dat$Date=="Year4",50])

chart.Correlation(soil.dat[,c(44,47,50,53)], histogram=TRUE, pch=19)
  # Activités enzymatiques corrélées

###################################################################

## Analyse de la distribution des niveaux trophiques ##

##################################################################

TLdf <- read.table("data_niveaux trophiques.txt",header=TRUE, sep='\t', dec=",")
TLdf <- Subset(TLdf,community!="T4A16")
# TLdf2 <- Subset(TLdf,Mod!="RR-PER")

# Histogrammes: hist() et ggplot ne donnent pas tout-à-fait la même chose
  # sans doute une histoire de "bindwidth" mais je n'ai pas réussi à les rendre identiques (?!)
x11(); par(mfrow=c(4,4))
for(i in levels(TLdf$Date)){
  for(j in levels(TLdf$Mod)){
    with(Subset(TLdf,Date==i&Mod==j),hist(TrL,main=paste(i,j)))
  }
}

ggplot(TLdf, aes(TrL, colour = Mod,fill=Mod)) +
  geom_histogram(binwidth =0.5,aes(y = ..density..)) + xlim(1, 4) +
  facet_wrap(Date~Mod,ncol=4)
ggplot(TLdf, aes(TrL, colour = Mod,fill=Mod)) +
  geom_histogram(binwidth =0.5)+ xlim(1, 4) +
  facet_wrap(Date~Mod,ncol=4)

# en transparence... faire varier "adjust"
range(TLdf$TrL)
ggplot(TLdf, aes(TrL, fill = Mod)) +
  geom_density(alpha = 0.2,adjust=0.8) + xlim(1, 4) +
  facet_wrap(~ Date)+
  theme_bw()
# dur d'interpréter...légér décalage vers les niveaux trophiques supérieurs dans les systèmes sans labour

# Autre représentation: densité cumulée

ggplot(TLdf, aes(TrL, colour = Mod)) +  stat_ecdf(geom="point") + facet_wrap(~Date)
ggplot(TLdf, aes(TrL, colour = Mod)) +  stat_ecdf(geom="step") + facet_wrap(~Date)

  # On voit une différence entre traitements mais on voit surtout que la grande majorité
    # des noeuds ont un niveau trophique autour de 2
      # Pourrait-t-on introduire une meilleure différentiation de la ressource de base (e.g., différentes sources de carbone)?
        # Ou faut-il revoir la description des préférences trophiques?










