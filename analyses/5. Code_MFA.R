
## Spécifier le directoire et charger les librairies ##
setwd("F:/5. RH/4. Post-doctorant/2014 SOFIA V Coudrain/dossier juin 2016/Scripts finaux")

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

## MFA pour mettre en lien les différentes variables ##
# Study the similarities between sites with respect to all the variables
  # and the linear relationships between variables
require(FactoMineR)
########################################################

# D'abord charger les données et faire un df des variables à utiliser #
dfBersierInd <- read.table("data_BersierInd.txt",header=TRUE, sep='\t', dec=",") # 
dfCohenInd <- read.table("data_CohenInd.txt",header=TRUE, sep='\t', dec=",") # 
dfBersierInd_wide <- dcast(dfBersierInd, community + Date+ Plot+ Block+Mod ~ BersierInd, value.var="Weighted")
dfCohenInd_wide <- dcast(dfCohenInd, community + Date+ Block+ Plot+Mod ~ CohenInd, value.var="value")
basedat <- read.table("data_basefile.txt",header=TRUE, sep='\t', dec=",")
basedat <- Subset(basedat,Date!="Year3")
# Changer le ratio C:N en N:C 
basedat$N.C_0.5 <- 1/(basedat$C.N_0.5)
IndicesET <- read.table("data_indicesET.txt",sep="\t",header=T, dec=",")
# Changer le ratio B:F en F:B
IndicesET$FBpath <- 1/(IndicesET$BFpath)
# Définir les groupes de variables à spécifier dans la MFA
# L'idée est qu'au sein d'un groupe de variables, chacune est représentante de l'idée générale du groupe
# Ne pas mettre ensemble des variables conflictuelles dans leur interprétation
  # 1) Variables d'état: basedat: "N.C_0.5", "C_soluble_0.5","C_total_0.5","N_total_0.5"
  # 2) Variables d'activité biologique: "quotient_resp_0.5","Arylsulfatase_0.5","Glucosidase_0.5","Phosphatase_0.5","Urease_0.5"
  # 3) Variables de fonction: "Nitrif_0.5","Mineraliz_0.5","Respiration_sol_0.5","Organisation_0.5"
  # 4) Variables de "communauté": dfBersierInd_wide: "LD","Con","FrTop"; dfCohenInd_wide: "S","L_cann" ; IndicesET:  "TLmean", "TLmax", "FBpath","NbL","DivL","EvenL" 

mfadf <- cbind(basedat[,c(1:5,92,23,38,41,14,44,47,50,53,28,30,35,32)],dfBersierInd_wide[,c(16,6,14)],dfCohenInd_wide[,c(29,10)],IndicesET[,c(6,7,9:12)])
str(mfadf)
## Replace the NA in nitrif by the mean of the corresponding modality
which(is.na(mfadf[,"Nitrif_0.5"])==T) # 59
mfadf[59,"Nitrif_0.5"] <- mean(Subset(mfadf,Date=="Year4"&Mod=="RT-RR")[1:3,"Nitrif_0.5"])
which(is.na(mfadf[,"Nitrif_0.5"])==T) # NULL
mfadf$Mod <- factor(mfadf$Mod, levels=c("CONV","RN","RT","RT-RR","RR-PER"))
names(mfadf)
# centrer réduire les données
mfadf_sc <- cbind(mfadf[,1:5],scale(mfadf[,-c(1:5)]))

library(PerformanceAnalytics)
chart.Correlation(mfadf_sc[,-c(1:5)], histogram=TRUE, pch=19)
# Uréase  pourrait être problématique
hist(mfadf_sc$Urease_0.5) # un peu bizarre...éventuellement enlever cette variable
par(mfrow=c(1,3)); hist(Subset(mfadf_sc,Date=="Year0")$Urease_0.5);hist(Subset(mfadf_sc,Date=="Year2")$Urease_0.5);hist(Subset(mfadf_sc,Date=="Year4")$Urease_0.5);par(mfrow=c(1,1))

#mfadf <- mfadf[,-c(14)] # enlever uréase 
#mfadf_sc <- mfadf_sc[,-c(14)] # enlever uréase 
names(mfadf_sc)
statevar <- mfadf[,6:9]
activvar <- mfadf[,10:14]
funcvar <- mfadf[,15:18]
commvar <- mfadf[,19:29]
statevar_sc <- scale(statevar)
activvar_sc <- scale(activvar)
funcvar_sc <- scale(funcvar)
commvar_sc <- scale(commvar)

coeffRV(statevar_sc, activvar_sc)$p.value # <0.01
coeffRV(statevar_sc, funcvar_sc)$p.value # <0.01
coeffRV(statevar_sc, commvar_sc)$p.value # <0.01
coeffRV(activvar_sc, funcvar_sc)$p.value # <0.01
coeffRV(activvar_sc, commvar_sc)$p.value # <0.01
coeffRV(funcvar_sc, commvar_sc)$p.value # <0.01

mfagroup <- c(ncol(statevar_sc), ncol(activvar_sc), ncol(funcvar_sc), ncol(commvar_sc))
allmfa <- MFA(mfadf_sc[,-c(1:5)], group=mfagroup , ncp=3, type=c("s","s","s","s"),
              name.group=c("statevar", "activvar", "funcvar", "commvar"))
barplot(allmfa$eig[,1]) # 33% explained by the first axis, little differences between the others
allmfa$group$Lg
allmfa$group$RV # standadized Lg
  # not much difference in contribution...
allmfa$group$contrib

allmfa$ind$within.inertia #  a site is all the more "homogeneous" that its superposed representations are close
  ## May allow to infer whereas sites from one modality are more homogeneous than other 
    ## --> groups can be predicted from each others...

palette(c("black","red","blue","green"))
plot(allmfa,choix="group",palette=palette())
x11();plot(allmfa,choix="var",invisible="none",hab="group",palette=palette())
x11();plot(allmfa,choix="ind",partial="all",habillage="group",palette=palette())
x11();plot(allmfa,choix="axes",habillage="group",palette=palette())

dimdesc(allmfa, axes = 1:2) # point out the variables and the categories that are the most characteristic according to each dimension 
allmfa$global.pca

df_dim12 <- cbind(basedat[,c(1:5)],allmfa$global.pca$ind$coord[,1:2])
s.class(allmfa$global.pca$ind$coord[,1:2], df_dim12$Mod, col = rainbow(7))
meandf_dim12 <- with(df_dim12,(aggregate(df_dim12[,6:7],by=list(Mod=Mod,Date=Date),mean))) # moyenne des scores par traitement
# Création des polygones convexes pour chaque année
grp.10.mfa <- df_dim12[df_dim12$Date == "Year0", ][chull(df_dim12[df_dim12$Date == "Year0", c("Dim.1", "Dim.2")]), ]  # hull values for Year0
grp.12.mfa  <- df_dim12[df_dim12$Date == "Year2", ][chull(df_dim12[df_dim12$Date == "Year2", c("Dim.1", "Dim.2")]), ]  # hull values for Year2
grp.14.mfa  <- df_dim12[df_dim12$Date == "Year4", ][chull(df_dim12[df_dim12$Date == "Year4", c("Dim.1", "Dim.2")]), ]  # hull values for Year4
hull.data.mfa <- rbind(grp.10.mfa, grp.12.mfa,grp.14.mfa)  #

x11(9,3)
# Représentation des polygones
p1 <- ggplot() + 
  geom_polygon(data=hull.data.mfa,aes(x=Dim.1,y=Dim.2,colour=Date,linetype=Date),alpha=0,size=1.2) + # ,colour="black"
  #guides(fill=F, colour=F, linetype=F) + # enlever la légende
  scale_linetype_discrete("",labels=c("T0","T2","T4"))  +
  scale_color_manual("",labels=c("T0","T2","T4"),values = c("green","blue","red"))+
  theme_classic(base_size=15)
# Ajout des scores
p2 <- p1 + geom_point(data=meandf_dim12,aes(x=Dim.1,y=Dim.2,label=Mod,shape=Mod,col=Date),size=3,alpha=1) +  
  geom_point(data=df_dim12,aes(x=Dim.1,y=Dim.2,shape=Mod),size=2,alpha=0.3) + # ajout des points moyens
  scale_shape_manual("",labels=c("CONV","RN","RT","RT-RR","RR-PER"),values = c(16,17,15,3,21))+
  coord_equal() 
p2


## A REFAIRE POUR CHAQUE ANNEE ##
mfadf_sc.T0 <- Subset(mfadf_sc,Date=="Year0")
mfadf_sc.T2 <- Subset(mfadf_sc,Date=="Year2")
mfadf_sc.T4 <- Subset(mfadf_sc,Date=="Year4")
chart.Correlation(mfadf_sc.T0[,-c(1:5)], histogram=TRUE, pch=19) # Respiration_sol (C minéralisation) est suspecte
chart.Correlation(mfadf_sc.T2[,-c(1:5)], histogram=TRUE, pch=19)
chart.Correlation(mfadf_sc.T4[,-c(1:5)], histogram=TRUE, pch=19)
  # il semble il y avoir quelques outliers dans certaines variables --> éventuellement à voir de plus près

mfa.T0 <- MFA(mfadf_sc.T0[,-c(1:5)], group=mfagroup , ncp=3, type=c("s","s","s","s"),
              name.group=c("statevar", "activvar", "funcvar", "commvar"))
barplot(mfa.T0$eig[,1])
mfa.T0$group$Lg
mfa.T0$group$RV # standadized Lg
mfa.T0$group$contrib

mfa.T0$ind$within.inertia #  

palette(c("black","red","blue","green"))
plot(mfa.T0,choix="group",palette=palette())
x11();plot(mfa.T0,choix="var",invisible="none",hab="group",palette=palette())
x11();plot(mfa.T0,choix="ind",partial="all",habillage="group",palette=palette())
x11();plot(mfa.T0,choix="axes",habillage="group",palette=palette())

dimdesc(mfa.T0, axes = 1:2) # point out the variables and the categories that are the most characteristic according to each dimension 

mfa.T0$global.pca
s.class(mfa.T0$global.pca$ind$coord[,1:2], mfadf_sc.T0$Mod, col = rainbow(7))

df_dimT0 <- cbind(mfadf_sc.T0[,c(1:5)],mfa.T0$global.pca$ind$coord[,1:2])
meandf_dimT0 <- with(df_dimT0,(aggregate(df_dimT0[,6:7],by=list(Mod=Mod),mean))) # moyenne des scores par traitement
# Création des polygones convexes pour chaque pratique
grp.mfaT0_CONV <- df_dimT0[df_dimT0$Mod == "CONV", ][chull(df_dimT0[df_dimT0$Mod == "CONV", c("Dim.1", "Dim.2")]), ]  # hull values for CONV
grp.mfaT0_RN <- df_dimT0[df_dimT0$Mod == "RN", ][chull(df_dimT0[df_dimT0$Mod == "RN", c("Dim.1", "Dim.2")]), ]  # hull values for RN
grp.mfaT0_RT <- df_dimT0[df_dimT0$Mod == "RT", ][chull(df_dimT0[df_dimT0$Mod == "RT", c("Dim.1", "Dim.2")]), ]  # hull values for RT
grp.mfaT0_RTRR <- df_dimT0[df_dimT0$Mod == "RT-RR", ][chull(df_dimT0[df_dimT0$Mod == "RT-RR", c("Dim.1", "Dim.2")]), ]  # hull values for RT-RR
grp.mfaT0_RRPER <- df_dimT0[df_dimT0$Mod == "RR-PER", ][chull(df_dimT0[df_dimT0$Mod == "RR-PER", c("Dim.1", "Dim.2")]), ]  # hull values for RR-PER
hull.data.mfaT0 <- rbind(grp.mfaT0_CONV, grp.mfaT0_RN,grp.mfaT0_RT,grp.mfaT0_RTRR,grp.mfaT0_RRPER)  #

p1T0 <- ggplot() + 
  geom_polygon(data=hull.data.mfaT0,aes(x=Dim.1,y=Dim.2,colour=Mod,linetype=Mod),alpha=0,size=1.2) + 
  scale_linetype_discrete("",labels=c("CONV","RN","RT","RT-RR","RR-PER"))  +
  scale_color_manual("",labels=c("CONV","RN","RT","RT-RR","RR-PER"),values = c("green","blue","red","orange","orchid"))+
  theme_classic(base_size=15)
# Ajout des scores
p2T0 <- p1T0 + geom_point(data=df_dimT0,aes(x=Dim.1,y=Dim.2,label=Mod,shape=Mod,colour=Mod),size=2,alpha=0.3) +  
  geom_point(data=meandf_dimT0,aes(x=Dim.1,y=Dim.2,shape=Mod,colour=Mod),size=3,alpha=1) + # ajout des points moyens
  scale_shape_manual("",labels=c("CONV","RN","RT","RT-RR","RR-PER"),values = c(16,17,15,3,21))+
 # scale_color_manual("",labels=c("CONV","RN","RT","RT-RR","RR-PER"),values = c("green","blue","red","orange","orchid"))+
  coord_equal() 
p2T0


#####

mfa.T2 <- MFA(mfadf_sc.T2[,-c(1:5)], group=mfagroup , ncp=3, type=c("s","s","s","s"),
              name.group=c("statevar", "activvar", "funcvar", "commvar"))
barplot(mfa.T2$eig[,1])
mfa.T2$group$Lg
mfa.T2$group$RV # standadized Lg
mfa.T2$group$contrib

mfa.T2$ind$within.inertia #  

palette(c("black","red","blue","green"))
plot(mfa.T2,choix="group",palette=palette())
x11();plot(mfa.T2,choix="var",invisible="none",hab="group",palette=palette())
x11();plot(mfa.T2,choix="ind",partial="all",habillage="group",palette=palette())
x11();plot(mfa.T2,choix="axes",habillage="group",palette=palette())

dimdesc(mfa.T2, axes = 1:2) # point out the variables and the categories that are the most characteristic according to each dimension 

s.class(mfa.T2$global.pca$ind$coord[,1:2], mfadf_sc.T2$Mod, col = rainbow(7))

df_dimT2 <- cbind(mfadf_sc.T2[,c(1:5)],mfa.T2$global.pca$ind$coord[,1:2])
meandf_dimT2 <- with(df_dimT2,(aggregate(df_dimT2[,6:7],by=list(Mod=Mod),mean))) # moyenne des scores par traitement
# Création des polygones convexes pour chaque pratique
grp.mfaT2_CONV <- df_dimT2[df_dimT2$Mod == "CONV", ][chull(df_dimT2[df_dimT2$Mod == "CONV", c("Dim.1", "Dim.2")]), ]  # hull values for CONV
grp.mfaT2_RN <- df_dimT2[df_dimT2$Mod == "RN", ][chull(df_dimT2[df_dimT2$Mod == "RN", c("Dim.1", "Dim.2")]), ]  # hull values for RN
grp.mfaT2_RT <- df_dimT2[df_dimT2$Mod == "RT", ][chull(df_dimT2[df_dimT2$Mod == "RT", c("Dim.1", "Dim.2")]), ]  # hull values for RT
grp.mfaT2_RTRR <- df_dimT2[df_dimT2$Mod == "RT-RR", ][chull(df_dimT2[df_dimT2$Mod == "RT-RR", c("Dim.1", "Dim.2")]), ]  # hull values for RT-RR
grp.mfaT2_RRPER <- df_dimT2[df_dimT2$Mod == "RR-PER", ][chull(df_dimT2[df_dimT2$Mod == "RR-PER", c("Dim.1", "Dim.2")]), ]  # hull values for RR-PER
hull.data.mfaT2 <- rbind(grp.mfaT2_CONV, grp.mfaT2_RN,grp.mfaT2_RT,grp.mfaT2_RTRR,grp.mfaT2_RRPER)  #

p1T2 <- ggplot() + 
  geom_polygon(data=hull.data.mfaT2,aes(x=Dim.1,y=Dim.2,colour=Mod,linetype=Mod),alpha=0,size=1.2) + 
  scale_linetype_discrete("",labels=c("CONV","RN","RT","RT-RR","RR-PER"))  +
  scale_color_manual("",labels=c("CONV","RN","RT","RT-RR","RR-PER"),values = c("green","blue","red","orange","orchid"))+
  theme_classic(base_size=15)
# Ajout des scores
p2T2 <- p1T2 + geom_point(data=df_dimT2,aes(x=Dim.1,y=Dim.2,label=Mod,shape=Mod,colour=Mod),size=2,alpha=0.3) +  
  geom_point(data=meandf_dimT2,aes(x=Dim.1,y=Dim.2,shape=Mod,colour=Mod),size=3,alpha=1) + # ajout des points moyens
  scale_shape_manual("",labels=c("CONV","RN","RT","RT-RR","RR-PER"),values = c(16,17,15,3,21))+
  # scale_color_manual("",labels=c("CONV","RN","RT","RT-RR","RR-PER"),values = c("green","blue","red","orange","orchid"))+
  coord_equal() 
p2T2

###

mfa.T4 <- MFA(mfadf_sc.T4[,-c(1:5)], group=mfagroup , ncp=3, type=c("s","s","s","s"),
              name.group=c("statevar", "activvar", "funcvar", "commvar"))
barplot(mfa.T4$eig[,1])
mfa.T4$group$Lg
mfa.T4$group$RV # standadized Lg
mfa.T4$group$contrib

mfa.T4$ind$within.inertia #  

palette(c("black","red","blue","green"))
plot(mfa.T4,choix="group",palette=palette())
x11();plot(mfa.T4,choix="var",invisible="none",hab="group",palette=palette())
x11();plot(mfa.T4,choix="ind",partial="all",habillage="group",palette=palette())
x11();plot(mfa.T4,choix="axes",habillage="group",palette=palette())

dimdesc(mfa.T4, axes = 1:2) # point out the variables and the categories that are the most characteristic according to each dimension 

mfa.T4$global.pca
s.class(mfa.T4$global.pca$ind$coord[,1:2], mfadf_sc.T4$Mod, col = rainbow(7))

df_dimT4 <- cbind(mfadf_sc.T4[,c(1:5)],mfa.T4$global.pca$ind$coord[,1:2])
meandf_dimT4 <- with(df_dimT4,(aggregate(df_dimT4[,6:7],by=list(Mod=Mod),mean))) # moyenne des scores par traitement
# Création des polygones convexes pour chaque pratique
grp.mfaT4_CONV <- df_dimT4[df_dimT4$Mod == "CONV", ][chull(df_dimT4[df_dimT4$Mod == "CONV", c("Dim.1", "Dim.2")]), ]  # hull values for CONV
grp.mfaT4_RN <- df_dimT4[df_dimT4$Mod == "RN", ][chull(df_dimT4[df_dimT4$Mod == "RN", c("Dim.1", "Dim.2")]), ]  # hull values for RN
grp.mfaT4_RT <- df_dimT4[df_dimT4$Mod == "RT", ][chull(df_dimT4[df_dimT4$Mod == "RT", c("Dim.1", "Dim.2")]), ]  # hull values for RT
grp.mfaT4_RTRR <- df_dimT4[df_dimT4$Mod == "RT-RR", ][chull(df_dimT4[df_dimT4$Mod == "RT-RR", c("Dim.1", "Dim.2")]), ]  # hull values for RT-RR
grp.mfaT4_RRPER <- df_dimT4[df_dimT4$Mod == "RR-PER", ][chull(df_dimT4[df_dimT4$Mod == "RR-PER", c("Dim.1", "Dim.2")]), ]  # hull values for RR-PER
hull.data.mfaT4 <- rbind(grp.mfaT4_CONV, grp.mfaT4_RN,grp.mfaT4_RT,grp.mfaT4_RTRR,grp.mfaT4_RRPER)  #

p1T4 <- ggplot() + 
  geom_polygon(data=hull.data.mfaT4,aes(x=Dim.1,y=Dim.2,colour=Mod,linetype=Mod),alpha=0,size=1.2) + 
  scale_linetype_discrete("",labels=c("CONV","RN","RT","RT-RR","RR-PER"))  +
  scale_color_manual("",labels=c("CONV","RN","RT","RT-RR","RR-PER"),values = c("green","blue","red","orange","orchid"))+
  theme_classic(base_size=15)
# Ajout des scores
p2T4 <- p1T4 + geom_point(data=df_dimT4,aes(x=Dim.1,y=Dim.2,label=Mod,shape=Mod,colour=Mod),size=2,alpha=0.3) +  
  geom_point(data=meandf_dimT4,aes(x=Dim.1,y=Dim.2,shape=Mod,colour=Mod),size=3,alpha=1) + # ajout des points moyens
  scale_shape_manual("",labels=c("CONV","RN","RT","RT-RR","RR-PER"),values = c(16,17,15,3,21))+
  # scale_color_manual("",labels=c("CONV","RN","RT","RT-RR","RR-PER"),values = c("green","blue","red","orange","orchid"))+
  coord_equal() 
p2T4


# En résumé #

g1 <- s.class(mfa.T0$global.pca$ind$coord[,1:2], mfadf_sc.T0$Mod, col = rainbow(7))
g2 <- s.class(mfa.T2$global.pca$ind$coord[,1:2], mfadf_sc.T2$Mod, col = rainbow(7))
g3 <- s.class(mfa.T4$global.pca$ind$coord[,1:2], mfadf_sc.T4$Mod, col = rainbow(7))
ADEgS(adeglist = list(g1,g2,g3))














