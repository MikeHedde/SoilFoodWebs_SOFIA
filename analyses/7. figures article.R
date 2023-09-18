# -----------------------------------------------------------------------------
# PROJECT:
#    Linking soil foodweb and soil functions dynamics under
#    changing agricultural practices managment
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# STEP:
#    1.0   Figures of the article
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# NOTES:
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# LIBRARIES:
library(ade4)
library(ggplot2)
library(Rmisc)
library(plyr)
library(cowplot)
library(FSA)
library(lme4)
library(nlme)
library(car)
library(dunn.test)
library(xlsx)
# -----------------------------------------------------------------------------

# directory
setwd("C:/1. Mike/5. RH/4. Post-doctorant/2014 SOFIA V Coudrain/dossier juin 2016/Scripts finaux")


df_totdiv <- read.xlsx("df_totdiv.xlsx", h = T, sheetIndex = 1, colIndex = c(2:11))


# -----------------------------------------------------------------------------
#Function
# -----------------------------------------------------------------------------
myfun <- function(mydata, measvar, lab_x, lab_y){
  pd <- position_dodge(0.1) # move them .05 to the left and right
  dt <- summarySE(mydata, measurevar="value", groupvars = c("Date", "Mod"))
  gdmean <- cbind.data.frame(Gmean = rep(mean(mydata$value),4), 
                             Gse = rep(sd(mydata$value)/sqrt(20),4), Date = rep("Year0", 4), Mod = c("CONV", "RN", "RT", "RT-RR"))
  ifelse(min(dt$value-dt$se)<0, MIN <- min(dt$value-dt$se), MIN <- 0)
  g <- ggplot(dt, aes(x = paste(Mod, Date), fill = Mod, alpha = Date))+
    geom_bar(aes(weight = value), position = pd, width = .5)+
    geom_errorbar(aes(ymin = value-se, ymax = value+se), width = .1, position = pd)+
    labs(x = "Years", y = lab_y)+
    scale_x_discrete(labels = rep(c(0,2,4), 4))+
    scale_alpha_manual(values = c(.3, .6, .9))+
    scale_fill_manual(values = c("steelblue4", "darkorchid4", "forestgreen", "olivedrab4"))+
    geom_vline(xintercept = c(3.5, 6.5, 9.5), linetype = 2)+
    geom_hline(yintercept = 0)+
    geom_text(aes(label = "CONV", x = 2, y = max(dt$value+dt$se)*1.2))+
    geom_text(aes(label = "RN", x = 5, y = max(dt$value+dt$se)*1.2))+
    geom_text(aes(label = "RT", x = 8, y = max(dt$value+dt$se)*1.2))+
    geom_text(aes(label = "RT-RR", x = 11, y = max(dt$value+dt$se)*1.2))+
    ylim(MIN*1.2, max(dt$value+dt$se)*1.2)+
    guides(alpha = F, fill = F)+
    theme (axis.title=element_text(size=10,face="bold"),
           axis.text=element_text(size=8))
  ls <- list(g, dt)
  ls
}

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Dataframe call
# -----------------------------------------------------------------------------
df_totdiv <- cbind(df_totdiv, indiceReseau[,-(1:4)])  
#indiceReseau vient du script "MH courbes cumulées TrL biom.R"
df_totdiv <- df_totdiv[df_totdiv$Mod != "RR-PER", ]
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# SCRIPT
# -----------------------------------------------------------------------------

soil.dat <- read.table("data_basefile.txt",header=TRUE, sep='\t', dec=",")
soil.dat <- Subset(soil.dat, Date!="Year3"& PlotID!="T4A16" & Mod != "RR-PER")

soil.fun5 <- data.frame(nitrif = soil.dat$Nitrif_0.5, 
                        Nmin = soil.dat$Mineraliz_0.5,
                        Nimm = soil.dat$Organisation_0.5, 
                        Cmin = soil.dat$Respiration_sol_0.5, 
                        glu = soil.dat$Glucosidase_0.5/soil.dat$C_total_0.5, 
                        ure = soil.dat$Urease_0.5/soil.dat$N_total_0.5)
soil.fun5$nitrif[47] <- mean(soil.fun5$nitrif[33:48], na.rm = T) # remplace une valeur manquante

soil.fun20 <- data.frame(nitrif = soil.dat$Nitrif_5.20, 
                        Nmin = soil.dat$Mineraliz_5.20,
                        Nimm = soil.dat$Organisation_5.20, 
                        Cmin = soil.dat$Respiration_sol_5.20, 
                        glu = soil.dat$Glucosidase_5.20/soil.dat$C_total_5.20, 
                        ure = soil.dat$Urease_5.20/soil.dat$N_total_5.20)

      #Figure 3: Carbon
      mdt <- cbind.data.frame(Date = soil.dat$Date, Mod = soil.dat$Mod, value = soil.fun5$Cmin)
      pCmin5 <- myfun(mydata = mdt, measvar = "value", lab_y = "C mineralisation") 

      mdt <- cbind.data.frame(Date = soil.dat$Date, Mod = soil.dat$Mod, value = soil.fun5$glu)
      pglu5 <- myfun(mydata = mdt, measvar = "value", lab_y = "Glucosidase\nactivity") 
      

      mdt <- cbind.data.frame(Date = soil.dat$Date, Mod = soil.dat$Mod, value = soil.fun5$nitrif)
      pNitrif5 <- myfun(mydata = mdt, measvar = "value", lab_y = "Nitrification") 
      
      
      mdt <- cbind.data.frame(Date = soil.dat$Date, Mod = soil.dat$Mod, value = soil.fun5$Nmin)
      pNmin5 <- myfun(mydata = mdt, measvar = "value", lab_y = "N\nmineralization") 

      
      mdt <- cbind.data.frame(Date = soil.dat$Date, Mod = soil.dat$Mod, value = soil.fun5$Nimm)
      pNimm5 <- myfun(mydata = mdt, measvar = "value", lab_y = "N immobilization") 
      

      mdt <- cbind.data.frame(Date = soil.dat$Date, Mod = soil.dat$Mod, value = soil.fun5$ure)
      pure5 <- myfun(mydata = mdt, measvar = "value", lab_y = "Urease\nactivity") 

      
      ggdraw()+
        draw_plot(pCmin5, x = 0, y = 0.66, width = 0.5, height = 0.33)+
        draw_plot(pglu5, x = 0.5, y = 0.66, width = 0.5, height = 0.33)+
        draw_plot(pNitrif5, x = 0, y = 0.33, width = 0.5, height = 0.33)+
        draw_plot(pNmin5, x = 0.5, y = 0.33, width = 0.5, height = 0.33)+
        draw_plot(pNimm5, x = 0, y = 0, width = 0.5, height = 0.33)+
        draw_plot(pure5, x = 0.5, y = 0, width = 0.5, height = 0.33)+
        draw_plot_label(c("A", "B", "C", "D", "E", "F"), x = c(0, 0.5, 0, 0.5, 0, 0.5), y = c(01, 1, .66, .66, 0.33, 0.33))
      

#--------------------------------------------

              #--------------------------------------------
              # analyses LMER
              #-------------------------------------------- 
              lmm.sol <- cbind.data.frame(soil.dat[,1:5], soil.fun5)
              
              Y = log(lmm.sol$ure)
              m1 <-lmer(Y~ Mod * Date  - 1 + (1|Block), data = lmm.sol, REML = T)
              
              #inspection des résidus
              E1 <- resid(m1)
              F1 <- fitted(m1)
              par(mfrow =c(3,2))
              hist(E1)  #normalité des résidus
              plot(E1~F1, xlab = "Fitted Values", ylab = "Normalized residuals") #indépendance des résidus
              abline(h = 0, lty = 2)
              plot(E1~lmm.sol$Mod, xlab = "Y variable",ylab = "Normalized residuals")
              abline(h = 0, lty = 2)
              boxplot(E1~lmm.sol$Date, xlab = "Y variable",ylab = "Normalized residuals")
              abline(h = 0, lty = 2)

              summary(m1)
              Anova(m1, type ="2")
              
              X = abbreviate(paste(lmm.sol$Mod, lmm.sol$Date))
              X = lmm.sol$Mod
              P = cbind.data.frame(dunn.test(Y, X)$comp,
                                   dunn.test(Y, X)$P.adjust)
              P <- P[P[,2]<0.05,]
              P
              #-------------------------------------------- 


#-------------------------------------------- 
for (i in c(7:ncol(df_totdiv))) df_totdiv[,i] <- as.numeric(as.character(df_totdiv[, i]))
              
# figure sur la redondance fonctionnelle
# redunMelted vient du script "Funct Redund.R"
mdt <- redunMelted[redunMelted$variable == "Total Network" & redunMelted$Mod != "RR-PER", ]
pFunRed <- myfun(mydata = mdt, measvar = "value", lab_y = "Functional\nredundancy")

# figure sur le rapport brown/green 
mdt <- cbind.data.frame(Date = df_totdiv$Date, Mod = df_totdiv$Mod, value = log(df_totdiv$BGpath))[df_totdiv$Mod != "RR-PER", ]
pBrwnGrn <- myfun(mydata = mdt, measvar = "value", lab_y =  "log Brown-to-\nGreen path ratio")

# figure sur le rapport bact/fung
mdt <- cbind.data.frame(Date = df_totdiv$Date, Mod = df_totdiv$Mod, value = log(df_totdiv$BFpath))[df_totdiv$Mod != "RR-PER", ]
pBactFung <- myfun(mydata = mdt, measvar = "value",  lab_y = "log Bacterial-to-\nFungal path ratio")

# figure sur la diversité des liens
mdt <- cbind.data.frame(Date = df_totdiv$Date, Mod = df_totdiv$Mod, value = df_totdiv$DivL)[df_totdiv$Mod != "RR-PER", ]
pLnkDiv <- myfun(mydata = mdt, measvar = "value",  lab_y = "Shannon's\nlinks diversity")

# figure sur le nb de liens
mdt <- cbind.data.frame(Date = df_totdiv$Date, Mod = df_totdiv$Mod, value = df_totdiv$NbL)[df_totdiv$Mod != "RR-PER", ]
pLnkNb <- myfun(mydata = mdt, measvar = "value", lab_y = "Number\nof links")

# figure sur l'AUC
mdt <- cbind.data.frame(Date = df_totdiv$Date, Mod = df_totdiv$Mod, value = df_totdiv$auc)[df_totdiv$Mod != "RR-PER", ]
pAUC <- myfun(mydata = mdt, measvar = "value", lab_y = "Area\nunder curve")

# figure sur la biomasse tot dans le réseau
mdt <- cbind.data.frame(Date = df_totdiv$Date, Mod = df_totdiv$Mod, value = df_totdiv$bmax/100000)[df_totdiv$Mod != "RR-PER", ]
pTotBiom <- myfun(mydata = mdt, measvar = "value", lab_y = "Total biomass\n(Mg per ha)")

# figure sur le niveau troph max
mdt <- cbind.data.frame(Date = df_totdiv$Date, Mod = df_totdiv$Mod, value = df_totdiv$TLmax)[df_totdiv$Mod != "RR-PER", ]
pmaxTL <- myfun(mydata = mdt, measvar = "value",  lab_y = "Max.\ntrophic level")

# figure sur le niveau troph moy
mdt <- cbind.data.frame(Date = df_totdiv$Date, Mod = df_totdiv$Mod, value = df_totdiv$TLmean)[df_totdiv$Mod != "RR-PER", ]
pmeanTL <- myfun(mydata = mdt, measvar = "value",  lab_y = "Mean\ntrophic level")

# figure sur le nombre de groupes trophiques
mdt <- cbind.data.frame(Date = df_totdiv$Date, Mod = df_totdiv$Mod, value = df_totdiv$nb_trophGroups)[df_totdiv$Mod != "RR-PER", ]
pnbTG <- myfun(mydata = mdt, measvar = "value",  lab_y = "Number of \ntrophic groups")

#X11()
plot_grid(pBrwnGrn, pBactFung, pmeanTL, pmaxTL, pTotBiom, pnbTG , pLnkNb, pAUC,   
          pLnkDiv, pFunRed, ncol = 2, nrow = 5, labels = "AUTO")
# -----------------------------------------------------------------------------
            #--------------------------------------------
            # analyses LMER
            #-------------------------------------------- 
            lmmodel <- cbind.data.frame(soil.dat[,1:5], foodweb, soil.fun5)
            Y = log(lmmodel$FunctRed)
            m1 <-lmer(Y~ Mod * Date - 1 + (1|Block), data = lmmodel, REML = T)
            
            #inspection des résidus
            E1 <- resid(m1)
            F1 <- fitted(m1)
            par(mfrow =c(3,2))
            hist(E1)  #normalité des résidus
            plot(E1~F1, xlab = "Fitted Values", ylab = "Normalized residuals") #indépendance des résidus
            abline(h = 0, lty = 2)
            plot(E1~lmmodel$Mod, xlab = "Y variable",ylab = "Normalized residuals")
            abline(h = 0, lty = 2)
            boxplot(E1~lmmodel$Date, xlab = "Y variable",ylab = "Normalized residuals")
            abline(h = 0, lty = 2)
            boxplot(E1~lmmodel$Block, xlab = "Y variable",ylab = "Normalized residuals")
            abline(h = 0, lty = 2)
            
            summary(m1)
            Anova(m1, type ="2")
            
            Y = log(lmmodel$Biom_tot)
            X = lmmodel$Mod
            P = cbind.data.frame(dunn.test(Y, X, method = "bonferroni")$comp,
                                 dunn.test(Y, X, method = "bonferroni")$P.adjust)
            P <- P[P[,2]<0.5,]
            P
            
            #-------------------------------------------- 
