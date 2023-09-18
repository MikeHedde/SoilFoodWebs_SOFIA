# -----------------------------------------------------------------------------
# PROJECT:
#    Linking soil foodweb and soil functions dynamics under
#    changing agricultural practices managment
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# STEP:
#    0.3   Comparison of TL-biomass curves for each year and each crop system
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# NOTES:
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# LIBRARIES:
library(zoo)
require(Bolstad2)
# -----------------------------------------------------------------------------

# directory
setwd("C:/Users/Mickael/Documents/5. RH/Post-doc/2014 SOFIA V Coudrain/dossier juin 2016/Scripts finaux")


# -----------------------------------------------------------------------------
# SCRIPT
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# FUNCTIONS
# -----------------------------------------------------------------------------
      myfun <- function(TL, BM){
      db <- data.frame(TrL = TL, bm = BM[BM>0])
      db <- db[order(db$TrL, decreasing = F),]
      db <- db[-c(1:2),]
      dbagg <- aggregate(db, list(db$TrL), sum)[, -2]
      TrL_prc <- dbagg[,1]/max(dbagg[,1])
      cumbm <- cumsum(dbagg$bm)
      db2 <- data.frame(TrL = dbagg[,1], cumbm = cumbm, cumprc = TrL_prc)
      db2
      }
      
      myfun_cumbm <- function(X){
        X_m <- aggregate(X$cumbm, list(X$TrL), mean)
        X_m$TrL<- as.numeric(as.character(X_m$Group.1))
        X_sd <- aggregate(X$cumbm, list(X$TrL), sd)
        X_m <- data.frame(TrL = X_m$TrL, mean = X_m$x, se = X_sd$x/sqrt(4))
        X_m
      }
      
      myfun_cumbmprc <- function(X){
        X_m <- aggregate(X$cumprc, list(X$TrL), mean)
        X_m$TrL<- as.numeric(as.character(X_m$Group.1))
        X_sd <- aggregate(X$cumprc, list(X$TrL), sd)
        X_m <- data.frame(TrL = X_m$TrL, mean = X_m$x, se = X_sd$x/sqrt(4))
        X_m
      }
# -----------------------------------------------------------------------------
 
     
# -----------------------------------------------------------------------------      
# liste des niveaux trophiques identifiés
# -----------------------------------------------------------------------------      
lev_TrL <- as.data.frame(levels(as.factor(unique(c(unlist(dfTLT4), unlist(dfTLT2), unlist(dfTLT0))))))
nb_trophGroups <- c(unlist(lapply(dfTLT0, nrow)), unlist(lapply(dfTLT2, nrow)), unlist(lapply(dfTLT4, nrow)))
      
colnames(lev_TrL) <- "TrL"

dbT0 <- list()
for (i in c(1:20)){
dbT0[[i]]        <- myfun(TL = dfTLT0[[i]]$TrL, BM = dataB.T0[[i]]$biomass)
dbT0[[i]]        <- merge(lev_TrL, dbT0[[i]], all.x = T)[-1,]
dbT0[[i]]$cumbm  <- na.locf(dbT0[[i]]$cumbm)
dbT0[[i]]$cumprc <- na.locf(dbT0[[i]]$cumprc)}

dbT2 <- list()
for (i in c(1:20)){
  dbT2[[i]] <-myfun(TL = dfTLT2[[i]]$TrL, BM = dataB.T2[[i]]$biomass)
  dbT2[[i]] <-merge(lev_TrL, dbT2[[i]], all.x = T)[-1,]
  dbT2[[i]]$cumbm <- na.locf(dbT2[[i]]$cumbm)
  dbT2[[i]]$cumprc <- na.locf(dbT2[[i]]$cumprc)}

dbT4 <- list()
for (i in c(1:20)){
  dbT4[[i]] <-myfun(TL = dfTLT4[[i]]$TrL, BM = dataB.T4[[i]]$biomass)
  dbT4[[i]] <-merge(lev_TrL, dbT4[[i]], all.x = T)[-1,]
  dbT4[[i]]$cumbm <- na.locf(dbT4[[i]]$cumbm)
  dbT4[[i]]$cumprc <- na.locf(dbT4[[i]]$cumprc)
  }
# -----------------------------------------------------------------------------     



# -----------------------------------------------------------------------------   
#Calcul des biomasses (brutes et pourcentage) cumulées 
# -----------------------------------------------------------------------------   
#CONV 
    convT0 <-rbind(dbT0[[3]], dbT0[[10]], dbT0[[11]], dbT0[[18]])
    convT2 <-rbind(dbT2[[3]], dbT2[[10]], dbT2[[11]], dbT2[[18]])
    convT4 <-rbind(dbT4[[3]], dbT4[[10]], dbT4[[11]], dbT4[[18]])
    
    #cumulated biomass
    convT0_m <- myfun_cumbm(convT0)
    convT2_m <- myfun_cumbm(convT2)
    convT4_m <- myfun_cumbm(convT4)

    #cumulated biomass percentage
    convT0_m_prc <- myfun_cumbmprc(convT0)
    convT2_m_prc <- myfun_cumbmprc(convT2)
    convT4_m_prc <- myfun_cumbmprc(convT4)

##RN
rnT0 <-rbind(dbT0[[1]], dbT0[[8]], dbT0[[12]], dbT0[[20]])
rnT2 <-rbind(dbT2[[1]], dbT2[[8]], dbT2[[12]], dbT2[[20]])
rnT4 <-rbind(dbT4[[1]], dbT4[[8]], dbT4[[12]], dbT4[[20]])

    #cumulated biomass
    rnT0_m <- myfun_cumbm(rnT0)
    rnT2_m <- myfun_cumbm(rnT2)
    rnT4_m <- myfun_cumbm(rnT4)
    
    #cumulated biomass percentage
    rnT0_m_prc <- myfun_cumbmprc(rnT0)
    rnT2_m_prc <- myfun_cumbmprc(rnT2)
    rnT4_m_prc <- myfun_cumbmprc(rnT4)
    
    
##RT
    rtT0 <-rbind(dbT0[[2]], dbT0[[9]], dbT0[[15]], dbT0[[16]])
    rtT2 <-rbind(dbT2[[2]], dbT2[[9]], dbT2[[15]], dbT2[[16]])
    rtT4 <-rbind(dbT4[[2]], dbT4[[9]], dbT4[[15]], dbT4[[16]])
    
    #cumulated biomass
    rtT0_m <- myfun_cumbm(rtT0)
    rtT2_m <- myfun_cumbm(rtT2)
    rtT4_m <- myfun_cumbm(rtT4)
    
    #cumulated biomass percentage
    rtT0_m_prc <- myfun_cumbmprc(rtT0)
    rtT2_m_prc <- myfun_cumbmprc(rtT2)
    rtT4_m_prc <- myfun_cumbmprc(rtT4)
    

##RT-RR
rtrrT0 <-rbind(dbT0[[5]], dbT0[[6]], dbT0[[14]], dbT0[[19]])
rtrrT2 <-rbind(dbT2[[5]], dbT2[[6]], dbT2[[14]], dbT2[[19]])
rtrrT4 <-rbind(dbT4[[5]], dbT4[[6]], dbT4[[14]], dbT4[[19]])

    #cumulated biomass
    rtrrT0_m <- myfun_cumbm(rtrrT0)
    rtrrT2_m <- myfun_cumbm(rtrrT2)
    rtrrT4_m <- myfun_cumbm(rtrrT4)
    
    #cumulated biomass percentage
    rtrrT0_m_prc <- myfun_cumbmprc(rtrrT0)
    rtrrT2_m_prc <- myfun_cumbmprc(rtrrT2)
    rtrrT4_m_prc <- myfun_cumbmprc(rtrrT4)
    


plot(rnT0_m$mean~rnT0_m$TrL, las =1, ylim = c(80, 100), type = "n",
     ylab = "Percentage of cumulated biomass", xlab = "Trophic level")
plot(rnT0_m$mean~rnT0_m$TrL, lwd = 3)
lines((rnT0_m$mean+rnT0_m$se)~rnT0_m$TrL, lty=2)
lines((rnT0_m$mean-rnT0_m$se)~rnT0_m$TrL, lty =2)
lines(rnT2_m$mean~rnT2_m$TrL, col = "blue", lwd = 3)
lines((rnT2_m$mean+rnT2_m$se)~rnT2_m$TrL, lty=2, col = "blue")
lines((rnT2_m$mean-rnT2_m$se)~rnT2_m$TrL, lty =2, col = "blue")
lines(rnT4_m$mean~rnT4_m$TrL, col = "red", lwd = 3)
lines((rnT4_m$mean+rnT4_m$se)~rnT4_m$TrL, lty=2, col = "red")
lines((rnT4_m$mean-rnT4_m$se)~rnT4_m$TrL, lty =2, col = "red")

plot(rtT2_m_prc$mean~rtT2_m_prc$TrL)
lines(rtrrT2_m_prc$mean~rtrrT2_m_prc$TrL)

plot(rtT0_m_prc$mean~rtT0_m_prc$TrL)
lines(rtrrT0_m_prc$mean~rtrrT0_m_prc$TrL)
lines(rnT0_m_prc$mean~rnT0_m_prc$TrL)
lines(convT0_m_prc$mean~convT0_m_prc$TrL)

plot(rtT4_m_prc$mean~rtT4_m_prc$TrL)
lines(rtrrT4_m_prc$mean~rtrrT4_m_prc$TrL)
# -----------------------------------------------------------------------------   


# -----------------------------------------------------------------------------   
# Comparaison des courbes cumulées de biomasse par année
# -----------------------------------------------------------------------------   
par(mfrow= c(1,3))

#T0
    plot(rnT0_m$mean~rnT0_m$TrL, las =1, type = "n", xlim = c(2, 3.5), ylim = c(6, 10),
         ylab = "Cumulated animal biomass (T ha-2)", xlab = "Trophic level", main = "T+0")
    polygon(x=c(rnT0_m$TrL, rev(rnT0_m$TrL)), y=c(rnT0_m$mean-rnT0_m$se, rev(rnT0_m$mean+rnT0_m$se))/100000, 
            col = rgb(1,0,0, 0.1), border=NA)
    lines(0.00001*rnT0_m$mean~rnT0_m$TrL, col = "red", lwd = 3)
    polygon(x=c(rtrrT0_m$TrL,rev(rtrrT0_m$TrL)), y=c(rtrrT0_m$mean-rtrrT0_m$se, rev(rtrrT0_m$mean+rtrrT0_m$se))/100000, 
            col = rgb(0,1,0, 0.1), border=NA)
    lines(0.00001*rtrrT0_m$mean~rtrrT0_m$TrL, col = "springgreen4", lwd = 3)
    polygon(x=c(convT0_m$TrL,rev(convT0_m$TrL)), y=c(convT0_m$mean-convT0_m$se, rev(convT0_m$mean+convT0_m$se))/100000, 
            col = rgb(0,0,1, 0.1), border=NA)
    lines(0.00001*convT0_m$mean~convT0_m$TrL, lwd = 3, col = "violet")
    polygon(x=c(rtT0_m$TrL,rev(rtT0_m$TrL)), y=c(rtT0_m$mean-rtT0_m$se, rev(rtT0_m$mean+rtT0_m$se))/100000, 
              col = rgb(0,0,0, 0.1), border=NA)
    lines(0.00001*rtT0_m$mean~rtT0_m$TrL, col = "grey", lwd = 3)
    legend("topleft", c("CONV", "RN", "RT", "RT-RR"), bty = "n", text.col = c("violet", "red", "grey", "springgreen4"))

#T2
    plot(rnT2_m$mean~rnT2_m$TrL, las =1, type = "n", xlim = c(2, 3.5), ylim = c(6, 10),
         ylab = "Cumulated animal biomass (T ha-2)", xlab = "Trophic level", main = "T+0")
    polygon(x=c(rnT2_m$TrL, rev(rnT2_m$TrL)), y=c(rnT2_m$mean-rnT2_m$se, rev(rnT2_m$mean+rnT2_m$se))/100000, 
            col = rgb(1,0,0, 0.1), border=NA)
    lines(0.00001*rnT2_m$mean~rnT2_m$TrL, col = "red", lwd = 3)
    polygon(x=c(rtrrT2_m$TrL,rev(rtrrT2_m$TrL)), y=c(rtrrT2_m$mean-rtrrT2_m$se, rev(rtrrT2_m$mean+rtrrT2_m$se))/100000, 
            col = rgb(0,1,0, 0.1), border=NA)
    lines(0.00001*rtrrT2_m$mean~rtrrT2_m$TrL, col = "springgreen4", lwd = 3)
    polygon(x=c(convT2_m$TrL,rev(convT2_m$TrL)), y=c(convT2_m$mean-convT2_m$se, rev(convT2_m$mean+convT2_m$se))/100000, 
            col = rgb(0,0,1, 0.1), border=NA)
    lines(0.00001*convT2_m$mean~convT2_m$TrL, lwd = 3, col = "violet")
    polygon(x=c(rtT2_m$TrL,rev(rtT2_m$TrL)), y=c(rtT2_m$mean-rtT2_m$se, rev(rtT2_m$mean+rtT2_m$se))/100000, 
            col = rgb(0,0,0, 0.1), border=NA)
    lines(0.00001*rtT2_m$mean~rtT2_m$TrL, col = "grey", lwd = 3)
    legend("topleft", c("CONV", "RN", "RT", "RT-RR"), bty = "n", text.col = c("violet", "red", "grey", "springgreen4"))

#T4
    plot(rnT4_m$mean~rnT4_m$TrL, las =1, type = "n", xlim = c(2, 3.5), ylim = c(6, 10),
         ylab = "Cumulated animal biomass (T ha-2)", xlab = "Trophic level", main = "T+0")
    polygon(x=c(rnT4_m$TrL, rev(rnT4_m$TrL)), y=c(rnT4_m$mean-rnT4_m$se, rev(rnT4_m$mean+rnT4_m$se))/100000, 
            col = rgb(1,0,0, 0.1), border=NA)
    lines(0.00001*rnT4_m$mean~rnT4_m$TrL, col = "red", lwd = 3)
    polygon(x=c(rtrrT4_m$TrL,rev(rtrrT4_m$TrL)), y=c(rtrrT4_m$mean-rtrrT4_m$se, rev(rtrrT4_m$mean+rtrrT4_m$se))/100000, 
            col = rgb(0,1,0, 0.1), border=NA)
    lines(0.00001*rtrrT4_m$mean~rtrrT4_m$TrL, col = "springgreen4", lwd = 3)
    polygon(x=c(convT4_m$TrL,rev(convT4_m$TrL)), y=c(convT4_m$mean-convT4_m$se, rev(convT4_m$mean+convT4_m$se))/100000, 
            col = rgb(0,0,1, 0.1), border=NA)
    lines(0.00001*convT4_m$mean~convT4_m$TrL, lwd = 3, col = "violet")
    polygon(x=c(rtT4_m$TrL,rev(rtT4_m$TrL)), y=c(rtT4_m$mean-rtT4_m$se, rev(rtT4_m$mean+rtT4_m$se))/100000, 
            col = rgb(0,0,0, 0.1), border=NA)
    lines(0.00001*rtT4_m$mean~rtT4_m$TrL, col = "grey", lwd = 3)
    legend("topleft", c("CONV", "RN", "RT", "RT-RR"), bty = "n", text.col = c("violet", "red", "grey", "springgreen4"))
# -----------------------------------------------------------------------------
 
    
# -----------------------------------------------------------------------------   
# Comparaison des courbes cumulées de pourcentage de biomasse
# -----------------------------------------------------------------------------   
par(mfrow= c(1,3))

plot(rnT0_m_prc$mean~rnT0_m_prc$TrL, las =1, type = "n", xlim = c(2, 4), ylim = c(0, 1),
     ylab = "Cumulated animal biomass (%)", xlab = "Trophic level", main = "T+0")
polygon(x=c(rnT0_m_prc$TrL, rev(rnT0_m_prc$TrL)), y=c(rnT0_m_prc$mean-rnT0_m_prc$se, rev(rnT0_m_prc$mean+rnT0_m_prc$se))/100000, 
        col = rgb(1,0,0, 0.1), border=NA)
lines(rnT0_m_prc$mean~rnT0_m_prc$TrL, col = "red", lwd = 3)
polygon(x=c(rtrrT0_m_prc$TrL,rev(rtrrT0_m_prc$TrL)), y=c(rtrrT0_m_prc$mean-rtrrT0_m_prc$se, rev(rtrrT0_m_prc$mean+rtrrT0_m_prc$se))/100000, 
        col = rgb(0,1,0, 0.1), border=NA)
lines(rtrrT0_m_prc$mean~rtrrT0_m_prc$TrL, col = "springgreen4", lwd = 3)
polygon(x=c(convT0_m_prc$TrL,rev(convT0_m_prc$TrL)), y=c(convT0_m_prc$mean-convT0_m_prc$se, rev(convT0_m_prc$mean+convT0_m_prc$se))/100000, 
        col = rgb(0,0,1, 0.1), border=NA)
lines(convT0_m_prc$mean~convT0_m_prc$TrL, lwd = 3, col = "violet")
polygon(x=c(rtT0_m_prc$TrL,rev(rtT0_m_prc$TrL)), y=c(rtT0_m_prc$mean-rtT0_m_prc$se, rev(rtT0_m_prc$mean+rtT0_m_prc$se))/100000, 
        col = rgb(0,0,0, 0.1), border=NA)
lines(rtT0_m_prc$mean~rtT0_m_prc$TrL, col = "grey", lwd = 3)
legend("topleft", c("CONV", "RN", "RT", "RT-RR"), bty = "n", text.col = c("violet", "red", "grey", "springgreen4"))

plot(rnT2_m_prc$mean~rnT2_m_prc$TrL, las =1, type = "n", xlim = c(2, 4), ylim = c(0, 1),
     ylab = "Cumulated animal biomass (%)", xlab = "Trophic level", main = "T+2")
polygon(x=c(rnT2_m_prc$TrL, rev(rnT2_m_prc$TrL)), y=c(rnT2_m_prc$mean-rnT2_m_prc$se, rev(rnT2_m_prc$mean+rnT2_m_prc$se))/100000, 
        col = rgb(1,0,0, 0.1), border=NA)
lines(rnT2_m_prc$mean~rnT2_m_prc$TrL, col = "red", lwd = 3)
polygon(x=c(rtrrT2_m_prc$TrL,rev(rtrrT2_m_prc$TrL)), y=c(rtrrT2_m_prc$mean-rtrrT2_m_prc$se, rev(rtrrT2_m_prc$mean+rtrrT2_m_prc$se))/100000, 
        col = rgb(0,1,0, 0.1), border=NA)
lines(rtrrT2_m_prc$mean~rtrrT2_m_prc$TrL, col = "springgreen4", lwd = 3)
polygon(x=c(convT2_m_prc$TrL,rev(convT2_m_prc$TrL)), y=c(convT2_m_prc$mean-convT2_m_prc$se, rev(convT2_m_prc$mean+convT2_m_prc$se))/100000, 
        col = rgb(0,0,1, 0.1), border=NA)
lines(convT2_m_prc$mean~convT2_m_prc$TrL, lwd = 3, col = "violet")
polygon(x=c(rtT2_m_prc$TrL,rev(rtT2_m_prc$TrL)), y=c(rtT2_m_prc$mean-rtT2_m_prc$se, rev(rtT2_m_prc$mean+rtT2_m_prc$se))/100000, 
        col = rgb(0,0,0, 0.1), border=NA)
lines(rtT2_m_prc$mean~rtT2_m_prc$TrL, col = "grey", lwd = 3)
legend("topleft", c("CONV", "RN", "RT", "RT-RR"), bty = "n", text.col = c("violet", "red", "grey", "springgreen4"))


plot(rnT0_m_prc$mean~rnT0_m_prc$TrL, las =1, type = "n", xlim = c(2, 4), ylim = c(0, 1),
     ylab = "Cumulated animal biomass (%)", xlab = "Trophic level", main = "T+4")
polygon(x=c(rnT4_m_prc$TrL, rev(rnT4_m_prc$TrL)), y=c(rnT4_m_prc$mean-rnT4_m_prc$se, rev(rnT4_m_prc$mean+rnT4_m_prc$se))/100000, 
        col = rgb(1,0,0, 0.1), border=NA)
lines(rnT4_m_prc$mean~rnT4_m_prc$TrL, col = "red", lwd = 3)
polygon(x=c(rtrrT4_m_prc$TrL,rev(rtrrT4_m_prc$TrL)), y=c(rtrrT4_m_prc$mean-rtrrT4_m_prc$se, rev(rtrrT4_m_prc$mean+rtrrT4_m_prc$se))/100000, 
        col = rgb(0,1,0, 0.1), border=NA)
lines(rtrrT4_m_prc$mean~rtrrT4_m_prc$TrL, col = "springgreen4", lwd = 3)
polygon(x=c(convT4_m_prc$TrL,rev(convT4_m_prc$TrL)), y=c(convT4_m_prc$mean-convT4_m_prc$se, rev(convT4_m_prc$mean+convT4_m_prc$se))/100000, 
        col = rgb(0,0,1, 0.1), border=NA)
lines(convT4_m_prc$mean~convT4_m_prc$TrL, lwd = 3, col = "violet")
polygon(x=c(rtT4_m_prc$TrL,rev(rtT4_m_prc$TrL)), y=c(rtT4_m_prc$mean-rtT4_m_prc$se, rev(rtT4_m_prc$mean+rtT4_m_prc$se))/100000, 
        col = rgb(0,0,0, 0.1), border=NA)
lines(rtT4_m_prc$mean~rtT4_m_prc$TrL, col = "grey", lwd = 3)
legend("topleft", c("CONV", "RN", "RT", "RT-RR"), bty = "n", text.col = c("violet", "red", "grey", "springgreen4"))
# -----------------------------------------------------------------------------   


# -----------------------------------------------------------------------------   
# Comparaison des courbes cumulées de pourcentage de biomasse par SdC
# -----------------------------------------------------------------------------   

par(mfrow= c(2,2))
COL = rgb(1,0,0, 0.1)
COL2 = rgb(0,1,0, 0.1)
soere_m_prc <- apply(cbind(rtT0_m_prc$mean, convT0_m_prc$mean, rnT0_m_prc$mean, rtrrT0_m_prc$mean), 1, mean)
soere_sd_prc <- apply(cbind(rtT0_m_prc$mean, convT0_m_prc$mean, rnT0_m_prc$mean, rtrrT0_m_prc$mean), 1, sd)/sqrt(4)

      plot(rnT0_m_prc$mean~rnT0_m_prc$TrL, las =1, type = "n", xlim = c(2, 4), ylim = c(50, 100),
           ylab = "Cumulated biomass (%)", xlab = "Trophic level", main = "CONV")
      #polygon(x=c(convT0_m_prc$TrL, rev(convT0_m_prc$TrL)), y=100*c(soere_m_prc-soere_sd_prc, rev(soere_m_prc+soere_sd_prc)), 
              col = COL2, border=NA)     
      #lines(100*soere_m_prc~convT0_m_prc$TrL,lty = 1)      
      polygon(x=c(convT0_m_prc$TrL, rev(convT0_m_prc$TrL)), y=100*c(convT0_m_prc$mean-convT0_m_prc$se, rev(convT0_m_prc$mean+convT0_m_prc$se)), 
              col = adjustcolor("steelblue4", alpha.f = 0.3), border=NA)
      lines(100*convT0_m_prc$mean~convT0_m_prc$TrL,lty = 2)
      polygon(x=c(convT2_m_prc$TrL, rev(convT2_m_prc$TrL)), y=100*c(convT2_m_prc$mean-convT2_m_prc$se, rev(convT2_m_prc$mean+convT2_m_prc$se)), 
              col = adjustcolor("steelblue4", alpha.f = 0.6), border=NA)
      lines(100*convT2_m_prc$mean~convT2_m_prc$TrL, lty = 3)
      polygon(x=c(convT4_m_prc$TrL, rev(convT4_m_prc$TrL)), y=100*c(convT4_m_prc$mean-convT4_m_prc$se, rev(convT4_m_prc$mean+convT4_m_prc$se)), 
              col = adjustcolor("steelblue4", alpha.f = 0.9), border=NA)
      lines(100*convT4_m_prc$mean~convT4_m_prc$TrL, lty = 4)
      
      legend("bottomright", c("T0", "T2", "T4"), bty = "n", lty = c(2,3, 4)) 
             col = c(adjustcolor("steelblue4", alpha.f = 0.3), adjustcolor("steelblue4", alpha.f = 0.6), 
                     adjustcolor("steelblue4", alpha.f = 0.9)))

          plot(rnT0_m_prc$mean~rnT0_m_prc$TrL, las =1, type = "n", xlim = c(2, 4), ylim = c(50, 100),
               ylab = "Cumulated biomass (%)", xlab = "Trophic level", main = "RN")
          #polygon(x=c(convT0_m_prc$TrL, rev(convT0_m_prc$TrL)), y=100*c(soere_m_prc-soere_sd_prc, rev(soere_m_prc+soere_sd_prc)), 
          #       col = adjustcolor("darkorchid4", alpha.f = 0.3), border=NA)     
          #lines(100*soere_m_prc~convT0_m_prc$TrL,lty = 1)  
          polygon(x=c(rnT0_m_prc$TrL, rev(rnT0_m_prc$TrL)), y=100*c(rnT0_m_prc$mean-rnT0_m_prc$se, rev(rnT0_m_prc$mean+rnT0_m_prc$se)), 
                  col = adjustcolor("darkorchid4", alpha.f = 0.3), border=NA)
          lines(100*rnT0_m_prc$mean~rnT0_m_prc$TrL,lty = 2)
          polygon(x=c(rnT2_m_prc$TrL, rev(rnT2_m_prc$TrL)), y=100*c(rnT2_m_prc$mean-rnT2_m_prc$se, rev(rnT2_m_prc$mean+rnT2_m_prc$se)), 
                  col = adjustcolor("darkorchid4", alpha.f = 0.6), border=NA)
          lines(100*rnT2_m_prc$mean~rnT2_m_prc$TrL, lty = 3)
          polygon(x=c(rnT4_m_prc$TrL, rev(rnT4_m_prc$TrL)), y=100*c(rnT4_m_prc$mean-rnT4_m_prc$se, rev(rnT4_m_prc$mean+rnT4_m_prc$se)), 
                  col = adjustcolor("darkorchid4", alpha.f = 0.9), border=NA)
          lines(100*rnT4_m_prc$mean~rnT4_m_prc$TrL, lty = 4)
          
          legend("bottomright", c( "T0", "T2", "T4"), bty = "n", lty = c(2,3, 4))

    plot(rnT0_m_prc$mean~rnT0_m_prc$TrL, las =1, type = "n", xlim = c(2, 4), ylim = c(50, 100),
         ylab = "Cumulated biomass (%)", xlab = "Trophic level", main = "RT")
    #polygon(x=c(convT0_m_prc$TrL, rev(convT0_m_prc$TrL)), y=100*c(soere_m_prc-soere_sd_prc, rev(soere_m_prc+soere_sd_prc)), 
    #       col = COL2, border=NA)     
    #lines(100*soere_m_prc~convT0_m_prc$TrL,lty = 1)  
    polygon(x=c(rtT0_m_prc$TrL, rev(rtT0_m_prc$TrL)), y=100*c(rtT0_m_prc$mean-rtT0_m_prc$se, rev(rtT0_m_prc$mean+rtT0_m_prc$se)), 
            col = adjustcolor("forestgreen", alpha.f = 0.3), border=NA)
    lines(100*rtT0_m_prc$mean~rtT0_m_prc$TrL,lty = 2)
    polygon(x=c(rtT2_m_prc$TrL, rev(rtT2_m_prc$TrL)), y=100*c(rtT2_m_prc$mean-rtT2_m_prc$se, rev(rtT2_m_prc$mean+rtT2_m_prc$se)), 
            col = adjustcolor("forestgreen", alpha.f = 0.6), border=NA)
    lines(100*rtT2_m_prc$mean~rtT2_m_prc$TrL, lty = 3)
    polygon(x=c(rtT4_m_prc$TrL, rev(rtT4_m_prc$TrL)), y=100*c(rtT4_m_prc$mean-rtT4_m_prc$se, rev(rtT4_m_prc$mean+rtT4_m_prc$se)), 
            col = adjustcolor("forestgreen", alpha.f = 0.9), border=NA)
    lines(100*rtT4_m_prc$mean~rtT4_m_prc$TrL, lty = 4)
    
    legend("bottomright", c("T0", "T2", "T4"), bty = "n", lty = c(2,3, 4))
    
    
        plot(rnT0_m_prc$mean~rnT0_m_prc$TrL, las =1, type = "n", xlim = c(2, 4), ylim = c(50, 100),
             ylab = "Cumulated biomass (%)", xlab = "Trophic level", main = "RT-RR")
        #polygon(x=c(convT0_m_prc$TrL, rev(convT0_m_prc$TrL)), y=100*c(soere_m_prc-soere_sd_prc, rev(soere_m_prc+soere_sd_prc)), 
        #        col = COL2, border=NA)     
        #lines(100*soere_m_prc~convT0_m_prc$TrL,lty = 1)  
        polygon(x=c(rtrrT0_m_prc$TrL, rev(rtrrT0_m_prc$TrL)), y=100*c(rtrrT0_m_prc$mean-rtrrT0_m_prc$se, rev(rtrrT0_m_prc$mean+rtrrT0_m_prc$se)), 
                col = adjustcolor("olivedrab4", alpha.f = 0.3), border=NA)
        lines(100*rtrrT0_m_prc$mean~rtrrT0_m_prc$TrL,lty = 2)
        polygon(x=c(rtrrT2_m_prc$TrL, rev(rtrrT2_m_prc$TrL)), y=100*c(rtrrT2_m_prc$mean-rtrrT2_m_prc$se, rev(rtrrT2_m_prc$mean+rtrrT2_m_prc$se)), 
                col = adjustcolor("olivedrab4", alpha.f = 0.6), border=NA)
        lines(100*rtrrT2_m_prc$mean~rtrrT2_m_prc$TrL, lty = 3)
        polygon(x=c(rtrrT4_m_prc$TrL, rev(rtrrT4_m_prc$TrL)), y=100*c(rtrrT4_m_prc$mean-rtrrT4_m_prc$se, rev(rtrrT4_m_prc$mean+rtrrT4_m_prc$se)), 
                col = adjustcolor("olivedrab4", alpha.f = 0.9), border=NA)
        lines(100*rtrrT4_m_prc$mean~rtrrT4_m_prc$TrL, lty = 4)
        
        legend("bottomright", c("T0", "T2", "T4"), bty = "n", lty = c(2,3, 4))
# -----------------------------------------------------------------------------   


# -----------------------------------------------------------------------------   
# Rapport entre les biomasses de 2 niveaux trophiques présentant un delta important
# -----------------------------------------------------------------------------   
convT4_thr1_i <- max(convT4_m$mean[convT4_m$TrL<2.4])
convT4_thr1_f <- max(convT4_m$mean[convT4_m$TrL<2.9])
convT4_thr1 <- (convT4_thr1_f-convT4_thr1_i)/convT4_thr1_f*100
convT4_thr1

rnT4_thr1_i <- max(rnT4_m$mean[rnT4_m$TrL<2.4])
rnT4_thr1_f <- max(rnT4_m$mean[rnT4_m$TrL<2.9])
rnT4_thr1 <- (rnT4_thr1_f-rnT4_thr1_i)/rnT4_thr1_f*100
rnT4_thr1

rtT4_thr1_i <- max(rtT4_m$mean[rtT4_m$TrL<2.4])
rtT4_thr1_f <- max(rtT4_m$mean[rtT4_m$TrL<2.9])
rtT4_thr1 <- (rtT4_thr1_f-rtT4_thr1_i)/rtT4_thr1_f*100
rtT4_thr1

rtrrT4_thr1_i <- max(rtrrT4_m$mean[rtrrT4_m$TrL<2.4])
rtrrT4_thr1_f <- max(rtrrT4_m$mean[rtrrT4_m$TrL<2.9])
rtrrT4_thr1 <- (rtrrT4_thr1_f-rtrrT4_thr1_i)/rtrrT4_thr1_f*100
rtrrT4_thr1
# -----------------------------------------------------------------------------   


# -----------------------------------------------------------------------------   
# Calcul des AUC
# -----------------------------------------------------------------------------   

indiceReseau <- unique(datatot[,1:4])
for (i in c(1:20)){
  indiceReseau$auc[i] <- sintegral(dbT0[[i]]$TrL, dbT0[[i]]$cumprc)$int
  indiceReseau$auc[i+20] <- sintegral(dbT2[[i]]$TrL, dbT2[[i]]$cumprc)$int
  indiceReseau$auc[i+40] <- sintegral(dbT4[[i]]$TrL, dbT4[[i]]$cumprc)$int
}

g <- ggplot(indiceReseau, aes (x = Date, y = auc))+
    facet_grid(~Mod)+
    theme_bw()+
    ylim(250, 350)+
    geom_boxplot()
g
# ----------------------------------------------------------------------------- 

# -----------------------------------------------------------------------------   
# Calcul des Biomasses totales
# -----------------------------------------------------------------------------   

for (i in c(1:20)){
  indiceReseau$bmax[i] <- max(dbT0[[i]]$cumbm)
  indiceReseau$bmax[i+20] <- max(dbT2[[i]]$cumbm)
  indiceReseau$bmax[i+40] <- max(dbT4[[i]]$cumbm)
}

g <- ggplot(indiceReseau, aes (x = Date, y = 0.00001*bmax))+
  facet_grid(~Mod)+
  ylab("Cumulated biomass (T ha-1)")+
  ylim (6, 12)+
  geom_boxplot()+
  theme_bw()
g
# ----------------------------------------------------------------------------- 

# -----------------------------------------------------------------------------   
# Calcul du niveau trophique max
# -----------------------------------------------------------------------------   

for (i in c(1:20)){
  indiceReseau$TLmax[i] <- max(dfTLT0[[i]]$TrL)
  indiceReseau$TLmax[i+20] <- max(dfTLT2[[i]]$TrL)
  indiceReseau$TLmax[i+40] <- max(dfTLT4[[i]]$TrL)
}

for (i in c(1:20)){
  indiceReseau$TLmean[i] <- mean(dfTLT0[[i]]$TrL)
  indiceReseau$TLmean[i+20] <- mean(dfTLT2[[i]]$TrL)
  indiceReseau$TLmean[i+40] <- mean(dfTLT4[[i]]$TrL)
}

indiceReseau$nb_trophGroups  <- nb_trophGroups

g <- ggplot(indiceReseau, aes (x = Date, y = TLmax))+
  facet_grid(~Mod)+
  ylab("maximum trophic level")+
  ylim (3, 4)+
  geom_boxplot()+
  theme_bw()
g
# ----------------------------------------------------------------------------- 

