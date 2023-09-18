# -----------------------------------------------------------------------------
# PROJECT:
#    Linking soil foodweb and soil functions dynamics under
#    changing agricultural practices managment
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# STEP:
#    0.?   Performing the multivariate temporal dynamics analysis
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# NOTES:
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# LIBRARIES:
library(ade4)
library(XLConnect)
library(FSA)
library(reshape2)
library(ggplot2)
library(gtools)
library(lme4)
library(nlme)
library(car)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# SCRIPT
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# FUNCTION
# -----------------------------------------------------------------------------
heatcor <- function(X){
  cormat <- round(cor(X),2)

        reorder_cormat <- function(cormat){ # Utiliser la corrélation entre les variables comme mésure de distance
        dd <- as.dist((1-cormat)/2)
        hc <- hclust(dd)
        cormat <-cormat[hc$order, hc$order]
      }
        get_upper_tri <- function(cormat){
          cormat[lower.tri(cormat)]<- NA
          return(cormat)
        }

      upper_tri <- get_upper_tri(cormat)
      melted_cormat <- melt(upper_tri, na.rm = TRUE)
      
      ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
        geom_tile(color = "white")+
        scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                             midpoint = 0, limit = c(-1,1), space = "Lab",
                             name="Pearson\nCorrelation") +
        theme_minimal()+ # minimal theme
        theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                         size = 12, hjust = 1))+
        coord_fixed()
      ggheatmap + 
        geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
        theme(
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          legend.justification = c(1, 0),
          legend.position = c(0.6, 0.7),
          legend.direction = "horizontal")+
        guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                     title.position = "top", title.hjust = 0.5))}
#--------------------------------------------


# directory
setwd("C:/Users/heddemic/Nextcloud2/Hedde M/5. RH/4. Post-doctorant/2014 SOFIA V Coudrain/dossier juin 2016/Scripts finaux")

#--------------------------------------------
# charger et mettre en forme les dataframe
#--------------------------------------------
 # indices des réseaux     

  df_totdiv <- read.xlsx2("df_totdiv.xlsx", sheetIndex = 1, h = T, rownames = 1) 
  df_totdiv <- cbind(df_totdiv, indiceReseau[,-(1:4)])  
                  #indiceReseau vient du script "MH courbes cumulées TrL biom.R"
                  # redunMelted vient du script "FUnct Redund.R"
  df_totdivA2 <- df_totdiv[df_totdiv$Mod != "RR-PER", ]
  fac <- df_totdivA2[,1:5]
  df_totdivA2 <- cbind(df_totdivA2[,-c(1:5)], FunctRed = redunMelted$value[redunMelted$variable == "Total Network" & redunMelted$Mod != "RR-PER"])
  
 # paramètres des sols     
  soil.dat <- read.table("data_basefile.txt",header=TRUE, sep='\t', dec=",")
  soil.dat <- Subset(soil.dat, Date!="Year3"& PlotID!="T4A16" & Mod != "RR-PER")
  
  soil.fun0 <- data.frame(nitrif = soil.dat$Nitrif_0.5, 
                         Nmin = soil.dat$Mineraliz_0.5,
                         Nimm = soil.dat$Organisation_0.5, 
                         Cmin = soil.dat$Respiration_sol_0.5, 
                         ure = soil.dat$Urease_0.5, 
                         glu = soil.dat$Glucosidase_0.5)
                         #glu_ure = log(soil.dat$Glucosidase_0.5)/log(soil.dat$Urease_0.5)
                         #effC_glu = soil.dat$Glucosidase_0.5/soil.dat$Respiration_sol_0.5, 
                         #effN_ure = soil.dat$Urease_0.5/soil.dat$Mineraliz_0.5)
  soil.fun0$nitrif[47] <- mean(soil.fun0$nitrif[33:48], na.rm = T) # remplace une valeur manquante
  
  
  soil.fun5 <- data.frame(nitrif = soil.dat$Nitrif_0.5, 
                          Nmin = soil.dat$Mineraliz_0.5,
                          Nimm = soil.dat$Organisation_0.5, 
                          Cmin = soil.dat$Respiration_sol_0.5, 
                          ure = soil.dat$Urease_0.5/soil.dat$N_total_0.5, 
                          glu = soil.dat$Glucosidase_0.5/soil.dat$C_total_0.5)
  soil.fun5$nitrif[47] <- mean(soil.fun5$nitrif[33:48], na.rm = T) # remplace une valeur manquante
  
  
  soil.fun20 <- data.frame(nitrif20 = soil.dat$Nitrif_5.20, 
                          Nmin20 = soil.dat$Mineraliz_5.20,
                          Nimm20 = soil.dat$Organisation_5.20, 
                          Cmin20 = soil.dat$Respiration_sol_5.20, 
                          ure20 = soil.dat$Urease_5.20, 
                          glu20 = soil.dat$Glucosidase_5.20)

  # Indices de réseau
  foodweb <- data.frame("Brown:Green" = as.numeric(as.character(df_totdivA2$BGpath)), 
                        "Bact:Fungi" = as.numeric(as.character(df_totdivA2$BFpath)), 
                        TL_mean = df_totdivA2$TLmean,
                        TL_max = df_totdivA2$TLmax,
                        TG_num = df_totdivA2$nb_trophGroups,
                        Lks_num = as.numeric(as.character(df_totdivA2$NbL)), 
                        Lks_div = as.numeric(as.character(df_totdivA2$DivL)),
                        FunctRed = df_totdivA2$FunctRed,
                        BottomHeavy = df_totdivA2$auc,
                        Biom_tot = df_totdivA2$bmax)
  #--------------------------------------------  
  
  
  #--------------------------------------------
  # Approche combinatoire pour isoler 2 x 6 variables
  #--------------------------------------------  
        # Test de toutes les CoStatis possibles pour identifier les meilleurs k-tableaux
        ind_fw <- as.matrix(colnames(foodweb))
        ind_sp <- as.matrix(colnames(soil.fun0))
        maxInd_fw = 6  #on  retient 6 indices de réseaux
        maxInd_sp = 6  #on retient 6 fonctions  
        comb_fw <-combinations(n = 10, r = maxInd_fw, v = ind_fw, repeats.allowed = FALSE) # on crée toutes les combinaisons 
        comb_sp <-combinations(n = 9, r = maxInd_sp, v = ind_sp, repeats.allowed = FALSE)# on crée toutes les combinaisons 
        # on crée une dataframe avec toutes les combinaisons d'indice de réseau x toutes les combinaisons de fonctions des sols
        comb_fw2 <- comb_fw[rep(1:nrow(comb_fw), each = nrow(comb_sp)),]
        comb_sp2 <- comb_sp[rep(1:nrow(comb_sp), times = nrow(comb_fw)),]
        comb <- cbind(comb_sp2, comb_fw2) #soit 17640 combinaisons
      
        # on extrait et on stocke dans un dataframe la pval et le RV du randtest de toutes les combinaisons de k-tableaux
        RDTEST<- matrix(NA, ncol = 2, nrow = nrow(comb))
        for ( i in 1:nrow(comb)){
          wit_s <- withinpca(soil.fun0[, colnames(soil.fun0) %in% comb[i, 1:maxInd_sp]], as.factor(fac$Date), scan = FALSE, scal = "partial") 
          kta_s <- ktab.within(wit_s, colnames = fac$Plot)
          
          wit_fw <- withinpca(foodweb[, colnames(foodweb) %in% comb[i, -c(1:maxInd_sp)]], as.factor(fac$Date), scan = FALSE, scal = "partial") 
          kta_fw <- ktab.within(wit_fw, colnames = fac$Plot)
          
          cost <- costatis(kta_s, kta_fw, scannf = F)
          RDTEST[i,1] <- costatis.randtest(kta_s, kta_fw)$pval
          RDTEST[i,2] <- costatis.randtest(kta_s, kta_fw)$obs
          RDTEST
        }
          
        comb2 <- cbind(comb, as.data.frame(RDTEST))
        colnames(comb2) <- c(paste("FW", 1:6, sep = ""), paste("Soil", 1:6, sep = ""), "pval", "RV")
        
        # on représente la distribution des pval 
        RVdist <- ggplot(comb2, aes(x = RV))+
          geom_histogram(binwidth = 0.015)+
          geom_freqpoly(binwidth = 0.015)+
          xlab("Distribution of the RV coefficient\nof all possible CoStatis (n = 999)")+
          geom_vline(xintercept = 0.61, color = "red", size = 1)
      
        pvaldist <- ggplot(comb2, aes(x = pval))+
          geom_histogram(binwidth = 0.02)+
          geom_freqpoly(binwidth = 0.02)+
          xlab("Distribution of the simulated p-value of \nMonte-Carlo test on all possible CoStatis (n = 999)")+
          geom_vline(xintercept = 0.001, color = "red", size = 1)
      
        plot_grid(RVdist, pvaldist, ncol = 1, nrow = 2, labels = "AUTO")
        
        comb3 <- comb2[comb2[,13] < 0.002, ]
        comb3[order(comb3[,14], decreasing = T),]
        
        
        #je préfère ne pas garder les AE brutes et garder le rapport B/F
        select_finale <- comb3[which(comb3[,1] == "Cmin" & comb3[,3] != "glu" & comb3[,2] != "glu" & comb3[,6] != "ure" & comb3[,14]>0.6),]
        select_finale <- select_finale[order(select_finale[,14], decreasing = T),]                           
        select_finale       # il reste 4 possibilités, je conserve la ligne 4
              sel = 4
  
          #Analyse finale
          wit_s <- withinpca(soil.fun0[, colnames(soil.fun0) %in% t(select_finale[4, 1:6])], as.factor(fac$Date), scan = FALSE, scal = "partial") 
          kta_s <- ktab.within(wit_s, colnames = fac$Plot)
          
          wit_fw <- withinpca(foodweb[, colnames(foodweb) %in% t(select_finale[4, 7:12])], as.factor(fac$Date), scan = FALSE, scal = "partial") 
          kta_fw <- ktab.within(wit_fw, colnames = fac$Plot)
          
          cost <- costatis(kta_s, kta_fw, scannf = F)
          costatis.randtest(kta_s, kta_fw)
          plot(cost)
        
        par(mfrow = c(1, 1))
        colnames(cost$lY) <- colnames(cost$lX) <- colnames(cost$l1) <- colnames(cost$c1) <- c("A1", "A2")
        ade4::s.arrow(rbind(cost$l1, cost$c1), boxes = F)
        s.match.class(cost$lY, cost$lX, soil.dat$Mod[1:16])
  #--------------------------------------------

        #--------------------------------------------
        # utilisation de toutes les variables
        #--------------------------------------------        
        wit_s5 <- withinpca(soil.fun5, as.factor(fac$Date), scan = FALSE, scal = "partial") 
        kta_s5 <- ktab.within(wit_s5, colnames = fac$Plot)
 
        wit_fw <- withinpca(foodweb, as.factor(fac$Date), scan = FALSE, scal = "partial") 
        kta_fw <- ktab.within(wit_fw, colnames = fac$Plot)
        
        cost5 <- costatis(kta_s5, kta_fw, scannf = F)
        costatis.randtest(kta_s5, kta_fw)
        
        colnames(cost5$lY) <- colnames(cost5$lX) <- colnames(cost5$l1) <- colnames(cost5$c1) <- c("A1", "A2")
        
        s.match.class(cost5$lY, cost5$lX, soil.dat$Mod[1:16])    
        proj <- rbind.data.frame(cost5$lX, cost5$lY)
        proj$data <- rep(c("soil", "foodweb"), each = 16)
        proj$mod <- rep(soil.dat$Mod[1:16], 2)
        gg <- merge(proj,aggregate(cbind(mean.x = A1,mean.y = A2)~data+mod,proj,mean), by = c("data", "mod"))
        gg2 <- aggregate(cbind(mean.x = A1,mean.y = A2)~data+mod,proj,mean)
        
        cc <- rbind.data.frame(cost5$c1, cost5$l1)
        cc_plot <- ggplot(cc, aes(A1, A2))+
          ylim(c(min(cc$A2)*1.005, max(cc$A2)*1.005))+
          xlim(c(min(cc$A1)*1.2, max(cc$A1)*1.25))+
          xlab("Component 1")+
          ylab("Component 2")+
          geom_vline(xintercept = 0, size = 2) + geom_hline(yintercept = 0, size = 2)+
          geom_segment(aes(x = rep(0, 16), xend = A1, y = rep(0,16), yend= A2), 
                       arrow = arrow(length=unit(0.3,"cm")), size = .5)+
          geom_label(label = rownames(cc), 
                     hjust=c(0, 1, 1,.5, 1, .5, .8, .5, 1, 0, 1, 1, 1, .5, .5,.5), 
                     vjust=c(0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.8, 0, 1, 0, 0, 0), 
                     color = c(rep("orange", 6), rep("maroon", 10)), size = 6, fontface = "bold")+
          theme_new()
        cc_plot
        
        proj_plot <- ggplot(gg, aes(x = A1, y = A2))+
          ylim(c(min(gg$A2)*1.005, max(gg$A2)*1.005))+
          xlim(c(min(gg$A1)*1.005, max(gg$A1)*1.005))+
          xlab("Component 1")+
          ylab("Component 2")+
          geom_vline(xintercept = 0, size = 2) + 
          geom_hline(yintercept = 0, size = 2)+
            geom_point(aes(shape = mod, col = data, fill = mod), size = 3, color = "black")+
            geom_segment(data = gg, aes(x= mean.x, y= mean.y, xend= A1, yend=A2), linetype = 3)+
            geom_label(data = gg2, aes(x = mean.x, y = mean.y, label = mod, fill = mod, col = data),
                       hjust=.5, vjust=.5, size = 6, fontface = "bold")+
          scale_shape_manual(values=c(21, 22, 23, 24), name = "")+
          scale_fill_manual(values = c("steelblue4", "darkorchid4", "forestgreen", "olivedrab4"))+
          scale_color_manual(values = c("black", "white"))+
          guides(size = FALSE, fill = F, shape = F, col = F)+
          theme_new()
           
        
        plot_grid(cc_plot, proj_plot, labels = "AUTO")
        #--------------------------------------------       



       

        
        
        
   
        
