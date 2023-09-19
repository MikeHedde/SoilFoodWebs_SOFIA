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
# FUNCTION
# -----------------------------------------------------------------------------
# 1. Compute network fluxes
Fluxfct <- function(DBtest, adjMat) {
  
  # Parameters definition
  nomsp   <- names(adjMat)[-1]      # List of nodes
  nbgp    <- length(nomsp)          # Total number of nodes 
  web     <- tibble(adjMat) %>%
            pivot_longer(cols = !ID, names_to = "variable") %>%
            dplyr::filter(value > 0) 
  DBtest <- as.data.frame(DBtest)
  web1 <- web %>%
    mutate(nbinter = 1:length(value)) %>%
    rename(prey = variable, pred = ID, pred_pref = value)
  G <- numeric(web1$nbinter)    
  B <- numeric(web1$nbinter)
  A <- diag(1, web1$nbinter)
  
  for (i in 1:nbgp) {
    predi <- subset(web1, pred == i)
    iprey <- predi$prey
    
    if (length(iprey) > 0) {
      nobiom <- iprey %in% c("BACT", "FUNG")                # if preys are bacteria or fungi
      Sum_nB <- sum(predi$pred_pref[nobiom])
      
    if (length(nobiom) > 0) {
      G[predi$numinter[nobiom]] <- predi$pred_pref[nobiom]
    }
      
      withbiom <- !nobiom
      if (any(withbiom)) {
        Sum_wB <- sum(predi$pred_pref[withbiom] * DBtest$biomass[DBtest$node == iprey[withbiom]])
        G[predi$numinter[withbiom]] <- (1 - Sum_nB) * predi$pred_pref[withbiom] * DBtest$biomass[DBtest$node == iprey[withbiom]] / Sum_wB
        
        # Utilisation de ifelse pour éviter les conditions répétées
        G[predi$numinter[withbiom]] <- ifelse(Sum_wB == 0, predi$pred_pref[withbiom] * DBtest$biomass[DBtest$node == iprey[withbiom]], G[predi$numinter[withbiom]])
        G[predi$numinter[nobiom]] <- ifelse(Sum_wB == 0, predi$pred_pref[nobiom] / Sum_nb, G[predi$numinter[nobiom]])
      }
      
      if (length(iprey) == 1) {
        G[predi$numinter] <- 1
      }
      
      preyi <- subset(web1, prey == i)
      ipred <- preyi$pred
      A[predi$numinter, preyi$numinter] <- A[predi$numinter, preyi$numinter] - G[predi$numinter] / DBtest$a[i]
    }
  }
  
  F <- solve(A, B)
  return(F)
}

# 2. integration des flux dans une matrice d'interactions
fun_web <- function(fl) {
  webmat <- replicate(20, matrix(0, nrow=nbgp, ncol=nbgp), simplify = F) # liste de matrices vides pour entrer les flux; le nombre de matrice correspond au nombre de parcelles (par ann?e)
  for (j in 1:20) {
    for (i in 1:(nbinter-1)) {
      webmat[[j]][web$prey[i], web$pred[i]] <- fl[[j]][i]    
      if (fl[[j]][i]<0) {print(fl[[j]][i])                # Les flux n?gatifs sont imprim?s 
        webmat[[j]][web$prey[i], web$pred[i]]<- 0       # ? sp?cifier si on veut ?liminer les flux n?gatifs, mais le probl?me est que ?a "cache" les probl?mes et donne une valeur artificielle de 0
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
