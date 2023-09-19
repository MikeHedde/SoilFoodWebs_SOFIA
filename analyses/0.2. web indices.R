# -----------------------------------------------------------------------------
# PROJECT:
#    Linking soil foodweb and soil functions dynamics under
#    changing agricultural practices managment
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# STEP:
#    0.2   Computing foodweb indices
# -----------------------------------------------------------------------------


myFoodWeb_indices <- function(web, biom){
# 1. Trophic levels
## Compute number, mean and max trophic level(s)
GetIND <- function(web){  
  tweb <- t(web)
  
  ## Line sums must equal 0
  rs <- rowSums(tweb)
  for(i in 1:length(tweb[,1]))
    tweb[i,tweb[i,]>0] = tweb[i,tweb[i,]>0]/rs[i]
  
  nb.TL <- try(solve(diag(length(tweb[,1])) - tweb), T)
  if(class(nb.TL) == "try-error")
    nbTL <- rep(NA, length(tweb[,1]))
  if(class(nb.TL) != "try-error")
    nbTL <- rowSums(nb.TL)
  return(TL <- list(nbTL, mean(nbTL), max(nbTL)))
}


# 2. Link number, diversity and evenness
Nblinks  <- function(web) {    
  flux   <- as.vector(web)
  flux   <- flux[flux!=0]
  nblink <- length(flux)
  nblink
}

Divlinks <- function(web) {
  flux   <- as.vector(web)
  flux   <- flux[flux!=0]
  Ftot   <- sum(flux)
  div    <- 0
  for (i in 1:length(flux))   
        div  <- div - flux[i]/Ftot*log2(flux[i]/Ftot)
  div
}

Evenlinks <- function(web) {
  even <- Divlinks(web)/Nblinks(web)
  even
}


# 3. Path length ratio
Path <- function(web, num){    
  tweb<- t(web)
  Vectpath <- rep(0, length(web[,1])) 
  Vectpath[num] <- 1
}

PropPath <- function(web, numP1, numP2, biom){  
  #biom : vector containing the biomass in the web dataframe
  path1 <- Path(web, numP1)
  path2 <- Path(web, numP2)
  proppath <-(sum(biom * path1)) / (sum(biom * path2))
  proppath
}

## Fungi-to-bacteria path length ratio
fb_plr <- PropPath(web = web, numP1 = "BACT", numP2 = "FUNG", biom = biom)

## Brown-to-green path length ratio
dp_plr <- PropPath(web = web, numP1 = "Plantparts", numP2 = "Detritus", biom = biom)

# Return results
return(res <- list(TL, nblink, div, even, fb_plr, dp_plr))
}
