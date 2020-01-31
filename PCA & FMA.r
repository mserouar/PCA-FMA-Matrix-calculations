# ---- Réduire ----
normalise <- function(X) {
  n <- nrow(X)
  P <- diag(rep(1/n, n))     # Matrice des poids
  E <- matrix(rep(1, n))     # Vecteur de 1
  G <- t(X) %*% P %*% E      # Centre de gravite
  Xc <- X- E %*% t(G)        # Matrice centree
  
  V <- t(Xc) %*% P %*% Xc    # Matrice variance covariance
  D <- diag(lapply(diag(V), sqrt)) #Matrice diagonale des ecart types (recuperee de V)
  Xr <- Xc %*% solve(D)      #Matrice centree reduite
}


# ---- Nipals ----
nipals <- function(Xr){
  P <- diag(rep(1/nrow(Xr), nrow(Xr))) # matrice des poids
  
  Xtmp <- Xr              # Xr la matrice des donnees, centrée réduite
  threshold <- 0.00000001   # limite pour la convergeance
  vectp <- matrix(nrow = dim(Xr)[2], ncol = dim(Xr)[2])    # creation matrice vide
  valp <- rep(NA, dim(Xr)[2])                              # creation vecteur vide
  new_axes <- matrix(nrow = dim(Xr)[1], ncol = dim(Xr)[2]) # creation matrice vide
  
  for (s in 1:dim(Xr)[2]){ # Pour chaque dimension calculable
    Fs <- Xtmp[,1]         # Initialisation simple
    repeat{
      us <- t(Xtmp) %*% Fs /as.numeric(t(Fs)%*% Fs) # Calcul du vecteur propre de rang s
      us <- us / sqrt(as.numeric(t(us) %*% us))     # Normalisation du vecteur propre (Us / norme(us))
      Fs_old <- Fs                                  # stockage de Fs
      Fs <- Xtmp %*% us / as.numeric(t(us) %*% us)  # Actualisation de Fs grâce au vecteur propre
      diff_Fs <- Fs_old - Fs                        # Pour estimer si ça converge ou non 
      if (sqrt(as.numeric(t(diff_Fs) %*% diff_Fs)) < threshold) break # Si la norme de diff_Fs est petite on est bon
    }
    Xtmp <- Xtmp - Fs %*% t(us) # Retire de la matrice l'information donnée par l'axe trouvé 
    vectp[,s] <- us                                # Stocke les vecteurs propres
    valp[s] <- t(us) %*% t(Xr) %*% P %*% Xr %*% us # Calcule la valeur propre correspondante
    new_axes[,s] <- Fs                             # Stocke les coordonnées des nouveaux axes
  }
  #
  # vp <- apply(vect_propres, 2, function(vectp) t(vectp) %*% t(Xr) %*% P %*% Xr %*% vectp)
  list(valp = valp,  vectp = vectp, Fu = new_axes)
}

# ---- Fonction ACP ----
acp_hardcode <- function(X, cr = T){
  P <- diag(rep(1/nrow(X), nrow(X))) # Matrice diagonale des poids
  
  # centrée réduite ou non selon le choix
  if(cr)
    Xr <- normalise(X)
  else
    Xr <- X
  
  # Coordonnees des individus, valeurs propres et vecteurs propres
  res.acp <- nipals(Xr)
  
  # coordonnees des variables
  valp <- sapply(res.acp$valp, function(x) 1/sqrt(x)) # 1/sqrt(val propres)
  tmp <- rbind(valp, res.acp$Fu) # tableau sur lequel utiliser apply
                                 # premier element de chaque col est 1/sqrt(val propres)
  Gu <- apply(tmp, 2, function(Fs) Fs[1] * t(Xr) %*% P %*% Fs[2:nrow(tmp)]) # formule de liaison
  res.acp$Gu <- Gu
  
  res.acp
}

# ---- Fonction AFM ----
afm_hardcode <- function(X, groups){
  X_pond <- X # Initialisation de la matrice ponderee
  i <- 1      # compteur pour selectioner les groupes
  for (k in groups){ # Iteration sur chaque groupe
    i2 <- k+i-1      # Autre compteur pour selectionner les groupes
    Xtmp <- X[,i:i2] # Matrice d'un groupe
    valp1 <- acp_hardcode(Xtmp)$valp[1] # premiere valeur propre de l'ACP d'un groupe
    X_pond[,i:i2] <- normalise(Xtmp) / sqrt(valp1) # Ponderation du groupe
    i <- i + k # incrementation du compteur...
  }
  acp_hardcode(X_pond, cr = FALSE) # ACP globale
}


# ---- Tests des fonction ACP et AFM----
require(FactoMineR)

data("decathlon")
X_decat <- as.matrix(decathlon[,-13])
res.acp <- acp_hardcode(X_decat)
res.acp.factoMineR <- PCA(X_decat, graph = F, ncp = 12)
# valeurs propres --> precis a 13 decimales pres
mean(round(res.acp$valp, 13) == round(res.acp.factoMineR$eig[,1], 13))
# individus --> precis a 5 decimales pres
mean(round(abs(res.acp$Fu), 5) == round(abs(res.acp.factoMineR$ind$coord), 5))
# variables --> precis a 6 decimales pres
mean(round(abs(res.acp$Gu), 6) == round(abs(res.acp.factoMineR$var$coord), 6))


data("wine")
X_wine <- as.matrix(wine[-c(1:2,30:31)]) # sans les illustratives
res.afm <- afm_hardcode(X_wine, c(5,3,10,9))
res.afm.factoMineR <- MFA(wine, group=c(2,5,3,10,9,2), type=c("n",rep("s",5)), ncp=27,
                          name.group=c("orig","olf","vis","olfag","gust","ens"),
                          num.group.sup=c(1,6), graph=FALSE)
# valeurs propres --> precis a 13 decimales pres
mean(round(res.afm$valp[1:20], 13) == round(res.afm.factoMineR$eig[,1], 13))
# individus moyens --> precis a 5 decimales pres
mean(round(abs(res.afm$Fu[,1:20]), 5) == round(abs(res.afm.factoMineR$ind$coord), 5))
