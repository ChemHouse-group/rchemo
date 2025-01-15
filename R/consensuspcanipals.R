#consensuspcanipals <- function(Xbl, nlv = NULL, bscal = "frob", scal = FALSE, weights = NULL) {
consensuspcanipals <- function(Xlist, blockscaling = TRUE, weights = NULL, nlv, Xscaling = c("none", "pareto", "sd")[1],
                               tol = .Machine$double.eps^0.5, maxit = 200) {
  # Xlist : Liste de matrices (ou blocs)
  # nlv : Nombre de variables latentes à calculer
  # bscal : Type de mise à l'échelle des blocs
  # scal : Booléen indiquant si on veut scaler les blocs
  # weights : Vecteur de poids des observations
  
  # Initialiser les poids si non fournis
  if (is.null(weights)) {
    weights <- rep(1, ncol(Xlist[[1]])) # Poids égaux pour chaque observation
  }
  
  xmeanslist <- lapply(1:length(Xlist), function(i) .colmeans(Xlist[[i]], weights = weights))
  
  if((length(Xscaling)==1) & (length(Xlist)>1)){Xscaling = rep(Xscaling, length(Xlist))}
  
  xscaleslist <- list()
  for(i in 1:length(Xlist)){
    if(Xscaling[i] == "none"){xscaleslist[[i]] <- rep(1, ncol(Xlist[[i]]))}
    if(Xscaling[i] == "pareto"){xscaleslist[[i]] <- sqrt(sqrt(.colvars(Xlist[[i]], weights = weights)))}
    if(Xscaling[i] == "sd"){xscaleslist[[i]] <- sqrt(.colvars(Xlist[[i]], weights = weights))}
  }
  Xlist <- lapply(1:length(Xlist), function(i) scale(Xlist[[i]], center = xmeanslist[[i]], scale = xscaleslist[[i]]))
  
  if(blockscaling==TRUE){
    Xblockscaled <- blockscal(Xtrain = Xlist, weights = weights)
    Xlist <- Xblockscaled$Xtrain
    Xnorms <- Xblockscaled$disp
  }else{
    Xnorms <- NA
  }
  
  
  # Initialiser les paramètres (éventuellement similaires à "recovkw")
  par <- list(nlv = nlv, blockscaling = blockscaling, Xscaling = Xscaling, tol = tol, maxit = maxit)
  
  # Calculs de base sur les blocs de données (similaire à "blockscal")
  # fitmbl <- blockscal(Xtrain = Xlist, weights = weights)
  
  # Initialisation des matrices et vecteurs nécessaires
  nbl <- length(Xlist)
  n <- nrow(Xlist[[1]])
  sqrtw <- sqrt(weights)
  
  # Créer les matrices vides pour stocker les résultats
  U <- matrix(0, n, nlv)
  W <- matrix(0, nbl, nlv)
  # Tbl <- vector("list", nbl)
  Tbl <- lapply(1:nbl, function(k) matrix(0, n, nlv))
  # Tb <- vector("list", nlv)
  Tb <- lapply(1:nlv, function(k) matrix(0, n, nbl))
  # Wbl <- vector("list", nbl)
  Wbl <- lapply(1:nbl, function(k) matrix(0, n, nlv))
  lb <- matrix(0, nbl, nlv)
  mu <- numeric(nlv)
  niter <- numeric(nlv)
  
  # Iterations Nipals (comme la méthode utilisée en PCA)
  for (a in 1:nlv) {
    X <- do.call(cbind, Xlist)  # Concatenation horizontale des blocs
    u <- nipals(X)$u           # Appliquer NIPALS pour calculer les scores
    iter <- 1
    cont <- TRUE
    
    while (cont) {
      u0 <- u  # Sauvegarde des scores pour vérifier la convergence
      for (k in 1:nbl) {
        wk <- t(Xlist[[k]]) %*% u
        dk <- sqrt(sum(wk^2))
        wk <- wk / dk
        tk <- Xlist[[k]] %*% wk
        Tb[[a]][, k] <- tk
        Tbl[[k]][, a] <- (1/sqrtw) * tk
        Wbl[[k]][, a] <- wk
        lb[k, a] <- dk^2
      }
      res <- nipals(Tb[[a]])
      u <- res$u
      w <- res$v
      dif <- sum((u - u0)^2)
      iter <- iter + 1
      if (dif < par$tol || iter > par$maxit) {
        cont <- FALSE
      }
    }
    
    niter[a] <- iter - 1
    U[, a] <- u
    W[, a] <- w 
    mu[a] <- sum(lb[, a])  # Calcul des valeurs propres
    for (k in 1:nbl) {
      Xlist[[k]] <- Xlist[[k]] - u %*% t(u) %*% Xlist[[k]]
    }
  }
  
  # Résultat final pour les scores globaux
  T <- (diag(1 / sqrt(weights))) %*% (sqrt(mu) * U)###### A VERIFIER
  
  # Retourner le modèle final
  result <- list(T = T, U = U, W = W, Tbl = Tbl, Tb = Tb, Wbl = Wbl, lb = lb, mu = mu, 
                 xmeans = xmeanslist, xscales = xscaleslist, weights = weights, blockscaling = blockscaling, Xnorms = Xnorms, #fitmbl = fitmbl, 
                 niter = niter, params = par)
  return(result)
}

# Exécution de la PCA multibloc
Xbl <- list(matrix(rnorm(100), 10, 10), matrix(rnorm(100), 10, 10))  # Exemple de blocs
result <- consensuspcanipals(Xbl, nlv = 3)


transform.Consensuspcanipals <- function(object, X, nlv = NULL) {
  # Calcul des variables latentes (scores) à partir du modèle ajusté
  nlv <- if (is.null(nlv)) ncol(object$T) else min(nlv, ncol(object$T))
  
  X <- lapply(1:length(X), function(i) .mat(X[[i]]))
  
  # Xinit <- X
  p <- dim(object$P)[1]
  n <- dim(object$T)[1]
  # TT <- object$weights * object$T * object$T
  # tt <- colSums(TT)
  
  X <- lapply(1:length(X), function(i) scale(.mat(X[[i]]), center = object$xmeans[[i]], scale = object$xscales[[i]]))
  if(object$blockscaling==TRUE){X <- blockscal(Xtrain = X, weights = object$weights)$Xtrain}
  # Xconc <- do.call("cbind",X)
  
  U <- matrix(0, nrow(X[[1]]), nlv)
  TB <- matrix(0, nrow(X[[1]]), length(X))
  Tbl <- lapply(1:length(X), function(k) matrix(0, n, nlv))
  u <- tk <- numeric(nrow(X[[1]]))
  
  for (a in 1:nlv) {
    # u <- numeric(nrow(X[[1]]))
    for (k in 1:length(X)) {
      TB[,k] <- Tbl[[k]][,a] <- tk <- X[[k]] %*% object$Wbl[[k]][, a]
      # u <- u + X[[k]] %*% object$Wbl[[k]][, a]
    }
    U[, a] <- u <- (1 / sqrt(object$mu[a])) * TB %*% object$W[, a] # A VERFIIER
    for (k in 1:length(X)) {
      Px <- sqrt(object$lb[k, a]) %*% t(object$Wbl[[k]][, a])
      X[[k]] <- X[[k]] - (u %*% Px)
    }
  }
  T <- sqrt((object$mu[1:nlv])) * U
  return(list(T = T, Tbl = Tbl))
}


summary.Consensuspcanipals <- function(object, X) {

  # Q <- typeof(X[[1]][1, 1])
  nbl <- length(X)
  nlv <- ncol(object$T)
  sqrtw <- sqrt(object$weights)
  
  X <- lapply(1:length(X), function(i) .mat(X[[i]]))
  
  p <- dim(object$P)[1]
  n <- dim(object$T)[1]
  # TT <- object$weights * object$T * object$T
  # tt <- colSums(TT)
  
  X <- lapply(1:length(X), function(i) scale(.mat(X[[i]]), center = object$xmeans[[i]], scale = object$xscales[[i]]))
  if(object$blockscaling==TRUE){X <- blockscal(Xtrain = X, weights = object$weights)$Xtrain}
  
  # for (k in 1:nbl) {
  #   zX[[k]] <- sqrtw * zX[[k]]
  # }

  Xconc <- do.call(cbind, X)
  
  # Explained_X
  sstot <- rep(0, nbl)
  for (k in 1:nbl) {
    sstot[k] <- sum(object$weights * X[[k]] * X[[k]], na.rm = TRUE)
  }

  tt <- colSums(object$lb)
  pvar <- tt / sum(sstot)
  cumpvar <- cumsum(pvar)
  explvarx <- data.frame(lv = 1:nlv, var = tt, pvar = pvar, cumpvar = cumpvar)
  
  # Contribution of the blocks to global scores
  # z = fscale(object.lb, colsum(object.lb))
  z <- scale(object$lb, center=FALSE, scale = colSums(object$lb))
  
  # contr_block = DataFrame(z, string.("lv", 1:nlv))
  contr_block <- data.frame(z)
  colnames(contr_block) <- paste0("lv", 1:nlv)
  
  # Explained inertia for each block
  # z = fscale((object.lb)', sstot)'
  z <- t(scale(t(object$lb), center=FALSE, scale = sstot))
  
  # explX = DataFrame(z, string.("lv", 1:nlv))
  explX <- data.frame(z)
  colnames(explX) <- paste0("lv", 1:nlv)
  
  # Correlation between original variables and global scores
  z <- cor(Xconc, object$U)
  
  # corx2t = DataFrame(z, string.("lv", 1:nlv))
  corx2t <- data.frame(z)
  colnames(corx2t) <- paste0("lv", 1:nlv)
  
  # Correlation between block scores and global scores
  z <- vector("list", nlv)
  for (a in 1:nlv) {
    z[[a]] <- cor(object$Tb[[a]], object$U[, a])
  }
  
  cortb2t <- data.frame(do.call(cbind, z))
  colnames(cortb2t) <- paste0("lv", 1:nlv)
  
  
  return(list(explvarx = explvarx, contr_block = contr_block, explX = explX, corx2t = corx2t, cortb2t = cortb2t))#, rv = zrv, lg = zlg))
}
