nipals <- function(X, tol = 1e-6, maxit = 100) {
  n <- nrow(X)
  p <- ncol(X)
  
  # Initialisation
  u <- X[, which.max(sapply(1:ncol(X), function(i) norm(X[,i, drop=FALSE], type = "1")))]  # Trouver le maximum de la norme des colonnes
  v <- rep(0, p)
  v0 <- rep(0, p)
  cont <- TRUE
  iter <- 1
  
  while (cont) {
    v0 <- v  # Sauvegarde de la derniere valeur de v
    v <- crossprod(X, u)  # v = X' * u
    v <- v / sqrt(sum(v^2))  # Normalisation de v
    u <- X %*% v  # u = X * v
    dif <- sum((v - v0)^2)  # Calcul de la difference entre v et v0
    
    iter <- iter + 1
    if (dif < tol || iter > maxit) {
      cont <- FALSE  # Condition d'arret
    }
  }
  
  sv <- sqrt(sum(u^2)) # Valeur singuliere
  u <- u / sv  # Normalisation de u
  list(u = u, v = v, sv = sv, niter = iter - 1)
}

