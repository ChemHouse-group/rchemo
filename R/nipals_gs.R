#nipals_gram_schmidt
nipals_gs <- function(X, UUt, VVt, tol = 1e-6, maxit = 100) {
  n <- nrow(X)
  p <- ncol(X)
  
  # Initialisation
  u <- X[, which.max(sapply(1:ncol(X), function(i) norm(X[, i], type = "2")))]  # Trouver le maximum de la norme des colonnes
  v <- rep(0, p)
  v0 <- rep(0, p)
  cont <- TRUE
  iter <- 1
  
  while (cont) {
    v0 <- v  # Sauvegarde de la dernière valeur de v
    v <- crossprod(X, u)  # v = X' * u
    v <- v - VVt %*% v  # Orthogonalisation de v
    v <- v / norm(v, type = "2")  # Normalisation de v
    u <- X %*% v  # u = X * v
    u <- u - UUt %*% u  # Orthogonalisation de u
    dif <- sum((v - v0)^2)  # Calcul de la différence entre v et v0
    
    iter <- iter + 1
    if (dif < tol || iter > maxit) {
      cont <- FALSE  # Condition d'arrêt
    }
  }
  
  sv <- norm(u, type = "2")  # Valeur singulière
  u <- u / sv  # Normalisation de u
  list(u = u, v = v, sv = sv, niter = iter - 1)
}

# Exemple avec une matrice X aléatoire
set.seed(123)
X <- matrix(rnorm(15), nrow = 5, ncol = 3)

# Appel de la fonction nipals sans Gram-Schmidt
res <- nipals(X)
print(res$niter)  # Nombre d'itérations
print(res$sv)     # Valeur singulière
print(res$v)      # Vecteur de charge (loadings)
print(res$u)      # Vecteur de score (scores)
