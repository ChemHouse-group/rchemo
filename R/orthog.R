orthog <- function(X, Y, weights = NULL) {
  
  # Y is orthogonalized to X
  
  X <- .mat(X)
  n <- dim(X)[1]
  
  Y <- .mat(Y)
  
  if(is.null(weights))
    weights <- rep(1 / n, n)
  else
    weights <- weights / sum(weights)
  
  fm <- lm(Y ~ X - 1, weights = weights)
  
  b <- coef(fm)
  
  if(ncol(Y) > 1)
    Yortho <- fm$residuals
  else
    Yortho <- matrix(fm$residuals, ncol = 1, 
                     dimnames = list(row.names(Y), colnames(Y)))
  
  list(Y = Yortho, b = b, weights = weights)
  
}
