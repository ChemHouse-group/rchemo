soplsrda <- function(Xlist, y, Xscaling = c("none", "pareto", "sd")[1], Yscaling = c("none", "pareto", "sd")[1], weights = NULL, nlv) {

  if(is.factor(y)){y <- as.character(y)}

  Xlist <- lapply(1:length(Xlist), function(x) .mat(Xlist[[x]]))
  
  n <- nrow(Xlist[[1]])
  p <- sapply(1:length(Xlist), function(x) ncol(Xlist[[x]]))
  nlv <- sapply(1:length(nlv), function(x) min(nlv[x], n, p[x]))
  if(is.null(weights))
    weights <- rep(1, n)
  weights <- .mweights(weights)
  z <- dummy(y)
  fm <- soplsr(Xlist, z$Y, Xscaling = Xscaling, Yscaling = Yscaling, weights = weights, nlv = nlv)
  structure(
    list(fm = fm, lev = z$lev, ni = z$ni),
    class = c("Soplsrda"))       
}

transform.Soplsrda <- function(object, X, ...){
  
  X <- lapply(1:length(X), function(i) .mat(X[[i]]))
  nbl <- length(object$fm$fm)
  if(object$fm$nlv[1]>0){
    T <- transform(object$fm$fm[[1]], X[[1]])
  }else{
    T <- NULL
  }
  
  if (nbl > 1){
    for (i in 2:nbl){
      if(object$fm$nlv[i]>0){
        Xo = X[[i]] - T %*% object$fm$b[[i]]
        T = cbind(T, transform(object$fm$fm[[i]], Xo))
      }
    }
  }
  T
}

predict.Soplsrda <- function(object, X, ...) {
  X <- lapply(1:length(X), function(i) .mat(X[[i]]))
  rownam <- row.names(X[[1]])
  colnam <- "y1"

  posterior <- predict(object$fm, X)
  dimnames(posterior) <- list(rownam, object$lev)
  
  z <- apply(posterior, FUN = .findmax, MARGIN = 1)
  pred <- matrix(.replace_bylev(z, object$lev), ncol = 1)
  dimnames(pred) <- list(rownam, colnam)
 
  list(pred = pred, posterior = posterior)
}
