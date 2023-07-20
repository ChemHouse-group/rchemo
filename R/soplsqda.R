soplsqda <- function(X_list, Y, nlv_vect, blocknames = NULL, weights = NULL, ...) {

  if(is.numeric(Y)==FALSE){Y <- dummy(as.character(Y))$Y}
  if(is.matrix(Y)==FALSE){Y <- as.matrix(Y)}
  if(sum(sapply(1:length(Y), function(x) Y[x]%in%c(0,1)==FALSE))>0){Y <- dummy(as.character(Y))$Y}
  
  # model
  
  soplsmodel <- sopls(X_list, Y, nlv_vect, blocknames, weights) 
  
  # classif
    
  # outputs
  fm <- list()
  fm$W <- soplsmodel$W
  fm$T <- soplsmodel$Tx
  fm$C <- soplsmodel$C
  fm$BCoef <- soplsmodel$BCoef
  fm$VIP   <- soplsmodel$VIP
  fm$Ymeans   <- soplsmodel$Ymeans
  #names(fm$W) <- names(fm$T) <- names(fm$C) <- names(fm$BCoef) <- names(fm$VIP) <- names(fm$Ymeans) <- blocknames
  fm$Fit <- soplsmodel$Fit
  fm$lev <- colnames(Y)
  fm$ni <- apply(Y, MARGIN = 2, FUN = sum)
  names(fm$ni) <- fm$lev
  fm$call  <- match.call()
  
  fm
  
}

