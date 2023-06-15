
vip <- function(object, X, Y, nlv = NULL){#VIP calculation as in the `mixOmics' package.
  
  if (class(object)%in%c("Plsrda","Plsprobda")){
    if(is.numeric(Y)==FALSE){Y <- dummy(as.character(Y))$Y}
    Tx <- object$fm[[1]]$T
    if (sum(rownames(Tx)!=rownames(X))>0){stop("rownames of Tx and X are different")}
    if(is.null(nlv)==TRUE){nlv <- ncol(Tx)}
    p <- ncol(X)
    q <- length(object$fm[[1]]$lev)
    W <- object$fm[[1]]$W[,1:nlv, drop = FALSE]
    if (sum(rownames(W)!=colnames(X))>0){stop("rownames of W and colnames of X are different")}
  }
  
  if (class(object)%in%c("Plsr","Pls")){
    Tx <- object$T
    if (sum(rownames(Tx)!=rownames(X))>0){stop("rownames Tx and X are different")}
    if(is.null(nlv)==TRUE){nlv <- ncol(Tx)}
    p <- ncol(X)
    q <- ncol(Y)
    W <- object$W[,1:nlv, drop = FALSE]
    if (sum(rownames(W)!=colnames(X))>0){stop("rownames of W and colnames of X are different")}
  }

    VIP <- matrix(0, nrow = p, ncol = nlv)
    
    cor2 <- cor(Y, Tx, use = "pairwise")^2
    cor2 <- as.matrix(cor2, nrow = q)
    
    VIP[, 1] <- W[, 1, drop = FALSE]^2
    
    if (nlv > 1)
    {
      for (h in 2:nlv)
      {
        if (q == 1)
        {
          Rd <- cor2[, 1:h, drop = FALSE] 
          VIP[, h] <- Rd %*% t(W[, 1:h, drop = FALSE]^2) / sum(Rd)
        } else {
          Rd <- apply(cor2[, 1:h, drop = FALSE], 2, sum)
          VIP[, h] <- Rd %*% t(W[, 1:h, drop = FALSE]^2) / sum(Rd)
        }
      }
    }
    
    # output
    VIP <- sqrt(p * VIP)
    rownames(VIP) <- rownames(W)
    colnames(VIP) <- paste("comp", 1:nlv)
    
    return(VIP)
}


#VIP.R: Implementation of VIP (variable importance in projection)(*) for the `pls' package.
# VIP <- function(object) {
#   if (object$method != "oscorespls")
#     stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
#   if (nrow(object$Yloadings) > 1)
#     stop("Only implemented for single-response models")
#   
#   SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
#   Wnorm2 <- colSums(object$loading.weights^2)
#   SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*")
#   sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
# }