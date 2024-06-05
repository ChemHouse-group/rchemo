
vip <- function(object, X, Y = NULL, nlv = NULL){
  # Y == NULL => VIP calculation based on the proportion of Y-variance explained by the components
  # Y != NULL => VIP calculation as proposed by Tenenhaus (1998), or implemented in the `mixOmics' package.
  #               the explained variance is replaced by the redundancy.
  
  if(is.null(Y)==FALSE){
    if (class(object)[[1]]%in%c("Plsrda")){
      if(is.numeric(Y)==FALSE){Y <- dummy(as.character(Y))$Y}
      if(is.matrix(Y)==FALSE){Y <- as.matrix(Y)}
      if(sum(sapply(1:length(Y), function(x) Y[x]%in%c(0,1)==FALSE))>0){Y <- dummy(as.character(Y))$Y}
      Tx <- object$fm$T
      if (sum(rownames(Tx)!=rownames(X))>0){stop("rownames of Tx and X are different")}
      if(is.null(nlv)==TRUE){nlv <- ncol(Tx)}
      p <- ncol(X)
      q <- length(object$fm$lev)
      W <- object$fm$W[,1:nlv, drop = FALSE] 
      if (sum(rownames(W)!=colnames(X))>0){stop("rownames of W and colnames of X are different")}
    }
    
    if (class(object)[[1]]%in%c("Plsprobda")){
      if(is.numeric(Y)==FALSE){Y <- dummy(as.character(Y))$Y}
      if(is.matrix(Y)==FALSE){Y <- as.matrix(Y)}
      if(sum(sapply(1:length(Y), function(x) Y[x]%in%c(0,1)==FALSE))>0){Y <- dummy(as.character(Y))$Y}
      Tx <- object$fm[[1]]$T
      if (sum(rownames(Tx)!=rownames(X))>0){stop("rownames of Tx and X are different")}
      if(is.null(nlv)==TRUE){nlv <- ncol(Tx)}
      p <- ncol(X)
      q <- length(object$fm[[1]]$lev)
      W <- object$fm[[1]]$W[,1:nlv, drop = FALSE]
      if (sum(rownames(W)!=colnames(X))>0){stop("rownames of W and colnames of X are different")}
    }
    
    if (class(object)[[1]]%in%c("Plsr","Pls")){
      Tx <- object$T
      if (sum(rownames(Tx)!=rownames(X))>0){stop("rownames Tx and X are different")}
      if(is.matrix(Y)==FALSE){Y <- as.matrix(Y)}
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
    
    if (nlv > 1){
      for (h in 2:nlv){
        if (q == 1){
          Rd <- cor2[, 1:h, drop = FALSE] 
          VIP[, h] <- Rd %*% t(W[, 1:h, drop = FALSE]^2) / sum(Rd)
        }else{
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
  
  if(is.null(Y)==TRUE){
    if (class(object)[[1]]%in%c("Plsrda")){
      Tx <- object$fm$T 
      if (sum(rownames(Tx)!=rownames(X))>0){stop("rownames of Tx and X are different")}
      if(is.null(nlv)==TRUE){nlv <- ncol(Tx)}
      p <- ncol(X)
      #q <- length(object$fm$lev) 
      W <- object$fm$W[,1:nlv, drop = FALSE] 
      if (sum(rownames(W)!=colnames(X))>0){stop("rownames of W and colnames of X are different")}
      Q <- object$fm$C[,1:nlv, drop = FALSE] 
    }
    
    if (class(object)[[1]]%in%c("Plsprobda")){
      Tx <- object$fm[[1]]$T
      if (sum(rownames(Tx)!=rownames(X))>0){stop("rownames of Tx and X are different")}
      if(is.null(nlv)==TRUE){nlv <- ncol(Tx)}
      p <- ncol(X)
      #q <- length(object$fm[[1]]$lev)
      W <- object$fm[[1]]$W[,1:nlv, drop = FALSE]
      if (sum(rownames(W)!=colnames(X))>0){stop("rownames of W and colnames of X are different")}
      Q <- object$fm[[1]]$C[,1:nlv, drop = FALSE]
    }
    
    if (class(object)[[1]]%in%c("Plsr","Pls")){
      Tx <- object$T # object scores
      if (sum(rownames(Tx)!=rownames(X))>0){stop("rownames Tx and X are different")}
      if(is.null(nlv)==TRUE){nlv <- ncol(Tx)}
      p <- ncol(X)
      W <- object$W[,1:nlv, drop = FALSE]# object loading.weights
      if (sum(rownames(W)!=colnames(X))>0){stop("rownames of W and colnames of X are different")}
      Q <- object$C[,1:nlv, drop = FALSE]# Y loading weights
    }

    WW <- W * W/apply(W, 2, function(x) sum(x * x))
    
    if(length(dim(Q)) == 0){
      Q2 <- as.numeric(Q) * as.numeric(Q)
    } else {
      Q2 <- rowSums(t(Q * Q))
    }
    
    VIP <- matrix(0, nrow = p, ncol = nlv)
    for(i in 1:nlv){
      Q2TT <- Q2[1:i] * diag(crossprod(Tx))[i]
      VIP[,i] <- sqrt(p * apply(sweep(WW[, 1:i, drop=FALSE],2,Q2TT,"*"), 1, sum)/sum(Q2TT))
    }
    rownames(VIP) <- rownames(W)
    colnames(VIP) <- paste("comp", 1:nlv)
    
    return(VIP)
  }
}
