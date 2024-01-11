.mbplsnipals <- function(Xlist, Y, scaling = c("none", "pareto", "sd")[1], blockscaling = TRUE, weights = NULL, nlv) {
    Xlist <- lapply(1:length(Xlist), function(X) .mat(Xlist[[X]]))
    Y <- .mat(Y, "y") 
    
    if(is.null(weights))
      weights <- rep(1, nrow(Y))
    weights <- .mweights(weights)
    
    xmeanslist <- lapply(1:length(Xlist), function(X) .colmeans(Xlist[[X]], weights = weights))
    ymeans     <- .colmeans(Y, weights = weights) 
    xsdslist   <- lapply(1:length(Xlist), function(X) sqrt(.colvars(Xlist[[X]], weights = weights)))#*nrow(Xlist[[X]])/(nrow(Xlist[[X]])-1)))
    ysds       <- sqrt(.colvars(Y, weights = weights))#*nrow(Y)/(nrow(Y)-1))
    
    if((length(Xscaling)=1) & (length(Xlist)>1)){Xscaling = rep(Xscaling, length(Xlist))}
    
    for(i in 1:length(Xlist)){
      if(Xscaling[i] == "none"){
        Xlist[[i]] <- .center(Xlist[[i]], xmeanslist[[i]])
      }
      if(Xscaling[i] == "pareto"){
        Xlist[[i]] <- .center(Xlist[[i]], xmeanslist[[i]])
        Xlist[[i]] <- scale(Xlist[[i]], center = FALSE, scale = sqrt(xsdslist[[i]]))
      }
      if(Xscaling[i] == "sd"){
        Xlist[[i]] <- .center(Xlist[[i]], xmeanslist[[i]])
        Xlist[[i]] <- scale(Xlist[[i]], center = FALSE, scale = xsdslist[[i]])
      }
    }
    
    if(blockscaling==TRUE){
      Xblockscaled <- blockscal(Xtrain = Xlist, weights = weights)
      Xlist <- Xblockscaled$Xtrain
      Xnorms <- Xblockscaled$disp
    }else{
      Xnorms <- NA
    }
    
    if(Yscaling == "none"){
      Y <- .center(Y, ymeans)
    }
    if(Yscaling == "pareto"){
      Y <- .center(Y, ymeans)
      Y <- scale(Y, center = FALSE, scale = sqrt(ysds))
    }
    if(Yscaling == "sd"){
      Y <- .center(Y, ymeans)
      Y <- scale(Y, center = FALSE, scale = ysds)
    }
    
    X <- do.call("cbind",Xlist)
    
    zdim <- dim(X)
    n <- zdim[1]
    p <- zdim[2]
    q <- dim(Y)[2]
    
    nam <- paste("lv", seq_len(nlv), sep = "")
    
    T <- matrix(nrow = n, ncol = nlv, dimnames = list(row.names(X), nam))                     
    R <- W <- P <- matrix(nrow = p, ncol = nlv, dimnames = list(colnames(X), nam)) 
    C <- matrix(nrow = q, ncol = nlv, dimnames = list(colnames(Y), nam))                     
    TT <- vector(length = nlv)
    for(a in seq_len(nlv)) {
        XtY <- crossprod(weights * X, Y)
        # = t(D %*% X) %*% Y = t(X) %*% D %*% Y
        if(q == 1) {
            w <- XtY
            w <- w / sqrt(sum(w * w))
        }
        else {
            w <- svd(XtY, nu = 1, nv = 0)$u
        }
        t <- X %*% w
        tt <- sum(weights * t * t)                                
        c <- crossprod(weights * Y, t)    / tt
        zp <- crossprod(weights * X, t) / tt
        X <- X - tcrossprod(t, zp)
        Y <- Y - tcrossprod(t, c)
        T[, a] <- t
        W[, a] <- w
        P[, a] <- zp
        C[, a] <- c
        TT[a] <- tt
    }
    R <- W %*% solve(crossprod(P, W))
    structure(
        list(T = T, P = P, R = R, W = W, C = C, TT = TT,
             xmeans = xmeanslist, ymeans = ymeans, xsds = xsdslist, ysds = ysds, weights = weights, scaling = scaling, blockscaling = blockscaling, Xnorms = Xnorms, U = NULL),
        class = c("Mbplsr", "Mbpls"))
}
