.mbplsrannar <- function(Xlist, Y, scaling = c("none", "pareto", "sd")[1], blockscaling = TRUE, weights = NULL, nlv) {
    Xlist <- lapply(1:length(Xlist), function(X) .mat(Xlist[[X]]))
    Y <- .mat(Y, "y")  
    
    if(is.null(weights))
        weights <- rep(1, n)
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
    #zdim <- dim(X)
    n <- dim(X)[1]
    q <- dim(Y)[2]

    nam <- paste("lv", seq_len(nlv), sep = "")
    U <- T <- Tclass <- matrix(nrow = n, ncol = nlv, 
                               dimnames = list(row.names(X), nam))                     
    TT <- vector(length = nlv)
    Xd <- sqrt(weights) * X
    Yd <- sqrt(weights) * Y
    XtX <- tcrossprod(Xd)
    YtY <- tcrossprod(Yd)
    XY <- XtX %*% YtY    
    I <- diag(n)
    for(a in seq_len(nlv)) {
        t <- .eigpow(XY)$v
        u <- YtY %*% t
        utemp <-    u / sum(u * t)
        wtw <- c(crossprod(utemp, XtX) %*% utemp)
        tclass <- t * sqrt(wtw) / sqrt(weights)
        tt <- sum(weights * tclass * tclass)    
        G <- I - tcrossprod(t)
        XtX <- G %*% (XtX) %*% G 
        YtY <- G %*% YtY %*% G
        XY <- XtX %*% YtY
        T[, a] <- t
        Tclass[, a] <- tclass
        U[, a] <- u
        TT[a] <- tt
        }
    W <- crossprod(Xd, U)
    W <- .scale(W, scale = .colnorms(W))
    Z <- solve(crossprod(T))
    Z <- .scale(Z, scale = sqrt(TT))
    P <- crossprod(Xd, T) %*% Z
    C <- crossprod(Yd, T) %*% Z
    R <- W %*% solve(crossprod(P, W))
    structure(
        list(T = Tclass, P = P, R = R, W = W, C = C, TT = TT,
            xmeans = xmeanslist, ymeans = ymeans, xsds = xsdslist, ysds = ysds, weights = weights, scaling = scaling, blockscaling = blockscaling, Xnorms = Xnorms, U = U),
        class = c("Mbplsr", "Mbpls"))    
    }
