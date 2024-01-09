mbplsr <- function(Xlist, Y, scaling = c("centered", "pareto", "ctreduced")[1], blockscaling = TRUE, weights = NULL, nlv) {
    Xlist <- lapply(1:length(Xlist), function(X) .mat(Xlist[[X]]))
    Y <- .mat(Y, "y")
    
    if(is.null(weights))
        weights <- rep(1, nrow(Y))
    weights <- .mweights(weights)
    
    xmeanslist <- lapply(1:length(Xlist), function(X) .colmeans(Xlist[[X]], weights = weights))
    ymeans     <- .colmeans(Y, weights = weights) 
    xsdslist   <- lapply(1:length(Xlist), function(X) sqrt(.colvars(Xlist[[X]], weights = weights)))#*nrow(Xlist[[X]])/(nrow(Xlist[[X]])-1)))
    ysds       <- sqrt(.colvars(Y, weights = weights))#*nrow(Y)/(nrow(Y)-1))
   
    Y <- .center(Y, ymeans)
    
    if((length(scaling)=1) & (length(Xlist)>1)){scaling= rep(scaling, length(Xlist))}
    
    for(i in 1:length(Xlist)){
      if(scaling[i] == "centered"){
        Xlist[[i]] <- .center(Xlist[[i]], xmeanslist[[i]])
        # Y <- .center(Y, ymeans)
      }
      if(scaling[i] == "pareto"){
        Xlist[[i]] <- .center(Xlist[[i]], xmeanslist[[i]])
        Xlist[[i]] <- scale(Xlist[[i]], center = FALSE, scale = sqrt(xsdslist[[i]]))
        # Y <- .center(Y, ymeans)
        # Y <- scale(Y, center = FALSE, scale = sqrt(ysds))
      }
      if(scaling[i] == "ctreduced"){
        Xlist[[i]] <- .center(Xlist[[i]], xmeanslist[[i]])
        Xlist[[i]] <- scale(Xlist[[i]], center = FALSE, scale = xsdslist[[i]])
        # Y <- .center(Y, ymeans)
        # Y <- scale(Y, center = FALSE, scale = ysds)
      }
    }
    
    if(blockscaling==TRUE){
      Xscaling <- blockscal(Xtrain = Xlist, weights = weights)
      Xlist <- Xscaling$Xtrain
      Xnorms <- Xscaling$disp
    }else{
      Xnorms <- NA
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
    Xd <- weights * X
    # = D %*% X = d * X = X * d
    tXY <- crossprod(Xd, Y)
    # = t(D %*% X) %*% Y = t(X) %*% D %*% Y
    for(a in seq_len(nlv)) {
        if(q == 1) w <- tXY
            else {
                u <- svd(t(tXY), nu = 1, nv = 0)$u
                ## Same as
                ## u <- svd(tXY, nu = 0, nv = 1)$v
                ## u <- eigen(crossprod(tXY), symmetric = TRUE)$vectors[, 1]
                w <- tXY %*% u
            } 
        w <- w / sqrt(sum(w * w))
        r <- w
        if(a > 1)
            for(j in seq_len(a - 1)) 
                    r <- r - sum(P[, j] * w) * R[, j]
        t <- X %*% r 
        tt <- sum(weights * t * t)         
        c <- crossprod(tXY, r) / tt
        zp <- crossprod(Xd, t) / tt 
        tXY <- tXY - tcrossprod(zp, c) * tt    
        T[, a] <- t
        P[, a] <- zp
        W[, a] <- w
        R[, a] <- r
        C[, a] <- c
        TT[a] <- tt
    }
    structure(
        list(T = T, P = P, R = R, W = W, C = C, TT = TT,
             xmeans = xmeanslist, ymeans = ymeans, xsds = xsdslist, ysds = ysds, weights = weights, scaling = scaling, blockscaling = blockscaling, Xnorms = Xnorms, U = NULL),
        class = c("Mbplsr"))
}

