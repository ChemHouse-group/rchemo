mbplsr <- function(Xlist, Y, blockscaling = TRUE, weights = NULL, nlv, Xscaling = c("none", "pareto", "sd")[1], Yscaling = c("none", "pareto", "sd")[1]) {
    Xlist <- lapply(1:length(Xlist), function(i) .mat(Xlist[[i]]))
    Y <- .mat(Y, "y")
    
    if(is.null(weights))
        weights <- rep(1, nrow(Y))
    weights <- .mweights(weights)
    
    xmeanslist <- lapply(1:length(Xlist), function(i) .colmeans(Xlist[[i]], weights = weights))
    ymeans     <- .colmeans(Y, weights = weights) 
    
    if((length(Xscaling)=1) & (length(Xlist)>1)){Xscaling = rep(Xscaling, length(Xlist))}
    
    xscaleslist <- list()
    for(i in 1:length(Xlist)){
      if(Xscaling[i] == "none"){xscaleslist[[i]] <- rep(1, ncol(Xlist[[i]]))}
      if(Xscaling[i] == "pareto"){xscaleslist[[i]] <- sqrt(sqrt(.colvars(Xlist[[i]], weights = weights)))}
      if(Xscaling[i] == "sd"){xscaleslist[[i]] <- sqrt(.colvars(Xlist[[i]], weights = weights))}
    }
    Xlist <- lapply(1:length(Xlist), function(i) scale(Xlist[[i]], center = xmeanslist[[i]], scale = xscaleslist[[i]]))
    
    if(blockscaling==TRUE){
      Xblockscaled <- blockscal(Xtrain = Xlist, weights = weights)
      Xlist <- Xblockscaled$Xtrain
      Xnorms <- Xblockscaled$disp
    }else{
      Xnorms <- NA
    }
    
    if(Yscaling == "none"){yscales <- rep(1, ncol(Y))}
    if(Yscaling == "pareto"){yscales <- sqrt(sqrt(.colvars(Y, weights = weights)))}
    if(Yscaling == "sd"){yscales <- sqrt(.colvars(Y, weights = weights))}
    Y <- scale(Y, center = ymeans, scale = yscales)
    
    Xconc <- do.call("cbind",Xlist)
    zdim <- dim(Xconc)
    n <- zdim[1]
    p <- zdim[2]
    q <- dim(Y)[2]
    nam <- paste("lv", seq_len(nlv), sep = "")
    T <- matrix(nrow = n, ncol = nlv, dimnames = list(row.names(Xconc), nam))                     
    R <- W <- P <- matrix(nrow = p, ncol = nlv, dimnames = list(colnames(Xconc), nam)) 
    C <- matrix(nrow = q, ncol = nlv, dimnames = list(colnames(Y), nam))                     
    TT <- vector(length = nlv)
    Xd <- weights * Xconc
    # = D %*% Xconc = d * Xconc = Xconc * d
    tXY <- crossprod(Xd, Y)
    # = t(D %*% Xconc) %*% Y = t(Xconc) %*% D %*% Y
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
        t <- Xconc %*% r 
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
             xmeans = xmeanslist, ymeans = ymeans, xscales = xscaleslist, yscales = yscales, weights = weights, blockscaling = blockscaling, Xnorms = Xnorms, U = NULL),
        class = c("Mbplsr"))
}

summary.Mbplsr <- function(object, X, ...) {
    zdim <- dim(object$T)
    n <- zdim[1]
    nlv <- zdim[2]
    
    X <- lapply(1:length(X), function(i) scale(.mat(X[[i]]), center = object$xmeans[[i]], scale = object$xscales[[i]]))
    
    if(object$blockscaling==TRUE){X <- blockscal(Xtrain = X, weights = object$weights)$Xtrain}

    Xconc <- do.call("cbind",X)
    
    sstot <- sum(object$weights * Xconc * Xconc, na.rm = TRUE)
    tt <- object$TT
    ## Only valid if scores T are orthogonal (or approximate)
    tt.adj <- colSums(object$P * object$P) * tt
    pvar <- tt.adj / sstot
    cumpvar <- cumsum(pvar)
    xvar <- tt.adj / n
    explvar <- data.frame(nlv = seq(nlv), var = xvar, pvar = pvar, cumpvar = cumpvar)
    row.names(explvar) <- seq(nlv)
    list(explvarx = explvar)
}

transform.Mbplsr <- function(object, X, ..., nlv = NULL) {
    a <- dim(object$T)[2]
    if(is.null(nlv)){
      nlv <- a
    }else{
      nlv <- min(a, nlv)
    }

    X <- lapply(1:length(X), function(i) scale(.mat(X[[i]]), center = object$xmeans[[i]], scale = object$xscales[[i]]))
    
    if(object$blockscaling==TRUE){X <- lapply(1:length(X), function(i) X[[i]]/object$Xnorms[i])}
    
    T <- do.call("cbind",X) %*% object$R[, seq_len(nlv), drop = FALSE]
    colnames(T) <- paste("lv", seq_len(dim(T)[2]), sep = "")
    T
}

coef.Mbplsr <- function(object, ..., nlv = NULL) {
    ## Works also for nlv = 0
    a <- dim(object$T)[2]
    if(is.null(nlv)){
      nlv <- a
    }else{
      nlv <- min(a, nlv)
    }
        
    p <- nrow(object$R)
    q <- nrow(object$C)
    beta <- t(object$C)[seq_len(nlv), , drop = FALSE]
    
    plist <- lapply(1:length(object$xmeans), function(P) length(object$xmeans[[P]]))
    cumsumPlist <- cumsum(plist)
    cumsumPlist1 <- c(1, cumsumPlist + 1)
    
    Blist <- list()
    for(i in 1:length(object$xmeans)){
      Blist[[i]] <- object$R[cumsumPlist1[i]:cumsumPlist[i], seq_len(nlv), drop = FALSE] %*% beta
      if(object$blockscaling == TRUE){
        Blist[[i]] <- Blist[[i]]/object$Xnorms[i]
      }
      Blist[[i]] <- Blist[[i]] * matrix(rep(object$yscales, each = plist[[i]]), ncol = q) / t(matrix(rep(object$xscales[[i]], each = q), ncol = plist[[i]]))
    }
    
    B <- do.call("rbind",Blist)
    int <- object$ymeans - t(do.call("c",object$xmeans)) %*% B
    list(int = int, B = B) 
}

predict.Mbplsr <- function(object, X, ..., nlv = NULL) {
    X <- lapply(1:length(X), function(i) .mat(X[[i]]))
    Xconc <- do.call("cbind",X)
    q <- length(object$ymeans)
    rownam <- row.names(Xconc)
    colnam <- paste("y", seq_len(q), sep = "")
    a <- dim(object$T)[2]
    if(is.null(nlv))
        nlv <- a 
    else 
        nlv <- seq(min(nlv), min(max(nlv), a))
    le_nlv <- length(nlv)
    pred <- vector(mode = "list", length = le_nlv)
    for(i in seq_len(le_nlv)) {
        z <- coef(object, nlv = nlv[i])
        zpred <- t(c(z$int) + t(Xconc %*% z$B))
        ## Same but faster than:
        ## zpred <- cbind(rep(1, m), Xconc) %*% rbind(z$int, z$B)
        dimnames(zpred) <- list(rownam, colnam)
        pred[[i]] <- zpred
        }
    names(pred) <- paste("lv", nlv, sep = "")
    if(le_nlv == 1)
        pred <- pred[[1]] 
    list(pred = pred)
}
