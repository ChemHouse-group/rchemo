plskern <- function(X, Y, weights = NULL, nlv, Xscaling = c("none", "pareto", "sd")[1], Yscaling = c("none", "pareto", "sd")[1]) {
    X <- .mat(X)
    Y <- .mat(Y, "y")     
    zdim <- dim(X)
    n <- zdim[1]
    p <- zdim[2]
    q <- dim(Y)[2]
    if(is.null(weights))
        weights <- rep(1, n)
    weights <- .mweights(weights)
    xmeans <- .colmeans(X, weights = weights)
    ymeans <- .colmeans(Y, weights = weights) 
    
    if(Xscaling == "none") {xscales <- rep(1, p)}
    if(Yscaling == "none") {yscales <- rep(1, q)}
    if(Xscaling == "sd") {xscales <- sqrt(.colvars(X, weights = weights))}
    if(Yscaling == "sd") {yscales <- sqrt(.colvars(Y, weights = weights))}
    if(Xscaling == "pareto") {xscales <- sqrt(sqrt(.colvars(X, weights = weights)))}
    if(Yscaling == "pareto") {yscales <- sqrt(sqrt(.colvars(Y, weights = weights)))}
    
    X <- scale(X, center = xmeans, scale = xscales)
    Y <- scale(Y, center = ymeans, scale = yscales)
    
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
        if(q == 1) 
            w <- tXY
        else {
            u <- svd(t(tXY), nu = 1, nv = 0)$u
            ## Same as
            ## u <- svd(tXY, nu = 0, nv = 1)$v
            ## u <- eigen(crossprod(tXY), symmetric = TRUE)$vectors[, 1]
            ## End
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
            xmeans = xmeans, ymeans = ymeans, xscales = xscales, yscales = yscales, 
            weights = weights, U = NULL),
        class = c("Plsr"))
}

transform.Plsr <- function(object, X, ..., nlv = NULL) {
    a <- dim(object$T)[2]
    if(is.null(nlv))
        nlv <- a
    else 
        nlv <- min(a, nlv)
    X <- scale(.mat(X), center = object$xmeans, scale = object$xscales)
    T <- X %*% object$R[, seq_len(nlv), drop = FALSE]
    colnames(T) <- paste("lv", seq_len(dim(T)[2]), sep = "")
    T
}

coef.Plsr <- function(object, ..., nlv = NULL) {
    ## Works also for nlv = 0
    a <- dim(object$T)[2]
    if(is.null(nlv))
        nlv <- a
    else 
        nlv <- min(a, nlv)
    znlv <- seq_len(nlv)
    p <- nrow(object$R)
    q <- nrow(object$C)
    beta <- t(object$C)[znlv, , drop = FALSE]
    
    K1 <- matrix(rep(object$yscales, each = p), ncol = q) 
    K2 <- t(matrix(rep(object$xscales, each = q), ncol = p))
    K <- K1 / K2
    B <- object$R[, znlv, drop = FALSE] %*% beta * K
    
    int <- object$ymeans - t(object$xmeans) %*% B
    list(int = int, B = B) 
}

predict.Plsr <- function(object, X, ..., nlv = NULL) {
    X <- .mat(X)
    q <- length(object$ymeans)
    rownam <- row.names(X)
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
        zpred <- t(c(z$int) + t(X %*% z$B))
        ## Same but faster than:
        ## zpred <- cbind(rep(1, m), X) %*% rbind(z$int, z$B)
        dimnames(zpred) <- list(rownam, colnam)
        pred[[i]] <- zpred
    }
    names(pred) <- paste("lv", nlv, sep = "")
    if(le_nlv == 1)
        pred <- pred[[1]] 
    list(pred = pred)
}



summary.Plsr <- function(object, X, ...) {
    zdim <- dim(object$T)
    n <- zdim[1]
    nlv <- zdim[2]
    X <- scale(.mat(X), center = object$xmeans, scale = object$xscales)
    sstot <- sum(object$weights * X * X, na.rm = TRUE)
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


