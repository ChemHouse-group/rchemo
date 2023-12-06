plsnipals <- function(X, Y, scaling = "Centered", weights = NULL, nlv) {
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
    xsds <- sqrt(.colvars(X, weights = weights)*nrow(X)/(nrow(X)-1))
    ysds <- sqrt(.colvars(Y, weights = weights)*nrow(Y)/(nrow(Y)-1))
    # xsds <- sqrt(apply(X, MARGIN = 2, FUN= var))
    # ysds <- sqrt(apply(Y, MARGIN = 2, FUN= var))
    
    if(scaling == "Centered"){
      X <- .center(X, xmeans)
      Y <- .center(Y, ymeans)
    }
    if(scaling == "Pareto"){
      X <- .center(X, xmeans)
      X <- scale(X, center = FALSE, scale = sqrt(xsds))
      Y <- .center(Y, ymeans)
      Y <- scale(Y, center = FALSE, scale = sqrt(ysds))
    }
    if(scaling == "CtReduced"){
      X <- .center(X, xmeans)
      X <- scale(X, center = FALSE, scale = xsds)
      Y <- .center(Y, ymeans)
      Y <- scale(Y, center = FALSE, scale = ysds)
    }
    
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
             xmeans = xmeans, ymeans = ymeans, xsds = xsds, ysds = ysds, weights = weights, scaling = scaling, U = NULL),
        class = c("Plsr", "Pls"))
}
