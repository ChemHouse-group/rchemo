plsnipals <- function(X, Y, weights = NULL, nlv, Xscaling = "none", Yscaling = "none") {
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
             xmeans = xmeans, ymeans = ymeans, xscales = xscales, yscales = yscales, weights = weights,  U = NULL),
        class = c("Plsr", "Pls"))
}
