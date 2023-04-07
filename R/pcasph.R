pcasph <- function(X, weights = NULL, nlv) {
    X <- .mat(X)
    zdim <- dim(X)
    n <- zdim[1]
    p <- zdim[2]
    nlv <- min(nlv, n, p)
    if(is.null(weights))
        weights <- rep(1, n)
    weights <- .mweights(weights)
    xmeans <- .colmeds_spa(X, delta = .001)
    X <- .center(X, xmeans)
    tX <- t(X)
    xnorms <- .colnorms(tX)
    tX <- .scale(tX, center = rep(0, n), xnorms)
    zX <- t(tX)
    z <- svd(sqrt(weights) * zX, nu = 0, nv = nlv)
    P <- z$v
    zT <- zX %*% P
    T <- X %*% P
    sv <- apply(zT, FUN = mad, MARGIN = 2) 
    u <- rev(order(sv))
    P <- P[, u, drop = FALSE]
    T <- T[, u, drop = FALSE]
    sv <- sv[u]    
    eig <- sv^2 
    row.names(T) <- row.names(X)
    row.names(P) <- colnames(X)
    colnames(T) <- colnames(P) <- paste("pc", seq_len(nlv), sep = "")
    structure(
        list(T = T, P = P, sv = sv, eig = eig,
             xmeans = xmeans, weights = weights, niter = NULL, conv = NULL),
        class = c("Pca"))
    }

