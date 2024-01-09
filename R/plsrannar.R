plsrannar <- function(X, Y, scaling = "centered", weights = NULL, nlv) {
    X <- .mat(X)
    Y <- .mat(Y, "y")     
    n <- dim(X)[1]
    q <- dim(Y)[2]
    if(is.null(weights))
        weights <- rep(1, n)
    weights <- .mweights(weights)
    xmeans <- .colmeans(X, weights = weights) 
    ymeans <- .colmeans(Y, weights = weights)
    xsds <- sqrt(.colvars(X, weights = weights))#*nrow(X)/(nrow(X)-1))
    ysds <- sqrt(.colvars(Y, weights = weights))#*nrow(Y)/(nrow(Y)-1))

    if(scaling == "centered"){
      X <- .center(X, xmeans)
      Y <- .center(Y, ymeans)
    }
    if(scaling == "pareto"){
      X <- .center(X, xmeans)
      X <- scale(X, center = FALSE, scale = sqrt(xsds))
      Y <- .center(Y, ymeans)
      Y <- scale(Y, center = FALSE, scale = sqrt(ysds))
    }
    if(scaling == "ctreduced"){
      X <- .center(X, xmeans)
      X <- scale(X, center = FALSE, scale = xsds)
      Y <- .center(Y, ymeans)
      Y <- scale(Y, center = FALSE, scale = ysds)
    }
    
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
            xmeans = xmeans, ymeans = ymeans, xsds = xsds, ysds = ysds, weights = weights, scaling = scaling, U = U),
        class = c("Plsr", "Pls"))    
    }
