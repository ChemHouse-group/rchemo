consensuspcanipals <- function(Xlist, blockscaling = TRUE, weights = NULL, nlv, Xscaling = c("none", "pareto", "sd")[1],
                               tol = .Machine$double.eps^0.5, maxit = 200) {

  if (is.null(weights)) {
    weights <- rep(1, nrow(Xlist[[1]])) 
  }
  
  xmeanslist <- lapply(1:length(Xlist), function(i) .colmeans(Xlist[[i]], weights = weights))
  
  if((length(Xscaling)==1) & (length(Xlist)>1)){Xscaling = rep(Xscaling, length(Xlist))}
  
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

  # Initialisation 
  nbl <- length(Xlist)
  n <- nrow(Xlist[[1]])
  sqrtw <- sqrt(weights)
  
  U <- matrix(0, n, nlv)
  WB <- matrix(0, nbl, nlv)
  Tbl <- lapply(1:nbl, function(k) matrix(0, n, nlv))
  Tb <- lapply(1:nlv, function(k) matrix(0, n, nbl))
  Wbl <- lapply(1:nbl, function(k) matrix(0, ncol(Xlist[[k]]), nlv))
  lb <- matrix(0, nbl, nlv)
  mu <- numeric(nlv)
  niter <- numeric(nlv)
  
  Xinit <- do.call(cbind, Xlist)
    
  # Iterations Nipals
  for (a in 1:nlv) {
    X <- do.call(cbind, Xlist)  
    u <- nipals(X, tol = tol, maxit = maxit)$u 
    iter <- 1
    cont <- TRUE
    
    while (cont) {
      u0 <- u  
      for (k in 1:nbl) {
        wk <- t(Xlist[[k]]) %*% u
        dk <- sqrt(sum(wk^2))
        wk <- wk / dk
        tk <- Xlist[[k]] %*% wk
        Tb[[a]][, k] <- tk
        Tbl[[k]][, a] <- (1/sqrtw) * tk
        Wbl[[k]][, a] <- wk
        lb[k, a] <- dk^2
      }
      res <- nipals(Tb[[a]], tol = tol, maxit = maxit)
      u <- res$u
      w <- res$v
      dif <- sum((u - u0)^2)
      iter <- iter + 1
      if (dif < tol || iter > maxit) {
        cont <- FALSE
      }
    }
    
    niter[a] <- iter - 1
    U[, a] <- u
    WB[, a] <- w 
    mu[a] <- res$sv^2# =sum(lb[, a])  # eigen values
    for (k in 1:nbl) {
      Xlist[[k]] <- Xlist[[k]] - u %*% (t(u) %*% Xlist[[k]])
    }
  }
  
  W <- crossprod(weights * Xinit, U) / sum(weights * U * U)
  W <- W / sqrt(sum(W * W))
  
  # global scores
  T <- (diag(1 / sqrt(weights))) %*% (matrix(rep(sqrt(mu),n),nrow=n, byrow=T) * U)
  
  P <- crossprod(weights * Xinit, T) / sum(weights * T * T)
  P <- P / sqrt(sum(P * P))
  
  structure(
    list(T = T, U = U, 
         P = P, W = W, 
         WB = WB, Tbl = Tbl, Tb = Tb, Wbl = Wbl, lb = lb, mu = mu, 
         xmeans = xmeanslist, xscales = xscaleslist, weights = weights, blockscaling = blockscaling, Xnorms = Xnorms, 
         niter = niter), 
    class = c("Consensuspcanipals"))
}


transform.Consensuspcanipals <- function(object, X, nlv = NULL) {
  nlv <- if (is.null(nlv)) ncol(object$T) else min(nlv, ncol(object$T))
  
  X <- lapply(1:length(X), function(i) .mat(X[[i]]))
  
  # Xinit <- X
  p <- dim(object$P)[1]
  n <- dim(object$T)[1]
  # TT <- object$weights * object$T * object$T
  # tt <- colSums(TT)
  
  X <- lapply(1:length(X), function(i) scale(.mat(X[[i]]), center = object$xmeans[[i]], scale = object$xscales[[i]]))
  if(object$blockscaling==TRUE){X <- blockscal(Xtrain = X, weights = object$weights)$Xtrain}
  # Xconc <- do.call("cbind",X)
  
  U <- matrix(0, nrow(X[[1]]), nlv)
  TB <- matrix(0, nrow(X[[1]]), length(X))
  Tbl <- lapply(1:length(X), function(k) matrix(0, n, nlv))
  u <- tk <- numeric(nrow(X[[1]]))
  
  for (a in 1:nlv) {
    # u <- numeric(nrow(X[[1]]))
    for (k in 1:length(X)) {
      TB[,k] <- Tbl[[k]][,a] <- tk <- X[[k]] %*% object$Wbl[[k]][, a]
      # u <- u + X[[k]] %*% object$Wbl[[k]][, a]
    }
    U[, a] <- u <- (1 / sqrt(object$mu[a])) * TB %*% object$WB[, a] # A VERFIIER
    for (k in 1:length(X)) {
      Px <- sqrt(object$lb[k, a]) %*% t(object$Wbl[[k]][, a])
      X[[k]] <- X[[k]] - (u %*% Px)
    }
  }
  T <- sqrt((object$mu[1:nlv])) * U
  return(list(T = T, Tbl = Tbl))
}


summary.Consensuspcanipals <- function(object, X) {

  X <- lapply(1:length(X), function(i) .mat(X[[i]]))
  Xinit <- X
  
  nbl <- length(X)
  nlv <- ncol(object$T)
  sqrtw <- sqrt(object$weights)
  p <- dim(object$P)[1]
  n <- dim(object$T)[1]
  
  X <- lapply(1:length(X), function(i) scale(.mat(X[[i]]), center = object$xmeans[[i]], scale = object$xscales[[i]]))
  if(object$blockscaling==TRUE){X <- blockscal(Xtrain = X, weights = object$weights)$Xtrain}
  
  # for (k in 1:nbl) {
  #   zX[[k]] <- sqrtw * zX[[k]]
  # }

  Xconc <- do.call(cbind, X)
  
  # Explained_X
  sstot <- rep(0, nbl)
  for (k in 1:nbl) {
    sstot[k] <- sum(object$weights * X[[k]] * X[[k]], na.rm = TRUE)
  }

  tt <- colSums(object$lb)
  pvar <- tt / sum(sstot)
  cumpvar <- cumsum(pvar)
  explvarx <- data.frame(lv = 1:nlv, var = tt, pvar = pvar, cumpvar = cumpvar)
  
  # Contribution of the blocks to global scores
  # z = fscale(object.lb, colsum(object.lb))
  z <- scale(object$lb, center=FALSE, scale = colSums(object$lb))
  
  # contr_block = DataFrame(z, string.("lv", 1:nlv))
  contr_block <- data.frame(z)
  colnames(contr_block) <- paste0("lv", 1:nlv)
  
  # Explained inertia for each block
  # z = fscale((object.lb)', sstot)'
  z <- t(scale(t(object$lb), center=FALSE, scale = sstot))
  
  # explX = DataFrame(z, string.("lv", 1:nlv))
  explX <- data.frame(z)
  colnames(explX) <- paste0("lv", 1:nlv)
  
  # Correlation between original variables and global scores
  z <- cor(Xconc, object$U)
  
  # corx2t = DataFrame(z, string.("lv", 1:nlv))
  corx2t <- data.frame(z)
  colnames(corx2t) <- paste0("lv", 1:nlv)
  
  # Correlation between block scores and global scores
  z <- vector("list", nlv)
  for (a in 1:nlv) {
    z[[a]] <- cor(object$Tb[[a]], object$U[, a])
  }
  
  cortb2t <- data.frame(do.call(cbind, z))
  colnames(cortb2t) <- paste0("lv", 1:nlv)
  
  # contr.ind
  TT <- object$weights * object$T * object$T
  # tt <- colSums(TT)
  contr_ind <- data.frame(.scale(TT, center = rep(0, nlv), scale = colSums(TT)))#tt))
  
  # contr_var, coord_var, cor_circle
  zT <- .scale(object$T, center = rep(0, nlv), scale = sqrt(colSums(TT)))#tt))
  cor_circle <- data.frame(t(object$weights * Xconc) %*% zT)#data.frame(t(object$weights * zX) %*% zT)
  Xinitconc <- as.matrix(do.call("cbind", Xinit))
  coord_var <- data.frame(crossprod(Xinitconc, object$weights * zT))#data.frame(crossprod(X, object$weights * zT))
  z <- coord_var^2
  contr_var <- data.frame(.scale(z, rep(0, nlv), colSums(z)))
  row.names(cor_circle) <- row.names(contr_var) <- row.names(coord_var) <- row.names(object$P)
  
  return(list(explvarx = explvarx, contr_block = contr_block, explX = explX, corx2t = corx2t, cortb2t = cortb2t, 
              contr_ind = contr_ind, contr_var = contr_var, coord_var = coord_var, cor_circle = cor_circle))#, rv = zrv, lg = zlg))
}
