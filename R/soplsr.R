
.soplsr <- function (Xlist, Y, Xscaling = c("none", "pareto", "sd")[1], Yscaling = c("none", "pareto", "sd")[1], weights = NULL, nlv){
  
  n <- nrow(Xlist[[1]])
  q <- ncol(Y)   
  nbl <- length(Xlist)
  
  if(length(nlv) == 1){nlv <- rep(nlv, nbl)}
  if(length(Xscaling) == 1){Xscaling <- rep(Xscaling, nbl)}
  
  if(is.null(weights))
    weights <- rep(1, nrow(Y))
  weights <- .mweights(weights)
  
  xmeanslist <- lapply(1:length(Xlist), function(X) .colmeans(Xlist[[X]], weights = weights))
  ymeans     <- .colmeans(Y, weights = weights) 
  
  xscaleslist <- list()
  for(i in 1:length(Xlist)){
    if(Xscaling[i] == "none"){xscaleslist[[i]] <- rep(1, ncol(Xlist[[i]]))}
    if(Xscaling[i] == "pareto"){xscaleslist[[i]] <- sqrt(sqrt(.colvars(Xlist[[i]], weights = weights)))}
    if(Xscaling[i] == "sd"){xscaleslist[[i]] <- sqrt(.colvars(Xlist[[i]], weights = weights))}
  }
  Xlist <- lapply(1:length(Xlist), function(X) scale(Xlist[[X]], center = xmeanslist[[X]], scale = xscaleslist[[X]]))
  
  if(Yscaling == "none"){yscales <- rep(1, ncol(Y))}
  if(Yscaling == "pareto"){yscales <- sqrt(sqrt(.colvars(Y, weights = weights)))}
  if(Yscaling == "sd"){yscales <- sqrt(.colvars(Y, weights = weights))}
  Y <- scale(Y, center = ymeans, scale = yscales)
  
  D <- diag(weights)
  
  fm <- list()
  b <- list()
  ymeanslist <- list()
  
  # First block
  if(nlv[1]>0){    
    fm[[1]] <- plskern(Xlist[[1]], Y, Xscaling = "none", Yscaling = "none",weights = weights, nlv = nlv[1])
    T <- fm[[1]]$T
    ymeanslist[[1]] <- fm[[1]]$ymeans
    pred <- predict(fm[[1]], Xlist[[1]])$pred
  }else{
    fm[[1]] <- NULL
    T <- NULL
    ymeanslist[[1]] <- apply(Y, FUN = mean, MARGIN = 2, na.rm = TRUE)
    pred <- matrix(rep(apply(Y, FUN = mean, MARGIN = 2, na.rm = TRUE), nrow(Y)), nrow=nrow(Y), byrow=TRUE)
    rownames(pred) <- rownames(Y)
  }
  b[[1]] <- NA
  
  # Other blocks
  if (nbl > 1){
    for (i in 2:nbl){
      if(nlv[i]>0){
        if(sum(nlv[1:(i-1)])>0){
          b[[i]] <- .ginv(t(T) %*% (D %*% T)) %*% t(T) %*% (D %*% Xlist[[i]])
          X <- Xlist[[i]] - T %*% b[[i]]
        }else{
          b[[i]] <- NA
          X <- Xlist[[i]]
        }
        fm[[i]] <- plskern(X, Y - pred, Xscaling = "none", Yscaling = "none", weights = weights, nlv = nlv[i])
        T <- cbind(T, fm[[i]]$T)
        ymeanslist[[i]] <- fm[[i]]$ymeans
        pred <- pred + predict(fm[[i]], X)$pred 
      }else{
        if(sum(nlv[1:(i-1)])>0){
          b[[i]] <- .ginv(t(T) %*% (D %*% T)) %*% t(T) %*% (D %*% Xlist[[i]])
          # X <- Xlist[[i]] - T %*% b[[i]]
        }else{
          b[[i]] <- NA
          # X <- Xlist[[i]]
        }
        fm[[i]] <- NULL
        #fm[[i]]$T <- NULL
        # T <- T
        ymeanslist[[i]] <- apply((Y - pred), FUN = mean, MARGIN = 2, na.rm = TRUE)
        pred <- pred 
        # pred <- pred + matrix(rep(apply((Y - pred), FUN = mean, MARGIN = 2, na.rm = TRUE), nrow(Y)), nrow=nrow(Y), byrow=TRUE)
      }
    }
  }
  
  pred <- pred * matrix(rep(yscales, nrow(Y)), nrow=nrow(Y), byrow=TRUE)+ matrix(rep(ymeans, nrow(Y)), nrow=nrow(Y), byrow=TRUE)

  structure(
    list(fm = fm, T = T, pred = pred, xmeans = xmeanslist, ymeans = ymeans, xscales = xscaleslist, yscales = yscales, b = b, weights = weights, nlv = nlv),
    class = c("Soplsr"))
}


soplsr <- function(Xlist, Y, Xscaling = c("none", "pareto", "sd")[1], Yscaling = c("none", "pareto", "sd")[1], weights = NULL, nlv){

  nbl <- length(Xlist)  
  zXlist <- lapply(1:nbl, function(k).mat(Xlist[[k]]))

  .soplsr(zXlist, .mat(Y), Xscaling = Xscaling, Yscaling = Yscaling, weights = weights, nlv = nlv)
}



transform.Soplsr <- function(object, X, ...){
  
  X <- lapply(1:length(X), function(i) scale(.mat(X[[i]]), center = object$xmeans[[i]], scale = object$xscales[[i]]))
  
  nbl <- length(object$fm)
  
  if(sum(object$nlv)==0){T <- NULL}
  
  if(object$nlv[1]>0){
    T <- transform(object$fm[[1]], X[[1]])
  }
  
  if (nbl > 1){
    for (i in 2:nbl){
      if(object$nlv[i]>0){
        if(sum(object$nlv[1:(i-1)])>0){
          Xo = X[[i]] - T %*% object$b[[i]]
          T = cbind(T, transform(object$fm[[i]], Xo))
        }else{
          T <- transform(object$fm[[i]], X[[i]])
        }
      }
    }
  }
  T
}


predict.Soplsr <- function (object, X, ...){
  
  X <- lapply(1:length(X), function(i) scale(.mat(X[[i]]), center = object$xmeans[[i]], scale = object$xscales[[i]]))
  
  nbl <- length(X)
  m <- nrow(X[[1]])
  
  if(object$nlv[1]>0){
    T <- transform(object$fm[[1]], X[[1]])
    pred <- predict(object$fm[[1]], X[[1]])$pred
    # pred <- matrix(rep(object$fm[[1]]$ymeans, m), nrow=m, byrow=TRUE) + T %*% t(object$fm[[1]]$C)
  }else{
    T <- NULL
    pred <- matrix(rep(object$ymeans, m), nrow=m, byrow=TRUE)*0
    rownames(pred) <- rownames(X[[1]])
    colnames(pred) <- colnames(object$pred)
  }

  if (nbl > 1){
    for (i in 2:nbl){
      if(object$nlv[i]>0){
        if(sum(object$nlv[1:(i-1)])>0){
          Xo <- X[[i]] - T %*% object$b[[i]]
          # zT <- transform(object$fm[[i]], Xo)
          pred <- pred + predict(object$fm[[i]], Xo)$pred
          # pred <- pred + matrix(rep(object$fm[[i]]$ymeans, m), nrow=m, byrow=TRUE) + zT %*% t(object$fm[[i]]$C)
          T = cbind(T, transform(object$fm[[i]], Xo))
        }else{
          Xo <- X[[i]]
          pred <- predict(object$fm[[i]], Xo)$pred
          T <- transform(object$fm[[i]], Xo)
        }
      # }else{
        # Xo <- X[[i]]
        # pred <- pred 
      # }
      }
    }
  }
  
  pred <- pred * matrix(rep(object$yscales, m), nrow=m, byrow=TRUE) + matrix(rep(object$ymeans, m), nrow=m, byrow=TRUE) 

  pred
}
           
         