
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
  xsdslist   <- lapply(1:length(Xlist), function(X) sqrt(.colvars(Xlist[[X]], weights = weights)))
  #xsdslist   <- lapply(1:length(Xlist), function(X) sqrt(.colvars(Xlist[[X]], weights = weights))*nrow(Xlist[[X]])/(nrow(Xlist[[X]])-1))
  ysds       <- sqrt(.colvars(Y, weights = weights))
  #ysds       <- sqrt(.colvars(Y, weights = weights)*nrow(Y)/(nrow(Y)-1))
  
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
  
  D <- diag(weights)
  
  fm <- list()
  b <- list()
  ymeanslist <- list()
  
  # First block
  if(nlv[1]>0){    
    fm[[1]] <- plskern(Xlist[[1]], Y, scaling = "none", weights = weights, nlv = nlv[1])
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
        fm[[i]] <- plskern(X, Y - pred, scaling = "none", weights = weights, nlv = nlv[i])
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
  
  if(Yscaling == "pareto"){
    pred <- pred * matrix(rep(sqrt(ysds), nrow(Y)), nrow=nrow(Y), byrow=TRUE)
  }
  if(Yscaling == "sd"){
    pred <- pred * matrix(rep(ysds, nrow(Y)), nrow=nrow(Y), byrow=TRUE)
  }
  pred <- pred + matrix(rep(ymeans, nrow(Y)), nrow=nrow(Y), byrow=TRUE)
  structure(
    list(fm = fm, T = T, pred = pred, xmeans = xmeanslist, ymeans = ymeans, xsds = xsdslist, ysds = ysds, b = b, weights = weights, Xscaling = Xscaling, Yscaling = Yscaling, nlv = nlv),
    class = c("Soplsr"))
}


soplsr <- function(Xlist, Y, Xscaling = c("none", "pareto", "sd")[1], Yscaling = c("none", "pareto", "sd")[1], weights = NULL, nlv){

  nbl <- length(Xlist)  
  zXlist <- lapply(1:nbl, function(k).mat(Xlist[[k]]))

  .soplsr(zXlist, .mat(Y), Xscaling = Xscaling, Yscaling = Yscaling, weights = weights, nlv = nlv)
}



transform.Soplsr <- function(object, Xlist){
  
  Xlist <- lapply(1:length(Xlist), function(X) .mat(Xlist[[X]]))
  
  for(i in 1:length(Xlist)){
    if(object$Xscaling[i] == "none"){
     Xlist[[i]] <- .center(Xlist[[i]], object$xmeans[[i]])
    }
    if(object$Xscaling[i] == "pareto"){
      Xlist[[i]] <- .center(Xlist[[i]], object$xmeans[[i]])
      Xlist[[i]] <- scale(Xlist[[i]], center = FALSE, scale = sqrt(object$xsds[[i]]))
    }
    if(object$Xscaling[i] == "sd"){
      Xlist[[i]] <- .center(Xlist[[i]], object$xmeans[[i]])
      Xlist[[i]] <- scale(Xlist[[i]], center = FALSE, scale = object$xsds[[i]])
    }
  }
  
  nbl <- length(object$fm)
  
  if(object$nlv[1]>0){
    T <- transform(object$fm[[1]], Xlist[[1]])
  }else{
    T <- NULL
  }
  
  if (nbl > 1){
    for (i in 2:nbl){
      if(object$nlv[i]>0){
        X = Xlist[[i]] - T %*% object$b[[i]]
        T = cbind(T, transform(object$fm[[i]], X))
      }
    }
  }
  T
}


predict.Soplsr <- function (object, Xlist){
  
  Xlist <- lapply(1:length(Xlist), function(X) .mat(Xlist[[X]]))
  
  for(i in 1:length(Xlist)){
    if(object$Xscaling[i] == "none"){
      Xlist[[i]] <- .center(Xlist[[i]], object$xmeans[[i]])
    }
    if(object$Xscaling[i] == "pareto"){
      Xlist[[i]] <- .center(Xlist[[i]], object$xmeans[[i]])
      Xlist[[i]] <- scale(Xlist[[i]], center = FALSE, scale = sqrt(object$xsds[[i]]))
    }
    if(object$Xscaling[i] == "sd"){
      Xlist[[i]] <- .center(Xlist[[i]], object$xmeans[[i]])
      Xlist[[i]] <- scale(Xlist[[i]], center = FALSE, scale = object$xsds[[i]])
    }
  }
  
  # nbl <- length(object$fm)
  nbl <- length(Xlist)
  m <- nrow(Xlist[[1]])
  
  if(object$nlv[1]>0){
    T <- transform(object$fm[[1]], Xlist[[1]])
    pred <- predict(object$fm[[1]], Xlist[[1]])$pred
    # pred <- matrix(rep(object$fm[[1]]$ymeans, m), nrow=m, byrow=TRUE) + T %*% t(object$fm[[1]]$C)
  }else{
    T <- NULL
    pred <- matrix(rep(object$ymeans, m), nrow=m, byrow=TRUE)*0
    rownames(pred) <- rownames(object$pred)
    colnames(pred) <- colnames(object$pred)
  }

  if (nbl > 1){
    for (i in 2:nbl){
      if(object$nlv[i]>0){
        X <- Xlist[[i]] - T %*% object$b[[i]]
        # zT <- transform(object$fm[[i]], X)
        pred <- pred + predict(object$fm[[i]], X)$pred
        # pred <- pred + matrix(rep(object$fm[[i]]$ymeans, m), nrow=m, byrow=TRUE) + zT %*% t(object$fm[[i]]$C)
        T = cbind(T, transform(object$fm[[i]], X))
      }else{
        # X <- Xlist[[i]]
        # pred <- pred + matrix(rep(object$fm[[i]]$ymeans, m), nrow=m, byrow=TRUE)
        pred <- pred # TO CHECK
      }
    }
  }
  
  if(object$Yscaling == "none"){
    pred <- pred
  }
  if(object$Yscaling == "pareto"){
    pred <- pred * matrix(rep(sqrt(object$ysds), m), nrow=m, byrow=TRUE) 
  }
  if(object$Yscaling == "sd"){
    pred <- pred * matrix(rep(object$ysds, m), nrow=m, byrow=TRUE) 
  }
  pred <- pred + matrix(rep(object$ymeans, m), nrow=m, byrow=TRUE) 
  pred
}
           
         