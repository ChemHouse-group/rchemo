
.soplsr <- function (Xlist, Y, scaling = c("Centered", "Pareto", "CtReduced")[1], weights = NULL, nlv){
  
  # Xlist <- lapply(1:length(Xlist), function(X) .mat(Xlist[[X]]))
  # Y <- .mat(Y)
  n <- nrow(Xlist[[1]])
  q <- ncol(Y)   
  nbl <- length(Xlist)
  
  if(length(nlv) == 1){nlv <- rep(nlv, nbl)}
  if(length(scaling) == 1){scaling <- rep(scaling, nbl)}
  
  if(is.null(weights))
    weights <- rep(1, nrow(Y))
  weights <- .mweights(weights)
  D <- diag(weights)
  
  fm <- list()
  b <- list()
  ymeans <- list()
  
  # First block
  if(nlv[1]>0){    
    fm[[1]] <- plskern(Xlist[[1]], Y, scaling = scaling[1], weights = weights, nlv = nlv[1])
    T <- fm[[1]]$T
    ymeans[[1]] <- fm[[1]]$ymeans
    pred <- predict(fm[[1]], Xlist[[1]])$pred
  }else{
    fm[[1]] <- NULL
    T <- NULL
    ymeans[[1]] <- apply(Y, FUN = mean, MARGIN = 2, na.rm = TRUE)
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
        fm[[i]] <- plskern(X, Y - pred, scaling = scaling[i], weights = weights, nlv = nlv[i])
        T <- cbind(T, fm[[i]]$T)
        ymeans[[i]] <- fm[[i]]$ymeans
        pred <- pred + predict(fm[[i]], X)$pred 
      }else{
        if(sum(nlv[1:(i-1)])>0){
          b[[i]] <- .ginv(t(T) %*% (D %*% T)) %*% t(T) %*% (D %*% Xlist[[i]])
          X <- Xlist[[i]] - T %*% b[[i]]
        }else{
          b[[i]] <- NA
          X <- Xlist[[i]]
        }
        fm[[i]] <- NULL
        #fm[[i]]$T <- NULL
        # T <- T
        ymeans[[i]] <- apply((Y - pred), FUN = mean, MARGIN = 2, na.rm = TRUE)
        pred <- pred + matrix(rep(apply((Y - pred), FUN = mean, MARGIN = 2, na.rm = TRUE), nrow(Y)), nrow=nrow(Y), byrow=TRUE)
      }
    }
  }
  structure(
    list(fm = fm, T = T, pred = pred, ymeans = ymeans, b = b, weights = weights, scaling = scaling, nlv = nlv),
    class = c("Soplsr"))
}


soplsr <- function(Xlist, Y, scaling = c("Centered", "Pareto", "CtReduced")[1], weights = NULL, nlv){
  
  # if(is.null(weights))
  #   weights <- rep(1, nrow(Y))
  # weights <- .mweights(weights)
  
  nbl <- length(Xlist)  
  zXlist <- lapply(1:nbl, function(k).mat(Xlist[[k]]))

  .soplsr(zXlist, .mat(Y), scaling = scaling, weights = weights, nlv = nlv)
}



           
# transform.Soplsr(object, Xlist)
# Compute latent variables (LVs = scores T) from a fitted model.
# %*% `object` : The fitted model.
# %*% `Xlist` : A list of blocks (matrices) for which LVs are computed.

transform.Soplsr <- function(object, Xlist){
  
  Xlist <- lapply(1:length(Xlist), function(X) .mat(Xlist[[X]]))
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

# predict.Soplsr(object, Xlist)
# Compute Y-predictions from a fitted model.
# %*% `object` : The fitted model.
# %*% `Xlist` : A list of X-data for which predictions are computed.

predict.Soplsr <- function (object, Xlist){
  
  Xlist <- lapply(1:length(Xlist), function(X) .mat(Xlist[[X]]))
  nbl <- length(object$fm)
  if(object$nlv[1]>0){
    T <- transform(object$fm[[1]], Xlist[[1]])
    #pred <- t(object$fm[[1]]$ymeans) + T %*% t(object$fm[[1]]$C)
    pred <- matrix(rep(object$ymeans[[1]], nrow(Y)), nrow=nrow(Y), byrow=TRUE) + T %*% t(object$fm[[1]]$C)
  }else{
    T <- NULL
    pred <- matrix(rep(object$ymeans[[1]], nrow(Y)), nrow=nrow(Y), byrow=TRUE)
    rownames(pred) <- rownames(Y)
  }

  if (nbl > 1){
    for (i in 2:nbl){
      if(object$nlv[i]>0){
        X <- Xlist[[i]] - T %*% object$b[[i]]
        zT <- transform(object$fm[[i]], X)
        #pred <- pred + t(object$fm[[i]]$ymeans) + zT %*% t(object$fm[[i]]$C)
        pred <- pred + matrix(rep(object$ymeans[[i]], nrow(Y)), nrow=nrow(Y), byrow=TRUE) + zT %*% t(object$fm[[i]]$C)
        T = cbind(T, transform(object$fm[[i]], X))
      }else{
        X <- Xlist[[i]]
        pred <- pred + matrix(rep(object$ymeans[[i]], nrow(Y)), nrow=nrow(Y), byrow=TRUE)
        #T <- T
      }
    }
  }
  pred
}
           
         