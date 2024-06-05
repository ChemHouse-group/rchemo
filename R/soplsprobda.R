.soplsprobda <- function(Xlist, y, Xscaling = c("none", "pareto", "sd")[1], Yscaling = c("none", "pareto", "sd")[1], weights = NULL, nlv, fun, prior = c("unif", "prop")) {
    prior <- match.arg(prior)
    if(is.factor(y))
        y <- as.character(y)
    Xlist <- lapply(1:length(Xlist), function(x) .mat(Xlist[[x]]))
    
    n <- nrow(Xlist[[1]]) 
    p <- sapply(1:length(Xlist), function(x) ncol(Xlist[[x]]))
    nlv <- sapply(1:length(nlv), function(x) min(nlv[x], n, p[x]))
    if(is.null(weights))
        weights <- rep(1, n)
    weights <- .mweights(weights)
    zd <- dummy(y)
    fm <- list()
    fm[[1]] <- soplsr(Xlist, zd$Y, Xscaling = Xscaling, Yscaling = Yscaling, weights = weights, nlv = nlv)
    ## Should be:
    ## z <- transform(fm[[1]], Xlist)
    ## But same as:
    z <- fm[[1]]$T
    
    if(is.null(z)){
      fm[[2]] <- list()
      nlev <- length(zd$lev)
      fm[[2]]$ct <- rep(0, nlev)
      fm[[2]]$W <- NULL
      fm[[2]]$wprior <- switch(prior,
                               "unif" = rep(1 / nlev, nlev),
                               "prop" = zd$ni / sum(zd$ni))
      fm[[2]]$lev <- zd$lev
      fm[[2]]$ni <- zd$ni
    }else{
      fm[[2]] <- fun(z, y, prior = prior)
    }
    structure(list(fm = fm, lev = zd$lev, ni = zd$ni), 
              class = "Soplsprobda")       
}

soplslda <- function(Xlist, y, Xscaling = c("none", "pareto", "sd")[1], Yscaling = c("none", "pareto", "sd")[1], weights = NULL, nlv, prior = c("unif", "prop"))
    .soplsprobda(Xlist, y, Xscaling = Xscaling, Yscaling = Yscaling, weights = weights, nlv = nlv, fun = rchemo::lda, prior = prior)

soplsqda <- function(Xlist, y, Xscaling = c("none", "pareto", "sd")[1], Yscaling = c("none", "pareto", "sd")[1], weights = NULL, nlv, prior = c("unif", "prop"))
    .soplsprobda(Xlist, y, Xscaling = Xscaling, Yscaling = Yscaling, weights = weights, nlv = nlv, fun = rchemo::qda, prior = prior)


transform.Soplsprobda <- function(object, X, ...){
  
  X <- lapply(1:length(X), function(i) .mat(X[[i]]))
  transform(object$fm[[1]], X)
  
}

predict.Soplsprobda <- function(object, X, ...) {
  X <- lapply(1:length(X), function(i) .mat(X[[i]]))
  rownam <- row.names(X[[1]])
  colnam <- "y1"
  
  if(sum(object$fm[[1]]$nlv)>0){
    z <- transform(object$fm[[1]], X)
    
    zres <- predict(object$fm[[2]], z)
    pred <- zres$pred
    posterior <- zres$posterior
  }else{
    posterior <- matrix(rep(object$fm[[2]]$wprior, each = nrow(X[[1]])), 
                        ncol = length(object$fm[[2]]$lev),
                        dimnames = list(rownam, object$fm[[2]]$lev))
    pred <- matrix(object$fm[[2]]$lev[which.max(object$fm[[2]]$wprior)[1]], 
                   nrow = nrow(X[[1]]),
                   dimnames = list(rownam, colnam))
  }

  list(pred = pred, posterior = posterior)
}
