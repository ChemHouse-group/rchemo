.soplsprobda <- function(Xlist, y, scaling = c("centered", "pareto", "ctreduced")[1], blockscaling = TRUE, weights = NULL, nlv, fun, prior = c("unif", "prop")) {
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
    fm[[1]] <- soplsr(Xlist, zd$Y, scaling = scaling, weights = weights, nlv = nlv)
    ## Should be:
    ## z <- transform(fm[[1]], X)
    ## But same as:
    z <- fm[[1]]$T
    
    fm[[2]] <- fun(z, y, prior = prior)
    structure(list(fm = fm, lev = zd$lev, ni = zd$ni), 
              class = "Soplsprobda")       
}

soplslda <- function(Xlist, y, scaling = c("centered", "pareto", "ctreduced")[1], weights = NULL, nlv, prior = c("unif", "prop"))
    .soplsprobda(Xlist, y, scaling = scaling, weights = weights, nlv = nlv, fun = rchemo::lda, prior = prior)

soplsqda <- function(Xlist, y, scaling = c("centered", "pareto", "ctreduced")[1], weights = NULL, nlv, prior = c("unif", "prop"))
    .soplsprobda(Xlist, y, scaling = scaling, weights = weights, nlv = nlv, fun = rchemo::qda, prior = prior)


transform.Soplsprobda <- function(object, Xlist){
  
  Xlist <- lapply(1:length(Xlist), function(X) .mat(Xlist[[X]]))
  transform(object$fm[[1]], Xlist)
  
}

predict.Soplsprobda <- function(object, Xlist) {
  Xlist <- lapply(1:length(Xlist), function(x) .mat(Xlist[[x]]))
  rownam <- row.names(Xlist[[1]])
  colnam <- "y1"
  
  if(sum(object$fm[[1]]$nlv)>0){
    z <- transform(object$fm[[1]], Xlist)
    
    zres <- predict(object$fm[[2]], z)
    pred <- zres$pred
    posterior <- zres$posterior
  }else{
    pred <- posterior <- NULL
  }

  list(pred = pred, posterior = posterior)
}
