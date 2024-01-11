.mbplsprobda <- function(Xlist, y, scaling = c("none", "pareto", "sd")[1], blockscaling = TRUE, weights = NULL, nlv, fun, prior = c("unif", "prop")) {
    prior <- match.arg(prior)
    if(is.factor(y))
        y <- as.character(y)
    Xlist <- lapply(1:length(Xlist), function(x) .mat(Xlist[[x]]))
    # zdim <- dim(X)
    n <- nrow(Xlist[[1]]) # n <- zdim[1]
    p <- sum(sapply(1:length(Xlist), function(x) ncol(Xlist[[x]]))) # p <- zdim[2]
    nlv <- min(nlv, n, p)
    if(is.null(weights))
        weights <- rep(1, n)
    weights <- .mweights(weights)
    zd <- dummy(y)
    fm <- list()
    fm[[1]] <- mbplsr(Xlist, zd$Y, scaling = scaling, blockscaling = blockscaling, weights = weights, nlv = nlv)
    ## Should be:
    ## z <- transform(fm[[1]], X)
    ## But same as:
    z <- fm[[1]]$T
    fm[[2]] <- vector(length = nlv, mode = "list")
    for(i in seq_len(nlv))
        fm[[2]][[i]] <- fun(z[, seq_len(i), drop = FALSE], y, prior = prior)
    structure(list(fm = fm, lev = zd$lev, ni = zd$ni), 
              class = "Mbplsprobda")       
}

mbplslda <- function(Xlist, y, scaling = c("none", "pareto", "sd")[1], blockscaling = TRUE, weights = NULL, nlv, prior = c("unif", "prop"))
    .mbplsprobda(Xlist, y, scaling, blockscaling, weights, nlv, fun = rchemo::lda, prior = prior)

mbplsqda <- function(Xlist, y, scaling = c("none", "pareto", "sd")[1], blockscaling = TRUE, weights = NULL, nlv, prior = c("unif", "prop"))
    .mbplsprobda(Xlist, y, scaling, blockscaling, weights, nlv, fun = rchemo::qda, prior = prior)


predict.Mbplsprobda <- function(object, Xlist, ..., nlv = NULL) {
    Xlist <- lapply(1:length(Xlist), function(x) .mat(Xlist[[x]]))
    A <- length(object$fm[[2]])
    if(is.null(nlv))
        nlv <- A
    else 
        nlv <- seq(max(1, min(nlv)), min(max(nlv), A))
    le_nlv <- length(nlv)
    posterior <- pred <- vector(mode = "list", length = le_nlv)
    for(i in seq_len(le_nlv)) {
        znlv <- nlv[i]
        z <- transform(object$fm[[1]], Xlist, nlv = znlv)
        zres <- predict(object$fm[[2]][[znlv]], z)
        pred[[i]] <- zres$pred
        posterior[[i]] <- zres$posterior
    }
    names(posterior) <- names(pred) <- paste("lv", nlv, sep = "")
    if(le_nlv == 1) {
        pred <- pred[[1]] 
        posterior <- posterior[[1]]
    }
    list(pred = pred, posterior = posterior)
}
