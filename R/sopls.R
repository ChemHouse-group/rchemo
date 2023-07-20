sopls <- function(X_list, Y, nlv_vect, blocknames = NULL, weights = NULL, ...) {

  nbl <- length(X_list)
  ## an object 'blocks' (= a list with the block indexes) is created 
  nbcolX <- sapply(1:nbl, function(X){ncol(X_list[[X]])})
  nbcolXcum <- c(0,cumsum(nbcolX))
  blocks <- lapply(1:nbl, function(X){(nbcolXcum[X]+1):nbcolXcum[X+1]})
  
  ## nb lv by block
  if(length(nlv_vect) == 1){
    for(i in 1:nbl){
      nlv_vect[i] <- min(nbcolX[i], nlv_vect)
    }
  }
  
  X <- do.call("cbind",X_list)
  obsnames <- rownames(X)
  Xnames <- colnames(X)
  Ynames <- colnames(Y)
  X <- .mat(X)
  Y <- .mat(Y)     
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- dim(Y)[2]
  
  ## Case "sum(nlv_vect) = 0"
  if(sum(nlv_vect) == 0) {
    
    fm <- list()
    fm$W <- matrix(rep(0, p), nrow = p, byrow = TRUE,
                   dimnames = list(Xnames, "lv0"))
    fm$T <- matrix(rep(0, n), nrow = n, byrow = TRUE,
                   dimnames = list(obsnames, "comp0"))
    fm$C <- matrix(rep(0, q), nrow = q, byrow = TRUE,
                   dimnames = list(Ynames, "comp0"))
    
    # predictions
    dots <- list(...)
    ymeans <- .xmean(Y, weights = weights) #list(...)$weights)
    fm$Fit <- matrix(rep(ymeans, n), nrow = n, byrow = TRUE,
                     dimnames = list(obsnames, Ynames))
    
    fm$blocks <- 0
    fm$nlv_vect <- nlv_vect
    
    fm$Ymeans <- matrix(rep(ymeans, n), nrow = n, byrow = TRUE)
    
    # coef by block
    fm$BCoef <- list()
    fm$Bcoef$int <- lapply(1:nbl, function(x) matrix(ymeans, nrow = 1, dimnames = list("intercept",Ynames)))
                           
    fm$Bcoef$B <- lapply(1:nbl, function(x) matrix(0, nrow = length(blocks[[x]]), ncol=ncol(Y), 
                                                   dimnames = list(colnames(X[,blocks[[x]]]),Ynames)))
    
    # vip by block
    fm$VIP <- lapply(1:nbl, function(x) matrix(1, nrow = length(blocks[[x]]), ncol=ncol(Y), 
                                               dimnames = list(colnames(X[,blocks[[x]]]),Ynames)))
    names(fm$BCoef) <- names(fm$VIP) <- blocknames
    fm$call  <- match.call()
    
  }
  
  ## Other cases
  else {
    
    ## Used for defining the block indexes of the calculated scores
    zblocks <- data.frame(
      numcol = seq_len(sum(nlv_vect)), 
      bl = rep(seq_len(nbl), nlv_vect)
    ) 
    
    ## Reorganisation of the data based on the blocks defined in argument 'blocks'
    newdat <- .blocksel(X, blocks)
    X <- newdat$X
    newblocks <- newdat$blocks
    
    W <- list()
    Tx <- list()
    C <- list()
    BCoef <- list()
    VIP <- list()
    Ymeans <- list()
    Fit <- matrix(rep(0, n*q), nrow = n, byrow = TRUE,
                  dimnames = list(obsnames, Ynames))
    Ymeans[[1]] <- matrix(rep(0, n*q), nrow = n, byrow = TRUE)
    
    for(i in 1:nbl){
      if(nlv_vect[i] == 0){
        dots <- list(...)
        ymeans <- .xmean(Y-Fit, weights = weights) #list(...)$weights)
        
        W[[i]] <- matrix(rep(0, length(blocks[[i]])), nrow = length(blocks[[i]]), byrow = TRUE,
                         dimnames = list(colnames(X[,blocks[[i]],drop=FALSE]), "lv0"))
        Tx[[i]] <- matrix(rep(0, n), nrow = n, byrow = TRUE,
                          dimnames = list(obsnames, "comp0"))
        C[[i]] <- matrix(rep(0, q), nrow = q, byrow = TRUE,
                          dimnames = list(Ynames, "comp0"))
        BCoef[[i]] <- list()
        BCoef[[i]]$int <- matrix(ymeans, nrow = 1)
        BCoef[[i]]$B <- matrix(0, nrow = length(blocks[[i]]), ncol=ncol(Y))
        rownames(BCoef[[i]]$B) <- colnames(X[,blocks[[i]],drop=FALSE])
        colnames(BCoef[[i]]$B) <- colnames(BCoef[[i]]$int) <- Ynames
        VIP[[i]] <- matrix(1, nrow = length(blocks[[i]]), ncol=ncol(Y), dimnames=list(colnames(X[,blocks[[i]],drop=FALSE]),Ynames))
        Ymeans[[i]] <- matrix(rep(ymeans, n), nrow = n, byrow = TRUE)
        Fit <- Fit + Ymeans[[i]]
      }else{
        if(which(nlv_vect>0)[1]!=i){
           z <- orthog(T, X[, newblocks[[i]], drop = FALSE], fmpls$weights) # newblocks is orthogonalized to Tr
        
          fmpls <- plskern(
            X = z$Y,
            Y = Y - Fit,
            weights = weights,
            nlv = nlv_vect[i])
        }else{
          fmpls <- plskern(
            X[, newblocks[[i]], drop = FALSE],
            Y - Fit,
            weights = weights,
            nlv = nlv_vect[i])
        }
        W[[i]] <- fmpls$W
        Tx[[i]] <- fmpls$T
        T <- do.call("cbind",Tx)
        Ymeans[[i]] <- matrix(rep(fmpls$ymeans, n), nrow = n, byrow = TRUE)
        C[[i]] <- fmpls$C
        Fit <- Fit + Ymeans[[i]] + fmpls$T %*% t(fmpls$C) # t(fmpls$C) = beta
        
        # Bcoef (not in the original function)
        BCoef[[i]] <- coef(fmpls)# equivalent to fm$R %*% t(fm$C) 
        colnames(BCoef[[i]]$B) <- colnames(BCoef[[i]]$int) <- Ynames
        
        # VIP (not in the rnirs function)
        #VIP[[i]] <- vip(object=fmpls,X=X[, newblocks[[i]], drop = FALSE],Y=Y,nlv=nlv_vect[i])
        VIP[[i]] <- vip(object=fmpls,X=X[, newblocks[[i]], drop = FALSE],Y=Y-Fit,nlv=nlv_vect[i])
        
        ## Block indexes for the scores
        blocks[[i]] <- zblocks$numcol[zblocks$bl == i]
      }
    }
    
    # outputs
    fm <- list()
    fm$W <- W
    fm$T <- Tx
    fm$C <- C
    fm$BCoef <- BCoef
    fm$VIP   <- VIP
    fm$Ymeans   <- Ymeans
    names(fm$W) <- names(fm$T) <- names(fm$C) <- names(fm$BCoef) <- names(fm$VIP) <- names(fm$Ymeans) <- blocknames
    fm$Fit <- Fit
    fm$call  <- match.call()
  }
  
  fm
  
}
