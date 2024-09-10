
plsr_plsda_allsteps <- function(X, Xname = NULL, Xscaling = c("none","pareto","sd")[1], 
                                Y, Yscaling = c("none","pareto","sd")[1], weights = NULL,
                                newX = NULL, newXname = NULL,
                                
                                method = c("plsr", "plsrda","plslda","plsqda")[1],
                                prior = c("unif", "prop")[1], 
                                
                                step = c("nlvtest","permutation","model","prediction")[1],
                                nlv, 
                                modeloutput = c("scores","loadings","coef","vip"), 
                                
                                cvmethod = c("kfolds","loo")[1], 
                                nbrep = 30, 
                                seed = 123, 
                                samplingk = NULL, 
                                nfolds = 10, 
                                npermut = 30, 
                                
                                criterion = c("err","rmse")[1], 
                                selection = c("localmin","globalmin","1std")[1],
                                
                                import = c("R","ChemFlow","W4M")[1],
                                outputfilename = NULL
                               
){
  
  # IMPORTATION DES DONNEES -------------------------------------------------------------------------------------------
  
  ### librairies
  library(rchemo)
  
  if(is.null(Xname)==TRUE){Xname <- "X"}
  if((is.null(newX)==FALSE)&(is.null(newXname)==TRUE)){newXname <- "Xnew"}
  
  if(length(method)>1){stop("length of method must be 1")}
  if(length(prior)>1){stop("length of prior must be 1")} 
  if(length(cvmethod)>1){stop("length of cvmethod must be 1")}
  if(length(step)>1){stop("length of step must be 1")} 
  
  if(length(criterion)>1){stop("length of criterion must be 1")}
  if(length(selection)>1){stop("length of selection must be 1")} 
  if(length(nbrep)>1){stop("length of nbrep must be 1")}
  if(length(seed)>1){stop("length of seed must be 1")} 
  if(length(nfolds)>1){stop("length of nfolds must be 1")}
  if(length(npermut)>1){stop("length of npermut must be 1")} 
 
  options(max.print=99999)

  if(import == "R"){# X matrix n x p; Y matrix n x q, in the Global Environment, with observation names in rownames
    Y           <- as.matrix(Y)
    n           <- nrow(Y)

    # X           <- X
  }
  
  if(import == "ChemFlow"){# X matrix n x p; Y matrix n x q, tabulated tables (.txt), with observation names in the first column
    Y           <- read.table(Y, sep="\t", header=TRUE, na.strings = c("","NA"), check.names=FALSE)
    rownames(Y) <- Y[,1]
    Y[,1]       <- NULL
    Y           <- as.matrix(Y)
    n           <- nrow(Y)

    X           <- read.table(X,sep="\t",dec=".",header=TRUE, na.strings=c("","NA"), check.names=FALSE)
    rownames(X) <- X[,1]
    X[,1]       <- NULL
  }
  
  if(import == "W4M"){# X matrix p x n; Y matrix n x q, tabulated tables (.txt), with observation names in the header of X and the first column of Y
    Y           <- read.table(Y, sep="\t", header=TRUE, na.strings = c("","NA"), check.names=FALSE)
    rownames(Y) <- Y[,1]
    Y[,1]       <- NULL
    Y           <- as.matrix(Y)
    n           <- nrow(Y)

    X           <- read.table(X,sep="\t",dec=".",header=TRUE, na.strings=c("","NA"), check.names=FALSE)
    X           <- t(X)
    rownames(X) <- X[,1]
    X[,1]       <- NULL
  }
  
  SamplesNotFound1 <- NULL
  SamplesNotFound2 <- NULL
  if(length(which(rownames(X)%in%rownames(Y)))<length(rownames(X))){
    SamplesNotFound1 <- rownames(X)[(rownames(X)%in%rownames(Y))==FALSE]
  } 
  if(length(which(rownames(Y)%in%rownames(X)))<length(rownames(Y))){
    SamplesNotFound2 <- rownames(Y)[(rownames(Y)%in%(rownames(X)))==FALSE]
  }
  if(is.null(SamplesNotFound1)==FALSE & is.null(SamplesNotFound2)==FALSE){
    stop("\n- - - - - - - - -\n",paste0("Some samples of ", Xname, " are not found in the Y: \n"),paste(SamplesNotFound1, collapse="\n"),
         "\n All samples must be in the Y. Please, check your data.",
         "\n\n",paste0("Some samples of the Y are not found in ", Xname, ": \n"), paste(SamplesNotFound2, collapse="\n"),
         "\n Please remove them or do missing value imputation.","\n- - - - - - - - -\n")
  }
  if(is.null(SamplesNotFound1)==FALSE & is.null(SamplesNotFound2)==TRUE){
    stop("\n- - - - - - - - -\n",paste0("Some samples of ", Xname, " are not found in the Y: \n"),paste(SamplesNotFound1, collapse="\n"),
         "\n All samples must be in the Y. Please, check your data.","\n- - - - - - - - -\n")
  }
  if(is.null(SamplesNotFound1)==TRUE & is.null(SamplesNotFound2)==FALSE){
    stop("\n- - - - - - - - -\n",paste0("Some samples of the Y are not found in ", Xname, ": \n"), paste(SamplesNotFound2, collapse="\n"),
         "\n Please remove them or do missing value imputation.","\n- - - - - - - - -\n")
  }
  X <- data.frame(X)

  if((step=="prediction")&(is.null(newX)==TRUE)){
    newX     <- X
    newXname <- Xname
  }
  if((step=="prediction")&(is.null(newX)==FALSE)){
    if(import == "R"){# X matrix n x p; Y matrix n x q, in the Global Environment, with observation names in rownames
      # newX       <- newX
    }
    if(import == "ChemFlow"){# X matrix n x p; Y matrix n x q, tabulated tables (.txt), with observation names in the first column
      newX       <- read.table(newX,sep="\t",dec=".",header=TRUE, na.strings=c("","NA"), check.names=FALSE)
      rownames(newX) <- newX[,1]
      newX[,1]       <- NULL
    }
    if(import == "W4M"){# X matrix p x n; Y matrix n x q, tabulated tables (.txt), with observation names in the header of X and the first column of Y
      newX       <- read.table(newX,sep="\t",dec=".",header=TRUE, na.strings=c("","NA"), check.names=FALSE)
      newX       <- t(newX)
      rownames(newX) <- newX[,1]
      newX[,1]       <- NULL
    }
    
    VariablesNotFound1 <- NULL
    VariablesNotFound2 <- NULL
    if(length(which(colnames(newX)%in%colnames(X)))<length(colnames(newX))){
      VariablesNotFound1 <- colnames(newX)[(colnames(newX)%in%colnames(X))==FALSE]
    } 
    if(length(which(colnames(X)%in%colnames(newX)))<length(colnames(X))){
      VariablesNotFound2 <- colnames(X)[(colnames(X)%in%(colnames(newX)))==FALSE]
    }
    if(is.null(VariablesNotFound1)==FALSE & is.null(VariablesNotFound2)==FALSE){
      stop("\n- - - - - - - - -\n",paste0("Some variables of ", newXname, " are not found in the X: \n"),paste(VariablesNotFound1, collapse="\n"),
           "\n All variables must be in the X. Please, check your data.",
           "\n\n",paste0("Some variables of the X are not found in ", newXname, ": \n"), paste(VariablesNotFound2, collapse="\n"),
           "\n Please remove them or do missing value imputation.","\n- - - - - - - - -\n")
    }
    if(is.null(VariablesNotFound1)==FALSE & is.null(VariablesNotFound2)==TRUE){
      stop("\n- - - - - - - - -\n",paste0("Some variables of ", newXname, " are not found in the X: \n"),paste(VariablesNotFound1, collapse="\n"),
           "\n All variables must be in the X. Please, check your data.","\n- - - - - - - - -\n")
    }
    if(is.null(VariablesNotFound1)==TRUE & is.null(VariablesNotFound2)==FALSE){
      stop("\n- - - - - - - - -\n",paste0("Some variables of the X are not found in ", newXname, ": \n"), paste(VariablesNotFound2, collapse="\n"),
           "\n Please remove them or do missing value imputation.","\n- - - - - - - - -\n")
    }
    newX <- data.frame(newX)
  }
    

  ### nlvtest
  
  if (step == "nlvtest"){

    # partition
    set.seed(seed = seed)
    if(cvmethod=="loo"){
      nbrep=1
      CVtype <- segmkf(n = nrow(X), K = nrow(X), type = "random", nrep = 1)
      }
    if((cvmethod=="kfolds") & (is.null(samplingk)==TRUE)){
      CVtype <- segmkf(n = nrow(X), K = nfolds, type = "random", nrep = nbrep)
      for(r in 1:nbrep){
        for(ss in 1:nfolds){
          CVtype[[r]][[ss]] <- sort(CVtype[[r]][[ss]])
        }
      }
    }
    if((cvmethod=="kfolds") & (is.null(samplingk)==FALSE)){
      samplingktab <- table(samplingk)
      partCVtype <- list()
      for(s in 1:length(names(samplingktab))){
        partCVtype[[s]] <- segmkf(n = samplingktab[s], K = nfolds, type = "random", nrep = nbrep)
        for(r in 1:nbrep){
          for(ss in 1:length(partCVtype[[s]][[r]])){
            partCVtype[[s]][[r]][[ss]] <- which(samplingk==names(samplingktab)[s])[partCVtype[[s]][[r]][[ss]]]
          }
        }
      }
      CVtype <- list()
      for(r in 1:nbrep){
        CVtype[[r]] <- list()
        for(ss in 1:nfolds){
          CVtype[[r]][[ss]] <- c(0)
          for(s in 1:length(names(samplingktab))){
            CVtype[[r]][[ss]] <- c(CVtype[[r]][[ss]], partCVtype[[s]][[r]][[ss]])
          }
          CVtype[[r]][[ss]] <- sort(CVtype[[r]][[ss]][-1])
        }
      }
    }
    
    if(method == c("plsr", "plsrda","plslda","plsqda")[1]){
      
      repnlvtest <- cbind.data.frame(expand.grid(nblv=0:nlv, rep=1:nbrep),rmse = rep(NA,(nbrep*(1+nlv))))
      for(i in 1:nbrep){
        reppredCV <- lapply(1:(nlv+1), function(l) matrix(NA, nrow=n, ncol = ncol(Y)))
        for(k in 1:length(CVtype[[i]])){
          Xlisttrain <- X[-CVtype[[i]][[k]],,drop=FALSE]
          Ytrain <- Y[-CVtype[[i]][[k]],,drop=FALSE]
          Xlisttest <- X[CVtype[[i]][[k]],,drop=FALSE]
          Ytest <- Y[CVtype[[i]][[k]],,drop=FALSE]
          repmodel <- plskern(X=Xlisttrain, Y=Ytrain, Xscaling = Xscaling, Yscaling = Yscaling, weights = weights[-CVtype[[i]][[k]]], nlv=nlv)
          predictions <- predict(repmodel, Xlisttest, nlv = 0:nlv)
          for(l in 1:(nlv+1)){
            reppredCV[[l]][CVtype[[i]][[k]],] <- predictions$pred[[l]][order(CVtype[[i]][[k]]),]
          }
        }
        for(l in 1:(nlv+1)){
          SqErrCV <- (Y-reppredCV[[l]])^2
          repnlvtest[(repnlvtest$rep==i)&(repnlvtest$nblv==(l-1)),"rmse"] <- sqrt(mean(SqErrCV))
        }
      }
      
      resnlvtesttable <- cbind.data.frame(expand.grid(nblv=0:nlv), rmse_mean = rep(NA,(1+nlv)), rmse_sd = rep(NA,(1+nlv)))
      for(j in 0:nlv){
        resnlvtesttable[resnlvtesttable$nblv==j,"rmse_mean"] <- mean(repnlvtest$rmse[repnlvtest$nblv==j], na.rm = TRUE)
        resnlvtesttable[resnlvtesttable$nblv==j,"rmse_sd"] <- sd(repnlvtest$rmse[repnlvtest$nblv==j], na.rm = TRUE)
      }
      
      optimum = rep(0, (nlv+1))
      if(selection == "globalmin"){
        optimum[which.min(resnlvtesttable$rmse_mean)[1]] <- 1
      }
      if(selection == "localmin"){
        diff<- sapply(1:nlv, function(i) resnlvtesttable$rmse_mean[i]<resnlvtesttable$rmse_mean[i+1])
        optimum[which(diff==TRUE)[1]] <- 1
      }
      if(selection == "1std"){
        seuil <- resnlvtesttable$rmse_mean[which.min(resnlvtesttable$rmse_mean)[1]]+resnlvtesttable$rmse_sd[which.min(resnlvtesttable$rmse_mean)[1]]
        # if(resnlvtesttable$rmse_mean[1] <= seuil) {
        #   optimum[1] <- 1
        # }else{
          intermed <- resnlvtesttable[1:which.min(resnlvtesttable$rmse_mean)[1],,drop=FALSE]
          optimum[min(which(intermed$rmse_mean<=seuil))] <- 1
        # }
      }
      
      resnlvtesttable <- cbind(resnlvtesttable, optimum = optimum)
      
      if(is.null(outputfilename)==FALSE){
        write.table(data.frame(rownames(resnlvtesttable),resnlvtesttable), 
                    file = paste0(outputfilename,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                    eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                    col.names = TRUE, qmethod = c("escape", "double"),
                    fileEncoding = "")
      }
    }
    
    
    if(method == c("plsr", "plsrda","plslda","plsqda")[2]){
      
      repnlvtest <- cbind.data.frame(expand.grid(nblv=0:nlv, rep=1:nbrep),rmse = rep(NA,(nbrep*(1+nlv))))
      repYnlvtest <- cbind.data.frame(expand.grid(nblv=0:nlv, rep=1:nbrep),err = rep(NA,(nbrep*(1+nlv))))
      for(i in 1:nbrep){
        reppredCV <- lapply(1:(nlv+1), function(l) matrix(NA, nrow=n, ncol = ncol(dummy(as.matrix(Y))$Y)))
        repYpredCV <- lapply(1:(nlv+1), function(l) matrix(NA, nrow=n, ncol = ncol(Y)))
        for(k in 1:length(CVtype[[i]])){
          Xlisttrain <- X[-CVtype[[i]][[k]],,drop=FALSE]
          Ytrain <- Y[-CVtype[[i]][[k]],,drop=FALSE]
          Xlisttest <- X[CVtype[[i]][[k]],,drop=FALSE]
          Ytest <- Y[CVtype[[i]][[k]],,drop=FALSE]
          repmodel <- plsrda(X=Xlisttrain, y=Ytrain, Xscaling = Xscaling, Yscaling = Yscaling, weights = weights[-CVtype[[i]][[k]]], nlv=nlv)
          predictions <- predict(repmodel, Xlisttest, nlv = 0:nlv)
          for(l in 1:(nlv+1)){
            reppredCV[[l]][CVtype[[i]][[k]],] <- predictions$posterior[[l]][order(CVtype[[i]][[k]]),]
            repYpredCV[[l]][CVtype[[i]][[k]],] <- predictions$pred[[l]][order(CVtype[[i]][[k]]),]
          }
        }
        for(l in 1:(nlv+1)){
          SqErrCV <- (dummy(as.matrix(Y))$Y-reppredCV[[l]])^2
          repnlvtest[(repnlvtest$rep==i)&(repnlvtest$nblv==(l-1)),"rmse"] <- sqrt(mean(SqErrCV))
          repYnlvtest[(repnlvtest$rep==i)&(repnlvtest$nblv==(l-1)),"err"] <- err(repYpredCV[[l]],Y)
        }
      }
      
      if(criterion == "rmse"){
        resnlvtesttable <- cbind.data.frame(expand.grid(nblv=0:nlv), rmse_mean = rep(NA,(1+nlv)), rmse_sd = rep(NA,(1+nlv)))
        for(j in 0:nlv){
          resnlvtesttable[resnlvtesttable$nblv==j,"rmse_mean"] <- mean(repnlvtest$rmse[repnlvtest$nblv==j], na.rm = TRUE)
          resnlvtesttable[resnlvtesttable$nblv==j,"rmse_sd"] <- sd(repnlvtest$rmse[repnlvtest$nblv==j], na.rm = TRUE)
        }
        
        optimum = rep(0, (nlv+1))
        if(selection == "globalmin"){
          optimum[which.min(resnlvtesttable$rmse_mean)[1]] <- 1
        }
        if(selection == "localmin"){
          diff<- sapply(1:nlv, function(i) resnlvtesttable$rmse_mean[i]<resnlvtesttable$rmse_mean[i+1])
          optimum[which(diff==TRUE)[1]] <- 1
        }
        if(selection == "1std"){
          seuil <- resnlvtesttable$rmse_mean[which.min(resnlvtesttable$rmse_mean)[1]]+resnlvtesttable$rmse_sd[which.min(resnlvtesttable$rmse_mean)[1]]
          # if(resnlvtesttable$rmse_mean[1] <= seuil) {
          #   optimum[1] <- 1
          # }else{
            intermed <- resnlvtesttable[1:which.min(resnlvtesttable$rmse_mean)[1],,drop=FALSE]
            optimum[min(which(intermed$rmse_mean<=seuil))] <- 1
          # }
        }
        
        resnlvtesttable <- cbind(resnlvtesttable, optimum = optimum)
      }
      
      if(criterion == "err"){
        resnlvtesttable <- cbind.data.frame(expand.grid(nblv=0:nlv), err_mean = rep(NA,(1+nlv)), err_sd = rep(NA,(1+nlv)))
        
        for(j in 0:nlv){
          resnlvtesttable[resnlvtesttable$nblv==j,"err_mean"] <- mean(repYnlvtest$err[repYnlvtest$nblv==j], na.rm = TRUE)
          resnlvtesttable[resnlvtesttable$nblv==j,"err_sd"] <- sd(repYnlvtest$err[repYnlvtest$nblv==j], na.rm = TRUE)
        }
        
        optimum = rep(0, (nlv+1))
        if(selection == "globalmin"){
          optimum[which.min(resnlvtesttable$err_mean)[1]] <- 1
        }
        if(selection == "localmin"){
          diff<- sapply(1:nlv, function(i) resnlvtesttable$err_mean[i]<resnlvtesttable$err_mean[i+1])
          optimum[which(diff==TRUE)[1]] <- 1
        }
        if(selection == "1std"){
          seuil <- resnlvtesttable$err_mean[which.min(resnlvtesttable$err_mean)[1]]+resnlvtesttable$err_sd[which.min(resnlvtesttable$err_mean)[1]]
          # if(resnlvtesttable$err_mean[1] <= seuil) {
          #   optimum[1] <- 1
          # }else{
            intermed <- resnlvtesttable[1:which.min(resnlvtesttable$err_mean)[1],,drop=FALSE]
            optimum[min(which(intermed$err_mean<=seuil))] <- 1
          # }
        }
        
        resnlvtesttable <- cbind(resnlvtesttable, optimum = optimum)
      }
      
      if(is.null(outputfilename)==FALSE){
        write.table(data.frame(rownames(resnlvtesttable),resnlvtesttable), 
                    file = paste0(outputfilename,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                    eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                    col.names = TRUE, qmethod = c("escape", "double"),
                    fileEncoding = "")
      }
    }
    
    
    if(method == c("plsr", "plsrda","plslda","plsqda")[3]){
      
      repnlvtest <- cbind.data.frame(expand.grid(nblv=0:nlv, rep=1:nbrep),rmse = rep(NA,(nbrep*(1+nlv))))
      repYnlvtest <- cbind.data.frame(expand.grid(nblv=0:nlv, rep=1:nbrep),err = rep(NA,(nbrep*(1+nlv))))
      for(i in 1:nbrep){
        reppredCV <- lapply(1:(nlv+1), function(l) matrix(NA, nrow=n, ncol = ncol(dummy(Y)$Y)))
        repYpredCV <- lapply(1:(nlv+1), function(l) matrix(NA, nrow=n, ncol = ncol(Y)))
        for(k in 1:length(CVtype[[i]])){
          Xlisttrain <- X[-CVtype[[i]][[k]],,drop=FALSE]
          Ytrain <- Y[-CVtype[[i]][[k]],,drop=FALSE]
          Xlisttest <- X[CVtype[[i]][[k]],,drop=FALSE]
          Ytest <- Y[CVtype[[i]][[k]],,drop=FALSE]
          #nlv = 0
          repmodel0lv <- plsrda(X=Xlisttrain, y=Ytrain, Xscaling = Xscaling , Yscaling = Yscaling , weights = weights[-CVtype[[i]][[k]]], nlv=1)
          predictions0lv <- predict(repmodel0lv, Xlisttest, nlv = 0)
          reppredCV[[1]][CVtype[[i]][[k]],] <- predictions0lv$posterior[order(CVtype[[i]][[k]]),]
          repYpredCV[[1]][CVtype[[i]][[k]],] <- predictions0lv$pred[order(CVtype[[i]][[k]]),]
          #nlv>0
          repmodel <- plslda(X=Xlisttrain, y=Ytrain, Xscaling = Xscaling , Yscaling = Yscaling, weights = weights[-CVtype[[i]][[k]]], nlv=nlv, prior = prior)
          predictions <- predict(repmodel, Xlisttest, nlv = 1:nlv)
          for(l in 2:(nlv+1)){
            reppredCV[[l]][CVtype[[i]][[k]],] <- predictions$posterior[[(l-1)]][order(CVtype[[i]][[k]]),]
            repYpredCV[[l]][CVtype[[i]][[k]],] <- predictions$pred[[(l-1)]][order(CVtype[[i]][[k]]),]
          }
        }
        for(l in 1:(nlv+1)){
          SqErrCV <- (dummy(as.matrix(Y))$Y-reppredCV[[l]])^2
          repnlvtest[(repnlvtest$rep==i)&(repnlvtest$nblv==(l-1)),"rmse"] <- sqrt(mean(SqErrCV))
          repYnlvtest[(repnlvtest$rep==i)&(repnlvtest$nblv==(l-1)),"err"] <- err(repYpredCV[[l]],Y)
        }
      }
      
      if(criterion == "rmse"){
        resnlvtesttable <- cbind.data.frame(expand.grid(nblv=0:nlv), rmse_mean = rep(NA,(1+nlv)), rmse_sd = rep(NA,(1+nlv)))
        for(j in 0:nlv){
          resnlvtesttable[resnlvtesttable$nblv==j,"rmse_mean"] <- mean(repnlvtest$rmse[repnlvtest$nblv==j], na.rm = TRUE)
          resnlvtesttable[resnlvtesttable$nblv==j,"rmse_sd"] <- sd(repnlvtest$rmse[repnlvtest$nblv==j], na.rm = TRUE)
        }
        
        optimum = rep(0, (nlv+1))
        if(selection == "globalmin"){
          optimum[which.min(resnlvtesttable$rmse_mean)[1]] <- 1
        }
        if(selection == "localmin"){
          diff<- sapply(1:nlv, function(i) resnlvtesttable$rmse_mean[i]<resnlvtesttable$rmse_mean[i+1])
          optimum[which(diff==TRUE)[1]] <- 1
        }
        if(selection == "1std"){
          seuil <- resnlvtesttable$rmse_mean[which.min(resnlvtesttable$rmse_mean)[1]]+resnlvtesttable$rmse_sd[which.min(resnlvtesttable$rmse_mean)[1]]
          # if(resnlvtesttable$rmse_mean[1] <= seuil) {
          #   optimum[1] <- 1
          # }else{
            intermed <- resnlvtesttable[1:which.min(resnlvtesttable$rmse_mean)[1],,drop=FALSE]
            optimum[min(which(intermed$rmse_mean<=seuil))] <- 1
          # }
        }
        
        resnlvtesttable <- cbind(resnlvtesttable, optimum = optimum)
      }
      
      if(criterion == "err"){
        resnlvtesttable <- cbind.data.frame(expand.grid(nblv=0:nlv), err_mean = rep(NA,(1+nlv)), err_sd = rep(NA,(1+nlv)))
        
        for(j in 0:nlv){
          resnlvtesttable[resnlvtesttable$nblv==j,"err_mean"] <- mean(repYnlvtest$err[repYnlvtest$nblv==j], na.rm = TRUE)
          resnlvtesttable[resnlvtesttable$nblv==j,"err_sd"] <- sd(repYnlvtest$err[repYnlvtest$nblv==j], na.rm = TRUE)
        }
        
        optimum = rep(0, (nlv+1))
        if(selection == "globalmin"){
          optimum[which.min(resnlvtesttable$err_mean)[1]] <- 1
        }
        if(selection == "localmin"){
          diff<- sapply(1:nlv, function(i) resnlvtesttable$err_mean[i]<resnlvtesttable$err_mean[i+1])
          optimum[which(diff==TRUE)[1]] <- 1
        }
        if(selection == "1std"){
          seuil <- resnlvtesttable$err_mean[which.min(resnlvtesttable$err_mean)[1]]+resnlvtesttable$err_sd[which.min(resnlvtesttable$err_mean)[1]]
          # if(resnlvtesttable$err_mean[1] <= seuil) {
          #   optimum[1] <- 1
          # }else{
            intermed <- resnlvtesttable[1:which.min(resnlvtesttable$err_mean)[1],,drop=FALSE]
            optimum[min(which(intermed$err_mean<=seuil))] <- 1
          # }
        }
        
        resnlvtesttable <- cbind(resnlvtesttable, optimum = optimum)
      }
      
      if(is.null(outputfilename)==FALSE){
        write.table(data.frame(rownames(resnlvtesttable),resnlvtesttable), 
                    file = paste0(outputfilename,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                    eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                    col.names = TRUE, qmethod = c("escape", "double"),
                    fileEncoding = "")
      }
    }
    
    if(method == c("plsr", "plsrda","plslda","plsqda")[4]){
      
      repnlvtest <- cbind.data.frame(expand.grid(nblv=0:nlv, rep=1:nbrep),rmse = rep(NA,(nbrep*(1+nlv))))
      repYnlvtest <- cbind.data.frame(expand.grid(nblv=0:nlv, rep=1:nbrep),err = rep(NA,(nbrep*(1+nlv))))
      for(i in 1:nbrep){
        reppredCV <- lapply(1:(nlv+1), function(l) matrix(NA, nrow=n, ncol = ncol(dummy(as.matrix(Y))$Y)))
        repYpredCV <- lapply(1:(nlv+1), function(l) matrix(NA, nrow=n, ncol = ncol(Y)))
        for(k in 1:length(CVtype[[i]])){
          Xlisttrain <- X[-CVtype[[i]][[k]],,drop=FALSE]
          Ytrain <- Y[-CVtype[[i]][[k]],,drop=FALSE]
          Xlisttest <- X[CVtype[[i]][[k]],,drop=FALSE]
          Ytest <- Y[CVtype[[i]][[k]],,drop=FALSE]
          #nlv = 0
          repmodel0lv <- plsrda(X=Xlisttrain, y=Ytrain, Xscaling = Xscaling, Yscaling = Yscaling, weights = weights[-CVtype[[i]][[k]]], nlv=1)
          predictions0lv <- predict(repmodel0lv, Xlisttest, nlv = 0)
          reppredCV[[1]][CVtype[[i]][[k]],] <- predictions0lv$posterior[order(CVtype[[i]][[k]]),]
          repYpredCV[[1]][CVtype[[i]][[k]],] <- predictions0lv$pred[order(CVtype[[i]][[k]]),]
          #nlv>0
          repmodel <- plsqda(X=Xlisttrain, y=Ytrain, Xscaling = Xscaling, Yscaling = Yscaling, weights = weights[-CVtype[[i]][[k]]], nlv=nlv, prior = prior)
          predictions <- predict(repmodel, Xlisttest, nlv = 1:nlv)
          for(l in 2:(nlv+1)){
            reppredCV[[l]][CVtype[[i]][[k]],] <- predictions$posterior[[(l-1)]][order(CVtype[[i]][[k]]),]
            repYpredCV[[l]][CVtype[[i]][[k]],] <- predictions$pred[[(l-1)]][order(CVtype[[i]][[k]]),]
          }
        }
        for(l in 1:(nlv+1)){
          SqErrCV <- (dummy(as.matrix(Y))$Y-reppredCV[[l]])^2
          repnlvtest[(repnlvtest$rep==i)&(repnlvtest$nblv==(l-1)),"rmse"] <- sqrt(mean(SqErrCV))
          repYnlvtest[(repnlvtest$rep==i)&(repnlvtest$nblv==(l-1)),"err"] <- err(repYpredCV[[l]],Y)
        }
      }
      
      if(criterion == "rmse"){
        resnlvtesttable <- cbind.data.frame(expand.grid(nblv=0:nlv), rmse_mean = rep(NA,(1+nlv)), rmse_sd = rep(NA,(1+nlv)))
        for(j in 0:nlv){
          resnlvtesttable[resnlvtesttable$nblv==j,"rmse_mean"] <- mean(repnlvtest$rmse[repnlvtest$nblv==j], na.rm = TRUE)
          resnlvtesttable[resnlvtesttable$nblv==j,"rmse_sd"] <- sd(repnlvtest$rmse[repnlvtest$nblv==j], na.rm = TRUE)
        }
        
        optimum = rep(0, (nlv+1))
        if(selection == "globalmin"){
          optimum[which.min(resnlvtesttable$rmse_mean)[1]] <- 1
        }
        if(selection == "localmin"){
          diff<- sapply(1:nlv, function(i) resnlvtesttable$rmse_mean[i]<resnlvtesttable$rmse_mean[i+1])
          optimum[which(diff==TRUE)[1]] <- 1
        }
        if(selection == "1std"){
          seuil <- resnlvtesttable$rmse_mean[which.min(resnlvtesttable$rmse_mean)[1]]+resnlvtesttable$rmse_sd[which.min(resnlvtesttable$rmse_mean)[1]]
          # if(resnlvtesttable$rmse_mean[1] <= seuil) {
          #   optimum[1] <- 1
          # }else{
            intermed <- resnlvtesttable[1:which.min(resnlvtesttable$rmse_mean)[1],,drop=FALSE]
            optimum[min(which(intermed$rmse_mean<=seuil))] <- 1
          # }
        }
        
        resnlvtesttable <- cbind(resnlvtesttable, optimum = optimum)
      }
      
      if(criterion == "err"){
        resnlvtesttable <- cbind.data.frame(expand.grid(nblv=0:nlv), err_mean = rep(NA,(1+nlv)), err_sd = rep(NA,(1+nlv)))
        
        for(j in 0:nlv){
          resnlvtesttable[resnlvtesttable$nblv==j,"err_mean"] <- mean(repYnlvtest$err[repYnlvtest$nblv==j], na.rm = TRUE)
          resnlvtesttable[resnlvtesttable$nblv==j,"err_sd"] <- sd(repYnlvtest$err[repYnlvtest$nblv==j], na.rm = TRUE)
        }
        
        optimum = rep(0, (nlv+1))
        if(selection == "globalmin"){
          optimum[which.min(resnlvtesttable$err_mean)[1]] <- 1
        }
        if(selection == "localmin"){
          diff<- sapply(1:nlv, function(i) resnlvtesttable$err_mean[i]<resnlvtesttable$err_mean[i+1])
          optimum[which(diff==TRUE)[1]] <- 1
        }
        if(selection == "1std"){
          seuil <- resnlvtesttable$err_mean[which.min(resnlvtesttable$err_mean)[1]]+resnlvtesttable$err_sd[which.min(resnlvtesttable$err_mean)[1]]
          # if(resnlvtesttable$err_mean[1] <= seuil) {
          #   optimum[1] <- 1
          # }else{
            intermed <- resnlvtesttable[1:which.min(resnlvtesttable$err_mean)[1],,drop=FALSE]
            optimum[min(which(intermed$err_mean<=seuil))] <- 1
          # }
        }
        
        resnlvtesttable <- cbind(resnlvtesttable, optimum = optimum)
      }
      
      if(is.null(outputfilename)==FALSE){
        write.table(data.frame(rownames(resnlvtesttable),resnlvtesttable), 
                    file = paste0(outputfilename,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                    eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                    col.names = TRUE, qmethod = c("escape", "double"),
                    fileEncoding = "")
      }
    }
    return(resnlvtesttable)
  }
  
  ### model
  
  if (step == "model"){
    if(method == c("plsr", "plsrda","plslda","plsqda")[1]){
      resmodel <- plskern(X = X, 
                         Y = Y, 
                         Xscaling = Xscaling, 
                         Yscaling = Yscaling, 
                         weights = weights, 
                         nlv = nlv)
      if("scores" %in% modeloutput){
        if(nlv>0){
          Tx <- resmodel$T
          if(is.null(outputfilename)==FALSE){
            write.table(data.frame(rownames(Tx),Tx), 
                        file = paste0(outputfilename,"_scores.txt"), append = FALSE, quote = TRUE, sep = "\t",
                        eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                        col.names = TRUE, qmethod = c("escape", "double"),
                        fileEncoding = "")
          }
        }
      }
      if("loadings" %in% modeloutput){
        if(nlv>0){
          P <- resmodel$P
          if(is.null(outputfilename)==FALSE){
            write.table(data.frame(rownames(P),P), 
                        file = paste0(outputfilename,"_loadings.txt"), append = FALSE, quote = TRUE, sep = "\t",
                        eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                        col.names = TRUE, qmethod = c("escape", "double"),
                        fileEncoding = "")
          }
          
        }
      }
      if("coef" %in% modeloutput){
        if(nlv>0){
          coefmatrix <- coef(resmodel, nlv = nlv)$B
          if(is.null(outputfilename)==FALSE){
            write.table(data.frame(rownames(coefmatrix),coefmatrix), 
                        file = paste0(outputfilename,"_coef.txt"), append = FALSE, quote = TRUE, sep = "\t",
                        eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                        col.names = TRUE, qmethod = c("escape", "double"),
                        fileEncoding = "")
          }
          
        }
      }
      if("vip" %in% modeloutput){
        if(nlv>0){
          WW <- resmodel$W * resmodel$W/apply(resmodel$W, 2, function(x) sum(x * x))
          
          Q <- resmodel$C
          if(length(dim(Q)) == 0){
            Q2 <- as.numeric(Q) * as.numeric(Q)
          } else {
            Q2 <- rowSums(t(Q * Q))
          }
          
          vipmatrix <- matrix(0, nrow = nrow(resmodel$W), ncol = nlv)
          for(i in 1:nlv){
            Q2TT <- Q2[1:i] * diag(crossprod(resmodel$T))[i]
            vipmatrix[,i] <- sqrt(nrow(resmodel$W) * apply(sweep(WW[, 1:i, drop=FALSE],2,Q2TT,"*"), 1, sum)/sum(Q2TT))
          }
          rownames(vipmatrix) <- rownames(resmodel$W)
          colnames(vipmatrix) <- paste("comp", 1:nlv)
          
          if(is.null(outputfilename)==FALSE){
            write.table(data.frame(rownames(vipmatrix),vipmatrix), 
                        file = paste0(outputfilename,"_vip.txt"), append = FALSE, quote = TRUE, sep = "\t",
                        eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                        col.names = TRUE, qmethod = c("escape", "double"),
                        fileEncoding = "")
          }
          
        }
      }
    }

    
    if(method == c("plsr", "plsrda","plslda","plsqda")[2]){
      resmodel <- plsrda(X = X, 
                           y = Y, 
                           Xscaling = Xscaling, 
                           Yscaling = Yscaling, 
                           weights = weights, 
                           nlv = nlv)
      if("scores" %in% modeloutput){
        if(nlv>0){
          Tx <- resmodel$fm$T
          if(is.null(outputfilename)==FALSE){
            write.table(data.frame(rownames(Tx),Tx), 
                        file = paste0(outputfilename,"_scores.txt"), append = FALSE, quote = TRUE, sep = "\t",
                        eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                        col.names = TRUE, qmethod = c("escape", "double"),
                        fileEncoding = "")
          }
          
        }
      }
      if("loadings" %in% modeloutput){
        if(nlv>0){
          P <- resmodel$fm$P
          if(is.null(outputfilename)==FALSE){
            write.table(data.frame(rownames(P),P), 
                        file = paste0(outputfilename,"_loadings.txt"), append = FALSE, quote = TRUE, sep = "\t",
                        eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                        col.names = TRUE, qmethod = c("escape", "double"),
                        fileEncoding = "")
          }
          
        }
      }
      if("coef" %in% modeloutput){
        if(nlv>0){
          coefmatrix <- coef(resmodel$fm, nlv = nlv)$B
          if(is.null(outputfilename)==FALSE){
            write.table(data.frame(rownames(coefmatrix),coefmatrix), 
                        file = paste0(outputfilename,"_coef.txt"), append = FALSE, quote = TRUE, sep = "\t",
                        eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                        col.names = TRUE, qmethod = c("escape", "double"),
                        fileEncoding = "")
          }
          
        }
      }
      if("vip" %in% modeloutput){
        if(nlv>0){
          WW <- resmodel$fm$W * resmodel$fm$W/apply(resmodel$fm$W, 2, function(x) sum(x * x))
          
          Q <- resmodel$fm$C
          if(length(dim(Q)) == 0){
            Q2 <- as.numeric(Q) * as.numeric(Q)
          } else {
            Q2 <- rowSums(t(Q * Q))
          }
          
          vipmatrix <- matrix(0, nrow = nrow(resmodel$fm$W), ncol = nlv)
          for(i in 1:nlv){
            Q2TT <- Q2[1:i] * diag(crossprod(resmodel$fm$T))[i]
            vipmatrix[,i] <- sqrt(nrow(resmodel$fm$W) * apply(sweep(WW[, 1:i, drop=FALSE],2,Q2TT,"*"), 1, sum)/sum(Q2TT))
          }
          rownames(vipmatrix) <- rownames(resmodel$fm$W)
          colnames(vipmatrix) <- paste("comp", 1:nlv)
          
          if(is.null(outputfilename)==FALSE){
            write.table(data.frame(rownames(vipmatrix),vipmatrix), 
                        file = paste0(outputfilename,"_vip.txt"), append = FALSE, quote = TRUE, sep = "\t",
                        eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                        col.names = TRUE, qmethod = c("escape", "double"),
                        fileEncoding = "")
          }
          
        }
      }
    }
    
    if(method == c("plsr", "plsrda","plslda","plsqda")[3]){
      resmodel <- plslda(X = X, 
                           y = Y, 
                           Xscaling = Xscaling, 
                           Yscaling = Yscaling, 
                           weights = weights, 
                           nlv = nlv,
                           prior = prior)
      if("scores" %in% modeloutput){
        if(nlv>0){
          Tx <- resmodel$fm[[1]]$T
          if(is.null(outputfilename)==FALSE){
            write.table(data.frame(rownames(Tx),Tx), 
                        file = paste0(outputfilename,"_scores.txt"), append = FALSE, quote = TRUE, sep = "\t",
                        eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                        col.names = TRUE, qmethod = c("escape", "double"),
                        fileEncoding = "")
          }
          
        }
      }
      if("loadings" %in% modeloutput){
        if(nlv>0){
          P <- resmodel$fm[[1]]$P
          if(is.null(outputfilename)==FALSE){
            write.table(data.frame(rownames(P),P), 
                        file = paste0(outputfilename,"_loadings.txt"), append = FALSE, quote = TRUE, sep = "\t",
                        eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                        col.names = TRUE, qmethod = c("escape", "double"),
                        fileEncoding = "")
          }
          
        }
      }
      if("coef" %in% modeloutput){
        if(nlv>0){
          coefmatrix <- coef(resmodel$fm[[1]], nlv = nlv)$B
          if(is.null(outputfilename)==FALSE){
            write.table(data.frame(rownames(coefmatrix),coefmatrix), 
                        file = paste0(outputfilename,"_coef.txt"), append = FALSE, quote = TRUE, sep = "\t",
                        eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                        col.names = TRUE, qmethod = c("escape", "double"),
                        fileEncoding = "")
          }
          
        }
      }
      if("vip" %in% modeloutput){
        if(nlv>0){
          WW <- resmodel$fm[[1]]$W * resmodel$fm[[1]]$W/apply(resmodel$fm[[1]]$W, 2, function(x) sum(x * x))
          
          Q <- resmodel$fm[[1]]$C
          if(length(dim(Q)) == 0){
            Q2 <- as.numeric(Q) * as.numeric(Q)
          } else {
            Q2 <- rowSums(t(Q * Q))
          }
          
          vipmatrix <- matrix(0, nrow = nrow(resmodel$fm[[1]]$W), ncol = nlv)
          for(i in 1:nlv){
            Q2TT <- Q2[1:i] * diag(crossprod(resmodel$fm[[1]]$T))[i]
            vipmatrix[,i] <- sqrt(nrow(resmodel$fm[[1]]$W) * apply(sweep(WW[, 1:i, drop=FALSE],2,Q2TT,"*"), 1, sum)/sum(Q2TT))
          }
          rownames(vipmatrix) <- rownames(resmodel$fm[[1]]$W)
          colnames(vipmatrix) <- paste("comp", 1:nlv)
          
          if(is.null(outputfilename)==FALSE){
            write.table(data.frame(rownames(vipmatrix),vipmatrix), 
                        file = paste0(outputfilename,"_vip.txt"), append = FALSE, quote = TRUE, sep = "\t",
                        eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                        col.names = TRUE, qmethod = c("escape", "double"),
                        fileEncoding = "")
          }
          
        }
      }
    }
    
    if(method == c("plsr", "plsrda","plslda","plsqda")[4]){
      resmodel <- plsqda(X = X, 
                           y = Y, 
                           Xscaling = Xscaling, 
                           Yscaling = Yscaling, 
                           weights = weights, 
                           nlv = nlv,
                           prior = prior)
      if("scores" %in% modeloutput){
        if(nlv>0){
          Tx <- resmodel$fm[[1]]$T
          if(is.null(outputfilename)==FALSE){
            write.table(data.frame(rownames(Tx),Tx), 
                        file = paste0(outputfilename,"_scores.txt"), append = FALSE, quote = TRUE, sep = "\t",
                        eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                        col.names = TRUE, qmethod = c("escape", "double"),
                        fileEncoding = "")
          }
          
        }
      }
      if("loadings" %in% modeloutput){
        if(nlv>0){
          P <- resmodel$fm[[1]]$P
          if(is.null(outputfilename)==FALSE){
            write.table(data.frame(rownames(P),P), 
                        file = paste0(outputfilename,"_loadings.txt"), append = FALSE, quote = TRUE, sep = "\t",
                        eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                        col.names = TRUE, qmethod = c("escape", "double"),
                        fileEncoding = "")
          }
          
        }
      }
      if("coef" %in% modeloutput){
        if(nlv>0){
          coefmatrix <- coef(resmodel$fm[[1]], nlv = nlv)$B
          if(is.null(outputfilename)==FALSE){
            write.table(data.frame(rownames(coefmatrix),coefmatrix), 
                        file = paste0(outputfilename,"_coef.txt"), append = FALSE, quote = TRUE, sep = "\t",
                        eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                        col.names = TRUE, qmethod = c("escape", "double"),
                        fileEncoding = "")
          }
          
        }
      }
      if("vip" %in% modeloutput){
        if(nlv>0){
          WW <- resmodel$fm[[1]]$W * resmodel$fm[[1]]$W/apply(resmodel$fm[[1]]$W, 2, function(x) sum(x * x))
          
          Q <- resmodel$fm[[1]]$C
          if(length(dim(Q)) == 0){
            Q2 <- as.numeric(Q) * as.numeric(Q)
          } else {
            Q2 <- rowSums(t(Q * Q))
          }
          
          vipmatrix <- matrix(0, nrow = nrow(resmodel$fm[[1]]$W), ncol = nlv)
          for(i in 1:nlv){
            Q2TT <- Q2[1:i] * diag(crossprod(resmodel$fm[[1]]$T))[i]
            vipmatrix[,i] <- sqrt(nrow(resmodel$fm[[1]]$W) * apply(sweep(WW[, 1:i, drop=FALSE],2,Q2TT,"*"), 1, sum)/sum(Q2TT))
          }
          rownames(vipmatrix) <- rownames(resmodel$fm[[1]]$W)
          colnames(vipmatrix) <- paste("comp", 1:nlv)
          
          if(is.null(outputfilename)==FALSE){
            write.table(data.frame(rownames(vipmatrix),vipmatrix), 
                        file = paste0(outputfilename,"_vip.txt"), append = FALSE, quote = TRUE, sep = "\t",
                        eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                        col.names = TRUE, qmethod = c("escape", "double"),
                        fileEncoding = "")
          }
          
        }
      }
    }
    outputs <- list()
    if("scores" %in% modeloutput){
      outputs$scores <- Tx
    }
    if("loadings" %in% modeloutput){
      outputs$loadings <- P
    }
    if("coef" %in% modeloutput){
      outputs$coef <- coefmatrix
    }
    if("vip" %in% modeloutput){
      outputs$vip <- vipmatrix
    }
    return(outputs)
  }
  
  
  ### permutation
  
  if (step == "permutation"){

    set.seed(seed = seed)
    permut_list <- lapply(1:npermut, function(i){sample(1:nrow(Y),nrow(Y))})
    permut_list[[(npermut+1)]] <- 1:nrow(Y)
    
    # partition
    set.seed(seed = seed)
    if(cvmethod=="loo"){
      nbrep=1
      CVtype <- segmkf(n = nrow(X), K = nrow(X), type = "random", nrep = 1)
    }
    if((cvmethod=="kfolds") & (is.null(samplingk)==TRUE)){
      CVtype <- segmkf(n = nrow(X), K = nfolds, type = "random", nrep = nbrep)
      for(r in 1:nbrep){
        for(ss in 1:nfolds){
          CVtype[[r]][[ss]] <- sort(CVtype[[r]][[ss]])
        }
      }
    }
    if((cvmethod=="kfolds") & (is.null(samplingk)==FALSE)){
      samplingktab <- table(samplingk)
      partCVtype <- list()
      for(s in 1:length(names(samplingktab))){
        partCVtype[[s]] <- segmkf(n = samplingktab[s], K = nfolds, type = "random", nrep = nbrep)
        for(r in 1:nbrep){
          for(ss in 1:length(partCVtype[[s]][[r]])){
            partCVtype[[s]][[r]][[ss]] <- which(samplingk==names(samplingktab)[s])[partCVtype[[s]][[r]][[ss]]]
          }
        }
      }
      CVtype <- list()
      for(r in 1:nbrep){
        CVtype[[r]] <- list()
        for(ss in 1:nfolds){
          CVtype[[r]][[ss]] <- c(0)
          for(s in 1:length(names(samplingktab))){
            CVtype[[r]][[ss]] <- c(CVtype[[r]][[ss]], partCVtype[[s]][[r]][[ss]])
          }
          CVtype[[r]][[ss]] <- sort(CVtype[[r]][[ss]][-1])
        }
      }
    }
    
    if(method == c("plsr", "plsrda","plslda","plsqda")[1]){
      
      coefficientRV <- function(X, Y) {
        tr <- function(Z) {
          sum(diag(Z))
        }
        Y <- scale(Y, scale = FALSE)
        X <- scale(X, scale = FALSE)
        
        W1 <- t(X) %*% X
        W2 <- t(Y) %*% Y
        W3 <- t(X) %*% Y
        W4 <- t(Y) %*% X
        rv <- tr(W3 %*% W4)/(tr(W1 %*% W1) * tr(W2 %*% W2))^0.5
        
        return(rv)
      }
      
      permut_dyssimilarity <- c()
      res_permut <- c()
      for(j in 1:(npermut+1)){
        Ypermut <- Y[permut_list[[j]],,drop=FALSE]
        permut_dyssimilarity[j] <- 1-coefficientRV(Ypermut,Y)
        
        repnlvpermut <- cbind.data.frame(expand.grid(permut=1:(npermut+1), rep=1:nbrep),rmse = rep(NA,(nbrep*(1+npermut))))
        for(i in 1:nbrep){
          reppredCV <- matrix(NA, nrow=n, ncol = ncol(Y))
          for(k in 1:length(CVtype[[i]])){
            Xlisttrain <- X[-CVtype[[i]][[k]],,drop=FALSE]
            Ytrain <- Ypermut[-CVtype[[i]][[k]],,drop=FALSE]
            Xlisttest <- X[CVtype[[i]][[k]],,drop=FALSE]
            Ytest <- Ypermut[CVtype[[i]][[k]],,drop=FALSE]
            repmodel <- plskern(X=Xlisttrain, Y=Ytrain, Xscaling = Xscaling, Yscaling = Yscaling, weights = weights[-CVtype[[i]][[k]]], nlv=nlv)
            predictions <- predict(repmodel, Xlisttest, nlv = nlv)
            reppredCV[CVtype[[i]][[k]],] <- predictions$pred[order(CVtype[[i]][[k]]),]
          }
          SqErrCV <- (Ypermut-reppredCV)^2
          repnlvpermut[(repnlvpermut$rep==i)&(repnlvpermut$permut==j),"rmse"] <- sqrt(mean(SqErrCV))
        }
        res_permut[j] <- apply(repnlvpermut[repnlvpermut$permut==j,"rmse", drop = FALSE], FUN = mean, MARGIN = 2, na.rm = TRUE)
      }
      
      #plot(permut_dyssimilarity,res_permut, ylim=c(0,1))
      #abline(h=res_permut[npermut+1])
      
      # Calcul de la p-valeur
      pval_permut <- mean(res_permut[1:npermut]<res_permut[npermut+1],na.rm=TRUE)
      
      #output
      respermuttesttable <- cbind.data.frame(permut_dyssimilarity,res_permut)
      if(is.null(outputfilename)==FALSE){
        write.table(data.frame(respermuttesttable), 
                    file = paste0(outputfilename,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                    eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                    col.names = TRUE, qmethod = c("escape", "double"),
                    fileEncoding = "")
      }
    }
    
    if(method == c("plsr", "plsrda","plslda","plsqda")[2]){
      
      permut_dyssimilarity <- c()
      res_permut <- c()
      for(j in 1:(npermut+1)){
        Ypermut <- Y[permut_list[[j]],,drop=FALSE]
        permut_dyssimilarity[j] <- err(Ypermut,Y)
        
        repnlvpermut <- cbind.data.frame(expand.grid(permut=1:(npermut+1), rep=1:nbrep),rmse = rep(NA,(nbrep*(1+npermut))))
        repnlvYpermut <- cbind.data.frame(expand.grid(permut=1:(npermut+1), rep=1:nbrep),err = rep(NA,(nbrep*(1+npermut))))
        
        for(i in 1:nbrep){
          reppredCV <- matrix(NA, nrow=n, ncol = ncol(dummy(as.matrix(Y))$Y))
          repYpredCV <- matrix(NA, nrow=n, ncol = ncol(Y))
          for(k in 1:length(CVtype[[i]])){
            Xlisttrain <- X[-CVtype[[i]][[k]],,drop=FALSE]
            Ytrain <- Ypermut[-CVtype[[i]][[k]],,drop=FALSE]
            Xlisttest <- X[CVtype[[i]][[k]],,drop=FALSE]
            Ytest <- Ypermut[CVtype[[i]][[k]],,drop=FALSE]
            repmodel <- plsrda(X=Xlisttrain, y=Ytrain, Xscaling = Xscaling, Yscaling = Yscaling, weights = weights[-CVtype[[i]][[k]]], nlv=nlv)
            predictions <- predict(repmodel, Xlisttest, nlv = nlv)
            reppredCV[CVtype[[i]][[k]],] <- predictions$posterior[order(CVtype[[i]][[k]]),]
            repYpredCV[CVtype[[i]][[k]],] <- predictions$pred[order(CVtype[[i]][[k]]),]
          }
          SqErrCV <- (dummy(Ypermut)$Y-reppredCV)^2
          repnlvpermut[(repnlvpermut$rep==i)&(repnlvpermut$permut==j),"rmse"] <- sqrt(mean(SqErrCV))
          repnlvYpermut[(repnlvYpermut$rep==i)&(repnlvYpermut$permut==j),"err"] <- err(repYpredCV,Ypermut)
        }
        if(criterion=="rmse"){
          res_permut[j] <- apply(repnlvpermut[repnlvpermut$permut==j,"rmse", drop = FALSE], FUN = mean, MARGIN = 2, na.rm = TRUE)
        }
        if(criterion=="err"){
          res_permut[j] <- apply(repnlvYpermut[repnlvYpermut$permut==j,"err", drop = FALSE], FUN = mean, MARGIN = 2, na.rm = TRUE)
        }
      }
    
      #plot(permut_dyssimilarity,res_permut, ylim=c(0,1))
      #abline(h=res_permut[npermut+1])
      
      # Calcul de la p-valeur
      pval_permut <- mean(res_permut[1:npermut]<res_permut[npermut+1],na.rm=TRUE)
      
      #output
      respermuttesttable <- cbind.data.frame(permut_dyssimilarity,res_permut)
      if(is.null(outputfilename)==FALSE){
        write.table(data.frame(respermuttesttable), 
                    file = paste0(outputfilename,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                    eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                    col.names = TRUE, qmethod = c("escape", "double"),
                    fileEncoding = "")
      }
    }
    
    if(method == c("plsr", "plsrda","plslda","plsqda")[3]){
      
      permut_dyssimilarity <- c()
      res_permut <- c()
      for(j in 1:(npermut+1)){
        Ypermut <- Y[permut_list[[j]],,drop=FALSE]
        permut_dyssimilarity[j] <- err(Ypermut,Y)
        
        repnlvpermut <- cbind.data.frame(expand.grid(permut=1:(npermut+1), rep=1:nbrep),rmse = rep(NA,(nbrep*(1+npermut))))
        repnlvYpermut <- cbind.data.frame(expand.grid(permut=1:(npermut+1), rep=1:nbrep),err = rep(NA,(nbrep*(1+npermut))))
        
        for(i in 1:nbrep){
          reppredCV <- matrix(NA, nrow=n, ncol = ncol(dummy(as.matrix(Y))$Y))
          repYpredCV <- matrix(NA, nrow=n, ncol = ncol(Y))
          for(k in 1:length(CVtype[[i]])){
            Xlisttrain <- X[-CVtype[[i]][[k]],,drop=FALSE]
            Ytrain <- Ypermut[-CVtype[[i]][[k]],,drop=FALSE]
            Xlisttest <- X[CVtype[[i]][[k]],,drop=FALSE]
            Ytest <- Ypermut[CVtype[[i]][[k]],,drop=FALSE]
            repmodel <- plslda(X=Xlisttrain, y=Ytrain, Xscaling = Xscaling, Yscaling = Yscaling, weights = weights[-CVtype[[i]][[k]]], nlv=nlv, prior=prior)
            predictions <- predict(repmodel, Xlisttest, nlv = nlv)
            reppredCV[CVtype[[i]][[k]],] <- predictions$posterior[order(CVtype[[i]][[k]]),]
            repYpredCV[CVtype[[i]][[k]],] <- predictions$pred[order(CVtype[[i]][[k]]),]
          }
          SqErrCV <- (dummy(Ypermut)$Y-reppredCV)^2
          repnlvpermut[(repnlvpermut$rep==i)&(repnlvpermut$permut==j),"rmse"] <- sqrt(mean(SqErrCV))
          repnlvYpermut[(repnlvYpermut$rep==i)&(repnlvYpermut$permut==j),"err"] <- err(repYpredCV,Ypermut)
        }
        if(criterion=="rmse"){
          res_permut[j] <- apply(repnlvpermut[repnlvpermut$permut==j,"rmse", drop = FALSE], FUN = mean, MARGIN = 2, na.rm = TRUE)
        }
        if(criterion=="err"){
          res_permut[j] <- apply(repnlvYpermut[repnlvYpermut$permut==j,"err", drop = FALSE], FUN = mean, MARGIN = 2, na.rm = TRUE)
        }
      }
  
      #plot(permut_dyssimilarity,res_permut, ylim=c(0,1))
      #abline(h=res_permut[npermut+1])
      
      # Calcul de la p-valeur
      pval_permut <- mean(res_permut[1:npermut]<res_permut[npermut+1],na.rm=TRUE)
      
      #output
      respermuttesttable <- cbind.data.frame(permut_dyssimilarity,res_permut)
      if(is.null(outputfilename)==FALSE){
        write.table(data.frame(respermuttesttable), 
                    file = paste0(outputfilename,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                    eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                    col.names = TRUE, qmethod = c("escape", "double"),
                    fileEncoding = "")
      }
    }
    
    if(method == c("plsr", "plsrda","plslda","plsqda")[4]){
      
      permut_dyssimilarity <- c()
      res_permut <- c()
      for(j in 1:(npermut+1)){
        Ypermut <- Y[permut_list[[j]],,drop=FALSE]
        permut_dyssimilarity[j] <- err(Ypermut,Y)
        
        repnlvpermut <- cbind.data.frame(expand.grid(permut=1:(npermut+1), rep=1:nbrep),rmse = rep(NA,(nbrep*(1+npermut))))
        repnlvYpermut <- cbind.data.frame(expand.grid(permut=1:(npermut+1), rep=1:nbrep),err = rep(NA,(nbrep*(1+npermut))))
        
        for(i in 1:nbrep){
          reppredCV <- matrix(NA, nrow=n, ncol = ncol(dummy(as.matrix(Y))$Y))
          repYpredCV <- matrix(NA, nrow=n, ncol = ncol(Y))
          for(k in 1:length(CVtype[[i]])){
            Xlisttrain <- X[-CVtype[[i]][[k]],,drop=FALSE]
            Ytrain <- Ypermut[-CVtype[[i]][[k]],,drop=FALSE]
            Xlisttest <- X[CVtype[[i]][[k]],,drop=FALSE]
            Ytest <- Ypermut[CVtype[[i]][[k]],,drop=FALSE]
            repmodel <- plsqda(X=Xlisttrain, y=Ytrain, Xscaling = Xscaling, Yscaling = Yscaling, weights = weights[-CVtype[[i]][[k]]], nlv=nlv, prior=prior)
            predictions <- predict(repmodel, Xlisttest, nlv = nlv)
            reppredCV[CVtype[[i]][[k]],] <- predictions$posterior[order(CVtype[[i]][[k]]),]
            repYpredCV[CVtype[[i]][[k]],] <- predictions$pred[order(CVtype[[i]][[k]]),]
          }
          SqErrCV <- (dummy(Ypermut)$Y-reppredCV)^2
          repnlvpermut[(repnlvpermut$rep==i)&(repnlvpermut$permut==j),"rmse"] <- sqrt(mean(SqErrCV))
          repnlvYpermut[(repnlvYpermut$rep==i)&(repnlvYpermut$permut==j),"err"] <- err(repYpredCV,Ypermut)
        }
        if(criterion=="rmse"){
          res_permut[j] <- apply(repnlvpermut[repnlvpermut$permut==j,"rmse", drop = FALSE], FUN = mean, MARGIN = 2, na.rm = TRUE)
        }
        if(criterion=="err"){
          res_permut[j] <- apply(repnlvYpermut[repnlvYpermut$permut==j,"err", drop = FALSE], FUN = mean, MARGIN = 2, na.rm = TRUE)
        }
      }
      
      #plot(permut_dyssimilarity,res_permut, ylim=c(0,1))
      #abline(h=res_permut[npermut+1])
      
      # Calcul de la p-valeur
      pval_permut <- mean(res_permut[1:npermut]<res_permut[npermut+1],na.rm=TRUE)
      
      #output
      respermuttesttable <- cbind.data.frame(permut_dyssimilarity,res_permut)
      
      if(is.null(outputfilename)==FALSE){
        write.table(data.frame(respermuttesttable), 
                    file = paste0(outputfilename,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                    eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                    col.names = TRUE, qmethod = c("escape", "double"),
                    fileEncoding = "")
      }
    }
    return(respermuttesttable)
  }
  
  ### prediction
  
  if (step == "prediction"){
    if(is.null(newX)){newX <- X}
    if(method == c("plsr", "plsrda","plslda","plsqda")[1]){
      resmodel <- plskern(X = X, 
                         Y = Y, 
                         Xscaling = Xscaling, 
                         Yscaling = Yscaling, 
                         weights = weights, 
                         nlv = nlv)
      respredict <- cbind(transform(resmodel, newX),pred = predict(resmodel, newX)$pred)
      if(is.null(outputfilename)==FALSE){
        write.table(data.frame(rownames(respredict),respredict), 
                    file = paste0(outputfilename,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                    eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                    col.names = TRUE, qmethod = c("escape", "double"),
                    fileEncoding = "")
      }
    }
    if(method == c("plsr", "plsrda","plslda","plsqda")[2]){
      resmodel <- plsrda(X = X, 
                           y = Y, 
                           Xscaling = Xscaling, 
                           Yscaling = Yscaling, 
                           weights = weights, 
                           nlv = nlv)
      respredict <- cbind.data.frame(transform(resmodel$fm, newX),pred=predict(resmodel, newX))
      if(is.null(outputfilename)==FALSE){
        write.table(data.frame(rownames(respredict),respredict), 
                    file = paste0(outputfilename,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                    eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                    col.names = TRUE, qmethod = c("escape", "double"),
                    fileEncoding = "")
      }
    }
    if(method == c("plsr", "plsrda","plslda","plsqda")[3]){
      resmodel <- plslda(X = X, 
                           y = Y, 
                           Xscaling = Xscaling, 
                           Yscaling = Yscaling, 
                           weights = weights, 
                           nlv = nlv,
                           prior = prior)
      respredict <- cbind.data.frame(transform(resmodel$fm[[1]], newX),pred=predict(resmodel, newX))
      if(is.null(outputfilename)==FALSE){
        write.table(data.frame(rownames(respredict),respredict), 
                    file = paste0(outputfilename,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                    eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                    col.names = TRUE, qmethod = c("escape", "double"),
                    fileEncoding = "")
      }
    }
    if(method == c("plsr", "plsrda","plslda","plsqda")[4]){
      resmodel <- plsqda(X = X, 
                           y = Y, 
                           Xscaling = Xscaling, 
                           Yscaling = Yscaling, 
                           weights = weights, 
                           nlv = nlv,
                           prior = prior)
      respredict <- cbind.data.frame(transform(resmodel$fm[[1]], newX),pred=predict(resmodel, newX))
      if(is.null(outputfilename)==FALSE){
        write.table(data.frame(rownames(respredict),respredict), 
                    file = paste0(outputfilename,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                    eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                    col.names = TRUE, qmethod = c("escape", "double"),
                    fileEncoding = "")
      }
    }
    return(respredict)
  }
}
