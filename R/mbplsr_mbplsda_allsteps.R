
mbplsr_mbplsda_allsteps <- function(import = c("R","ChemFlow","W4M")[1],
                                  Xlist, Xnames = NULL, Xscaling = c("none","pareto","sd")[1], 
                                  Y, Yscaling = c("none","pareto","sd")[1], weights = NULL,
                                  newXlist = NULL, newXnames = NULL,
                                  
                                  method = c("mbplsr", "mbplsrda","mbplslda","mbplsqda")[1],
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
                                  
                                  outputfilename = NULL
){
  
  # IMPORTATION DES DONNEES -------------------------------------------------------------------------------------------
  
  ### librairies
  library(rchemo)
  
  if(is.null(Xnames)==TRUE){Xnames <- paste0(rep("X",length(Xlist)),1:length(Xlist))}
  if((is.null(newXlist)==FALSE)&(is.null(newXnames)==TRUE)){newXnames <- paste0(rep("Xnew",length(newXlist)),1:length(newXlist))}
  
  if(length(Xscaling) == 1){Xscaling = rep(Xscaling, length(Xlist))}
  
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
  
      ### Y
  if(import == "R"){# Y matrix n x q, in the Global Environment, with observation names in rownames
    Y           <- as.matrix(Y)
  }
  if(import %in% c("ChemFlow","W4M")){# Y matrix n x q, tabulated table (.txt), with observation names in the first column
    Y           <- read.table(Y, sep="\t", header=TRUE, na.strings = c("","NA"), check.names=FALSE)
    rownames(Y) <- Y[,1]
    Y[,1]       <- NULL
    Y           <- as.matrix(Y)
  }
  
  n  <- nrow(Y)
  
  ### X 
  XBlockNb   <- length(Xlist)

  for(i in 1:XBlockNb){
    if(import == "R"){
      Xlist[[i]] <- Xlist[[i]]
    }
    if(import == "ChemFlow"){# X matrix n x p, tabulated tables (.txt), with observation names in the first column
      Xlist[[i]]           <- read.table(Xlist[[i]],sep="\t",dec=".",header=TRUE, na.strings=c("","NA"), check.names=FALSE)
      rownames(Xlist[[i]]) <- Xlist[[i]][,1]
      Xlist[[i]][,1]       <- NULL
    }
    if(import == "W4M"){# X matrix p x n tabulated tables (.txt), with observation names in the header of X and the first column of Y
      Xlist[[i]]           <- read.table(Xlist[[i]],sep="\t",dec=".",header=TRUE, na.strings=c("","NA"), check.names=FALSE)
      Xlist[[i]]           <- t(Xlist[[i]])
      rownames(Xlist[[i]]) <- Xlist[[i]][,1]
      Xlist[[i]][,1]       <- NULL
    }
    
    SamplesNotFound1 <- NULL
    SamplesNotFound2 <- NULL
    if(length(which(rownames(Xlist[[i]])%in%rownames(Y)))<length(rownames(Xlist[[i]]))){
      SamplesNotFound1 <- rownames(Xlist[[i]])[(rownames(Xlist[[i]])%in%rownames(Y))==FALSE]
    } 
    if(length(which(rownames(Y)%in%rownames(Xlist[[i]])))<length(rownames(Y))){
      SamplesNotFound2 <- rownames(Y)[(rownames(Y)%in%(rownames(Xlist[[i]])))==FALSE]
    }
    if(is.null(SamplesNotFound1)==FALSE & is.null(SamplesNotFound2)==FALSE){
      stop("\n- - - - - - - - -\n",paste0("Some samples of ", Xnames[i], " are not found in the Y: \n"),paste(SamplesNotFound1, collapse="\n"),
           "\n All samples must be in the Y. Please, check your data.",
           "\n\n",paste0("Some samples of the Y are not found in ", Xnames[i], ": \n"), paste(SamplesNotFound2, collapse="\n"),
           "\n Please remove them or do missing value imputation.","\n- - - - - - - - -\n")
    }
    if(is.null(SamplesNotFound1)==FALSE & is.null(SamplesNotFound2)==TRUE){
      stop("\n- - - - - - - - -\n",paste0("Some samples of ", Xnames[i], " are not found in the Y: \n"),paste(SamplesNotFound1, collapse="\n"),
           "\n All samples must be in the Y. Please, check your data.","\n- - - - - - - - -\n")
    }
    if(is.null(SamplesNotFound1)==TRUE & is.null(SamplesNotFound2)==FALSE){
      stop("\n- - - - - - - - -\n",paste0("Some samples of the Y are not found in ", Xnames[i], ": \n"), paste(SamplesNotFound2, collapse="\n"),
           "\n Please remove them or do missing value imputation.","\n- - - - - - - - -\n")
    }
    Xlist[[i]] <- data.frame(Xlist[[i]])
  }

  if((step=="prediction")&(is.null(newXlist)==TRUE)){
    newXlist <- Xlist
    newXnames <- Xnames
  }
  if((step=="prediction")&(is.null(newXlist)==FALSE)){
    ### newXlist
    newXBlockNb   <- length(newXlist)
    
    for(i in 1:newXBlockNb){
      if(import == "R"){
        newXlist[[i]] <- newXlist[[i]]
      }
      if(import == "ChemFlow"){
        newXlist[[i]]           <- read.table(newXlist[[i]],sep="\t",dec=".",header=TRUE, na.strings=c("","NA"), check.names=FALSE)
        rownames(newXlist[[i]]) <- newXlist[[i]][,1]
        newXlist[[i]][,1]       <- NULL
      }
      if(import == "W4M"){
        newXlist[[i]]           <- read.table(newXlist[[i]],sep="\t",dec=".",header=TRUE, na.strings=c("","NA"), check.names=FALSE)
        newXlist[[i]]           <- t(newXlist[[i]])
        rownames(newXlist[[i]]) <- newXlist[[i]][,1]
        newXlist[[i]][,1]       <- NULL
      }
      
      VariablesNotFound1 <- NULL
      VariablesNotFound2 <- NULL
      if(length(which(colnames(newXlist[[i]])%in%colnames(Xlist[[i]])))<length(colnames(newXlist[[i]]))){
        VariablesNotFound1 <- colnames(newXlist[[i]])[(colnames(newXlist[[i]])%in%colnames(Xlist[[i]]))==FALSE]
      } 
      if(length(which(colnames(Xlist[[i]])%in%colnames(newXlist[[i]])))<length(colnames(Xlist[[i]]))){
        VariablesNotFound2 <- colnames(Xlist[[i]])[(colnames(Xlist[[i]])%in%(colnames(newXlist[[i]])))==FALSE]
      }
      if(is.null(VariablesNotFound1)==FALSE & is.null(VariablesNotFound2)==FALSE){
        stop("\n- - - - - - - - -\n",paste0("Some variables of ", newXnames[i], " are not found in the XBlock: \n"),paste(VariablesNotFound1, collapse="\n"),
             "\n All variables must be in the XBlock. Please, check your data.",
             "\n\n",paste0("Some variables of the XBlock are not found in ", newXnames[i], ": \n"), paste(VariablesNotFound2, collapse="\n"),
             "\n Please remove them or do missing value imputation.","\n- - - - - - - - -\n")
      }
      if(is.null(VariablesNotFound1)==FALSE & is.null(VariablesNotFound2)==TRUE){
        stop("\n- - - - - - - - -\n",paste0("Some variables of ", newXnames[i], " are not found in the XBlock: \n"),paste(VariablesNotFound1, collapse="\n"),
             "\n All variables must be in the XBlock. Please, check your data.","\n- - - - - - - - -\n")
      }
      if(is.null(VariablesNotFound1)==TRUE & is.null(VariablesNotFound2)==FALSE){
        stop("\n- - - - - - - - -\n",paste0("Some variables of the XBlock are not found in ", newXnames[i], ": \n"), paste(VariablesNotFound2, collapse="\n"),
             "\n Please remove them or do missing value imputation.","\n- - - - - - - - -\n")
      }
      newXlist[[i]]      <- data.frame(newXlist[[i]])
    }
  }

  
  ### nlvtest
  
  if (step == "nlvtest"){

    # partition
    set.seed(seed = seed)
    if(cvmethod=="loo"){
      nbrep=1
      CVtype <- segmkf(n = nrow(Xlist[[1]]), K = nrow(Xlist[[1]]), type = "random", nrep = nbrep)
    }
    if((cvmethod=="kfolds") & (is.null(samplingk)==TRUE)){
      CVtype <- segmkf(n = nrow(Xlist[[1]]), K = nfolds, type = "random", nrep = nbrep)
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
    
    if(method == c("mbplsr", "mbplsrda","mbplslda","mbplsqda")[1]){
      
      repnlvtest <- cbind.data.frame(expand.grid(nblv=0:nlv, rep=1:nbrep),rmse = rep(NA,(nbrep*(1+nlv))))
      for(i in 1:nbrep){
        reppredCV <- lapply(1:(nlv+1), function(l) matrix(NA, nrow=n, ncol = ncol(Y)))
        for(k in 1:length(CVtype[[i]])){
          Xlisttrain <- lapply(1:length(Xlist),function(x) Xlist[[x]][-CVtype[[i]][[k]],,drop=FALSE])
          Ytrain <- Y[-CVtype[[i]][[k]],,drop=FALSE]
          Xlisttest <- lapply(1:length(Xlist),function(x) Xlist[[x]][CVtype[[i]][[k]],,drop=FALSE])
          Ytest <- Y[CVtype[[i]][[k]],,drop=FALSE]
          repmodel <- mbplsr(Xlist=Xlisttrain, Y=Ytrain, Xscaling = Xscaling, Yscaling = Yscaling, weights = weights[-CVtype[[i]][[k]]], nlv=nlv)
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
        diff<- sapply(1:nlv, function(i) resnlvtesttable$rmse_mean[i]<=resnlvtesttable$rmse_mean[i+1])
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
    
    
    if(method == c("mbplsr", "mbplsrda","mbplslda","mbplsqda")[2]){
      
      repnlvtest <- cbind.data.frame(expand.grid(nblv=0:nlv, rep=1:nbrep),rmse = rep(NA,(nbrep*(1+nlv))))
      repYnlvtest <- cbind.data.frame(expand.grid(nblv=0:nlv, rep=1:nbrep),err = rep(NA,(nbrep*(1+nlv))))
      for(i in 1:nbrep){
        reppredCV <- lapply(1:(nlv+1), function(l) matrix(NA, nrow=n, ncol = ncol(dummy(Y)$Y)))
        repYpredCV <- lapply(1:(nlv+1), function(l) matrix(NA, nrow=n, ncol = ncol(Y)))
        for(k in 1:length(CVtype[[i]])){
          Xlisttrain <- lapply(1:length(Xlist),function(x) Xlist[[x]][-CVtype[[i]][[k]],,drop=FALSE])
          Ytrain <- Y[-CVtype[[i]][[k]],,drop=FALSE]
          Xlisttest <- lapply(1:length(Xlist),function(x) Xlist[[x]][CVtype[[i]][[k]],,drop=FALSE])
          Ytest <- Y[CVtype[[i]][[k]],,drop=FALSE]
          repmodel <- mbplsrda(Xlist=Xlisttrain, y=Ytrain, Xscaling = Xscaling, Yscaling = Yscaling, weights = weights[-CVtype[[i]][[k]]], nlv=nlv)
          predictions <- predict(repmodel, Xlisttest, nlv = 0:nlv)
          for(l in 1:(nlv+1)){
            reppredCV[[l]][CVtype[[i]][[k]],] <- predictions$posterior[[l]][order(CVtype[[i]][[k]]),]
            repYpredCV[[l]][CVtype[[i]][[k]],] <- predictions$pred[[l]][order(CVtype[[i]][[k]]),]
          }
        }
        for(l in 1:(nlv+1)){
          SqErrCV <- (dummy(Y)$Y-reppredCV[[l]])^2
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
          diff<- sapply(1:nlv, function(i) resnlvtesttable$rmse_mean[i]<=resnlvtesttable$rmse_mean[i+1])
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
          diff<- sapply(1:nlv, function(i) resnlvtesttable$err_mean[i]<=resnlvtesttable$err_mean[i+1])
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
    
    
    if(method == c("mbplsr", "mbplsrda","mbplslda","mbplsqda")[3]){
      
      repnlvtest <- cbind.data.frame(expand.grid(nblv=0:nlv, rep=1:nbrep),rmse = rep(NA,(nbrep*(1+nlv))))
      repYnlvtest <- cbind.data.frame(expand.grid(nblv=0:nlv, rep=1:nbrep),err = rep(NA,(nbrep*(1+nlv))))
      for(i in 1:nbrep){
        reppredCV <- lapply(1:(nlv+1), function(l) matrix(NA, nrow=n, ncol = ncol(dummy(Y)$Y)))
        repYpredCV <- lapply(1:(nlv+1), function(l) matrix(NA, nrow=n, ncol = ncol(Y)))
        for(k in 1:length(CVtype[[i]])){
          Xlisttrain <- lapply(1:length(Xlist),function(x) Xlist[[x]][-CVtype[[i]][[k]],,drop=FALSE])
          Ytrain <- Y[-CVtype[[i]][[k]],,drop=FALSE]
          Xlisttest <- lapply(1:length(Xlist),function(x) Xlist[[x]][CVtype[[i]][[k]],,drop=FALSE])
          Ytest <- Y[CVtype[[i]][[k]],,drop=FALSE]
          #nlv = 0
          repmodel0lv <- mbplsrda(Xlist=Xlisttrain, y=Ytrain, Xscaling = Xscaling, Yscaling = Yscaling, weights = weights[-CVtype[[i]][[k]]], nlv=1)
          predictions0lv <- predict(repmodel0lv, Xlisttest, nlv = 0)
          reppredCV[[1]][CVtype[[i]][[k]],] <- predictions0lv$posterior[order(CVtype[[i]][[k]]),]
          repYpredCV[[1]][CVtype[[i]][[k]],] <- predictions0lv$pred[order(CVtype[[i]][[k]]),]
          #nlv>0
          repmodel <- mbplslda(Xlist=Xlisttrain, y=Ytrain, Xscaling = Xscaling, Yscaling = Yscaling, weights = weights[-CVtype[[i]][[k]]], nlv=nlv, prior = prior)
          predictions <- predict(repmodel, Xlisttest, nlv = 1:nlv)
          for(l in 2:(nlv+1)){
            reppredCV[[l]][CVtype[[i]][[k]],] <- predictions$posterior[[(l-1)]][order(CVtype[[i]][[k]]),]
            repYpredCV[[l]][CVtype[[i]][[k]],] <- predictions$pred[[(l-1)]][order(CVtype[[i]][[k]]),]
          }
        }
        for(l in 1:(nlv+1)){
          SqErrCV <- (dummy(Y)$Y-reppredCV[[l]])^2
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
          diff<- sapply(1:nlv, function(i) resnlvtesttable$rmse_mean[i]<=resnlvtesttable$rmse_mean[i+1])
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
          diff<- sapply(1:nlv, function(i) resnlvtesttable$err_mean[i]<=resnlvtesttable$err_mean[i+1])
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
    
    if(method == c("mbplsr", "mbplsrda","mbplslda","mbplsqda")[4]){
      
      repnlvtest <- cbind.data.frame(expand.grid(nblv=0:nlv, rep=1:nbrep),rmse = rep(NA,(nbrep*(1+nlv))))
      repYnlvtest <- cbind.data.frame(expand.grid(nblv=0:nlv, rep=1:nbrep),err = rep(NA,(nbrep*(1+nlv))))
      for(i in 1:nbrep){
        reppredCV <- lapply(1:(nlv+1), function(l) matrix(NA, nrow=n, ncol = ncol(dummy(Y)$Y)))
        repYpredCV <- lapply(1:(nlv+1), function(l) matrix(NA, nrow=n, ncol = ncol(Y)))
        for(k in 1:length(CVtype[[i]])){
          Xlisttrain <- lapply(1:length(Xlist),function(x) Xlist[[x]][-CVtype[[i]][[k]],,drop=FALSE])
          Ytrain <- Y[-CVtype[[i]][[k]],,drop=FALSE]
          Xlisttest <- lapply(1:length(Xlist),function(x) Xlist[[x]][CVtype[[i]][[k]],,drop=FALSE])
          Ytest <- Y[CVtype[[i]][[k]],,drop=FALSE]
          #nlv = 0
          repmodel0lv <- mbplsrda(Xlist=Xlisttrain, y=Ytrain, Xscaling = Xscaling, Yscaling = Yscaling, weights = weights[-CVtype[[i]][[k]]], nlv=1)
          predictions0lv <- predict(repmodel0lv, Xlisttest, nlv = 0)
          reppredCV[[1]][CVtype[[i]][[k]],] <- predictions0lv$posterior[order(CVtype[[i]][[k]]),]
          repYpredCV[[1]][CVtype[[i]][[k]],] <- predictions0lv$pred[order(CVtype[[i]][[k]]),]
          #nlv>0
          repmodel <- mbplsqda(Xlist=Xlisttrain, y=Ytrain, Xscaling = Xscaling, Yscaling = Yscaling, weights = weights[-CVtype[[i]][[k]]], nlv=nlv, prior = prior)
          predictions <- predict(repmodel, Xlisttest, nlv = 1:nlv)
          for(l in 2:(nlv+1)){
            reppredCV[[l]][CVtype[[i]][[k]],] <- predictions$posterior[[(l-1)]][order(CVtype[[i]][[k]]),]
            repYpredCV[[l]][CVtype[[i]][[k]],] <- predictions$pred[[(l-1)]][order(CVtype[[i]][[k]]),]
          }
        }
        for(l in 1:(nlv+1)){
          SqErrCV <- (dummy(Y)$Y-reppredCV[[l]])^2
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
          diff<- sapply(1:nlv, function(i) resnlvtesttable$rmse_mean[i]<=resnlvtesttable$rmse_mean[i+1])
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
          diff<- sapply(1:nlv, function(i) resnlvtesttable$err_mean[i]<=resnlvtesttable$err_mean[i+1])
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
    if(method == c("mbplsr", "mbplsrda","mbplslda","mbplsqda")[1]){
      resmodel <- mbplsr(Xlist = Xlist, 
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

    
    if(method == c("mbplsr", "mbplsrda","mbplslda","mbplsqda")[2]){
      resmodel <- mbplsrda(Xlist = Xlist, 
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
    
    if(method == c("mbplsr", "mbplsrda","mbplslda","mbplsqda")[3]){
      resmodel <- mbplslda(Xlist = Xlist, 
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
    
    if(method == c("mbplsr", "mbplsrda","mbplslda","mbplsqda")[4]){
      resmodel <- mbplsqda(Xlist = Xlist, 
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
      CVtype <- segmkf(n = nrow(Xlist[[1]]), K = nrow(Xlist[[1]]), type = "random", nrep = nbrep)
    }
    if((cvmethod=="kfolds") & (is.null(samplingk)==TRUE)){
      CVtype <- segmkf(n = nrow(Xlist[[1]]), K = nfolds, type = "random", nrep = nbrep)
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
    
    if(method == c("mbplsr", "mbplsrda","mbplslda","mbplsqda")[1]){
      
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
            Xlisttrain <- lapply(1:length(Xlist),function(x) Xlist[[x]][-CVtype[[i]][[k]],,drop=FALSE])
            Ytrain <- Ypermut[-CVtype[[i]][[k]],,drop=FALSE]
            Xlisttest <- lapply(1:length(Xlist),function(x) Xlist[[x]][CVtype[[i]][[k]],,drop=FALSE])
            Ytest <- Ypermut[CVtype[[i]][[k]],,drop=FALSE]
            repmodel <- mbplsr(Xlist=Xlisttrain, Y=Ytrain, Xscaling = Xscaling, Yscaling = Yscaling, weights = weights[-CVtype[[i]][[k]]], nlv=nlv)
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
    
    if(method == c("mbplsr", "mbplsrda","mbplslda","mbplsqda")[2]){
      
      permut_dyssimilarity <- c()
      res_permut <- c()
      for(j in 1:(npermut+1)){
        Ypermut <- Y[permut_list[[j]],,drop=FALSE]
        permut_dyssimilarity[j] <- err(Ypermut,Y)
        
        repnlvpermut <- cbind.data.frame(expand.grid(permut=1:(npermut+1), rep=1:nbrep),rmse = rep(NA,(nbrep*(1+npermut))))
        repnlvYpermut <- cbind.data.frame(expand.grid(permut=1:(npermut+1), rep=1:nbrep),err = rep(NA,(nbrep*(1+npermut))))
        
        for(i in 1:nbrep){
          reppredCV <- matrix(NA, nrow=n, ncol = ncol(dummy(Y)$Y))
          repYpredCV <- matrix(NA, nrow=n, ncol = ncol(Y))
          for(k in 1:length(CVtype[[i]])){
            Xlisttrain <- lapply(1:length(Xlist),function(x) Xlist[[x]][-CVtype[[i]][[k]],,drop=FALSE])
            Ytrain <- Ypermut[-CVtype[[i]][[k]],,drop=FALSE]
            Xlisttest <- lapply(1:length(Xlist),function(x) Xlist[[x]][CVtype[[i]][[k]],,drop=FALSE])
            Ytest <- Ypermut[CVtype[[i]][[k]],,drop=FALSE]
            repmodel <- mbplsrda(Xlist=Xlisttrain, y=Ytrain, Xscaling = Xscaling, Yscaling = Yscaling, weights = weights[-CVtype[[i]][[k]]], nlv=nlv)
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
    
    if(method == c("mbplsr", "mbplsrda","mbplslda","mbplsqda")[3]){
      
      permut_dyssimilarity <- c()
      res_permut <- c()
      for(j in 1:(npermut+1)){
        Ypermut <- Y[permut_list[[j]],,drop=FALSE]
        permut_dyssimilarity[j] <- err(Ypermut,Y)
        
        repnlvpermut <- cbind.data.frame(expand.grid(permut=1:(npermut+1), rep=1:nbrep),rmse = rep(NA,(nbrep*(1+npermut))))
        repnlvYpermut <- cbind.data.frame(expand.grid(permut=1:(npermut+1), rep=1:nbrep),err = rep(NA,(nbrep*(1+npermut))))
        
        for(i in 1:nbrep){
          reppredCV <- matrix(NA, nrow=n, ncol = ncol(dummy(Y)$Y))
          repYpredCV <- matrix(NA, nrow=n, ncol = ncol(Y))
          for(k in 1:length(CVtype[[i]])){
            Xlisttrain <- lapply(1:length(Xlist),function(x) Xlist[[x]][-CVtype[[i]][[k]],,drop=FALSE])
            Ytrain <- Ypermut[-CVtype[[i]][[k]],,drop=FALSE]
            Xlisttest <- lapply(1:length(Xlist),function(x) Xlist[[x]][CVtype[[i]][[k]],,drop=FALSE])
            Ytest <- Ypermut[CVtype[[i]][[k]],,drop=FALSE]
            repmodel <- mbplslda(Xlist=Xlisttrain, y=Ytrain, Xscaling = Xscaling, Yscaling = Yscaling, weights = weights[-CVtype[[i]][[k]]], nlv=nlv, prior=prior)
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
    
    if(method == c("mbplsr", "mbplsrda","mbplslda","mbplsqda")[4]){
      
      permut_dyssimilarity <- c()
      res_permut <- c()
      for(j in 1:(npermut+1)){
        Ypermut <- Y[permut_list[[j]],,drop=FALSE]
        permut_dyssimilarity[j] <- err(Ypermut,Y)
        
        repnlvpermut <- cbind.data.frame(expand.grid(permut=1:(npermut+1), rep=1:nbrep),rmse = rep(NA,(nbrep*(1+npermut))))
        repnlvYpermut <- cbind.data.frame(expand.grid(permut=1:(npermut+1), rep=1:nbrep),err = rep(NA,(nbrep*(1+npermut))))
        
        for(i in 1:nbrep){
          reppredCV <- matrix(NA, nrow=n, ncol = ncol(dummy(Y)$Y))
          repYpredCV <- matrix(NA, nrow=n, ncol = ncol(Y))
          for(k in 1:length(CVtype[[i]])){
            Xlisttrain <- lapply(1:length(Xlist),function(x) Xlist[[x]][-CVtype[[i]][[k]],,drop=FALSE])
            Ytrain <- Ypermut[-CVtype[[i]][[k]],,drop=FALSE]
            Xlisttest <- lapply(1:length(Xlist),function(x) Xlist[[x]][CVtype[[i]][[k]],,drop=FALSE])
            Ytest <- Ypermut[CVtype[[i]][[k]],,drop=FALSE]
            repmodel <- mbplsqda(Xlist=Xlisttrain, y=Ytrain, Xscaling = Xscaling, Yscaling = Yscaling, weights = weights[-CVtype[[i]][[k]]], nlv=nlv, prior=prior)
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
    if(is.null(newXlist)){newXlist <- Xlist}
    if(method == c("mbplsr", "mbplsrda","mbplslda","mbplsqda")[1]){
      resmodel <- mbplsr(Xlist = Xlist, 
                         Y = Y, 
                         Xscaling = Xscaling, 
                         Yscaling = Yscaling, 
                         weights = weights, 
                         nlv = nlv)
      respredict <- cbind(transform(resmodel, newXlist),pred = predict(resmodel, newXlist)$pred)
      if(is.null(outputfilename)==FALSE){
        write.table(data.frame(rownames(respredict),respredict), 
                    file = paste0(outputfilename,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                    eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                    col.names = TRUE, qmethod = c("escape", "double"),
                    fileEncoding = "")
      }
      
    }
    if(method == c("mbplsr", "mbplsrda","mbplslda","mbplsqda")[2]){
      resmodel <- mbplsrda(Xlist = Xlist, 
                           y = Y, 
                           Xscaling = Xscaling, 
                           Yscaling = Yscaling, 
                           weights = weights, 
                           nlv = nlv)
      respredict <- cbind.data.frame(transform(resmodel$fm, newXlist),pred=predict(resmodel, newXlist))
      if(is.null(outputfilename)==FALSE){
        write.table(data.frame(rownames(respredict),respredict), 
                    file = paste0(outputfilename,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                    eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                    col.names = TRUE, qmethod = c("escape", "double"),
                    fileEncoding = "")
      }
      
    }
    if(method == c("mbplsr", "mbplsrda","mbplslda","mbplsqda")[3]){
      resmodel <- mbplslda(Xlist = Xlist, 
                           y = Y, 
                           Xscaling = Xscaling, 
                           Yscaling = Yscaling, 
                           weights = weights, 
                           nlv = nlv,
                           prior = prior)
      respredict <- cbind.data.frame(transform(resmodel$fm[[1]], newXlist),pred=predict(resmodel, newXlist))
      if(is.null(outputfilename)==FALSE){
        write.table(data.frame(rownames(respredict),respredict), 
                    file = paste0(outputfilename,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                    eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                    col.names = TRUE, qmethod = c("escape", "double"),
                    fileEncoding = "")
      }
      
    }
    if(method == c("mbplsr", "mbplsrda","mbplslda","mbplsqda")[4]){
      resmodel <- mbplsqda(Xlist = Xlist, 
                           y = Y, 
                           Xscaling = Xscaling, 
                           Yscaling = Yscaling, 
                           weights = weights, 
                           nlv = nlv,
                           prior = prior)
      respredict <- cbind.data.frame(transform(resmodel$fm[[1]], newXlist),pred=predict(resmodel, newXlist))
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
