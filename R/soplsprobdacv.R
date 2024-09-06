
.soplsprobdacv <- function(funda, Xlist, y, Xscaling = c("none", "pareto", "sd")[1], Yscaling = c("none", "pareto", "sd")[1], weights = NULL, nlvlist = list(), prior = c("unif", "prop"), nbrep = 30, cvmethod = "kfolds", seed = 123, samplingk = NULL, nfolds = 7, optimisation = c("global","sequential")[1], criterion = c("err", "rmse")[1], selection = c("localmin","globalmin","1std")[1]){
  
  Y <- y

  # verifications
  
  if(length(Xlist) != length(nlvlist)){stop("the number of nbcolXlist or nlvlist elements is not correct")}
  if((is.vector(Y)==FALSE) & (sum(rownames(Xlist[[1]])!= rownames(Y))>0)){stop("the rownames of Xlist and Y are different")}
  if((is.vector(Y)==TRUE) & (nrow(Xlist[[1]])!= length(Y))){stop("the row numbers of Xlist and Y are different")}
  if((is.null(samplingk)==FALSE) & (nrow(Xlist[[1]])!= length(samplingk))){stop("the length of samplingk is not correct")}
  if((cvmethod!="loo") & (is.null(nfolds)==FALSE)){
    if(nrow(Xlist[[1]])<nfolds){stop("the value of nfolds is not correct")}
  }
  if((nbrep==1)&(selection=="1std")){stop("nbrep must be >1 when selection is '1std'")}
  if((cvmethod=="loo") & (selection=="1std")){stop("selection must not be '1std' when cvmethod is 'loo'")}
  if((cvmethod=="loo") & (nbrep>1)){stop("nbrep must not be set to 1 when cvmethod is 'loo'")}
  if((cvmethod=="kfolds") & (nfolds<=1)){stop("nfolds must not be >1 when cvmethod is 'kfolds'")}
  
  # argument computations
  Y <- .mat(Y)
  n <- nrow(Y)
  nXblocks <- length(Xlist)
  
  Yvariables    <- colnames(Y)

  inertieY <- sum(diag(t(as.matrix(scale(dummy(Y)$Y, scale=FALSE)))%*%(as.matrix(scale(dummy(Y)$Y, scale=FALSE)))))

  # partition
  set.seed(seed = seed)
  if(cvmethod=="loo"){CVtype<- segmkf(n = nrow(Xlist[[1]]), K = nrow(Xlist[[1]]), type = "random", nrep = 1)}#=as.list(1:n)}
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
  
  ############################################################################################################################
  ## GLOBAL optimisation
  if((optimisation=="global")|(nXblocks==1)){
    
    lvcombi <- as.matrix(expand.grid(nlvlist))
    colnames(lvcombi) <- paste0(rep("Xlist",nXblocks),1:nXblocks)
    rownames(lvcombi) <- paste0("lvcombi",1:nrow(lvcombi))
    
    # output tables
    res_errCV <- res_ExplVarCV <- res_rmseCV <- matrix(NA, nrow=nrow(lvcombi), ncol=2, dimnames=list(rownames(lvcombi),c("mean","sd")))
    Rep_errCV <- Rep_rmseCV <- Rep_ExplVarCV <- array(0, dim=c(nrow(lvcombi),nbrep, 1), dimnames = list(paste0("lvcombi",1:nrow(lvcombi)), paste0("rep",1:nbrep), "value"))

    #optimYpredCV <- list()

    for(i in 1:nbrep){
      
      for(j in 1:nrow(lvcombi)){
        reppredCV <- matrix(NA, nrow=n, ncol = length(unique(Y)), dimnames=list(rownames(Xlist),sort(unique(Y))))
        reppredYCV <- matrix(NA, nrow=n, ncol = length(Yvariables), dimnames=list(rownames(Xlist),Yvariables))
        for(k in 1:length(CVtype[[i]])){
          Xlisttrain <-lapply(1:length(Xlist),function(x) Xlist[[x]][-CVtype[[i]][[k]],,drop=FALSE])
          Ytrain <- Y[-CVtype[[i]][[k]],,drop=FALSE]
          Xlisttest <- lapply(1:length(Xlist),function(x) Xlist[[x]][CVtype[[i]][[k]],,drop=FALSE])
          Ytest <- Y[CVtype[[i]][[k]],,drop=FALSE]
          repmodel <- funda(Xlist=Xlisttrain, y=Ytrain, Xscaling = Xscaling, Yscaling = Yscaling, weights = weights[-CVtype[[i]][[k]]], nlv=lvcombi[j,], prior = prior)
          reppredCV[CVtype[[i]][[k]],] <- predict(repmodel, Xlisttest)$posterior[order(CVtype[[i]][[k]]),]
          reppredYCV[CVtype[[i]][[k]],] <- predict(repmodel, Xlisttest)$pred[order(CVtype[[i]][[k]]),]
        }
        SqErrCV <- (dummy(Y)$Y-reppredCV)^2

        Rep_rmseCV[j,i,]    <- sqrt(mean(SqErrCV))
        Rep_ExplVarCV[j,i,] <- mean(1-(apply(SqErrCV,2,mean)/(apply(dummy(Y)$Y,2,var)*(n-1)/n)))
        Rep_errCV[j,i,] <- mean(reppredYCV != Y)*100
      }
    }

    res_rmseCV[,"mean"] <- apply(Rep_rmseCV, c(1,3),mean, na.rm=T)
    res_rmseCV[,"sd"]   <- apply(Rep_rmseCV, c(1,3),sd, na.rm=T)

    res_ExplVarCV[,"mean"] <- apply(Rep_ExplVarCV, c(1,3),mean, na.rm=T)
    res_ExplVarCV[,"sd"]   <- apply(Rep_ExplVarCV, c(1,3),sd, na.rm=T)

    res_errCV[,"mean"] <- apply(Rep_errCV, c(1,3),mean, na.rm=T)# A VERIFIER QD MULTI Y
    res_errCV[,"sd"]   <- apply(Rep_errCV, c(1,3),sd, na.rm=T)# A VERIFIER QD MULTI Y
    
    if(criterion == "rmse"){
      # optim combination for each total number of components // res_rmse_Ysel[,"mean"]
   
      nlvsum=apply(lvcombi,1,sum)
      res_nlvsum_rmseCV_Ysel <- data.frame(
        index = sapply((unique(nlvsum)), function(i) which(rownames(res_rmseCV)==rownames(res_rmseCV[(nlvsum==i),,drop=FALSE][which.min(res_rmseCV[(nlvsum==i),"mean"]),,drop=FALSE]))),
        combnum = sapply((unique(nlvsum)), function(i) rownames(res_rmseCV[(nlvsum==i),,drop=FALSE][which.min(res_rmseCV[(nlvsum==i),"mean"]),,drop=FALSE])),
        t((sapply((unique(nlvsum)), function(i) unlist(data.frame(lvcombi[(nlvsum==i),,drop=FALSE][which.min(res_rmseCV[(nlvsum==i),"mean"]),,drop=FALSE]))))),
        totalnlv = unique(nlvsum),
        t((sapply((unique(nlvsum)), function(i) unlist(data.frame(res_rmseCV[(nlvsum==i),c("mean","sd"),drop=FALSE][which.min(res_rmseCV[(nlvsum==i),"mean"]),,drop=FALSE])))))
      )
      res_nlvsum_rmseCV_Ysel <- res_nlvsum_rmseCV_Ysel[order(res_nlvsum_rmseCV_Ysel[,"totalnlv"], decreasing=FALSE),]
      
      if(selection=="localmin"){
        if(nrow(res_nlvsum_rmseCV_Ysel)>1){
          # sign of the difference of accuracies to select the optim combination with the lower total number of components
          rtsdiff <- c(NA,sapply(2:nrow(res_nlvsum_rmseCV_Ysel), function(i)((res_nlvsum_rmseCV_Ysel$mean[i]-res_nlvsum_rmseCV_Ysel$mean[i-1])<=0)))
          kchoix <- min(res_nlvsum_rmseCV_Ysel[(which(rtsdiff==FALSE)-1)[1],"index"],res_nlvsum_rmseCV_Ysel[nrow(res_nlvsum_rmseCV_Ysel),"index"],na.rm=TRUE)
        }else{
          kchoix <- 1
        }
      }
      if(selection=="globalmin"){
        kchoix <- res_nlvsum_rmseCV_Ysel[which.min(res_nlvsum_rmseCV_Ysel$mean)[1],"index"]
      }
      if(selection=="1std"){
        if(nrow(res_nlvsum_rmseCV_Ysel)>1){
          # one standard error rule to select the optim number of components
          minmean    <- which.min(res_nlvsum_rmseCV_Ysel$mean)[1]
          threshmean <- res_nlvsum_rmseCV_Ysel$mean[minmean] + res_nlvsum_rmseCV_Ysel$sd[minmean]
          if((minmean == 1) | (sum(res_nlvsum_rmseCV_Ysel$mean[1:minmean]<=threshmean)==0)){
            kchoix <- 1
          }else{
            kchoix <- min(res_nlvsum_rmseCV_Ysel[which(res_nlvsum_rmseCV_Ysel$mean[1:minmean]<=threshmean),"index"])
          }
        }else{
          kchoix <- 1
        }
      }
    }
    if(criterion == "err"){
      # optim combination for each total number of components // res_err_Ysel[,"mean"]

      nlvsum=apply(lvcombi,1,sum)
      res_nlvsum_errCV_Ysel <- data.frame(
        index = sapply((unique(nlvsum)), function(i) which(rownames(res_errCV)==rownames(res_errCV[(nlvsum==i),,drop=FALSE][which.min(res_errCV[(nlvsum==i),"mean"]),,drop=FALSE]))),
        combnum = sapply((unique(nlvsum)), function(i) rownames(res_errCV[(nlvsum==i),,drop=FALSE][which.min(res_errCV[(nlvsum==i),"mean"]),,drop=FALSE])),
        t((sapply((unique(nlvsum)), function(i) unlist(data.frame(lvcombi[(nlvsum==i),,drop=FALSE][which.min(res_errCV[(nlvsum==i),"mean"]),,drop=FALSE]))))),
        totalnlv = unique(nlvsum),
        t((sapply((unique(nlvsum)), function(i) unlist(data.frame(res_errCV[(nlvsum==i),c("mean","sd"),drop=FALSE][which.min(res_errCV[(nlvsum==i),"mean"]),,drop=FALSE])))))
      )
      res_nlvsum_errCV_Ysel <- res_nlvsum_errCV_Ysel[order(res_nlvsum_errCV_Ysel[,"totalnlv"], decreasing=FALSE),]
      
      if(selection=="localmin"){
        if(nrow(res_nlvsum_errCV_Ysel)>1){
          # sign of the difference of accuracies to select the optim combination with the lower total number of components
          rtsdiff <- c(NA,sapply(2:nrow(res_nlvsum_errCV_Ysel), function(i)((res_nlvsum_errCV_Ysel$mean[i]-res_nlvsum_errCV_Ysel$mean[i-1])<=0)))
          kchoix <- min(res_nlvsum_errCV_Ysel[(which(rtsdiff==FALSE)-1)[1],"index"],res_nlvsum_errCV_Ysel[nrow(res_nlvsum_errCV_Ysel),"index"],na.rm=TRUE)
        }else{
          kchoix <- 1
        }
      }
      if(selection=="globalmin"){
        kchoix <- res_nlvsum_errCV_Ysel[which.min(res_nlvsum_errCV_Ysel$mean)[1],"index"]
      }
      if(selection=="1std"){
        if(nrow(res_nlvsum_errCV_Ysel)>1){
          # one standard error rule to select the optim number of components
          minmean    <- which.min(res_nlvsum_errCV_Ysel$mean)[1]
          threshmean <- res_nlvsum_errCV_Ysel$mean[minmean] + res_nlvsum_errCV_Ysel$sd[minmean]
          if((minmean == 1) | (sum(res_nlvsum_errCV_Ysel$mean[1:minmean]<=threshmean)==0)){
            kchoix <- 1
          }else{
            kchoix <- min(res_nlvsum_errCV_Ysel[which(res_nlvsum_errCV_Ysel$mean[1:minmean]<=threshmean),"index"])
          }
        }else{
          kchoix <- 1
        }
      }
    }

    # outputs
    rts <- list(lvcombi=lvcombi, 
                optimcombi=unlist(lvcombi[kchoix,,drop=FALSE]), 
                optimExplVarCV=res_ExplVarCV[kchoix,,drop=FALSE], 
                res_rmseCV=res_rmseCV, 
                res_ExplVarCV=res_ExplVarCV, 
                res_errCV=res_errCV)
    
  }
  
  ############################################################################################################################
  
  ## SEQUENTIAL optimisation
  if((optimisation=="sequential") & (nXblocks>1)){
    
    ## list initialisations
    lvcombi <- list()

    optimYpredCV <- list()
    res_errCV <- res_rmseCV <- res_ExplVarCV <- list()
    Rep_errCV <- Rep_rmseCV <- Rep_ExplVarCV <- list()

    for (m in 1:nXblocks) {
      # combinations
      if(m==1){
        lvcombi[[1]] <- as.matrix(expand.grid(nlvlist[[1]]))
        colnames(lvcombi[[1]]) <- "Xlist1"
        rownames(lvcombi[[1]]) <- paste0("lvcombi",1:nrow(lvcombi[[1]]))
      }
      if(m>1){
        lvcombi[[m]] <- data.frame(matrix(rep(unlist(optimcombi),length(nlvlist[[m]])),nrow=length(nlvlist[[m]]),byrow=TRUE),as.matrix(expand.grid(nlvlist[[m]])))
        colnames(lvcombi[[m]]) <- paste0(rep("Xlist",m),1:m)
        rownames(lvcombi[[m]]) <- paste0("lvcombi",1:nrow(lvcombi[[m]]))
      }
      
      # output initialisation

      res_errCV[[m]] <- res_ExplVarCV[[m]] <- res_rmseCV[[m]] <- matrix(NA, nrow=nrow(lvcombi[[m]]), ncol=2, dimnames=list(rownames(lvcombi[[m]]),c("mean","sd")))
      Rep_errCV[[m]] <- Rep_rmseCV[[m]] <- Rep_ExplVarCV[[m]] <- array(0, dim=c(nrow(lvcombi[[m]]),nbrep, 1), dimnames = list(paste0("lvcombi",1:nrow(lvcombi[[m]])), paste0("rep",1:nbrep), "value"))

      for(i in 1:nbrep){
        for(j in 1:nrow(lvcombi[[m]])){
          reppredCV <- matrix(NA, nrow=n, ncol = length(unique(Y)), dimnames=list(rownames(Xlist),sort(unique(Y))))
          reppredYCV <- matrix(NA, nrow=n, ncol = length(Yvariables), dimnames=list(rownames(Xlist),Yvariables))
          
          for(k in 1:length(CVtype[[i]])){
            Xlisttrain <-lapply(1:m,function(x) Xlist[[x]][-CVtype[[i]][[k]],,drop=FALSE])
            Ytrain <- Y[-CVtype[[i]][[k]],,drop=FALSE]
            Xlisttest <- lapply(1:m,function(x) Xlist[[x]][CVtype[[i]][[k]],,drop=FALSE])
            Ytest <- Y[CVtype[[i]][[k]],,drop=FALSE]
            repmodel <- funda(Xlist=Xlisttrain, y=Ytrain, Xscaling = Xscaling, Yscaling = Yscaling, weights = weights[-CVtype[[i]][[k]]], nlv=unlist(lvcombi[[m]][j,]), prior = prior)
            reppredCV[CVtype[[k]],] <- predict(repmodel, Xlisttest)$posterior[order(CVtype[[i]][[k]]),]
            reppredYCV[CVtype[[k]],] <- predict(repmodel, Xlisttest)$pred[order(CVtype[[i]][[k]]),]
          }
          SqErrCV <- (dummy(Y)$Y-reppredCV)^2

          Rep_rmseCV[[m]][j,i,]    <- sqrt(mean(SqErrCV))
          Rep_ExplVarCV[[m]][j,i,] <- mean(1-(apply(SqErrCV,2,mean)/(apply(dummy(Y)$Y,2,var)*(n-1)/n)))

          Rep_errCV[[m]][j,i,]    <- mean(reppredYCV != Y)*100
        }
      }
      
      res_rmseCV[[m]][,"mean"] <- apply(Rep_rmseCV[[m]], c(1,3),mean, na.rm=T)
      res_rmseCV[[m]][,"sd"]   <- apply(Rep_rmseCV[[m]], c(1,3),sd, na.rm=T)
      
      res_ExplVarCV[[m]][,"mean"] <- apply(Rep_ExplVarCV[[m]], c(1,3),mean, na.rm=T)
      res_ExplVarCV[[m]][,"sd"]   <- apply(Rep_ExplVarCV[[m]], c(1,3),sd, na.rm=T)
      
      res_errCV[[m]][,"mean"] <- apply(Rep_errCV[[m]], c(1,3),mean, na.rm=T)
      res_errCV[[m]][,"sd"]   <- apply(Rep_errCV[[m]], c(1,3),sd, na.rm=T)
      
      if(criterion == "rmse"){
        # optim combination for each total number of components // res_rmse_Ysel[,"mean"]
        
        nlvsum=apply(lvcombi[[m]],1,sum)
        if(m==1){
          res_nlvsum_rmseCV_Ysel <- data.frame(
            index = sapply((unique(nlvsum)), function(i) which(rownames(res_rmseCV[[m]])==rownames(res_rmseCV[[m]][(nlvsum==i),,drop=FALSE][which.max(res_rmseCV[[m]][(nlvsum==i),"mean"]),,drop=FALSE]))),
            combnum = sapply((unique(nlvsum)), function(i) rownames(res_rmseCV[[m]][(nlvsum==i),,drop=FALSE][which.max(res_rmseCV[[m]][(nlvsum==i),"mean"]),,drop=FALSE])),
            lvcombi[[m]][,,drop=FALSE],
            totalnlv = unique(nlvsum),
            t((sapply((unique(nlvsum)), function(i) unlist(data.frame(res_rmseCV[[m]][(nlvsum==i),c("mean","sd"),drop=FALSE][which.max(res_rmseCV[[m]][(nlvsum==i),"mean"]),,drop=FALSE])))))
          )
        }
        if(m>1){
          res_nlvsum_rmseCV_Ysel <- data.frame(
            index = sapply((unique(nlvsum)), function(i) which(rownames(res_rmseCV[[m]])==rownames(res_rmseCV[[m]][(nlvsum==i),,drop=FALSE][which.max(res_rmseCV[[m]][(nlvsum==i),"mean"]),,drop=FALSE]))),
            combnum = sapply((unique(nlvsum)), function(i) rownames(res_rmseCV[[m]][(nlvsum==i),,drop=FALSE][which.max(res_rmseCV[[m]][(nlvsum==i),"mean"]),,drop=FALSE])),
            t((sapply((unique(nlvsum)), function(i) unlist(data.frame(lvcombi[[m]][(nlvsum==i),,drop=FALSE][which.max(res_rmseCV[[m]][(nlvsum==i),"mean"]),,drop=FALSE]))))),
            totalnlv = unique(nlvsum),
            t((sapply((unique(nlvsum)), function(i) unlist(data.frame(res_rmseCV[[m]][(nlvsum==i),c("mean","sd"),drop=FALSE][which.max(res_rmseCV[[m]][(nlvsum==i),"mean"]),,drop=FALSE])))))
          )
        }
        
        if(selection=="localmin"){
          if(nrow(res_nlvsum_rmseCV_Ysel)>1){
            # sign of the difference of accuracies to select the optim combination with the lower total number of components
            rtsdiff <- c(NA,sapply(2:nrow(res_nlvsum_rmseCV_Ysel), function(i)((res_nlvsum_rmseCV_Ysel$mean[i]-res_nlvsum_rmseCV_Ysel$mean[i-1])<=0)))
            kchoix <- min(res_nlvsum_rmseCV_Ysel[(which(rtsdiff==FALSE)-1)[1],"index"],res_nlvsum_rmseCV_Ysel[nrow(res_nlvsum_rmseCV_Ysel),"index"],na.rm=TRUE)
          }else{
            kchoix <- 1
          }
        }
        if(selection=="globalmin"){
          kchoix <- res_nlvsum_rmseCV_Ysel[which.min(res_nlvsum_rmseCV_Ysel$mean)[1],"index"]
        }
        if(selection=="1std"){
          if(nrow(res_nlvsum_rmseCV_Ysel)>1){
            # one standard error rule to select the optim number of components
            minmean    <- which.min(res_nlvsum_rmseCV_Ysel$mean)[1]
            threshmean <- res_nlvsum_rmseCV_Ysel$mean[minmean] + res_nlvsum_rmseCV_Ysel$sd[minmean]
            if((minmean == 1) | (sum(res_nlvsum_rmseCV_Ysel$mean[1:minmean]<=threshmean)==0)){
              kchoix <- 1
            }else{
              kchoix <- min(res_nlvsum_rmseCV_Ysel[which(res_nlvsum_rmseCV_Ysel$mean[1:minmean]<=threshmean),"index"])
            }
          }else{
            kchoix <- 1
          }
        }
      }
      
      if(criterion == "err"){
        # optim combination for each total number of components // res_err_Ysel[,"mean"]
        
        nlvsum=apply(lvcombi[[m]],1,sum)
        if(m==1){
          res_nlvsum_errCV_Ysel <- data.frame(
            index = sapply((unique(nlvsum)), function(i) which(rownames(res_errCV[[m]])==rownames(res_errCV[[m]][(nlvsum==i),,drop=FALSE][which.max(res_errCV[[m]][(nlvsum==i),"mean"]),,drop=FALSE]))),
            combnum = sapply((unique(nlvsum)), function(i) rownames(res_errCV[[m]][(nlvsum==i),,drop=FALSE][which.max(res_errCV[[m]][(nlvsum==i),"mean"]),,drop=FALSE])),
            lvcombi[[m]][,,drop=FALSE],
            totalnlv = unique(nlvsum),
            t((sapply((unique(nlvsum)), function(i) unlist(data.frame(res_errCV[[m]][(nlvsum==i),c("mean","sd"),drop=FALSE][which.max(res_errCV[[m]][(nlvsum==i),"mean"]),,drop=FALSE])))))
          )
        }
        if(m>1){
          res_nlvsum_errCV_Ysel <- data.frame(
            index = sapply((unique(nlvsum)), function(i) which(rownames(res_errCV[[m]])==rownames(res_errCV[[m]][(nlvsum==i),,drop=FALSE][which.max(res_errCV[[m]][(nlvsum==i),"mean"]),,drop=FALSE]))),
            combnum = sapply((unique(nlvsum)), function(i) rownames(res_errCV[[m]][(nlvsum==i),,drop=FALSE][which.max(res_errCV[[m]][(nlvsum==i),"mean"]),,drop=FALSE])),
            t((sapply((unique(nlvsum)), function(i) unlist(data.frame(lvcombi[[m]][(nlvsum==i),,drop=FALSE][which.max(res_errCV[[m]][(nlvsum==i),"mean"]),,drop=FALSE]))))),
            totalnlv = unique(nlvsum),
            t((sapply((unique(nlvsum)), function(i) unlist(data.frame(res_errCV[[m]][(nlvsum==i),c("mean","sd"),drop=FALSE][which.max(res_errCV[[m]][(nlvsum==i),"mean"]),,drop=FALSE])))))
          )
        }
        
        if(selection=="localmin"){
          if(nrow(res_nlvsum_errCV_Ysel)>1){
            # sign of the difference of accuracies to select the optim combination with the lower total number of components
            rtsdiff <- c(NA,sapply(2:nrow(res_nlvsum_errCV_Ysel), function(i)((res_nlvsum_errCV_Ysel$mean[i]-res_nlvsum_errCV_Ysel$mean[i-1])<=0)))
            kchoix <- min(res_nlvsum_errCV_Ysel[(which(rtsdiff==FALSE)-1)[1],"index"],res_nlvsum_errCV_Ysel[nrow(res_nlvsum_errCV_Ysel),"index"],na.rm=TRUE)
          }else{
            kchoix <- 1
          }
        }
        if(selection=="globalmin"){
          kchoix <- res_nlvsum_errCV_Ysel[which.min(res_nlvsum_errCV_Ysel$mean)[1],"index"]
        }
        if(selection=="1std"){
          if(nrow(res_nlvsum_errCV_Ysel)>1){
            # one standard error rule to select the optim number of components
            minmean    <- which.min(res_nlvsum_errCV_Ysel$mean)[1]
            threshmean <- res_nlvsum_errCV_Ysel$mean[minmean] + res_nlvsum_errCV_Ysel$sd[minmean]
            if((minmean == 1)| (sum(res_nlvsum_errCV_Ysel$mean[1:minmean]<=threshmean)==0)){
              kchoix <- 1
            }else{
              kchoix <- min(res_nlvsum_errCV_Ysel[which(res_nlvsum_errCV_Ysel$mean[1:minmean]<=threshmean),"index"])
            }
          }else{
            kchoix <- 1
          }
        }
      }
      
      # optim combination
      optimcombi <- as.vector(unlist(lvcombi[[m]][kchoix,]))
      names(optimcombi) <- paste0(rep("Xlist",m),1:m)
      nlvlist[[m]] <- lvcombi[[m]][kchoix,m]
    }# end loop on m
    

    ## output names
    names(lvcombi) <- names(res_errCV) <- names(res_rmseCV) <- names(res_ExplVarCV) <- paste0("nbBlocks",1:nXblocks)

    # outputs
    rts <- list(lvcombi=lvcombi,
                optimcombi=unlist(lvcombi[[nXblocks]][kchoix,,drop=FALSE]), 
                optimExplVarCV=res_ExplVarCV[[nXblocks]][kchoix,,drop=FALSE], 
                res_rmseCV=res_rmseCV, 
                res_ExplVarCV=res_ExplVarCV, 
                res_errCV=res_errCV)

  }
  rts$call   <- match.call()
  class(rts) <- c("Soplscv")
  rts
}


soplsldacv <- function(Xlist, y, Xscaling = c("none", "pareto", "sd")[1], Yscaling = c("none", "pareto", "sd")[1], weights = NULL, nlvlist = list(), prior = c("unif", "prop"), nbrep = 30, cvmethod = "kfolds", seed = 123, samplingk = NULL, nfolds = 7, optimisation = c("global","sequential")[1], criterion = c("err", "rmse")[1], selection = c("localmin","globalmin","1std")[1]){
  .soplsprobdacv(funda = soplslda, Xlist = Xlist, y = y, Xscaling = Xscaling, Yscaling = Yscaling, weights = weights, nlvlist = nlvlist, prior = prior, nbrep = nbrep, cvmethod = cvmethod, seed = seed, samplingk = samplingk, nfolds = nfolds, optimisation = optimisation, criterion = criterion, selection = selection)
}


soplsqdacv <- function(Xlist, y, Xscaling = c("none", "pareto", "sd")[1], Yscaling = c("none", "pareto", "sd")[1], weights = NULL, nlvlist = list(), prior = c("unif", "prop"), nbrep = 30, cvmethod = "kfolds", seed = 123, samplingk = NULL, nfolds = 7, optimisation = c("global","sequential")[1], criterion = c("err", "rmse")[1], selection = c("localmin","globalmin","1std")[1]){
  .soplsprobdacv(funda = soplsqda, Xlist = Xlist, y = y, Xscaling = Xscaling, Yscaling = Yscaling, weights = weights, nlvlist = nlvlist, prior = prior, nbrep = nbrep, cvmethod = cvmethod, seed = seed, samplingk = samplingk, nfolds = nfolds, optimisation = optimisation, criterion = criterion, selection = selection)
}
  
