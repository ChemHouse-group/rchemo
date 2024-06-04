
.soplsrcv <- function(Xlist, Y, Xscaling = c("none", "pareto", "sd")[1], Yscaling = c("none", "pareto", "sd")[1], weights = NULL, nlvlist=list(), nbrep=30, cvmethod="kfolds", seed = 123, samplingk=NULL, nfolds=7, optimisation="global", selection="1std", majorityvote=FALSE){
  
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
  if((cvmethod=="kfolds") & (nfolds<=1)){stop("nfolds must not be >1 when cvmethod is 'kfolds'")}
  
  # argument computations
  n <- nrow(Y)
  nXblocks <- length(Xlist)
  
  Yvariables    <- colnames(Y)
  #nYvariables   <- length(Yvariables)
  
  inertieY <- sum(diag(t(as.matrix(scale(Y, scale=FALSE)))%*%(as.matrix(scale(Y, scale=FALSE)))))

  # datasets
  set.seed(seed = seed)
  
  obsList     <- list()
  initialYrep <- list()
  Yrep        <- list()
  Xlistrep    <- list()
  
  initialY     <- Y
  
  if(nbrep==1){
    initialYrep[[1]] <- initialY
    Yrep[[1]]        <- Y
    Xlistrep[[1]]    <- Xlist
  }else{
    for(i in 1:nbrep){
      obsList[[i]]      <- sample(x=n, size=n, replace = FALSE, prob = NULL)
      initialYrep[[i]]  <- initialY[obsList[[i]],, drop=FALSE]
      Yrep[[i]]         <- Y[obsList[[i]],, drop=FALSE]
      Xlistrep[[i]]     <- lapply(1:length(Xlist), function (j) Xlist[[j]][obsList[[i]],,drop=FALSE])
    }
  }

  
  ############################################################################################################################
  ## GLOBAL optimisation
  if((optimisation=="global")|(nXblocks==1)){
    
    lvcombi <- as.matrix(expand.grid(nlvlist))
    colnames(lvcombi) <- paste0(rep("Xlist",nXblocks),1:nXblocks)
    rownames(lvcombi) <- paste0("lvcombi",1:nrow(lvcombi))
    
    # output tables
    res_ExplVarCV_byY <- res_rmseCV_byY <- matrix(NA, nrow=nrow(lvcombi), ncol=(ncol(Y)*2), 
                         dimnames=list(rownames(lvcombi),c(paste0("mean_",Yvariables),paste0("sd_",Yvariables))))
    
    res_ExplVarCV <- res_rmseCV <- matrix(NA, nrow=nrow(lvcombi), ncol=2, dimnames=list(rownames(lvcombi),c("mean","sd")))
    
    Rep_ExplVarCV_byY <- Rep_rmseCV_byY <- array(0, dim=c(nrow(lvcombi),nbrep, ncol(Y)), dimnames = list(rownames(lvcombi), paste0("rep",1:nbrep), Yvariables))
    
    Rep_rmseCV <- Rep_ExplVarCV <- array(0, dim=c(nrow(lvcombi),nbrep, 1), dimnames = list(paste0("lvcombi",1:nrow(lvcombi)), paste0("rep",1:nbrep), "value"))
    
  	#optimYpredCV <- list()
  	

    for(i in 1:nbrep){
      
      # partition
      set.seed(seed = seed)
      if(cvmethod=="loo"){CVtype=as.list(1:n)}
      if((cvmethod=="kfolds") & (is.null(samplingk)==TRUE)){
        CVtype <- segmkf(n = n, K = nfolds, type = "interleaved")$rep1
      }
      if((cvmethod=="kfolds") & (is.null(samplingk)==FALSE)){
        samplingktab <- table(samplingk)
        for(s in 1:length(names(samplingktab))){
          partCVtype[[s]] <- segmkf(n = samplingktab[s], K = nfolds, type = "interleaved")$rep1
          for(ss in 1:length(partCVtype[[s]])){
            partCVtype[[s]][[ss]] <- which(samplingk==names(samplingktab)[s])[partCVtype[[s]][[ss]]]
          }
        }
        CVtype <- list()
        for(sss in 1:nfolds){
          CVtype[[sss]] <- partCVtype[[1]][[sss]]
          for(s in 2:length(names(samplingktab))){
            CVtype[[sss]] <- c(CVtype[[sss]],partCVtype[[s]][[sss]])
          }
        }
      }
      
      
      for(j in 1:nrow(lvcombi)){
        reppredCV <- matrix(NA, nrow=n, ncol = length(Yvariables))
        for(k in 1:length(CVtype)){
          Xlisttrain <-lapply(1:length(Xlistrep[[i]]),function(x) Xlistrep[[i]][[x]][-CVtype[[k]],,drop=FALSE])
          Ytrain <- Yrep[[i]][-CVtype[[k]],,drop=FALSE]
          Xlisttest <- lapply(1:length(Xlistrep[[i]]),function(x) Xlistrep[[i]][[x]][CVtype[[k]],,drop=FALSE])
          Ytest <- Yrep[[i]][CVtype[[k]],,drop=FALSE]
          repmodel <- soplsr(Xlist=Xlisttrain, Y=Ytrain, Xscaling = Xscaling, Yscaling = Yscaling, weights = weights, nlv=lvcombi[j,])
          reppredCV[CVtype[[k]],] <- predict(repmodel, Xlisttest)
        }
        SqErrCV <- (Yrep[[i]]-reppredCV)^2
        
        Rep_rmseCV_byY[j,i,]    <- sqrt(apply(SqErrCV,2,mean))
        Rep_ExplVarCV_byY[j,i,] <- 1- ((Rep_rmseCV_byY[j,i,]^2)/matrix((apply(Yrep[[i]],2,var)*(n-1)/n), ncol=ncol(Yrep[[i]])))
        Rep_rmseCV[j,i,]    <- sqrt(apply(Rep_rmseCV_byY[j,i,,drop=FALSE]^2, FUN = mean, MARGIN = 1))
        Rep_ExplVarCV[j,i,] <- apply(Rep_ExplVarCV_byY[j,i,,drop=FALSE], FUN = mean, MARGIN = 1)
      }
    }

  	res_rmseCV_byY[,1:ncol(Y)]               <- apply(Rep_rmseCV_byY, c(1,3),mean, na.rm=T)
    res_rmseCV_byY[,(ncol(Y)+1):(ncol(Y)*2)] <- apply(Rep_rmseCV_byY, c(1,3),sd, na.rm=T)

    res_ExplVarCV_byY[,1:ncol(Y)]               <- apply(Rep_ExplVarCV_byY, c(1,3),mean, na.rm=T)
    res_ExplVarCV_byY[,(ncol(Y)+1):(ncol(Y)*2)] <- apply(Rep_ExplVarCV_byY, c(1,3),sd, na.rm=T)

    res_rmseCV[,"mean"] <- apply(Rep_rmseCV, c(1,3),mean, na.rm=T)
    res_rmseCV[,"sd"]   <- apply(Rep_rmseCV, c(1,3),sd, na.rm=T)
    
    res_ExplVarCV[,"mean"] <- apply(Rep_ExplVarCV, c(1,3),mean, na.rm=T)
    res_ExplVarCV[,"sd"]   <- apply(Rep_ExplVarCV, c(1,3),sd, na.rm=T)
    
    # optim combination for each total number of components // res_rmse_Ysel[,"mean"]
    if(majorityvote==TRUE){
      choiceYH <- rep(NA,length(Yvariables))
      nlvsum <- apply(lvcombi,1,sum)
      
      for(yy in 1:length(Yvariables)){
        res_rmseCV_Ysel <- cbind.data.frame(res_rmseCV_byY[,yy,drop=FALSE], res_rmseCV_byY[,(ncol(Y)+yy),drop=FALSE]) 
        colnames(res_rmseCV_Ysel) <- c("mean","sd")
        rownames(res_rmseCV_Ysel) <- rownames(res_rmseCV_byY)
        
        res_nlvsum_rmseCV_Ysel <- data.frame(
          index = sapply((unique(nlvsum)), function(i) which(rownames(res_rmseCV_Ysel)==rownames(res_rmseCV_Ysel[(nlvsum==i),,drop=FALSE][which.min(res_rmseCV_Ysel[(nlvsum==i),"mean"]),,drop=FALSE]))),
          combnum = sapply((unique(nlvsum)), function(i) rownames(res_rmseCV_Ysel[(nlvsum==i),,drop=FALSE][which.min(res_rmseCV_Ysel[(nlvsum==i),"mean"]),,drop=FALSE])),
          t((sapply((unique(nlvsum)), function(i) unlist(data.frame(lvcombi[(nlvsum==i),,drop=FALSE][which.min(res_rmseCV_Ysel[(nlvsum==i),"mean"]),,drop=FALSE]))))),
          totalnlv = unique(nlvsum),
          t((sapply((unique(nlvsum)), function(i) unlist(data.frame(res_rmseCV_Ysel[(nlvsum==i),c("mean","sd"),drop=FALSE][which.min(res_rmseCV_Ysel[(nlvsum==i),"mean"]),,drop=FALSE])))))
        )
        
        res_nlvsum_rmseCV_Ysel <- res_nlvsum_rmseCV_Ysel[order(res_nlvsum_rmseCV_Ysel[,"totalnlv"], decreasing=FALSE),]
        
        if(selection=="localmin"){
          if(nrow(res_nlvsum_rmseCV_Ysel)>1){
            # sign of the difference of accuracies to select the optim combination with the lower total number of components
            rtsdiff <- c(NA,sapply(2:nrow(res_nlvsum_rmseCV_Ysel), function(i)((res_nlvsum_rmseCV_Ysel$mean[i]-res_nlvsum_rmseCV_Ysel$mean[i-1])<=0)))
            choiceYH[yy] <- min(res_nlvsum_rmseCV_Ysel[(which(rtsdiff==FALSE)-1)[1],"index"],res_nlvsum_rmseCV_Ysel[nrow(res_nlvsum_rmseCV_Ysel),"index"],na.rm=TRUE)
          }else{
            choiceYH[yy]<- 1
          }
        }
        if(selection=="globalmin"){
          choiceYH[yy] <- res_nlvsum_rmseCV_Ysel[which.min(res_nlvsum_rmseCV_Ysel$mean)[1],"index"]
        }
        if(selection=="1std"){
          if(nrow(res_nlvsum_rmseCV_Ysel)>1){
            # one standard error rule to select the optim number of components
            minmean    <- which.min(res_nlvsum_rmseCV_Ysel$mean)[1]
            threshmean <- res_nlvsum_rmseCV_Ysel$mean[minmean] + res_nlvsum_rmseCV_Ysel$sd[minmean]
            if((minmean == 1) | (sum(res_nlvsum_rmseCV_Ysel$mean[1:minmean]>=threshmean)==0)){
              choiceYH[yy] <- 1
            }else{
              choiceYH[yy] <- max(res_nlvsum_rmseCV_Ysel[which(res_nlvsum_rmseCV_Ysel$mean[1:minmean]>=threshmean),"index"])
            }
          }else{
            choiceYH[yy] <- 1
          }
        }
      }
      tabfreq <- sapply(1:nrow(lvcombi), function(l) sum(choiceYH==l))
      kchoix <- which.max(tabfreq)[1] 
    }

    if(majorityvote==FALSE){
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
          if((minmean == 1) | (sum(res_nlvsum_rmseCV_Ysel$mean[1:minmean]>=threshmean)==0)){
            kchoix <- 1
          }else{
            kchoix <- max(res_nlvsum_rmseCV_Ysel[which(res_nlvsum_rmseCV_Ysel$mean[1:minmean]>=threshmean),"index"])
          }
        }else{
          kchoix <- 1
        }
      }
    }
    

    # outputs
    rts <- list(lvcombi=lvcombi, 
                #optimCombiLine=kchoix, 
                optimcombi=unlist(lvcombi[kchoix,,drop=FALSE]), 
                #optimExplVarCV=res_ExplVarCV[kchoix,,drop=FALSE], 
                rmseCV_byY=res_rmseCV_byY, 
                ExplVarCV_byY=res_ExplVarCV_byY, 
                rmseCV=res_rmseCV, 
                ExplVarCV=res_ExplVarCV)
  }
  
  ############################################################################################################################
  
  ## SEQUENTIAL optimisation
  if((optimisation=="sequential") & (nXblocks>1)){

    ## list initialisations
    lvcombi        <- list()

    #optimYpredCV   <- list()
    res_rmseCV_byY <- res_ExplVarCV_byY <- list()
    res_rmseCV     <- res_ExplVarCV     <- list()
    Rep_rmseCV_byY <- Rep_ExplVarCV_byY <- list()
    Rep_rmseCV     <- Rep_ExplVarCV     <- list()
    
    for(i in 1:nbrep){
      # partition
      set.seed(seed = seed)
      if(cvmethod=="loo"){CVtype=as.list(1:n)}
      if((cvmethod=="kfolds") & (is.null(samplingk)==TRUE)){
        CVtype <- segmkf(n = n, K = nfolds, type = "interleaved")$rep1
      }
      if((cvmethod=="kfolds") & (is.null(samplingk)==FALSE)){
        samplingktab <- table(samplingk)
        for(s in 1:length(names(samplingktab))){
          partCVtype[[s]] <- segmkf(n = samplingktab[s], K = nfolds, type = "interleaved")$rep1
          for(ss in 1:length(partCVtype[[s]])){
            partCVtype[[s]][[ss]] <- which(samplingk==names(samplingktab)[s])[partCVtype[[s]][[ss]]]
          }
        }
        CVtype <- list()
        for(sss in 1:nfolds){
          CVtype[[sss]] <- partCVtype[[1]][[sss]]
          for(s in 2:length(names(samplingktab))){
            CVtype[[sss]] <- c(CVtype[[sss]],partCVtype[[s]][[sss]])
          }
        }
      }
    }
    
    for (m in 1:nXblocks) {
      # combinations
      if(m==1){
        lvcombi[[1]] <- as.matrix(expand.grid(nlvlist[[1]]))
        colnames(lvcombi[[1]]) <- "Xlist1"
        rownames(lvcombi[[1]]) <- paste0("lvcombi",1:nrow(lvcombi[[1]]))
      }
      if(m>1){
        lvcombi[[m]] <- data.frame(matrix(rep(unlist(optimcombi),length(nlvlist[[m]])),nrow=length(nlvlist[[m]]),byrow=TRUE),as.matrix(expand.grid(nlvlist[[m]])))
        colnames(lvcombi[[m]]) <- paste0("Xlist",1:m)
        rownames(lvcombi[[m]]) <- paste0("lvcombi",1:nrow(lvcombi[[m]]))
      }
      
      # output initialisation
      res_ExplVarCV_byY[[m]] <- res_rmseCV_byY[[m]] <- matrix(NA, nrow=nrow(lvcombi[[m]]), ncol=(ncol(Y)*2), 
                                    dimnames=list(rownames(lvcombi[[m]]),c(paste0("mean_",Yvariables),paste0("sd_",Yvariables))))

      res_ExplVarCV[[m]] <- res_rmseCV[[m]] <- matrix(NA, nrow=nrow(lvcombi[[m]]), ncol=2, dimnames=list(rownames(lvcombi[[m]]),c("mean","sd")))

      Rep_ExplVarCV_byY[[m]] <- Rep_rmseCV_byY[[m]] <- array(0, dim=c(nrow(lvcombi[[m]]),nbrep, ncol(Y)), dimnames = list(paste0("lvcombi",1:nrow(lvcombi[[m]])), paste0("rep",1:nbrep), Yvariables))
      
      Rep_rmseCV[[m]] <- Rep_ExplVarCV[[m]] <- array(0, dim=c(nrow(lvcombi[[m]]),nbrep, 1), dimnames = list(paste0("lvcombi",1:nrow(lvcombi[[m]])), paste0("rep",1:nbrep), "value"))
      
 
      for(i in 1:nbrep){
        for(j in 1:nrow(lvcombi[[m]])){
          reppredCV <- matrix(NA, nrow=n, ncol = length(Yvariables))
          for(k in 1:length(CVtype)){
            Xlisttrain <-lapply(1:m,function(x) Xlistrep[[i]][[x]][-CVtype[[k]],,drop=FALSE])
            Ytrain <- Yrep[[i]][-CVtype[[k]],,drop=FALSE]
            Xlisttest <- lapply(1:m,function(x) Xlistrep[[i]][[x]][CVtype[[k]],,drop=FALSE])
            Ytest <- Yrep[[i]][CVtype[[k]],,drop=FALSE]
            repmodel <- soplsr(Xlist=Xlisttrain, Y=Ytrain, Xscaling = Xscaling, Yscaling = Yscaling, weights = weights, nlv=unlist(lvcombi[[m]][j,]))
            reppredCV[CVtype[[k]],] <- predict(repmodel, Xlisttest)
          }
          SqErrCV <- (Yrep[[i]]-reppredCV)^2

          Rep_rmseCV_byY[[m]][j,i,]    <- sqrt(apply(SqErrCV,2,mean))
          Rep_ExplVarCV_byY[[m]][j,i,] <- 1- ((Rep_rmseCV_byY[[m]][j,i,]^2)/matrix((apply(Yrep[[i]],2,var)*(n-1)/n), ncol=ncol(Yrep[[i]])))

          Rep_rmseCV[[m]][j,i,]    <- sqrt(apply(Rep_rmseCV_byY[[m]][j,i,,drop=FALSE]^2, FUN = mean, MARGIN = 1))
          Rep_ExplVarCV[[m]][j,i,] <- apply(Rep_ExplVarCV_byY[[m]][j,i,,drop=FALSE], FUN = mean, MARGIN = 1)
        }
      }
      
      res_rmseCV_byY[[m]][,1:ncol(Y)]                <- apply(Rep_rmseCV_byY[[m]], c(1,3),mean, na.rm=T)
      res_rmseCV_byY[[m]][,(ncol(Y)+1):(ncol(Y)*2)] <- apply(Rep_rmseCV_byY[[m]], c(1,3),sd, na.rm=T)

      res_ExplVarCV_byY[[m]][,1:ncol(Y)]                <- apply(Rep_ExplVarCV_byY[[m]], c(1,3),mean, na.rm=T)
      res_ExplVarCV_byY[[m]][,(ncol(Y)+1):(ncol(Y)*2)] <- apply(Rep_ExplVarCV_byY[[m]], c(1,3),sd, na.rm=T)
      
      res_rmseCV[[m]][,"mean"] <- apply(Rep_rmseCV[[m]], c(1,3),mean, na.rm=T)
      res_rmseCV[[m]][,"sd"]   <- apply(Rep_rmseCV[[m]], c(1,3),sd, na.rm=T)

      res_ExplVarCV[[m]][,"mean"] <- apply(Rep_ExplVarCV[[m]], c(1,3),mean, na.rm=T)
      res_ExplVarCV[[m]][,"sd"]   <- apply(Rep_ExplVarCV[[m]], c(1,3),sd, na.rm=T)

      # optim combination for each total number of components // res_rmse_Ysel[,"mean"]
      if(majorityvote==TRUE){
        choiceYH <- rep(NA,length(Yvariables))
        nlvsum <- apply(lvcombi[[m]],1,sum)
        
        for(yy in 1:length(Yvariables)){
          res_rmseCV_Ysel <- cbind.data.frame(res_rmseCV_byY[[m]][,yy,drop=FALSE], res_rmseCV_byY[[m]][,(ncol(Y)+yy),drop=FALSE]) 
          colnames(res_rmseCV_Ysel) <- c("mean","sd")
          rownames(res_rmseCV_Ysel) <- rownames(res_rmseCV_byY[[m]])
          
          if(m==1){
            res_nlvsum_rmseCV_Ysel <- data.frame(
              index = sapply((unique(nlvsum)), function(i) which(rownames(res_rmseCV_Ysel)==rownames(res_rmseCV_Ysel[(nlvsum==i),,drop=FALSE][which.min(res_rmseCV_Ysel[(nlvsum==i),"mean"]),,drop=FALSE]))),
              combnum = sapply((unique(nlvsum)), function(i) rownames(res_rmseCV_Ysel[(nlvsum==i),,drop=FALSE][which.min(res_rmseCV_Ysel[(nlvsum==i),"mean"]),,drop=FALSE])),
              lvcombi[[m]][,,drop=FALSE],
              totalnlv = unique(nlvsum),
              t((sapply((unique(nlvsum)), function(i) unlist(data.frame(res_rmseCV_Ysel[(nlvsum==i),c("mean","sd"),drop=FALSE][which.min(res_rmseCV_Ysel[(nlvsum==i),"mean"]),,drop=FALSE])))))
            )
          }
          if(m>=1){
            res_nlvsum_rmseCV_Ysel <- data.frame(
              index = sapply((unique(nlvsum)), function(i) which(rownames(res_rmseCV_Ysel)==rownames(res_rmseCV_Ysel[(nlvsum==i),,drop=FALSE][which.min(res_rmseCV_Ysel[(nlvsum==i),"mean"]),,drop=FALSE]))),
              combnum = sapply((unique(nlvsum)), function(i) rownames(res_rmseCV_Ysel[(nlvsum==i),,drop=FALSE][which.min(res_rmseCV_Ysel[(nlvsum==i),"mean"]),,drop=FALSE])),
              t((sapply((unique(nlvsum)), function(i) unlist(data.frame(lvcombi[[m]][(nlvsum==i),,drop=FALSE][which.min(res_rmseCV_Ysel[(nlvsum==i),"mean"]),,drop=FALSE]))))),
              totalnlv = unique(nlvsum),
              t((sapply((unique(nlvsum)), function(i) unlist(data.frame(res_rmseCV_Ysel[(nlvsum==i),c("mean","sd"),drop=FALSE][which.min(res_rmseCV_Ysel[(nlvsum==i),"mean"]),,drop=FALSE])))))
            )
          }
          res_nlvsum_rmseCV_Ysel <- res_nlvsum_rmseCV_Ysel[order(res_nlvsum_rmseCV_Ysel[,"totalnlv"], decreasing=FALSE),]
          
          if(selection=="localmin"){
            if(nrow(res_nlvsum_rmseCV_Ysel)>1){
              rtsdiff <- c(NA,sapply(2:nrow(res_nlvsum_rmseCV_Ysel), function(i)((res_nlvsum_rmseCV_Ysel$mean[i]-res_nlvsum_rmseCV_Ysel$mean[i-1])<=0)))
              choiceYH[yy] <- min(res_nlvsum_rmseCV_Ysel[(which(rtsdiff==FALSE)-1)[1],"index"],res_nlvsum_rmseCV_Ysel[nrow(res_nlvsum_rmseCV_Ysel),"index"],na.rm=TRUE)
            }else{
              choiceYH[yy]<- 1
            }
          }
          if(selection=="globalmin"){
            choiceYH[yy] <- res_nlvsum_rmseCV_Ysel[which.min(res_nlvsum_rmseCV_Ysel$mean)[1],"index"]
          }
          if(selection=="1std"){
            if(nrow(res_nlvsum_rmseCV_Ysel)>1){
              # one standard error rule
              minmean    <- which.min(res_nlvsum_rmseCV_Ysel$mean)[1]
              threshmean <- res_nlvsum_rmseCV_Ysel$mean[minmean] + res_nlvsum_rmseCV_Ysel$sd[minmean]
              if((minmean == 1) | (sum(res_nlvsum_rmseCV_Ysel$mean[1:minmean]>=threshmean)==0)){
                choiceYH[yy] <- 1
              }else{
                choiceYH[yy] <- max(res_nlvsum_rmseCV_Ysel[which(res_nlvsum_rmseCV_Ysel$mean[1:minmean]>=threshmean),"index"])
              }
            }else{
              choiceYH[yy] <- 1
            }
          }
        }
        tabfreq <- sapply(1:nrow(lvcombi[[m]]), function(l) sum(choiceYH==l))
        kchoix <- which.max(tabfreq)[1] 
      }

      if(majorityvote==FALSE){
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
        
        res_nlvsum_rmseCV_Ysel <- res_nlvsum_rmseCV_Ysel[order(res_nlvsum_rmseCV_Ysel[,"totalnlv"], decreasing=FALSE),]
        
        if(selection=="localmin"){
          if(nrow(res_nlvsum_rmseCV_Ysel)>1){
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
            # one standard error rule 
            minmean    <- which.min(res_nlvsum_rmseCV_Ysel$mean)[1]
            threshmean <- res_nlvsum_rmseCV_Ysel$mean[minmean] + res_nlvsum_rmseCV_Ysel$sd[minmean]
            if((minmean == 1) | (sum(res_nlvsum_rmseCV_Ysel$mean[1:minmean]>=threshmean)==0)){
              kchoix <- 1
            }else{
              kchoix <- max(res_nlvsum_rmseCV_Ysel[which(res_nlvsum_rmseCV_Ysel$mean[1:minmean]>=threshmean),"index"])
            }
          }else{
            kchoix <- 1
          }
        }
      }
      
      # optim combination
      optimcombi <- as.vector(unlist(lvcombi[[m]][kchoix,]))
      #names(optimcombi) <- paste0("Xlist",1:m)
      nlvlist[[m]] <- lvcombi[[m]][kchoix,m]
    }# end loop on m
 
    ## output names
    names(optimcombi) <- paste0("Xlist",1:nXblocks)
    names(lvcombi) <- names(res_rmseCV_byY) <- names(res_ExplVarCV_byY) <- names(res_rmseCV) <- names(res_ExplVarCV) <- paste0("nbBlocks",1:nXblocks)

    # outputs
    rts <- list(lvcombi=lvcombi,
                #optimCombiLine=kchoix, 
                optimcombi=unlist(lvcombi[[nXblocks]][kchoix,,drop=FALSE]), 
                #optimExplVarCV=res_ExplVarCV[[nXblocks]][kchoix,,drop=FALSE], 
                rmseCV_byY=res_rmseCV_byY, 
                ExplVarCV_byY=res_ExplVarCV_byY, 
                rmseCV=res_rmseCV, 
                ExplVarCV=res_ExplVarCV)

  }
  rts$call   <- match.call()
  class(rts) <- c("Soplscv")
  rts
}


soplsrcv <- function(Xlist, Y, Xscaling = c("none", "pareto", "sd")[1], Yscaling = c("none", "pareto", "sd")[1], weights = NULL, nlvlist=list(), nbrep=30, cvmethod="kfolds", seed = 123, samplingk=NULL, nfolds=7, optimisation="global", selection="1std", majorityvote=FALSE){
  .soplsrcv(Xlist = Xlist, Y = Y, Xscaling = Xscaling, Yscaling = Yscaling, weights = weights, nlvlist=nlvlist, nbrep = nbrep, cvmethod = cvmethod, seed = seed, samplingk = samplingk, nfolds = nfolds, optimisation = optimisation, selection = selection, majorityvote = majorityvote)
}

