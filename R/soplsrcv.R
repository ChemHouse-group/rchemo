
# SOPLScv function
#####################

# ARGUMENTS
## Xlist : A list of matrices or data frames of reference (= training) observations.
## Y : A n x q matrix or data frame, or a vector of length n, of responses.
## Xscaling: vector (of length Xlist) of variable scaling for each datablock, among "none" (mean-centering only), "pareto" (mean-centering and pareto scaling), "sd" (mean-centering and unit variance scaling). 
##          If "pareto" or "sd", uncorrected standard deviation is used.
## Yscaling: variable scaling for the Y-block, among "none" (mean-centering only), "pareto" (mean-centering and pareto scaling), "sd" (mean-centering and unit variance scaling). 
##          If "pareto" or "sd", uncorrected standard deviation is used.
## weights: a priori weights to the rows of the reference matrix in the calculations.
## nbrep : An integer, setting the number of CV repetitions. Default value is 30.
## cvmethod : "kfolds" for k-folds cross-validation
## seed : a numeric. Seed used for the repeated resampling, and if cvmethod is "kfolds" and samplingk is not NULL
## samplingk : A vector of length n. The elements are the values of a qualitative variable used for stratified partition creation
##            if NULL, the first observation is set in the first fold, the second observation in the second fold, etc...
## nfolds : An integer, setting the number of partitions to create. Default value is 7.
## optimisation : "global" or "sequential" optimisation of the number of components.
## nlvXlist : A list of same length as the number of X-blocks. Each component of the list gives the number of PLS components of the corresponding X-block to test.
## selection : a character indicating the selection method to use to choose the optimal combination of components, among "localmin","globalmin","NoSignifDecrease1"
##                "localmin": the optimal combination corresponds to the first local maximum of the mean CV global accuracy
##                "globalmin" : the optimal combination corresponds to the maximum mean CV global accuracy
##                "1std" : it corresponds to the first combination after which the mean cross-validated accuracy does not increase significantly
## majorityvote : only if optimisation is "global" or one X-block.
##             if majorityvote is TRUE, the optimal combination is chosen for each Y variable, with the chosen selection, before a majority vote.
##             if majorityvote is "FALSE, the optimal combination is simply chosen with the chosen selection.


# OUTPUTS
## lvcombi : matrix or list of matrices, of tested component combinations
## optimCombiLine : number of the combination line corresponding to the optimal one. 
##        In the case of a sequential optimisation, it is the number of the combination line in the model with all the X-blocks
## optimcombi : the number of PLS components of each X-block allowing the optimisation of the mean rmseCV
## res_rmseCV_byY : matrix or list of matrices of mean and sd of cross-validated RMSE in the model for each combination and each response variable
## res_rmseC_byY :matrix or list of matrices of mean and sd of RMSE in the model for each combination and each response variable
## res_ExplVarCV_byY : matrix or list of matrices of mean and sd of cross-validated explained variances in the model for each combination and each response variable
## res_ExplVarC_byY : matrix or list of matrices of mean and sd of explained variances in the model for each combination and each response variable
## optimExplVarCV : cross-validated explained variance for the optimal sopls model
## optimExplVarC : calibration explained variance for the optimal sopls model
## res_rmseCV : matrix or list of matrices of mean and sd of cross-validated RMSE in the model for each combination and response variables included in Yselection
## res_rmseC : matrix or list of matrices of mean and sd of RMSE in the model for each combination and response variables included in Yselection
## res_ExplVarCV : matrix or list of matrices of mean and sd of cross-validated explained variances in the model for each combination and response variables included in Yselection
## res_ExplVarC : matrix or list of matrices of mean and sd of explained variances in the model for each combination and response variables included in Yselection

#####################

# DEPENDENCIES

#source("C:\\Users\\mbrandolini\\Documents\\sauvegardeSept17\\DEVELOPPEMENTS\\main\\SOPLS(PM)\\TestDimSOPLS_function.R")


soplsrcv <- function(Xlist, Y, Xscaling = c("none", "pareto", "sd")[1], Yscaling = c("none", "pareto", "sd")[1], weights = NULL, nbrep=30, cvmethod="kfolds", seed = 123, samplingk=NULL, nfolds=7, optimisation="global", nlvXlist=list(), selection="1std", majorityvote=FALSE){
  
  # verifications
  
  if(length(Xlist) != length(nlvXlist)){stop("the number of nbcolXlist or nlvXlist elements is not correct")}
  if((is.vector(Y)==FALSE) & (sum(rownames(Xlist[[1]])!= rownames(Y))>0)){stop("the rownames of Xlist and Y are different")}
  if((is.vector(Y)==TRUE) & (nrow(Xlist[[1]])!= length(Y))){stop("the row numbers of Xlist and Y are different")}
  if((is.null(samplingk)==FALSE) & (nrow(Xlist[[1]])!= length(samplingk))){stop("the length of samplingk is not correct")}
  if((cvmethod!="loo") & (is.null(nfolds)==FALSE)){
    if(nrow(Xlist[[1]])<nfolds){stop("the value of nfolds is not correct")}
  }
  if((nbrep==1)&(selection=="1std")){stop("nbrep must be >1 when selection is '1std'")}
  
  
  # additional funtions
  
  inertie <-function(tab) { # function computing the total variance of a dataset
    tab<- scale(tab, scale=FALSE)
    tab<-as.matrix(tab)
    V<-t(tab)%*%tab
    sum(diag(V))
  }
  
  # argument computations
  n <- nrow(Y)
  nXblocks <- length(Xlist)
  
  Yvariables    <- colnames(Y)
  #nYvariables   <- length(Yvariables)
  
  inertieY <- inertie(Y) #sum(diag(t(as.matrix(scale(Y, scale=FALSE)))%*%(as.matrix(scale(Y, scale=FALSE)))))

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
    
    lvcombi <- as.matrix(expand.grid(nlvXlist))
    colnames(lvcombi) <- paste0(rep("Xlist",nXblocks),1:nXblocks)
    rownames(lvcombi) <- paste0("lvcombi",1:nrow(lvcombi))
    
    # output tables
    res_rmseCV_byY <- matrix(NA, nrow=nrow(lvcombi), ncol=(ncol(Y)*2), 
                         dimnames=list(rownames(lvcombi),c(paste0("mean_",Yvariables),paste0("sd_",Yvariables))))
    res_rmseC_byY  <- res_ExplVarCV_byY <- res_ExplVarC_byY <- res_rmseCV_byY
      
    res_rmseCV <- matrix(NA, nrow=nrow(lvcombi), ncol=2, dimnames=list(rownames(lvcombi),c("mean","sd")))
    res_rmseC  <- res_ExplVarCV <- res_ExplVarC <- res_rmseCV

    Rep_rmseCV_byY <- array(0, dim=c(nrow(lvcombi),nbrep, ncol(Y)), dimnames = list(rownames(lvcombi), paste0("rep",1:nbrep), Yvariables))
    Rep_rmseC_byY  <- Rep_ExplVarCV_byY <- Rep_ExplVarC_byY  <- Rep_rmseCV_byY
    
    Rep_rmseCV <- Rep_rmseC <- Rep_ExplVarCV <- Rep_ExplVarC  <- array(0, dim=c(nrow(lvcombi),nbrep, 1), dimnames = list(paste0("lvcombi",1:nrow(lvcombi)), paste0("rep",1:nbrep), "value"))
    
  	optimYpredCV <- list()
  	
    #####################################################################
    
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
        #reppred   <- matrix(NA, nrow=n, ncol = length(Yvariables))
        for(k in 1:length(CVtype)){
          Xlisttrain <-lapply(1:length(Xlistrep[[i]]),function(x) Xlistrep[[i]][[x]][-CVtype[[k]],,drop=FALSE])
          Ytrain <- Yrep[[i]][-CVtype[[k]],,drop=FALSE]
          Xlisttest <- lapply(1:length(Xlistrep[[i]]),function(x) Xlistrep[[i]][[x]][CVtype[[k]],,drop=FALSE])
          Ytest <- Yrep[[i]][CVtype[[k]],,drop=FALSE]
          repmodel <- soplsr(Xlist=Xlisttrain, Y=Ytrain, Xscaling = Xscaling, Yscaling = Yscaling, weights = weights, nlv=lvcombi[j,])
          reppredCV[CVtype[[k]],] <- predict(repmodel, Xlisttest)
          #reppred[-CVtype[[k]],] <- (repmodel$pred)/(length(CVtype)-1)
        }
        SqErrCV <- (Yrep[[i]]-reppredCV)^2
        #SqErrC <- (Yrep[[i]]-reppred)^2

        Rep_rmseCV_byY[j,i,]    <- sqrt(apply(SqErrCV,2,mean))#repmodel$rmseCV_byY[j,,drop=FALSE]
        #Rep_rmseC_byY[j,i,]     <- sqrt(apply(SqErrC,2,mean))#repmodel$rmseC_byY[j,,drop=FALSE]
        Rep_ExplVarCV_byY[j,i,] <- 1- ((Rep_rmseCV_byY[j,i,]^2)/matrix((apply(Yrep[[i]],2,var)*(n-1)/n), ncol=ncol(Yrep[[i]])))#repmodel$ExplVarCV_byY[j,,drop=FALSE]
        #Rep_ExplVarC_byY[j,i,]  <- 1- ((Rep_rmseC_byY[j,i,]^2)/matrix((apply(Yrep[[i]],2,var)*(n-1)/n), ncol=ncol(Yrep[[i]])))#repmodel$ExplVarC_byY[j,,drop=FALSE]
        
        Rep_rmseCV[j,i,]    <- sqrt(apply(Rep_rmseCV_byY[j,i,,drop=FALSE]^2, FUN = mean, MARGIN = 1))#repmodel$rmseCV[j]
        #Rep_rmseC[j,i,]     <- sqrt(apply(Rep_rmseC_byY[j,i,,drop=FALSE]^2, FUN = mean, MARGIN = 1))#repmodel$rmseC[j]
        Rep_ExplVarCV[j,i,] <- apply(Rep_ExplVarCV_byY[j,i,,drop=FALSE], FUN = mean, MARGIN = 1)#repmodel$ExplVarCV[j]
        #Rep_ExplVarC[j,i,]  <- apply(Rep_ExplVarC_byY[j,i,,drop=FALSE], FUN = mean, MARGIN = 1)#repmodel$ExplVarC[j]
      }
    }

  	res_rmseCV_byY[,1:ncol(Y)]                <- apply(Rep_rmseCV_byY, c(1,3),mean, na.rm=T)
    res_rmseCV_byY[,(ncol(Y)+1):(ncol(Y)*2)] <- apply(Rep_rmseCV_byY, c(1,3),sd, na.rm=T)
    
    # res_rmseC_byY[,1:ncol(Y)]                 <- apply(Rep_rmseC_byY, c(1,3),mean, na.rm=T)
    # res_rmseC_byY[,(ncol(Y)+1):(ncol(Y)*2)]  <- apply(Rep_rmseC_byY, c(1,3),sd, na.rm=T)
    
    res_ExplVarCV_byY[,1:ncol(Y)]                <- apply(Rep_ExplVarCV_byY, c(1,3),mean, na.rm=T)
    res_ExplVarCV_byY[,(ncol(Y)+1):(ncol(Y)*2)] <- apply(Rep_ExplVarCV_byY, c(1,3),sd, na.rm=T)
    
    # res_ExplVarC_byY[,1:ncol(Y)]                 <- apply(Rep_ExplVarC_byY, c(1,3),mean, na.rm=T)
    # res_ExplVarC_byY[,(ncol(Y)+1):(ncol(Y)*2)]  <- apply(Rep_ExplVarC_byY, c(1,3),sd, na.rm=T)
    
    res_rmseCV[,"mean"] <- apply(Rep_rmseCV, c(1,3),mean, na.rm=T)# A VERIFIER QD MULTI Y
    res_rmseCV[,"sd"]   <- apply(Rep_rmseCV, c(1,3),sd, na.rm=T)# A VERIFIER QD MULTI Y
    
    # res_rmseC[,"mean"] <- apply(Rep_rmseC, c(1,3),mean, na.rm=T)# A VERIFIER QD MULTI Y
    # res_rmseC[,"sd"]   <- apply(Rep_rmseC, c(1,3),sd, na.rm=T)# A VERIFIER QD MULTI Y
    
    res_ExplVarCV[,"mean"] <- apply(Rep_ExplVarCV, c(1,3),mean, na.rm=T)# A VERIFIER QD MULTI Y
    res_ExplVarCV[,"sd"]   <- apply(Rep_ExplVarCV, c(1,3),sd, na.rm=T)# A VERIFIER QD MULTI Y
    
    # res_ExplVarC[,"mean"] <- apply(Rep_ExplVarC, c(1,3),mean, na.rm=T)# A VERIFIER QD MULTI Y
    # res_ExplVarC[,"sd"]   <- apply(Rep_ExplVarC, c(1,3),sd, na.rm=T)# A VERIFIER QD MULTI Y
    
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
          t((sapply((unique(nlvsum)), function(i) unlist(data.frame(res_rmseCV_Ysel[(nlvsum==i),c("mean","sd"),drop=FALSE][which.min(res_rmseCV_Ysel[(nlvsum==i),"mean"]),,drop=FALSE])))))
        )
        
        if(selection=="localmin"){
          if(nrow(res_nlvsum_rmseCV_Ysel)>1){
            # sign of the difference of accuracies to select the optim combination with the lower total number of components
            rtsdiff <- c(NA,sapply(2:nrow(res_nlvsum_rmseCV_Ysel), function(i)((res_nlvsum_rmseCV_Ysel$mean[i]-res_nlvsum_rmseCV_Ysel$mean[i-1])<0)))
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
            choiceYH[yy]    <- res_nlvsum_rmseCV_Ysel[which(res_nlvsum_rmseCV_Ysel$mean<threshmean)[1],"index"]
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
        t((sapply((unique(nlvsum)), function(i) unlist(data.frame(res_rmseCV[(nlvsum==i),c("mean","sd"),drop=FALSE][which.min(res_rmseCV[(nlvsum==i),"mean"]),,drop=FALSE])))))
      )
      
      if(selection=="localmin"){
        if(nrow(res_nlvsum_rmseCV_Ysel)>1){
          # sign of the difference of accuracies to select the optim combination with the lower total number of components
          rtsdiff <- c(NA,sapply(2:nrow(res_nlvsum_rmseCV_Ysel), function(i)((res_nlvsum_rmseCV_Ysel$mean[i]-res_nlvsum_rmseCV_Ysel$mean[i-1])<0)))
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
          kchoix     <- res_nlvsum_rmseCV_Ysel[which(res_nlvsum_rmseCV_Ysel$mean<threshmean)[1],"index"]
        }else{
          kchoix <- 1
        }
      }
    }
    

    # outputs
    rts <- list(lvcombi=lvcombi, 
                optimCombiLine=kchoix, 
                optimcombi=unlist(lvcombi[kchoix,,drop=FALSE]), 
                res_rmseCV_byY=res_rmseCV_byY, 
                #res_rmseC_byY=res_rmseC_byY, 
                res_ExplVarCV_byY=res_ExplVarCV_byY, 
                #res_ExplVarC_byY=res_ExplVarC_byY, 
                optimExplVarCV=res_ExplVarCV[kchoix,,drop=FALSE], 
                #optimExplVarC=res_ExplVarC[kchoix,,drop=FALSE],
                res_rmseCV=res_rmseCV, 
                #res_rmseC=res_rmseC, 
                res_ExplVarCV=res_ExplVarCV#, 
                #res_ExplVarC=res_ExplVarC
                )
    
  }
  
  ############################################################################################################################
  
  ## SEQUENTIAL optimisation
  if((optimisation=="sequential") & (nXblocks>1)){
    
    #####################################################################
    
    ## list initialisations
    lvcombi            <- list()

    optimYpredCV      <- list()
    res_rmseCV_byY <- res_rmseC_byY  <- res_ExplVarCV_byY <- res_ExplVarC_byY <- list()
    res_rmseCV     <- res_rmseC      <- res_ExplVarCV     <- res_ExplVarC     <- list()
    Rep_rmseCV_byY   <- Rep_rmseC_byY    <- Rep_ExplVarCV_byY   <- Rep_ExplVarC_byY   <- list()
    Rep_rmseCV       <- Rep_rmseC        <- Rep_ExplVarCV       <- Rep_ExplVarC       <- list()
    
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
        lvcombi[[1]] <- as.matrix(expand.grid(nlvXlist[[1]]))
        colnames(lvcombi[[1]]) <- "Xlist1"
        rownames(lvcombi[[1]]) <- paste0("lvcombi",1:nrow(lvcombi[[1]]))
      }
      if(m>1){
        lvcombi[[m]] <- data.frame(matrix(rep(unlist(optimcombi),length(nlvXlist[[m]])),nrow=length(nlvXlist[[m]]),byrow=TRUE),as.matrix(expand.grid(nlvXlist[[m]])))
        colnames(lvcombi[[m]]) <- paste0(rep("Xlist",m),1:m)
        rownames(lvcombi[[m]]) <- paste0("lvcombi",1:nrow(lvcombi[[m]]))
      }
      
      # output initialisation
      res_rmseCV_byY[[m]] <- matrix(NA, nrow=nrow(lvcombi[[m]]), ncol=(ncol(Y)*2), 
                                    dimnames=list(rownames(lvcombi[[m]]),c(paste0("mean_",Yvariables),paste0("sd_",Yvariables))))
      res_rmseC_byY[[m]]  <- res_ExplVarCV_byY[[m]] <- res_ExplVarC_byY[[m]] <- res_rmseCV_byY[[m]]
      
      res_rmseCV[[m]] <- matrix(NA, nrow=nrow(lvcombi[[m]]), ncol=2, dimnames=list(rownames(lvcombi[[m]]),c("mean","sd")))
      res_rmseC[[m]]  <- res_ExplVarCV[[m]] <- res_ExplVarC[[m]] <- res_rmseCV[[m]]
      
      Rep_rmseCV_byY[[m]] <- array(0, dim=c(nrow(lvcombi[[m]]),nbrep, ncol(Y)), dimnames = list(paste0("lvcombi",1:nrow(lvcombi[[m]])), paste0("rep",1:nbrep), Yvariables))
      Rep_rmseC_byY[[m]]  <- Rep_ExplVarCV_byY[[m]] <- Rep_ExplVarC_byY[[m]]  <- Rep_rmseCV_byY[[m]]
      
      Rep_rmseCV[[m]]  <- Rep_rmseC[[m]]  <- Rep_ExplVarCV[[m]] <- Rep_ExplVarC[[m]]  <- array(0, dim=c(nrow(lvcombi[[m]]),nbrep, 1), dimnames = list(paste0("lvcombi",1:nrow(lvcombi[[m]])), paste0("rep",1:nbrep), "value"))
      
 
      for(i in 1:nbrep){
        for(j in 1:nrow(lvcombi[[m]])){
          reppredCV <- matrix(NA, nrow=n, ncol = length(Yvariables))
          #reppred <- matrix(NA, nrow=n, ncol = length(Yvariables))
          for(k in 1:length(CVtype)){
            Xlisttrain <-lapply(1:m,function(x) Xlistrep[[i]][[x]][-CVtype[[k]],,drop=FALSE])
            Ytrain <- Yrep[[i]][-CVtype[[k]],,drop=FALSE]
            Xlisttest <- lapply(1:m,function(x) Xlistrep[[i]][[x]][CVtype[[k]],,drop=FALSE])
            Ytest <- Yrep[[i]][CVtype[[k]],,drop=FALSE]
            repmodel <- soplsr(Xlist=Xlisttrain, Y=Ytrain, Xscaling = Xscaling, Yscaling = Yscaling, weights = weights, nlv=unlist(lvcombi[[m]][j,]))
            reppredCV[CVtype[[k]],] <- predict(repmodel, Xlisttest)
            #reppred[-CVtype[[k]],] <- repmodel$pred/(length(CVtype)-1)
          }
          SqErrCV <- (Yrep[[i]]-reppredCV)^2
          #SqErrC <- (Yrep[[i]]-reppred)^2
          
          Rep_rmseCV_byY[[m]][j,i,]    <- sqrt(apply(SqErrCV,2,mean))#repmodel$rmseCV_byY[j,,drop=FALSE]
          #Rep_rmseC_byY[j,i,]     <- sqrt(apply(SqErrC,2,mean))#repmodel$rmseC_byY[j,,drop=FALSE]
          Rep_ExplVarCV_byY[[m]][j,i,] <- 1- ((Rep_rmseCV_byY[[m]][j,i,]^2)/matrix((apply(Yrep[[i]],2,var)*(n-1)/n), ncol=ncol(Yrep[[i]])))#repmodel$ExplVarCV_byY[j,,drop=FALSE]
          #Rep_ExplVarC_byY[j,i,]  <- 1- ((Rep_rmseC_byY[[m]][j,i,]^2)/matrix((apply(Yrep[[i]],2,var)*(n-1)/n), ncol=ncol(Yrep[[i]])))#repmodel$ExplVarC_byY[j,,drop=FALSE]
          
          Rep_rmseCV[[m]][j,i,]    <- sqrt(apply(Rep_rmseCV_byY[[m]][j,i,,drop=FALSE]^2, FUN = mean, MARGIN = 1))#repmodel$rmseCV[j]
          #Rep_rmseC[j,i,]     <- sqrt(apply(Rep_rmseC_byY[[m]][j,i,,drop=FALSE]^2, FUN = mean, MARGIN = 1))#repmodel$rmseC[j]
          Rep_ExplVarCV[[m]][j,i,] <- apply(Rep_ExplVarCV_byY[[m]][j,i,,drop=FALSE], FUN = mean, MARGIN = 1)#repmodel$ExplVarCV[j]
          #Rep_ExplVarC[j,i,]  <- apply(Rep_ExplVarC_byY[[m]][j,i,,drop=FALSE], FUN = mean, MARGIN = 1)#repmodel$ExplVarC[j]
        }
      }
      
      res_rmseCV_byY[[m]][,1:ncol(Y)]                <- apply(Rep_rmseCV_byY[[m]], c(1,3),mean, na.rm=T)
      res_rmseCV_byY[[m]][,(ncol(Y)+1):(ncol(Y)*2)] <- apply(Rep_rmseCV_byY[[m]], c(1,3),sd, na.rm=T)
      
      # res_rmseC_byY[[m]][,1:ncol(Y)]                 <- apply(Rep_rmseC_byY[[m]], c(1,3),mean, na.rm=T)
      # res_rmseC_byY[[m]][,(ncol(Y)+1):(ncol(Y)*2)]  <- apply(Rep_rmseC_byY[[m]], c(1,3),sd, na.rm=T)
      
      res_ExplVarCV_byY[[m]][,1:ncol(Y)]                <- apply(Rep_ExplVarCV_byY[[m]], c(1,3),mean, na.rm=T)
      res_ExplVarCV_byY[[m]][,(ncol(Y)+1):(ncol(Y)*2)] <- apply(Rep_ExplVarCV_byY[[m]], c(1,3),sd, na.rm=T)
      
      # res_ExplVarC_byY[[m]][,1:ncol(Y)]                 <- apply(Rep_ExplVarC_byY[[m]], c(1,3),mean, na.rm=T)
      # res_ExplVarC_byY[[m]][,(ncol(Y)+1):(ncol(Y)*2)]  <- apply(Rep_ExplVarC_byY[[m]], c(1,3),sd, na.rm=T)
      
      res_rmseCV[[m]][,"mean"] <- apply(Rep_rmseCV[[m]], c(1,3),mean, na.rm=T)
      res_rmseCV[[m]][,"sd"]   <- apply(Rep_rmseCV[[m]], c(1,3),sd, na.rm=T)
      
      # res_rmseC[[m]][,"mean"] <- apply(Rep_rmseC[[m]], c(1,3),mean, na.rm=T)
      # res_rmseC[[m]][,"sd"]   <- apply(Rep_rmseC[[m]], c(1,3),sd, na.rm=T)
      
      res_ExplVarCV[[m]][,"mean"] <- apply(Rep_ExplVarCV[[m]], c(1,3),mean, na.rm=T)
      res_ExplVarCV[[m]][,"sd"]   <- apply(Rep_ExplVarCV[[m]], c(1,3),sd, na.rm=T)
      
      # res_ExplVarC[[m]][,"mean"] <- apply(Rep_ExplVarC[[m]], c(1,3),mean, na.rm=T)
      # res_ExplVarC[[m]][,"sd"]   <- apply(Rep_ExplVarC[[m]], c(1,3),sd, na.rm=T)
      

      # optim combination for each total number of components // res_rmse_Ysel[,"mean"]
      if(majorityvote==TRUE){
        choiceYH <- rep(NA,length(Yvariables))
        nlvsum <- apply(lvcombi[[m]],1,sum)
        
        for(yy in 1:length(Yvariables)){
          res_rmseCV_Ysel <- cbind.data.frame(res_rmseCV_byY[[m]][,yy,drop=FALSE], res_rmseCV_byY[[m]][,(ncol(Y)+yy),drop=FALSE]) 
          colnames(res_rmseCV_Ysel) <- c("mean","sd")
          rownames(res_rmseCV_Ysel) <- rownames(res_rmseCV_byY[[m]])
          
          res_nlvsum_rmseCV_Ysel <- data.frame(
            index = sapply((unique(nlvsum)), function(i) which(rownames(res_rmseCV_Ysel)==rownames(res_rmseCV_Ysel[(nlvsum==i),,drop=FALSE][which.min(res_rmseCV_Ysel[(nlvsum==i),"mean"]),,drop=FALSE]))),
            combnum = sapply((unique(nlvsum)), function(i) rownames(res_rmseCV_Ysel[(nlvsum==i),,drop=FALSE][which.min(res_rmseCV_Ysel[(nlvsum==i),"mean"]),,drop=FALSE])),
            t((sapply((unique(nlvsum)), function(i) unlist(data.frame(lvcombi[[m]][(nlvsum==i),,drop=FALSE][which.min(res_rmseCV_Ysel[(nlvsum==i),"mean"]),,drop=FALSE]))))),
            t((sapply((unique(nlvsum)), function(i) unlist(data.frame(res_rmseCV_Ysel[(nlvsum==i),c("mean","sd"),drop=FALSE][which.min(res_rmseCV_Ysel[(nlvsum==i),"mean"]),,drop=FALSE])))))
          )
          
          if(selection=="localmin"){
            if(nrow(res_nlvsum_rmseCV_Ysel)>1){
              # sign of the difference of accuracies to select the optim combination with the lower total number of components
              rtsdiff <- c(NA,sapply(2:nrow(res_nlvsum_rmseCV_Ysel), function(i)((res_nlvsum_rmseCV_Ysel$mean[i]-res_nlvsum_rmseCV_Ysel$mean[i-1])<0)))
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
              choiceYH[yy]    <- res_nlvsum_rmseCV_Ysel[which(res_nlvsum_rmseCV_Ysel$mean<threshmean)[1],"index"]
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
            t((sapply((unique(nlvsum)), function(i) unlist(data.frame(res_rmseCV[[m]][(nlvsum==i),c("mean","sd"),drop=FALSE][which.max(res_rmseCV[[m]][(nlvsum==i),"mean"]),,drop=FALSE])))))
          )
        }
        if(m>1){
          res_nlvsum_rmseCV_Ysel <- data.frame(
            index = sapply((unique(nlvsum)), function(i) which(rownames(res_rmseCV[[m]])==rownames(res_rmseCV[[m]][(nlvsum==i),,drop=FALSE][which.max(res_rmseCV[[m]][(nlvsum==i),"mean"]),,drop=FALSE]))),
            combnum = sapply((unique(nlvsum)), function(i) rownames(res_rmseCV[[m]][(nlvsum==i),,drop=FALSE][which.max(res_rmseCV[[m]][(nlvsum==i),"mean"]),,drop=FALSE])),
            t((sapply((unique(nlvsum)), function(i) unlist(data.frame(lvcombi[[m]][(nlvsum==i),,drop=FALSE][which.max(res_rmseCV[[m]][(nlvsum==i),"mean"]),,drop=FALSE]))))),
            t((sapply((unique(nlvsum)), function(i) unlist(data.frame(res_rmseCV[[m]][(nlvsum==i),c("mean","sd"),drop=FALSE][which.max(res_rmseCV[[m]][(nlvsum==i),"mean"]),,drop=FALSE])))))
          )
        }
        
        if(selection=="localmin"){
          if(nrow(res_nlvsum_rmseCV_Ysel)>1){
            # sign of the difference of accuracies to select the optim combination with the lower total number of components
            rtsdiff <- c(NA,sapply(2:nrow(res_nlvsum_rmseCV_Ysel), function(i)((res_nlvsum_rmseCV_Ysel$mean[i]-res_nlvsum_rmseCV_Ysel$mean[i-1])<0)))
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
            kchoix     <- res_nlvsum_rmseCV_Ysel[which(res_nlvsum_rmseCV_Ysel$mean<threshmean)[1],"index"]
          }else{
            kchoix <- 1
          }
        }
      }
      
      # optim combination
      optimcombi <- as.vector(unlist(lvcombi[[m]][kchoix,]))
      names(optimcombi) <- paste0(rep("Xlist",m),1:m)
      nlvXlist[[m]] <- lvcombi[[m]][kchoix,m]
    }# end loop on m
    

    ## output names
    names(lvcombi) <- names(res_rmseCV_byY) <- names(res_rmseC_byY) <- names(res_ExplVarCV_byY) <- names(res_ExplVarC_byY) <- names(res_rmseCV) <- names(res_rmseC) <- names(res_ExplVarCV) <- names(res_ExplVarC) <- paste0("nbBlocks",1:nXblocks)

    # outputs
    rts <- list(lvcombi=lvcombi,
                optimCombiLine=kchoix, 
                optimcombi=unlist(lvcombi[[nXblocks]][kchoix,,drop=FALSE]), 
                res_rmseCV_byY=res_rmseCV_byY, 
                #res_rmseC_byY=res_rmseC_byY, 
                res_ExplVarCV_byY=res_ExplVarCV_byY, 
                #res_ExplVarC_byY=res_ExplVarC_byY, 
                optimExplVarCV=res_ExplVarCV[[nXblocks]][kchoix,,drop=FALSE], 
                #optimExplVarC=res_ExplVarC[[nXblocks]][kchoix,,drop=FALSE],
                res_rmseCV=res_rmseCV, 
                #res_rmseC=res_rmseC, 
                res_ExplVarCV=res_ExplVarCV#, 
                #res_ExplVarC=res_ExplVarC               
                )

  }
  rts$call   <- match.call()
  class(rts) <- c("Soplscv")
  rts
}
############################################################################################################################

if(FALSE){
  rm(list=ls())
  library(FactoMineR)
  library(rchemo)
  
  liste_functions <- list.files("C:/Users/mbrandolini/Documents/DEVELOPPEMENTS/main/rchemo/R", full.names = TRUE)
  # import des données :
  lapply(liste_functions, source)
  
  data(wine)
  #n=nrow(wine)
  
  A=wine[,3:7]
  B=wine[,8:10]
  C=wine[,11:20]
  D=wine[,21:29]  
  E=data.frame(Overall.quality = wine[,30], row.names = rownames(wine)) 
  
  Xlist = list(A,B,C,E)
  Y = D
  nlvXlist=list(0:2,1:ncol(B),0:3,0:1)
    
  nbrep=3
  cvmethod="kfolds"
  seed = 123
  samplingk=NULL
  nfolds=7
  optimisation=c("global","sequential")[1]
  selection=c("1std","localmin","globalmin")[3]
  majorityvote=c(TRUE,FALSE)[1]
  Xscaling = c("none","pareto","sd")[3]
  Yscaling = c("none","pareto","sd")[3]
  weights = NULL
  
  test <- soplsrcv(Xlist, Y, Xscaling = Xscaling, Yscaling=Yscaling, weights = weights, nbrep=nbrep, cvmethod=cvmethod, seed = 123, samplingk=NULL, nfolds=nfolds, optimisation=optimisation, nlvXlist=nlvXlist, selection=selection, majorityvote=FALSE)
    
}
