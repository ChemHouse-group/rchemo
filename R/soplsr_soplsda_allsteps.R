

soplsr_soplsda_allsteps <- function(Xlist, Xnames = NULL, Xscaling = c("none","pareto","sd")[1], 
                                    Y, Yscaling = c("none","pareto","sd")[1], weights = NULL,
                                    newXlist = NULL, newXnames = NULL,
                                    
                                    method = c("soplsr", "soplsrda","soplslda","soplsqda")[1],
                                    prior = c("unif", "prop")[1], 
                                    
                                    step = c("nlvtest","permutation","model","prediction")[1],
                                    nlv = c(),
                                    nlvlist=list(), 
                                    modeloutput = c("scores","loadings","coef","vip"), 
                                    
                                    cvmethod = c("kfolds","loo")[1], 
                                    nbrep = 30, 
                                    seed = 123, 
                                    samplingk = NULL, 
                                    nfolds = 10, 
                                    npermut = 30, 
                                    
                                    optimisation = c("global","sequential")[1], 
                                    # majorityvote = FALSE, 
                                    criterion = c("err","rmse")[1], 
                                    selection = c("localmin","globalmin","1std")[1],
                                    
                                    import = c("R","ChemFlow","W4M")[1],
                                    outputfilename = NULL
                                  
){
  
  # IMPORT -------------------------------------------------------------------------------------------

  majorityvote = FALSE
  
  if(is.null(Xnames)==TRUE){Xnames <- paste0(rep("X",length(Xlist)),1:length(Xlist))}
  if((is.null(newXlist)==FALSE)&(is.null(newXnames)==TRUE)){newXnames <- paste0(rep("Xnew",length(newXlist)),1:length(newXlist))}
  
  if(length(Xscaling) == 1){Xscaling = rep(Xscaling, length(Xlist))}
  
  if(length(method)>1){stop("length of method must be 1")}
  if(length(prior)>1){stop("length of prior must be 1")} 
  if(length(cvmethod)>1){stop("length of cvmethod must be 1")}
  if(length(step)>1){stop("length of step must be 1")} 
  
  if(length(optimisation)>1){stop("length of optimisation must be 1")} 
  if(length(majorityvote)>1){stop("length of majorityvote must be 1")}
  if(length(criterion)>1){stop("length of criterion must be 1")}
  if(length(selection)>1){stop("length of selection must be 1")} 
  if(length(nbrep)>1){stop("length of nbrep must be 1")}
  if(length(seed)>1){stop("length of seed must be 1")} 
  if(length(nfolds)>1){stop("length of nfolds must be 1")}
  if(length(npermut)>1){stop("length of npermut must be 1")} 
  
  options(max.print=99999)
  
  ### Y
  if(import == "R"){# Y matrice n x p, avec identifiants en rownames, au format Rdata
    Y           <- as.matrix(Y)
  }
  if(import %in% c("ChemFlow","W4M")){# Y matrice n x p, avec identifiants en 1ere colonne, au format .txt
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
    if(import == "ChemFlow"){
      Xlist[[i]]           <- read.table(Xlist[[i]],sep="\t",dec=".",header=TRUE, na.strings=c("","NA"), check.names=FALSE)
      rownames(Xlist[[i]]) <- Xlist[[i]][,1]
      Xlist[[i]][,1]       <- NULL
    }
    if(import == "W4M"){
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
    ### newX
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
    if(method == c("soplsr", "soplsrda","soplslda","soplsqda")[1]){
      resnlvtest <- soplsrcv(Xlist = Xlist, 
                             Y = Y, 
                             Xscaling = Xscaling, 
                             Yscaling = Yscaling, 
                             weights = weights, 
                             nlvlist = nlvlist, 
                             nbrep = nbrep, 
                             cvmethod = cvmethod, 
                             seed = seed, 
                             samplingk = samplingk, 
                             nfolds = nfolds, 
                             optimisation = optimisation, 
                             selection = selection, 
                             majorityvote = majorityvote)
      if(optimisation == "sequential"){
        bcombi <- sapply(1:length(resnlvtest$lvcombi), function(nr) nrow(resnlvtest$lvcombi[[nr]]))
        ncombi <- sum(bcombi)
        allcombi <- matrix(0, ncol = length(Xlist), nrow = ncombi,
                           dimnames = list(c(sapply(1:length(resnlvtest$lvcombi), function(i) paste0(names(resnlvtest$lvcombi)[i], "_",rownames(resnlvtest$lvcombi[[i]])))), paste0("XBlock",1:length(Xlist))))
        for(i in 1 : length(Xlist)){
          allcombi[(c(0,cumsum(bcombi))+1)[i]:(c(0,cumsum(bcombi)))[i+1],1:i] <- unlist(resnlvtest$lvcombi[[i]][c(1,cumsum(bcombi))[1]:c(1,cumsum(bcombi))[2],1:i])
        }
        optimum <- rep(0, nrow(allcombi))
        for(i in 1:nrow(allcombi)){
          if(sum(resnlvtest$optimcombi != allcombi[i,])==0){
            optimum[i] <- 1
          }
        }
        resnlvtesttable <- cbind(allcombi, nlvsum = apply(allcombi, MARGIN = 1, FUN = sum, na.rm = TRUE), optimum = optimum, do.call("rbind",resnlvtest$rmseCV))
        
      }else{
        optimum <- rep(0, nrow(resnlvtest$lvcombi))
        for(i in 1:nrow(resnlvtest$lvcombi)){
          if(sum(resnlvtest$optimcombi != resnlvtest$lvcombi[i,])==0){
            optimum[i] <- 1
          }
        }
        resnlvtesttable <- cbind(resnlvtest$lvcombi, nlvsum = apply(resnlvtest$lvcombi, MARGIN = 1, FUN = sum, na.rm = TRUE), optimum = optimum, resnlvtest$rmseCV)
      }
      if(is.null(outputfilename)==FALSE){
        write.table(data.frame(rownames(resnlvtesttable),resnlvtesttable), 
                    file = paste0(outputfilename,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                    eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                    col.names = TRUE, qmethod = c("escape", "double"),
                    fileEncoding = "")
      }
      
    }
    
    
    if(method == c("soplsr", "soplsrda","soplslda","soplsqda")[2]){
      resnlvtest <- soplsrdacv(Xlist = Xlist, 
                             y = Y, 
                             Xscaling = Xscaling, 
                             Yscaling = Yscaling, 
                             weights = weights, 
                             nlvlist = nlvlist, 
                             nbrep = nbrep, 
                             cvmethod = cvmethod, 
                             seed = seed, 
                             samplingk = samplingk, 
                             nfolds = nfolds, 
                             optimisation = optimisation, 
                             criterion = criterion,
                             selection = selection)
      if(optimisation == "sequential"){
        bcombi <- sapply(1:length(resnlvtest$lvcombi), function(nr) nrow(resnlvtest$lvcombi[[nr]]))
        ncombi <- sum(bcombi)
        allcombi <- matrix(0, ncol = length(Xlist), nrow = ncombi,
                           dimnames = list(c(sapply(1:length(resnlvtest$lvcombi), function(i) paste0(names(resnlvtest$lvcombi)[i], "_",rownames(resnlvtest$lvcombi[[i]])))), paste0("XBlock",1:length(Xlist))))
        for(i in 1 : length(Xlist)){
          allcombi[(c(0,cumsum(bcombi))+1)[i]:(c(0,cumsum(bcombi)))[i+1],1:i] <- unlist(resnlvtest$lvcombi[[i]][c(1,cumsum(bcombi))[1]:c(1,cumsum(bcombi))[2],1:i])
        }
        optimum <- rep(0, nrow(allcombi))
        for(i in 1:nrow(allcombi)){
          if(sum(resnlvtest$optimcombi != allcombi[i,])==0){
            optimum[i] <- 1
          }
        }
        if(criterion == "err"){
          resnlvtesttable <- cbind(allcombi, nlvsum = apply(allcombi, MARGIN = 1, FUN = sum, na.rm = TRUE), optimum = optimum, do.call("rbind",resnlvtest$res_errCV))
        }else{
          resnlvtesttable <- cbind(allcombi, nlvsum = apply(allcombi, MARGIN = 1, FUN = sum, na.rm = TRUE), optimum = optimum, do.call("rbind",resnlvtest$res_rmseCV))
        }
      }
      if(optimisation == "global"){
        optimum <- rep(0, nrow(resnlvtest$lvcombi))
        for(i in 1:nrow(resnlvtest$lvcombi)){
          if(sum(resnlvtest$optimcombi != resnlvtest$lvcombi[i,])==0){
            optimum[i] <- 1
          }
        }
        if(criterion == "err"){
          resnlvtesttable <- cbind(resnlvtest$lvcombi, nlvsum = apply(resnlvtest$lvcombi, MARGIN = 1, FUN = sum, na.rm = TRUE), optimum = optimum, resnlvtest$res_errCV)
        }else{
          resnlvtesttable <- cbind(resnlvtest$lvcombi, nlvsum = apply(resnlvtest$lvcombi, MARGIN = 1, FUN = sum, na.rm = TRUE), optimum = optimum, resnlvtest$res_rmseCV)
        }
      }
      if(is.null(outputfilename)==FALSE){
        write.table(data.frame(rownames(resnlvtesttable),resnlvtesttable), 
                    file = paste0(outputfilename,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                    eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                    col.names = TRUE, qmethod = c("escape", "double"),
                    fileEncoding = "")
      }
      
    }
    
    
    if(method == c("soplsr", "soplsrda","soplslda","soplsqda")[3]){
      resnlvtest <- soplsldacv(Xlist = Xlist, 
                               y = Y, 
                               Xscaling = Xscaling, 
                               Yscaling = Yscaling, 
                               weights = weights, 
                               nlvlist = nlvlist, 
                               prior = prior,
                               nbrep = nbrep, 
                               cvmethod = cvmethod, 
                               seed = seed, 
                               samplingk = samplingk, 
                               nfolds = nfolds, 
                               optimisation = optimisation, 
                               criterion = criterion,
                               selection = selection)
      if(optimisation == "sequential"){
        bcombi <- sapply(1:length(resnlvtest$lvcombi), function(nr) nrow(resnlvtest$lvcombi[[nr]]))
        ncombi <- sum(bcombi)
        allcombi <- matrix(0, ncol = length(Xlist), nrow = ncombi,
                           dimnames = list(c(sapply(1:length(resnlvtest$lvcombi), function(i) paste0(names(resnlvtest$lvcombi)[i], "_",rownames(resnlvtest$lvcombi[[i]])))), paste0("XBlock",1:length(Xlist))))
        for(i in 1 : length(Xlist)){
          allcombi[(c(0,cumsum(bcombi))+1)[i]:(c(0,cumsum(bcombi)))[i+1],1:i] <- unlist(resnlvtest$lvcombi[[i]][c(1,cumsum(bcombi))[1]:c(1,cumsum(bcombi))[2],1:i])
        }
        optimum <- rep(0, nrow(allcombi))
        for(i in 1:nrow(allcombi)){
          if(sum(resnlvtest$optimcombi != allcombi[i,])==0){
            optimum[i] <- 1
          }
        }
        if(criterion == "err"){
          resnlvtesttable <- cbind(allcombi, nlvsum = apply(allcombi, MARGIN = 1, FUN = sum, na.rm = TRUE), optimum = optimum, do.call("rbind",resnlvtest$res_errCV))
        }else{
          resnlvtesttable <- cbind(allcombi, nlvsum = apply(allcombi, MARGIN = 1, FUN = sum, na.rm = TRUE), optimum = optimum, do.call("rbind",resnlvtest$res_rmseCV))
        }
      }
      if(optimisation == "global"){
        optimum <- rep(0, nrow(resnlvtest$lvcombi))
        for(i in 1:nrow(resnlvtest$lvcombi)){
          if(sum(resnlvtest$optimcombi != resnlvtest$lvcombi[i,])==0){
            optimum[i] <- 1
          }
        }
        if(criterion == "err"){
          resnlvtesttable <- cbind(resnlvtest$lvcombi, nlvsum = apply(resnlvtest$lvcombi, MARGIN = 1, FUN = sum, na.rm = TRUE), optimum = optimum, resnlvtest$res_errCV)
        }else{
          resnlvtesttable <- cbind(resnlvtest$lvcombi, nlvsum = apply(resnlvtest$lvcombi, MARGIN = 1, FUN = sum, na.rm = TRUE), optimum = optimum, resnlvtest$res_rmseCV)
        }
      }
      if(is.null(outputfilename)==FALSE){
        write.table(data.frame(rownames(resnlvtesttable),resnlvtesttable), 
                    file = paste0(outputfilename,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                    eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                    col.names = TRUE, qmethod = c("escape", "double"),
                    fileEncoding = "")
      }
      
    }
    
    if(method == c("soplsr", "soplsrda","soplslda","soplsqda")[4]){
      resnlvtest <- soplsqdacv(Xlist = Xlist, 
                               y = Y, 
                               Xscaling = Xscaling, 
                               Yscaling = Yscaling, 
                               weights = weights, 
                               nlvlist = nlvlist, 
                               prior = prior,
                               nbrep = nbrep, 
                               cvmethod = cvmethod, 
                               seed = seed, 
                               samplingk = samplingk, 
                               nfolds = nfolds, 
                               optimisation = optimisation, 
                               criterion = criterion,
                               selection = selection)
      if(optimisation == "sequential"){
        bcombi <- sapply(1:length(resnlvtest$lvcombi), function(nr) nrow(resnlvtest$lvcombi[[nr]]))
        ncombi <- sum(bcombi)
        allcombi <- matrix(0, ncol = length(Xlist), nrow = ncombi,
                           dimnames = list(c(sapply(1:length(resnlvtest$lvcombi), function(i) paste0(names(resnlvtest$lvcombi)[i], "_",rownames(resnlvtest$lvcombi[[i]])))), paste0("XBlock",1:length(Xlist))))
        for(i in 1 : length(Xlist)){
          allcombi[(c(0,cumsum(bcombi))+1)[i]:(c(0,cumsum(bcombi)))[i+1],1:i] <- unlist(resnlvtest$lvcombi[[i]][c(1,cumsum(bcombi))[1]:c(1,cumsum(bcombi))[2],1:i])
        }
        optimum <- rep(0, nrow(allcombi))
        for(i in 1:nrow(allcombi)){
          if(sum(resnlvtest$optimcombi != allcombi[i,])==0){
            optimum[i] <- 1
          }
        }
        if(criterion == "err"){
          resnlvtesttable <- cbind(allcombi, nlvsum = apply(allcombi, MARGIN = 1, FUN = sum, na.rm = TRUE), optimum = optimum, do.call("rbind",resnlvtest$res_errCV))
        }else{
          resnlvtesttable <- cbind(allcombi, nlvsum = apply(allcombi, MARGIN = 1, FUN = sum, na.rm = TRUE), optimum = optimum, do.call("rbind",resnlvtest$res_rmseCV))
        }
      }
      if(optimisation == "global"){
        optimum <- rep(0, nrow(resnlvtest$lvcombi))
        for(i in 1:nrow(resnlvtest$lvcombi)){
          if(sum(resnlvtest$optimcombi != resnlvtest$lvcombi[i,])==0){
            optimum[i] <- 1
          }
        }
        if(criterion == "err"){
          resnlvtesttable <- cbind(resnlvtest$lvcombi, nlvsum = apply(resnlvtest$lvcombi, MARGIN = 1, FUN = sum, na.rm = TRUE), optimum = optimum, resnlvtest$res_errCV)
        }else{
          resnlvtesttable <- cbind(resnlvtest$lvcombi, nlvsum = apply(resnlvtest$lvcombi, MARGIN = 1, FUN = sum, na.rm = TRUE), optimum = optimum, resnlvtest$res_rmseCV)
        }
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
    if(method == c("soplsr", "soplsrda","soplslda","soplsqda")[1]){
      resmodel <- soplsr(Xlist = Xlist, 
                             Y = Y, 
                             Xscaling = Xscaling, 
                             Yscaling = Yscaling, 
                             weights = weights, 
                             nlv = nlv)
      if("scores" %in% modeloutput){
        Tx <- list()
        for(i in 1:length(Xlist)){
          if(resmodel$nlv[i]>0){
            Tx[[i]] <- resmodel$fm[[i]]$T
            if(is.null(outputfilename)==FALSE){
              write.table(data.frame(rownames(Tx[[i]]),Tx[[i]]), 
                          file = paste0(outputfilename,"_scores_Block",i,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                          eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                          col.names = TRUE, qmethod = c("escape", "double"),
                          fileEncoding = "")
            }
            
          }
        }
      }
      if("loadings" %in% modeloutput){
        P <- list()
        for(i in 1:length(Xlist)){
          if(resmodel$nlv[i]>0){
            P[[i]] <- resmodel$fm[[i]]$P
            if(is.null(outputfilename)==FALSE){
              write.table(data.frame(rownames(P[[i]]),P[[i]]), 
                          file = paste0(outputfilename,"_loadings_Block",i,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                          eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                          col.names = TRUE, qmethod = c("escape", "double"),
                          fileEncoding = "")
            }
            
          }
        }
      }
      if("coef" %in% modeloutput){
        coeflist <- list()
        for(i in 1:length(Xlist)){
          if(resmodel$nlv[i]>0){
            coeflist[[i]] <- coef(resmodel$fm[[i]], nlv = resmodel$nlv[i])$B
            if(is.null(outputfilename)==FALSE){
              write.table(data.frame(rownames(coeflist[[i]]),coeflist[[i]]), 
                          file = paste0(outputfilename,"_coef_Block",i,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                          eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                          col.names = TRUE, qmethod = c("escape", "double"),
                          fileEncoding = "")
            }
            
          }
        }
      }
      if("vip" %in% modeloutput){
        viplist <- list()
        for(i in 1:length(Xlist)){
          if(resmodel$nlv[i]>0){
            viplist[[i]] <- vip(resmodel$fm[[i]], Xlist[[i]], Y=NULL, nlv = resmodel$nlv[i])
            if(is.null(outputfilename)==FALSE){
              write.table(data.frame(rownames(viplist[[i]]),viplist[[i]]), 
                          file = paste0(outputfilename,"_vip_Block",i,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                          eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                          col.names = TRUE, qmethod = c("escape", "double"),
                          fileEncoding = "")
            }
            
          }
        }
      }
    }
    
    if(method == c("soplsr", "soplsrda","soplslda","soplsqda")[2]){
      resmodel <- soplsrda(Xlist = Xlist, 
                           y = Y, 
                           Xscaling = Xscaling, 
                           Yscaling = Yscaling, 
                           weights = weights, 
                           nlv = nlv)
      if("scores" %in% modeloutput){
        Tx <- list()
        for(i in 1:length(Xlist)){
          if(resmodel$fm$nlv[i]>0){
            Tx[[i]] <- resmodel$fm$fm[[i]]$T
            if(is.null(outputfilename)==FALSE){
              write.table(data.frame(rownames(Tx[[i]]),Tx[[i]]), 
                          file = paste0(outputfilename,"_scores_Block",i,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                          eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                          col.names = TRUE, qmethod = c("escape", "double"),
                          fileEncoding = "")
            }
            
          }
        }
      }
      if("loadings" %in% modeloutput){
        P <- list()
        for(i in 1:length(Xlist)){
          if(resmodel$fm$nlv[i]>0){
            P[[i]] <- resmodel$fm$fm[[i]]$P
            if(is.null(outputfilename)==FALSE){
              write.table(data.frame(rownames(P[[i]]),P[[i]]), 
                          file = paste0(outputfilename,"_loadings_Block",i,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                          eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                          col.names = TRUE, qmethod = c("escape", "double"),
                          fileEncoding = "")
            }
            
          }
        }
      }
      if("coef" %in% modeloutput){
        coeflist <- list()
        for(i in 1:length(Xlist)){
          if(resmodel$fm$nlv[i]>0){
            coeflist[[i]] <- coef(resmodel$fm$fm[[i]], nlv = resmodel$fm$nlv[i])$B
            if(is.null(outputfilename)==FALSE){
              write.table(data.frame(rownames(coeflist[[i]]),coeflist[[i]]), 
                          file = paste0(outputfilename,"_coef_Block",i,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                          eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                          col.names = TRUE, qmethod = c("escape", "double"),
                          fileEncoding = "")
            }
            
          }
        }
      }
      if("vip" %in% modeloutput){
        viplist <- list()
        for(i in 1:length(Xlist)){
          if(resmodel$fm$nlv[i]>0){
            viplist[[i]] <- vip(resmodel$fm$fm[[i]], Xlist[[i]], Y=NULL, nlv = resmodel$fm$nlv[i])
            if(is.null(outputfilename)==FALSE){
              write.table(data.frame(rownames(viplist[[i]]),viplist[[i]]), 
                          file = paste0(outputfilename,"_vip_Block",i,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                          eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                          col.names = TRUE, qmethod = c("escape", "double"),
                          fileEncoding = "")
            }
            
          }
        }
      }
    }
    
    if(method == c("soplsr", "soplsrda","soplslda","soplsqda")[3]){
      resmodel <- soplslda(Xlist = Xlist, 
                           y = Y, 
                           Xscaling = Xscaling, 
                           Yscaling = Yscaling, 
                           weights = weights, 
                           nlv = nlv,
                           prior = prior)
      if("scores" %in% modeloutput){
        Tx <- list()
        for(i in 1:length(Xlist)){
          if(resmodel$fm[[1]]$nlv[i]>0){
            Tx[[i]] <- resmodel$fm[[1]]$fm[[i]]$T
            if(is.null(outputfilename)==FALSE){
              write.table(data.frame(rownames(Tx[[i]]),Tx[[i]]), 
                          file = paste0(outputfilename,"_scores_Block",i,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                          eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                          col.names = TRUE, qmethod = c("escape", "double"),
                          fileEncoding = "")
            }
            
          }
        }
      }
      if("loadings" %in% modeloutput){
        P <- list()
        for(i in 1:length(Xlist)){
          if(resmodel$fm[[1]]$nlv[i]>0){
            P[[i]] <- resmodel$fm[[1]]$fm[[i]]$P
            if(is.null(outputfilename)==FALSE){
              write.table(data.frame(rownames(P[[i]]),P[[i]]), 
                          file = paste0(outputfilename,"_loadings_Block",i,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                          eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                          col.names = TRUE, qmethod = c("escape", "double"),
                          fileEncoding = "")
            }
            
          }
        }
      }
      if("coef" %in% modeloutput){
        coeflist <- list()
        for(i in 1:length(Xlist)){
          if(resmodel$fm[[1]]$nlv[i]>0){
            coeflist[[i]] <- coef(resmodel$fm[[1]]$fm[[i]], nlv = resmodel$fm[[1]]$nlv[i])$B
            if(is.null(outputfilename)==FALSE){
              write.table(data.frame(rownames(coeflist[[i]]),coeflist[[i]]), 
                          file = paste0(outputfilename,"_coef_Block",i,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                          eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                          col.names = TRUE, qmethod = c("escape", "double"),
                          fileEncoding = "")
            }
            
          }
        }
      }
      if("vip" %in% modeloutput){
        viplist <- list()
        for(i in 1:length(Xlist)){
          if(resmodel$fm[[1]]$nlv[i]>0){
            viplist[[i]] <- vip(resmodel$fm[[1]]$fm[[i]], Xlist[[i]], Y=NULL, nlv = resmodel$fm[[1]]$nlv[i])
            if(is.null(outputfilename)==FALSE){
              write.table(data.frame(rownames(viplist[[i]]),viplist[[i]]), 
                          file = paste0(outputfilename,"_vip_Block",i,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                          eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                          col.names = TRUE, qmethod = c("escape", "double"),
                          fileEncoding = "")
            }
            
          }
        }
      }
    }
    
    if(method == c("soplsr", "soplsrda","soplslda","soplsqda")[4]){
      resmodel <- soplsqda(Xlist = Xlist, 
                           y = Y, 
                           Xscaling = Xscaling, 
                           Yscaling = Yscaling, 
                           weights = weights, 
                           nlv = nlv,
                           prior = prior)
      if("scores" %in% modeloutput){
        Tx <- list()
        for(i in 1:length(Xlist)){
          if(resmodel$fm[[1]]$nlv[i]>0){
            Tx[[i]] <- resmodel$fm[[1]]$fm[[i]]$T
            if(is.null(outputfilename)==FALSE){
              write.table(data.frame(rownames(Tx[[i]]),Tx[[i]]), 
                          file = paste0(outputfilename,"_scores_Block",i,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                          eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                          col.names = TRUE, qmethod = c("escape", "double"),
                          fileEncoding = "")
            }
            
          }
        }
      }
      if("loadings" %in% modeloutput){
        P <- list()
        for(i in 1:length(Xlist)){
          if(resmodel$fm[[1]]$nlv[i]>0){
            P[[i]] <- resmodel$fm[[1]]$fm[[i]]$P
            if(is.null(outputfilename)==FALSE){
              write.table(data.frame(rownames(P[[i]]),P[[i]]), 
                          file = paste0(outputfilename,"_loadings_Block",i,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                          eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                          col.names = TRUE, qmethod = c("escape", "double"),
                          fileEncoding = "")
            }
            
          }
        }
      }
      if("coef" %in% modeloutput){
        coeflist <- list()
        for(i in 1:length(Xlist)){
          if(resmodel$fm[[1]]$nlv[i]>0){
            coeflist[[i]] <- coef(resmodel$fm[[1]]$fm[[i]], nlv = resmodel$fm[[1]]$nlv[i])$B
            if(is.null(outputfilename)==FALSE){
              write.table(data.frame(rownames(coeflist[[i]]),coeflist[[i]]), 
                          file = paste0(outputfilename,"_coef_Block",i,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                          eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                          col.names = TRUE, qmethod = c("escape", "double"),
                          fileEncoding = "")
            }
            
          }
        }
      }
      if("vip" %in% modeloutput){
        viplist <- list()
        for(i in 1:length(Xlist)){
          if(resmodel$fm[[1]]$nlv[i]>0){
            viplist[[i]] <- vip(resmodel$fm[[1]]$fm[[i]], Xlist[[i]], Y=NULL, nlv = resmodel$fm[[1]]$nlv[i])
            if(is.null(outputfilename)==FALSE){
              write.table(data.frame(rownames(viplist[[i]]),viplist[[i]]), 
                          file = paste0(outputfilename,"_vip_Block",i,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                          eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                          col.names = TRUE, qmethod = c("escape", "double"),
                          fileEncoding = "")
            }
            
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
      outputs$coef <- coeflist
    }
    if("vip" %in% modeloutput){
      outputs$vip <- viplist
    }
    return(outputs)
  }
  
  
  ### permutation
  
  if (step == "permutation"){
    if(method == c("soplsr", "soplsrda","soplslda","soplsqda")[1]){
      
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
      
      set.seed(seed = seed)
      
      permut_list <- lapply(1:npermut, function(i){sample(1:nrow(Y),nrow(Y))})
      permut_list[[(npermut+1)]] <- 1:nrow(Y)
      
      permut_dyssimilarity <- c()

      res_permut <- c()
      for(i in 1:(npermut+1)){
        Ypermut <- Y[permut_list[[i]],,drop=FALSE]
        rownames(Ypermut) <- rownames(Y)
        permut_dyssimilarity[i] <- 1-coefficientRV(Ypermut,Y)
        res_permut[i] <- soplsrcv(Xlist = Xlist, 
                 Y = Ypermut, 
                 Xscaling = Xscaling, 
                 Yscaling = Yscaling, 
                 weights = weights, 
                 nlvlist = as.list(nlv), 
                 nbrep = nbrep, 
                 cvmethod = cvmethod, 
                 seed = seed, 
                 samplingk = samplingk, 
                 nfolds = nfolds, 
                 optimisation = c("global","sequential")[1], 
                 selection = c("localmin","globalmin","1std")[2], 
                 majorityvote = FALSE)$rmsecv
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
    
    if(method == c("soplsr", "soplsrda","soplslda","soplsqda")[2]){
      set.seed(seed = seed)
      
      permut_list <- lapply(1:npermut, function(i){sample(1:nrow(Y),nrow(Y))})
      permut_list[[(npermut+1)]] <- 1:nrow(Y)
      
      permut_dyssimilarity <- c()

      res_permut <- c()
      for(i in 1:(npermut+1)){
        Ypermut <- Y[permut_list[[i]],,drop=FALSE]
        rownames(Ypermut) <- rownames(Y)
        permut_dyssimilarity[i] <- err(Ypermut,Y)
        if(criterion == "err"){
          res_permut[i] <- soplsrdacv(Xlist = Xlist, 
                                      y = Ypermut, 
                                      Xscaling = Xscaling, 
                                      Yscaling = Yscaling, 
                                      weights = weights, 
                                      nlvlist = as.list(nlv), 
                                      nbrep = nbrep, 
                                      cvmethod = cvmethod, 
                                      seed = seed, 
                                      samplingk = samplingk, 
                                      nfolds = nfolds, 
                                      optimisation = c("global","sequential")[1], 
                                      criterion = criterion,
                                      selection = c("localmin","globalmin","1std")[2])$res_errCV[,"mean"]
        }else{
          res_permut[i] <- soplsrdacv(Xlist = Xlist, 
                                      y = Ypermut, 
                                      Xscaling = Xscaling, 
                                      Yscaling = Yscaling, 
                                      weights = weights, 
                                      nlvlist = as.list(nlv), 
                                      nbrep = nbrep, 
                                      cvmethod = cvmethod, 
                                      seed = seed, 
                                      samplingk = samplingk, 
                                      nfolds = nfolds, 
                                      optimisation = c("global","sequential")[1], 
                                      criterion = criterion,
                                      selection = c("localmin","globalmin","1std")[2])$resnlvtest$res_rmseCV[,"mean"]
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
    
    if(method == c("soplsr", "soplsrda","soplslda","soplsqda")[3]){
      set.seed(seed = seed)
      
      permut_list <- lapply(1:npermut, function(i){sample(1:nrow(Y),nrow(Y))})
      permut_list[[(npermut+1)]] <- 1:nrow(Y)
      
      permut_dyssimilarity <- c()

      res_permut <- c()
      for(i in 1:(npermut+1)){
        Ypermut <- Y[permut_list[[i]],,drop=FALSE]
        rownames(Ypermut) <- rownames(Y)
        permut_dyssimilarity[i] <- err(Ypermut,Y)
        if(criterion == "err"){
          res_permut[i] <- soplsldacv(Xlist = Xlist, 
                                      y = Ypermut, 
                                      Xscaling = Xscaling, 
                                      Yscaling = Yscaling, 
                                      weights = weights, 
                                      nlvlist = as.list(nlv), 
                                      prior = prior,
                                      nbrep = nbrep, 
                                      cvmethod = cvmethod, 
                                      seed = seed, 
                                      samplingk = samplingk, 
                                      nfolds = nfolds, 
                                      optimisation = c("global","sequential")[1], 
                                      criterion = criterion,
                                      selection = c("localmin","globalmin","1std")[2])$res_errCV[,"mean"]
        }else{
          res_permut[i] <- soplsldacv(Xlist = Xlist, 
                                      y = Ypermut, 
                                      Xscaling = Xscaling, 
                                      Yscaling = Yscaling, 
                                      weights = weights, 
                                      nlvlist = as.list(nlv), 
                                      prior = prior,
                                      nbrep = nbrep, 
                                      cvmethod = cvmethod, 
                                      seed = seed, 
                                      samplingk = samplingk, 
                                      nfolds = nfolds, 
                                      optimisation = c("global","sequential")[1], 
                                      criterion = criterion,
                                      selection = c("localmin","globalmin","1std")[2])$resnlvtest$res_rmseCV[,"mean"]
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
    if(method == c("soplsr", "soplsrda","soplslda","soplsqda")[4]){
      set.seed(seed = seed)
      
      permut_list <- lapply(1:npermut, function(i){sample(1:nrow(Y),nrow(Y))})
      permut_list[[(npermut+1)]] <- 1:nrow(Y)
      
      permut_dyssimilarity <- c()

      res_permut <- c()
      for(i in 1:(npermut+1)){
        Ypermut <- Y[permut_list[[i]],,drop=FALSE]
        rownames(Ypermut) <- rownames(Y)
        permut_dyssimilarity[i] <- err(Ypermut,Y)
        if(criterion == "err"){
          res_permut[i] <- soplsqdacv(Xlist = Xlist, 
                                      y = Ypermut, 
                                      Xscaling = Xscaling, 
                                      Yscaling = Yscaling, 
                                      weights = weights, 
                                      nlvlist = as.list(nlv), 
                                      prior = prior,
                                      nbrep = nbrep, 
                                      cvmethod = cvmethod, 
                                      seed = seed, 
                                      samplingk = samplingk, 
                                      nfolds = nfolds, 
                                      optimisation = c("global","sequential")[1], 
                                      criterion = criterion,
                                      selection = c("localmin","globalmin","1std")[2])$res_errCV[,"mean"]
        }else{
          res_permut[i] <- soplsqdacv(Xlist = Xlist, 
                                      y = Ypermut, 
                                      Xscaling = Xscaling, 
                                      Yscaling = Yscaling, 
                                      weights = weights, 
                                      nlvlist = as.list(nlv),
                                      prior = prior,
                                      nbrep = nbrep, 
                                      cvmethod = cvmethod, 
                                      seed = seed, 
                                      samplingk = samplingk, 
                                      nfolds = nfolds, 
                                      optimisation = c("global","sequential")[1], 
                                      criterion = criterion,
                                      selection = c("localmin","globalmin","1std")[2])$resnlvtest$res_rmseCV[,"mean"]
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
    if(method == c("soplsr", "soplsrda","soplslda","soplsqda")[1]){
      resmodel <- soplsr(Xlist = Xlist, 
                         Y = Y, 
                         Xscaling = Xscaling, 
                         Yscaling = Yscaling, 
                         weights = weights, 
                         nlv = nlv)
      respredict <- cbind(transform(resmodel, newXlist),predict(resmodel, newXlist))
      if(is.null(outputfilename)==FALSE){
        write.table(data.frame(rownames(respredict),respredict), 
                    file = paste0(outputfilename,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                    eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                    col.names = TRUE, qmethod = c("escape", "double"),
                    fileEncoding = "")
      }
      
    }
    if(method == c("soplsr", "soplsrda","soplslda","soplsqda")[2]){
      resmodel <- soplsrda(Xlist = Xlist, 
                           y = Y, 
                           Xscaling = Xscaling, 
                           Yscaling = Yscaling, 
                           weights = weights, 
                           nlv = nlv)
      respredict <- cbind.data.frame(transform(resmodel, newXlist),pred=predict(resmodel, newXlist))
      if(is.null(outputfilename)==FALSE){
        write.table(data.frame(rownames(respredict),respredict), 
                    file = paste0(outputfilename,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                    eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                    col.names = TRUE, qmethod = c("escape", "double"),
                    fileEncoding = "")
      }
      
    }
    if(method == c("soplsr", "soplsrda","soplslda","soplsqda")[3]){
      resmodel <- soplslda(Xlist = Xlist, 
                           y = Y, 
                           Xscaling = Xscaling, 
                           Yscaling = Yscaling, 
                           weights = weights, 
                           nlv = nlv,
                           prior = prior)
      respredict <- cbind.data.frame(transform(resmodel, newXlist),pred=predict(resmodel, newXlist))
      if(is.null(outputfilename)==FALSE){
        write.table(data.frame(rownames(respredict),respredict), 
                    file = paste0(outputfilename,".txt"), append = FALSE, quote = TRUE, sep = "\t",
                    eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                    col.names = TRUE, qmethod = c("escape", "double"),
                    fileEncoding = "")
      }
      
    }
    if(method == c("soplsr", "soplsrda","soplslda","soplsqda")[4]){
      resmodel <- soplsqda(Xlist = Xlist, 
                           y = Y, 
                           Xscaling = Xscaling, 
                           Yscaling = Yscaling, 
                           weights = weights, 
                           nlv = nlv,
                           prior = prior)
      respredict <- cbind.data.frame(transform(resmodel, newXlist),pred=predict(resmodel, newXlist))
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
