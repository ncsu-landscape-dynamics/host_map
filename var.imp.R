var.select <- function(){
  require(biomod2)
  require(plyr)
  
  algos <- list(c('SRE'), c('GLM'), c('RF'), c('MAXENT.Phillips'))
  names(algos) <- c('SRE', 'GLM', 'RF', 'MAXENT.Phillips')
  #algos <- list(c('SRE')); names(algos) <- c('SRE')
  
  i <- 1
  k <- 1; k.list <- list()
  j <- 1; j.list <- list(); j.data <- data.frame()
  t.list <- list()
  
  for(i in 1:length(algos)){
    i.algo <- algos[[i]]
    n.algo <- names(algos)[i]
    i.t1 <- Sys.time()
    for(j in 1:1){#reps of same algo
      for(k in 1:max(envi.cv$cluster)){#reps through clusters
        if(k==1){k.names <- names(envi); k.stack <-NULL}
        
        k.test <- ldply(.data=k.names,
                        .fun=function(X){
                          myExpl <- stack(envi[[X]], k.stack)
                          PA.df <- as.data.frame(myResp); PA.df[is.na(PA.df)] <- FALSE
                          PA.fact <- sum(PA.df==F) / sum(myResp, na.rm=T)
                          PA.df$myResp[which(PA.df==F)[round((1:sum(myResp, na.rm=T))*PA.fact)]] <- TRUE
                          PA.df$myResp <- as.logical(PA.df$myResp)
                          
                          myOptions <- BIOMOD_ModelingOptions('GLM'=list(test='none'))
                          myData <- BIOMOD_FormatingData(resp.var = myResp,
                                                         expl.var = myExpl,
                                                         resp.xy = myRespXY,
                                                         resp.name = myRespName,
                                                         PA.nb.rep = 1,
                                                         PA.strategy = 'user.defined',
                                                         PA.table=PA.df)
                          
                          # 3. Computing the models
                          myModels <- biomod2::BIOMOD_Modeling(data=myData, 
                                                               models.options = myOptions,
                                                               models = i.algo, 
                                                               NbRunEval = 1, #number of runs
                                                               DataSplit = 100, #50,
                                                               Prevalence = NULL,
                                                               VarImport = 0,
                                                               models.eval.meth = c('ROC', 'TSS'),
                                                               SaveObj = T, # recommended to leave true
                                                               rescal.all.models = F, #experimental don't use
                                                               do.full.models = F,
                                                               modeling.id=paste(myRespName,"FirstModeling",sep=""))
                          
                          # Evaluation
                          myEval <- get_evaluations(myModels) # get all models evaluation
                          eval.v <- c(myEval[c('ROC','TSS'),"Testing.data",,,]); return(mean(eval.v))
                          # myEval <- t(myEval[c('ROC', 'TSS'), "Testing.data",,,]); myEval <- cbind(myEval, data.frame(rowMeans(myEval))); colnames(myEval)[3] <- 'score'; return(myEval)
                        })
        row.names(k.test) <- k.names
        
        k.out <- data.frame(cbind(envi.cv$cluster[row.names(k.test)], k.test$V1))
        colnames(k.out) <- c('cluster', 'score')
        k.out <- k.out[order(k.out$cluster),]
        k.out[,3] <- as.integer(rank(-k.out$score))
        k.var <- row.names(k.out)[which.min(k.out$V3)]
        k.score <- k.out$score[which.max(k.out$score)]
        k.cluster <- k.out[k.var, 'cluster']
        k.envi <- envi[[k.var]]
        k.stack <- stack(k.envi, k.stack)
        k.list[[k]] <- list('score'=k.score, 'stack'=k.stack)
        
        if(k==1){redux <- 3; k.redux <- data.frame()
        for(i in 1:length(unique(k.out$cluster))){
          i.clust <- k.out[k.out$cluster==i,]
          i.redux <- i.clust[order(i.clust$score, decreasing=T),][1:redux,]
          k.redux <- rbind(k.redux, i.redux)}
        k.out <- k.redux}
        k.names <-  row.names(k.out)[which(k.out$cluster!=k.cluster)]
      }
      
      j.length <- max(envi.cv$cluster)
      j.list[[j]] <- k.list[[which.max(lapply(X=c(1:j.length), FUN=function(X){k.list[[X]]$score}))]]$stack
      j.stack <- j.list[[j]]
      j.score <- max(unlist(lapply(X=c(1:j.length), FUN=function(X){k.list[[X]]$score})))
      j.vars <- rev(names(j.stack))
      j.v1 <- data.frame(t(data.frame(1:j.length))[0,])
      j.v2 <- data.frame(t(data.frame(j.vars)))
      j.v3 <- rbind.fill(j.v1, j.v2); colnames(j.v3) <- paste('var', 1:j.length, sep = '')
      j.v4 <- data.frame('algo'=n.algo, 'score'=j.score, j.v3)
      j.data <- rbind(j.data, j.v4)
    }
    i.t2 <- Sys.time(); i.time <- i.t2-i.t1
    t.list[[i]] <- i.time
    #unique(i.data$var)
    #i.data <- rbind(i.data, j.data)
  }
  i.data <- ldply(.data=c(1:length(algos)),
                  .fun=function(X){
                    X.df <- data.frame('algo'=j.data$algo[X], 'score'=j.data$score[X],
                                       'var'=data.frame(t(j.data[X, which(!colnames(j.data)%in%c('algo', 'score'))]))[,1])
                    if(any(is.na(X.df$var))){X.df$var[which(is.na(X.df$var))] <- paste('NA', which(is.na(X.df$var)), sep='_')}
                    return(X.df)})  
  i.d2 <- ldply(.data=unique(i.data$var), .fun=function(X){data.frame('var'=X, 'score'=sum(i.data$score[which(i.data$var==X)]))})
  i.d2 <- i.d2[order(i.d2$score, decreasing=T),]
  i.d2 <- i.d2[1:max(envi.cv$cluster),]
  i.stack <- envi[[i.d2[1:max(envi.cv$cluster),'var'][!grepl('NA_',i.d2[1:max(envi.cv$cluster),'var'])]]]
  
  return(i.stack)
}