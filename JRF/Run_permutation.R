#' Run JRF for permuted data.
#'
#' MAIN FUNCTION -- > JRF_permutation
#' 
#' INPUT
#' 
#' X            list object containing expression data for each class
#' ntree        number of trees
#' mtry         number of variables to be sampled at each node
#' genes.name   vector containing gene names 
#' M            number of permutations
#'  
#' OUTPUT: importance score of gene-gene interactions for permuted data.
#'
#'
#' OTHER FUNCTIONS -- > importance  and  JRF_onetarget
#' 
#' importance     compute importance score for an object of class JRF 
#' (this file is a modified version of file importance contained in package randomForest, A. Liaw and M. Wiener (2002))
#' 
#' JRF_onetarget  for each class, model the expression of a target gene as a function of the expression of other genes via random forest. 
#'                class specific tree ensemble are designed to borrow information across them. 
#' (this file is a modified version of file randomForest contained in package randomForest, A. Liaw and M. Wiener (2002))
#'   
#'
#' @export 
#"Run_permutation" <-  function(X, ...)UseMethod("JRF")


"JRF_permutation" <-
  function(X, ntree,mtry,genes.name,perm) {
    
    nclasses<-length(X)
    sampsize<-rep(0,nclasses)
    
    for (j in 1:nclasses) sampsize[j]<-dim(X[[j]])[2]
    
    totsize<-tot<-max(sampsize);
    p<-dim(X[[1]])[1];
    imp<-array(0,c(p,length(genes.name),nclasses))
    
    imp.final<-matrix(0,p*(p-1)/2,nclasses);
    vec1<-matrix(rep(genes.name,p),p,p)
    vec2<-t(vec1)
    vec1<-vec1[lower.tri(vec1,diag=FALSE)]
    vec2<-vec2[lower.tri(vec2,diag=FALSE)]
    
    index<-seq(1,p)
    
    for (j in 1:length(genes.name)){
      
      covar<-matrix(0,(p-1)*nclasses,tot)             
      y<-matrix(0,nclasses,tot)             
      
      set.seed((perm-1)*nclasses+1)
      for (c in 1:nclasses)  {
        y[c,seq(1,sampsize[c])]<-as.matrix(X[[c]][j,sample(sampsize[c])])
        covar[seq((c-1)*(p-1)+1,c*(p-1)),seq(1,sampsize[c])]<-X[[c]][-j,]
      }
      
      jrf.out<-JRF_onetarget(x=covar,y=y,mtry=mtry,importance=TRUE,sampsize=sampsize,nclasses=nclasses,ntree=ntree)
      
      for (s in 1:nclasses) imp[-j,j,s]<-importance(jrf.out,scale=FALSE)[seq((p-1)*(s-1)+1,(p-1)*(s-1)+p-1)]  #- save importance score for net1
      
    }
    
    # --- Derive importance score for each interaction 
    for (s in 1:nclasses){ 
      imp.s<-imp[,,s]; t.imp<-t(imp.s)
      imp.final[,s]<-(imp.s[lower.tri(imp.s,diag=FALSE)]+t.imp[lower.tri(t.imp,diag=FALSE)])/2        
    }
    out<-as.data.frame(imp.final)
    colnames(out)<-paste0('importance',1:nclasses)
    return(out)
    
  }

"Run_permutation" <-
  function(X, ntree,mtry,genes.name,M) {
    p<-length(genes.name)
    nclasses<-length(X)
    perm<-array(0,c((p^2-p)/2,M,nclasses))
    
    for (j in 1:M)  perm[,j,]<-as.matrix(JRF_permutation(X, ntree,mtry,genes.name,j))  
    
    return(perm)
  }
    