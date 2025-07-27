#' Derive threshold for importance scores and return final networks.
#'
#' INPUT
#' 
#' out.jrf      output of JRF function
#' out.perm     output of Run_permutation function
#' TH          fdr threshold
#'
#' @export 
#"JRF_network" <-  function(out, ...)UseMethod("JRF")


JRF_network <- function(out.jrf,out.perm,TH) {
 
  nclasses<-dim(out.perm)[3]
  M<-dim(out.perm)[2]; 

  for (net in 1:nclasses) { 
    j.np<-sort(out.jrf[,2+net],decreasing=TRUE)
    FDR<-matrix(0,dim(out.perm)[1],1); 
    for (s in 1:length(j.np)){ 
      FP<-sum(sum(out.perm[,,net]>=j.np[s]))/M
      FDR[s]<-FP/s;
      if (FDR[s]>TH) {th<-j.np[s];  break;}
    }
  if (net==1) out<-list(out.jrf[out.jrf[,2+net]>th,seq(1,2)])
  if (net>1) out<-list(out, out.jrf[out.jrf[,2+net]>th,seq(1,2)])
  }
 
  return(out)
}