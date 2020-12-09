#' @title Create artificial 1D NMR by Jres feature projection
#' @description Function uses cauchy distribution (aka Lorentz distr) parameterised with max feature position on f2+f1 to determine location and 1D fwhm estimate to determine scale. Cauchy is scaled in intensity according to observed intensity of the resp. feature.
#' @param plist Dataframe of Jres peak summaries
#' @param jr Jres spectrum
#' @param sf Spectrometer frequency (MHz)
#' @param lwi estimated full widht at half max (fwhm) of 1D (usually TSP signal is fine, see spec.qc function)
#' @author torben.kimhofer@murdoch.edu.au
#' @export


projFeat_f1Tof2=function(plist, jr, sf, lw.ppm){
  
  f2.ppm=as.numeric(colnames(jr))
  #f1.hz=as.numeric(rownames(jr))
  
  projM=lapply(1:nrow(plist), function(i, ppl=plist, js=jr, lw=lw.ppm, f2p=f2.ppm){
    feat=ppl[i,]
    center=feat$cent.f2+(feat$cent.f1/sf)
    
    
    out=dcauchy(f2p, loc=center, scale=lw)
    #idx=order(out, decreasing = T)[1:4]
    # scale caucy based on feature inensity
    out=out*(feat$Int/max(out))
    out
    
    
  })
  
  
  
  if(is.list(projM)){
    projM=do.call(rbind, projM[sapply(projM, length)>0])
    projM[is.na(projM)]=0
    proj=apply(projM, 2, sum)
  }else{
    projM[is.na(projM)]=0
    proj=apply(projM, 1, sum)
  }
  
  return(proj)
  
  
}
