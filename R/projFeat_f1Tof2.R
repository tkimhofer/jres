#' @title Create artificial 1D NMR by Jres feature projection
#' @description Function uses cauchy distribution (aka Lorentz distr) parameterised with max feature position on f2+f1 to determine location and 1D fwhm estimate to determine scale. Cauchy is scaled in intensity according to observed intensity of the resp. feature.
#' @param plist Dataframe of Jres peak summaries
#' @param jr Jres spectrum
#' @param sf Spectrometer frequency (MHz)
#' @param lwi estimated full widht at half max (fwhm) of 1D (usually TSP signal is fine, see spec.qc function)
#' @author [Torben Kimhofer](https://tkimhofer.com)
#' @export



projFeat_f1Tof2=function(plist, jr, sf, lwi){
  
  f2.ppm=as.numeric(colnames(jr))
  #f1.hz=as.numeric(rownames(jr))
  
  
  projM=apply(plist, 1, function(feat, js=jr, lw=lwi){
    center=feat[2]+((feat[3])/sf)
    idx.row=feat[8]:feat[9]
    idx.col=feat[10]:feat[11]
    
    if(length(idx.row)>2 & length(idx.col)>4){
      sub=js[idx.row, idx.col]
      
      out=dcauchy(f2.ppm, loc=center, scale=lw)
      #idx=order(out, decreasing = T)[1:4]
      # scale caucy based on feature inensity
      out=out*(feat[12]/max(out))
      out
    }
    
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
