#' 1D spectra equidistant binning
#' @param spec NMR matrix, spectra in rows (matrix)
#' @param ppm Chemcial shift position of NMR matrix columns (vector)
#' @param bins Desired bin size (ppm)
#' @author [Torben Kimhofer](https://tkimhofer.com)
#' @export


bin.specs=function(spec, ppm, bins=0.005){
  if(ppm[1]<ppm[2]){ppm=rev(ppm); spec=spec[,ncol(spec):1]}
  
  # determine idx length for respectiv bin size
  d=ppm[1]-ppm[2]
  d.idx=abs(floor(bins/d))
  cat('Exact bin size is', abs(ppm[1]-ppm[d.idx]), '\n')
  # how many new variables (bins) and which ppm belong to one bin
  idx.new=floor(ncol(spec)/d.idx)
  idx.bin=rep(1:idx.new, each=d.idx)
  
  # perform binning (X)
  X.bin=sapply(1:idx.new, function(i, idx=idx.bin, X=spec){
    res=apply(X[,which(idx==i)], 1, sum)
  })
  
  # binning of ppm variable (mean of each bin)
  ppm.bin=sapply(1:idx.new, function(i, idx=idx.bin, pp=ppm){
    mean(pp[which(idx==i)])
  })
  
  return(list(X.bin, ppm.bin))
}



