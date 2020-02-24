#' @title Noise threshold estimation for 2D NMR spectrum
#' @description
#' Noise in JRES experiments is believed to be Rayleigh (RL) -like distributed (RL looks similar to a right-skewed normal dist). A good  noise cut-off is three times noise standard deviation (SD), however, due to right-skewness there will be few noise spikes  (Ref PhD thesis HM Parson, U Birmingham).
#' For each t1
#' 1. bin spectum
#' 2. divide binned spectrum in 32 equidistant sections
#' 3. noise level is three times smallest SD found in all sections
#' 4. take lower 10% value as noise cutoff for entire spectrum
#'
#' Noise is likely to be different in different t1 areas, therefore a future application might focus on de-noising JRES/2D with location-specific noise values.
#'
#' @param jres J-coupling resolved NMR spectrum (df or matrix, t1 x t2)
#' @param binsize Bin size in f2 dimension (ppm)
#' @param sections window size f1, f2
#' @param sdlim Multiplier for SD to determine noise cutoff in each bin (usually value of three)
#' @param probs Quantile probability determining final noise cutoff for 2D
#' @author [Torben Kimhofer](https://tkimhofer.com)
#' @importFrom stats sd quantile
#' @export

noise2D=function(jres, binsize=0.02, sections=32, sdlim=3, probs=0.1){
  
  ppm=as.numeric(colnames(jres))
  sbin=bin.specs(jres, ppm, bins=binsize)[[1]]
  
  # assign section to each bin
  idx=rep(1:sections, each=ceiling(ncol(sbin)/sections))[1:ncol(sbin)]
  
  out=apply(sbin, 1, function(t1, sec=unique(idx)){
    sec.sd=sapply(sec, function(s, t1.cut=t1, iid=idx){
      sd(t1.cut[iid==s])
    })
    min(sec.sd)*sdlim
  })
  
  noise=quantile(out, probs, na.rm=T)
  
  return(noise)
  
}
