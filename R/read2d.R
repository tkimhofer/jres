#' Read-in J-coupling resolved NMR spectra generated with a Bruker spectrometer
#' @param path File path to J-res NMR experiment(s)
#' @param procs_exp Bruker processing experiment (usually 1)
#' @param n_max Maximum number of experiments to read
#' @param inter bilinear interpolation
#' @export
#' @author torben.kimhofer@murdoch.edu.au
#' @section
#path='/Users/TKimhofer/Downloads/'

read2d <- function( path,  procs_exp=1, n_max=1000, interp=T) {
  
  # check folder structure
  f_list <- .checkFiles(datapath = path,  procs_exp, n_max)
  
  # extract parameters from acqus and procs (for f1 and f2)
  pars <- t(.extract_pars ( f_list ))
  
  # convert to numeric where possible
  dtype_num<-apply(pars, 2, function(x){
    !any(is.na(suppressWarnings(as.numeric(x))))
  })
  
  pars<-data.frame(pars)
  pars[,dtype_num]<-apply(pars[,dtype_num], 2, as.numeric)
  rownames(pars)<-f_list[[2]]
  
  # read in binary file and 
  out <- lapply(1:length(f_list[[1]]), function(s){
    
    # chem shift
    csF2_ppm <- .chem_shift(swidth=pars$af2_SW[s], offset=pars$pf2_OFFSET[s], si=pars$pf2_SI[s])
    csF1_hz <- .chem_shift(swidth=pars$af1_SW_h[s], offset=pars$pf1_OFFSET[s]*pars$pf1_SF[s], si=pars$pf1_SI[s])
    
    byteorda=c('little'=0, 'big'=1)
    names(byteorda)[match(pars$af2_BYTORDA[s], byteorda)]
    
    
    spec <- readBin(paste0(f_list[[1]][s], .Platform$file.sep, 'pdata', .Platform$file.sep, procs_exp, .Platform$file.sep, '2rr'), 
                    what = "int", 
                    n = pars$pf1_FTSIZE[s]*pars$pf2_FTSIZE[s],
                    size = 4, 
                    signed = T, 
                    endian = names(byteorda)[match(pars$af2_BYTORDA[s], byteorda)]
    ) # this is spectra
    spec <- ( spec * (2^pars$af2_NC[s]) )
    nspec <- length(spec)
    
    spec2d=reshape_submatToMat(pars$pf2_SI[s], pars$pf1_SI[s], pars$pf2_XDIM[s], pars$pf1_XDIM[s], spec)
    colnames(spec2d) <- csF2_ppm
    rownames(spec2d) <- csF1_hz
    
    return(spec2d)
    
  })
  
  # if(interp==TRUE) {
  #   
  #   # not implented yet
  #   # vec_span=lapply(out, function(x) {
  #   #   browser()
  #   #   ppm=as.numeric(colnames(x))
  #   #   hz=as.numeric(rownames(x))
  #   #   data.frame(f2=c(min(ppm), max(ppm), diff(ppm[1:2])), 
  #   #              f1=c(min(hz), max(hz), diff(hz[1:2])))
  #   # 
  #   # })
  # }
  
  return(out)
}

#' Select all 2d NMR spectral files with intact file directories
#' @param datapath File path to J-res NMR experiment(s)
#' @param procs_exp Bruker processing experiment (usually 1)
#' @param n_max Maximum number of experiments to read
#' @details Files needed: acqus, procs, proc2s, 2rr
#' @author torben.kimhofer@murdoch.edu.au
#' @section

.checkFiles <- function(datapath, procs_exp=1, n_max=10) {
  
  datapath=gsub(paste0(.Platform$file.sep, '$'), '', datapath)
  
  # searches for f2 processing parameter files
  f_procs <- list.files(path = datapath, pattern ="procs",
                        all.files = FALSE, full.names = TRUE, recursive = TRUE,
                        ignore.case = TRUE)
  f_procs=grep(paste0('pdata', .Platform$file.sep, procs_exp), f_procs, value=T)
  
  # searches for f1 processing parameter files
  f_proc2s <- list.files(path = datapath, pattern = "^proc2s$",
                         all.files = FALSE, full.names = TRUE, recursive = TRUE,
                         ignore.case = TRUE)
  f_proc2s=grep(paste0('pdata', .Platform$file.sep, procs_exp), f_proc2s, value=T)
  
  # searches for 2rr files
  f_2rr <- list.files(path = datapath, pattern = "^2rr$", all.files = FALSE,
                      full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  f_2rr=grep(paste0('pdata', .Platform$file.sep, procs_exp), f_2rr, value=T)
  
  #browser()
  if(n_max < length(f_proc2s) ) { f_proc2s=f_proc2s[1:n_max]; message('Reached n_max - not all spectra read-in.') }
  if(length(f_procs) == 0 | length(f_proc2s) == 0 | length(f_2rr) == 0) { message('No spectrum found'); return(NULL) }
  
  lev=strsplit(f_proc2s[1], .Platform$file.sep)[[1]] # three last are descriptors
  dpath=strsplit(datapath, .Platform$file.sep)[[1]]
  dpath_desc=grep(dpath[length(dpath)], lev)
  
  #browser()
  if( (length(lev) - dpath_desc) == 3 ) { # single file read-in
    exp_no=dpath[dpath_desc]
    # no filter needed as length statement filter for completeness
    p_intact=datapath
  }else{
    p_procs=gsub('/pdata/1/procs$', '', f_procs)
    p_proc2s=gsub('/pdata/1/proc2s', '', f_proc2s)
    p_2rr=gsub('/pdata/1/2rr$', '', f_2rr)
    
    p_intact=intersect( intersect(p_procs, p_proc2s), p_2rr)
    exp_no=as.numeric(gsub(paste0(datapath, .Platform$file.sep), '', p_intact))
    
    if(any(is.na(exp_no))) { stop('Experiment number(s) is/are not numeric!') }
    
    idx=order(exp_no)
    p_intact=p_intact[idx]
    exp_no=exp_no[idx]
    
  }
  
  return(list(path=p_intact, exp_no=exp_no))
}



#' Extract metadata from procs and proc2s files (JCAMP-JDX)
#' @param f_list List of directories and exp IDs (output checkFiles function)
#' @author torben.kimhofer@murdoch.edu.au
#' @section
.extract_pars <- function( f_list ) {
  
  
  out=sapply(f_list[[1]], function(fil){
    
    #browser()
    # procs
    f_procs=paste0(fil, .Platform$file.sep, 'pdata/1/procs')
    # extract procs information for t2
    fhand <- file(f_procs, open = "r")
    f_procs <- readLines(fhand, n = -1, warn = FALSE)
    close(fhand)
    
    out=strsplit(gsub('^##\\$', '',  grep('^##\\$', f_procs, value=T, fixed = F), fixed = F), '=')
    d_procs_val=gsub('^ ', '', sapply(out, '[[', 2))
    names(d_procs_val) = paste0('pf2_', sapply(out, '[[', 1))
    
    # proc2s
    f_proc2s=paste0(fil, .Platform$file.sep, 'pdata/1/proc2s')
    # extract procs information for t2
    fhand <- file(f_proc2s, open = "r")
    f_proc2s <- readLines(fhand, n = -1, warn = FALSE)
    close(fhand)
    
    out=strsplit(gsub('^##\\$', '',  grep('^##\\$', f_proc2s, value=T, fixed = F), fixed = F), '=')
    d_proc2s_val=gsub('^ ', '', sapply(out, '[[', 2))
    names(d_proc2s_val) = paste0('pf1_', sapply(out, '[[', 1))
    
    
    # acqus
    f_acqu=paste0(fil, .Platform$file.sep, 'acqus')
    # extract procs information for t2
    fhand <- file(f_acqu, open = "r")
    f_acqu <- readLines(fhand, n = -1, warn = FALSE)
    close(fhand)
    
    out=strsplit(gsub('^##\\$', '',  grep('^##\\$', f_acqu, value=T, fixed = F), fixed = F), '=')
    d_acqu_val=gsub('^ ', '', sapply(out, '[[', 2))
    names(d_acqu_val) = paste0('af2_', sapply(out, '[[', 1))
    
    # change date
    idx=grep('date', names(d_acqu_val), ignore.case = T)
    d_acqu_val[idx]=as.character(as.POSIXct(x = '01/01/1970 00:00:00', format='%d/%m/%Y %H:%M:%S')+(as.numeric(d_acqu_val[idx])))
    
    # acqu2s
    f_acqu2=paste0(fil, .Platform$file.sep, 'acqu2s')
    # extract procs information for t2
    fhand <- file(f_acqu2, open = "r")
    f_acqu2 <- readLines(fhand, n = -1, warn = FALSE)
    close(fhand)
    
    out=strsplit(gsub('^##\\$', '',  grep('^##\\$', f_acqu2, value=T, fixed = F), fixed = F), '=')
    d_acqu2_val=gsub('^ ', '', sapply(out, '[[', 2))
    names(d_acqu2_val) = paste0('af1_', sapply(out, '[[', 1))
    
    pars = c(d_acqu_val, d_acqu2_val, d_procs_val, d_proc2s_val)
    
    return(pars)
    
  })
  
  
  return(out)
  
}

#' Calculate chemical shift vector for a NMR dimension
#' @param swidth Sweep width / spectral width
#' @param offset Offset
#' @param si Number of data points of the respective NMR dimension
#' @author torben.kimhofer@murdoch.edu.au
#' @section
.chem_shift <- function(swidth, offset, si){
  
  dppm <- swidth/(si - 1) # ho
  cshift <- seq(offset, (offset - swidth), by = -dppm)
  
  return(cshift)
  
}




#' Calibrate chemical shift vector for a NMR dimension
#' @param spec2d list of 2D matrices with row and column names representing chemical shift values
#' @param f2_ra F2 range of expected peak of calibration signal
#' @param f1_ra F1 range of expected peak of calibration signal
#' @param f2_fix F2 axis position of calibrated signal
#' @param f2_fix F1 axis position of calibrated signal
#' @export
#' @author torben.kimhofer@murdoch.edu.au
#' @section
.calibrate2d <- function( spec2d, f2_ra=c(-0.1, 0.1), f1_ra=c(260, 3020), f2_fix=0, f1_fix=0) {
  
  out=lapply(spec2d, function(sp){
    
    f2=as.numeric(colnames(sp))
    f1=as.numeric(rownames(sp))
    
    idx_f2=get.idx(f2_ra, f2)
    idx_f1=get.idx(f1_ra, f1)
    
    idx_max=which(sp[idx_f1, idx_f2] == max(sp[idx_f1, idx_f2]), arr.ind=T)
    
    f1_max=f1[idx_f1][idx_max[1,1]]
    f2_max=f2[idx_f2][idx_max[1,2]]
    
    f2_new<-(f2-f2_max)+f2_fix
    f1_new<-(f1-f1_max)+f1_fix
    
    colnames(sp)<-f2_new
    rownames(sp)<-f1_new
    
    sp
  })
  
  return(out)
}
























