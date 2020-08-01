# select all 2d NMR spectral files with intact file directories

# files needed: acqus, procs, proc2s, 2rr

# datapath='/Volumes/ANPC_ext1/COVID_urine_highRes_processes/'
# procs_exp=1

checkFiles <- function(datapath, procs_exp=1, n_max=10) {
  
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







# extract metadata from procs and proc2s files (CAMP-JDX)
# f_list - list of directories and exp IDs (output checkFiles function)

extract_pars <- function( f_list ) {
  
  
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
    

chem_shift <- function(swidth, offset, si){
  
  dppm <- swidth/(si - 1) # ho
  cshift <- seq(offset, (offset - swidth), by = -dppm)
  
  return(cshift)
  
}







