#' Plotting of Jres including projections, inclusion of resp. 1D
#' @param ds Jres with f1 in rows, f2 in cols
#' @param t2.lim Plotting range in f2 dimension
#' @param t1.lim Plotting range in f1 dimension
#' @param spec.1d 1D NMR spectrum 
#' @param ppm.1d Matching ppm scale for 1D NMR spec
#' @param z.probs Probability value for quantile function (to exlude lower intensities, e.g. noise)
#' @param SF Spec frequency
#' @param bbox dataframe of detected features for bounding box drawing (optional, can be set to NULL)
#' @param title Figure title
#' @param pPeak Probability threshold for drawing a feature's bounding box, part of bbox df and determined by SVM classification
#' @param tile Should spec be tilted by 45 degree (depreciated)
#' @param addProjpp Vector of feature projections created by function *projFeat_f1_Tof2*
#' @author tkimhofer@gmail.com
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @importFrom reshape2 melt
#' @importFrom scales pretty_breaks
#' @importFrom colorRamps matlab.like2
#' @importFrom MetaboMate get.idx

plotjres_overlay1D_bboxNull=function(ds, t2.lim=c(3,3.1), t1.lim=c(-10,10), spec.1d, ppm.1d, z.probs=0,  SF=600 , bbox,  title='', pPeak= 0.6, tilt=F, addProjpp=NULL){

  linMap <- function(x, from, to)
    (x - min(x)) / max(x - min(x)) * (to - from) + from
  
  ds[is.na(ds)]=0
  
  
  # subset data to desired ranges (f1, f2 and ppm for 1D)
  t2.ppm=as.numeric(colnames(ds))
  t1.Hz=as.numeric(rownames(ds))
  
  t2.idx=get.idx(range = t2.lim, t2.ppm)  # 2D
  t1.idx=get.idx(range = t1.lim, t1.Hz)  # 2D
  x.1d.idx=get.idx(t2.lim, ppm.1d)                          # 1D)
  
  ds=as.data.frame(ds[t1.idx, t2.idx], stringsAsFactors = F)
  
  hz2d=t1.Hz[t1.idx]
  ppm2d=t2.ppm[t2.idx]
  
  
  # produce projections and 1D
  t1.sum=data.frame(sum=apply(ds, 1, sum, na.rm=T), ppm=hz2d, group='Sum')
  t1.sky=data.frame(sum=apply(ds, 1, max, na.rm=T), ppm=hz2d, group='Skyline')
  t1.sum$sum=linMap(t1.sum$sum, from=t2.lim[1], to=t2.lim[2])
  t1.sky$sum=linMap(t1.sky$sum, from=t2.lim[1], to=t2.lim[2])
  t1.sum=rbind(t1.sum, t1.sky)
  
  t2.sum=data.frame(sum=apply(ds, 2, sum, na.rm=T), ppm=ppm2d, group='Sum')
  t2.sky=data.frame(sum=apply(ds, 2, max, na.rm=T), ppm=ppm2d, group='Skyline')
  t2.1d=data.frame(sum=spec.1d[x.1d.idx], ppm=ppm.1d[x.1d.idx], group='1D')
  
  
  
  # # produce p45 deg tilted JRES
  #  out=t(sapply(1:length(t1.Hz),function(d, dss=ds, xh=t2.ppm*SF, yh=t1.Hz){
  #   cord=dss[d,]
  #   xnew=xh-yh[d]
  #   sInter = approxfun(xnew,cord)
  #   sInter(xh)
  # } ))
  # colnames(out)=t2.ppm
  # rownames(out)=t1.Hz
  
  
  # produce projections and 1D
  # t2.tilt=data.frame(sum=apply(out, 2, max, na.rm=T), ppm=ppm2d, group='T45 Skyline')
  # t2.tilt1=data.frame(sum=apply(out, 2, sum, na.rm=T), ppm=ppm2d, group='T45 Sum')
  
  # map projections and 1D to same scale
  t2.sum$sum=linMap(t2.sum$sum, from=t1.lim[1], to=t1.lim[2])
  t2.sky$sum=linMap(t2.sky$sum, from=t1.lim[1], to=t1.lim[2])
  t2.1d$sum=linMap(t2.1d$sum, from=t1.lim[1], to=t1.lim[2])
  # t2.tilt$sum=linMap(t2.tilt$sum, from=t1.lim[1], to=t1.lim[2])
  # t2.tilt1$sum=linMap(t2.tilt1$sum, from=t1.lim[1], to=t1.lim[2])
  t2.sum=rbind(t2.sum, t2.1d, t2.sky) #, t2.tilt, t2.tilt1)
  
  if(!is.null(addProjpp)){
    #browser()
    t2.projpp=data.frame(sum=addProjpp, ppm=t2.ppm, group='Feat proj')
    t2.projpp$sum=linMap(t2.projpp$sum, from=t1.lim[1], to=t1.lim[2])
    t2.sum=rbind(t2.sum, t2.projpp)
    
  }
  
  
  t2.sum$group=factor(t2.sum$group, levels=c('1D', 'Sum', 'Skyline', 'Feat proj')) #, 'T45 Sum', 'T45 Skyline'))
  
  # define nontilted or tilted data as main df for plotting
  # if(tilt==T){
  #   dsm=melt(out)
  #   colnames(dsm)=c('y', 'x' ,'z')
  # }else{
  ds$ppm2=hz2d
  dsm=melt(ds, id.vars = 'ppm2')
  dsm$ppm2=as.numeric(dsm$ppm2)
  dsm$variable=as.numeric(as.character(dsm$variable))
  colnames(dsm)=c('y', 'x' ,'z')
  #}
  
  #dsm$z[dsm$z<0]=0
  #dsm$z=dsm$z # scale intensities
  z.lim=quantile(dsm$z, probs=z.probs, na.rm = T) # define intensity cutoff
  dsm=dsm[which(dsm$z>z.lim),]
  
  if(!is.null(bbox)){
    bbox$P.peakOri[bbox$P.peakOri>1]=1
    
    # bbox$cent.f2=bbox$f2.ppm
    # bbox$cent.f1=bbox$f1.Hz
    # bbox$bb.width=
    
    bbox$bb.f2h=((bbox$cent.f2*SF)-bbox$bb.width.f2)/SF
    bbox$bb.f2l=((bbox$cent.f2*SF)+bbox$bb.width.f2)/SF
    bbox$bb.f1h=(bbox$cent.f1)-bbox$bb.width.f1
    bbox$bb.f1l=(bbox$cent.f1)+bbox$bb.width.f1
    
    #browser()
    bbox$idtext=paste0('p=',round(bbox$P.peakOri,2))
  }
  
  dsm$z=minmax(dsm$z)
  # create plots
  d2d=ggplot()+
    #geom_contour(data=subset(dsm1, z>z.lim), aes(x=x, y=y, z=z, colour=..level..),bins=50)+
    geom_tile(data=dsm, aes_string(x='x', y='y', fill='z'))+
    #geom_rect(data=bbox[order(bbox$P.peakOri),], aes_string(xmax='bb.f2l', xmin='bb.f2h', ymax='bb.f1h', ymin='bb.f1l', colour='P.peakOri', size='P.peakOri'), fill=NA)+
    scale_size(range=c(0.1,0.9), guide=F)+
    #scale_linetype_manual(values=c('TRUE'=1, 'FALSE'=9), guide=F)+
    #geom_text(data=bbox, aes_string(x='bb.f2h', y='bb.f1l', label='idtext'), hjust=0, size=3,  colour='black')+
    scale_x_reverse(breaks=pretty_breaks(), limits=rev(t2.lim))+
    scale_fill_gradientn(colours=matlab.like2(10))+
    scale_colour_gradient(low = "black", high = "white", limits=c(0,1))+
    #scale_y_continuous(limits=rev(t1.lim), breaks=pretty_breaks(), trans = 'reverse')+
    theme_bw()+
    theme(
      text= element_text(family='Helvetica'),
      legend.position=c(0.01,0.1),
      legend.justification=c(0,1),
      #legend.key.size = unit(0.7, units='cm'),
      legend.direction = 'horizontal',
      legend.box = 'horizontal',
      legend.box.background = element_rect(colour = "black"),
      legend.key = element_blank(),
      legend.margin=margin(0.1,0.1,0.1,0.1, unit= "lines"),
      legend.title=NULL,
      legend.text=element_text(size=8))+
    labs(x=expression(f[2]~~'(ppm)'),
         y=expression(f[1]~~'(Hz)'),
         fill=expression(Int['sc']),
         colour=expression('p'['SVM']))
  
  if(!is.null(bbox)){
    d2d=d2d+
      geom_rect(data=bbox[order(bbox$P.peakOri),], aes_string(xmax='bb.f2l', xmin='bb.f2h', ymax='bb.f1h', ymin='bb.f1l', colour='P.peakOri', size='P.peakOri'), fill=NA)+
      geom_text(data=bbox, aes_string(x='bb.f2h', y='bb.f1l', label='idtext'), hjust=0, size=3,  colour='black')
  }
  
  # t1 projections plot
  dv=ggplot(t1.sum, aes_string('ppm', 'sum', linetype='group', colour='group'))+
    geom_line(size=0.4)+
    scale_linetype_manual(values=c('1D'=1, 'Skyline'=8, 'Sum'=8), name=NULL)+
    scale_colour_manual(values=c('1D'='black', 'Skyline'='#D85140', 'Sum'='#5383EC'), name=NULL)+
    theme_bw()+
    theme(text= element_text(family='Helvetica'),
          axis.text = element_text(colour='transparent'),
          axis.ticks = element_line(colour='transparent'),
          legend.position=c(0.98,0.99),
          legend.justification=c(1,1),
          panel.grid.minor.x = element_blank(),
          legend.margin=margin(-1,0.1,0.1,0.1, unit= "lines"),
          legend.title=NULL,
          legend.spacing.y = unit(0, 'cm'),
          legend.text=element_text(size=7),
          legend.key.size = unit(1, 'lines'))+
    labs(x='', y=expression('minmax'))+
    scale_x_continuous(limits=t1.lim, breaks=pretty_breaks())+
    coord_flip()
  
  # t2 projections
  x1=ggplot(t2.sum, aes_string('ppm', 'sum', linetype='group', colour='group'))+
    geom_line(size=0.4)+
    scale_x_continuous(limits=rev(t2.lim), trans = 'reverse', breaks=pretty_breaks())+
    scale_linetype_manual(values=c('1D'=1, 'Skyline'=8, 'Sum'=8, 'T45 Skyline'=1, 'T45 Sum'=1, 'Feat proj'=1), name=NULL)+
    scale_colour_manual(values=c('1D'='black', 'Skyline'='#D85140', 'Sum'='#5383EC', 'T45 Skyline'='#F1BF42', 'T45 Sum'='#58A65C', 'Feat proj'='red'), name=NULL)+
    theme_bw()+
    theme(text= element_text(family='Helvetica'),
          axis.text = element_text(colour='transparent'),
          axis.ticks = element_line(colour='transparent'),
          legend.position=c(0.12,0.99),
          legend.justification=c(1,1),
          panel.grid.minor.y = element_blank(),
          legend.margin=margin(0, unit="lines"),
          legend.title=NULL,
          legend.text=element_text(size=7),
          legend.spacing.y = unit(0, 'cm'),
          legend.key.size = unit(1, 'lines'))+
    labs(x='', y=expression('minmax'))
  
  if(class(title)=='list'){
    title=paste(unlist(title), collapse ='\n')
  }
  
  
  empty <- ggplot()+
    geom_point(aes(0, 1), colour="white") +
    scale_x_continuous(expand=c(0,0), limits=c(0,1))+
    scale_y_continuous(expand=c(0,0), limits=c(0,1))+
    # theme(axis.line=element_blank(),
    #       axis.text.x=element_blank(),
    #       axis.text.y=element_blank(),
    #       axis.ticks=element_blank(),
    #       axis.title.x=element_blank(),
    #       axis.title.y=element_blank(),
    #       legend.position="none",
    #       panel.background=element_blank(),
    #       panel.border=element_blank(),
    #       panel.grid.major=element_blank(),
    #       panel.grid.minor=element_blank(),
  #       plot.background=element_blank())+
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )+
    #labs(title=title)+
    geom_text(aes(x=0, y=1, label=title), size=3, hjust=0, vjust=1)
  #arrange the plots together, with appropriate height and width for each row and column
  grid.arrange(x1,empty, d2d,dv, ncol=2, nrow=2, widths=c(4,1), heights=c(1.3, 4))
  
}
