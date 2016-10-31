# NEJM p-values
NEJM_pval <-function(x){ ifelse(x >= .001, signif(x,2) , 'P<0.001') }
#  In general, P values larger than 0.01 should be reported to two decimal places, 
# those between 0.01 and 0.001 to three decimal places; 
# P values smaller than 0.001 should be reported as P<0.001.

######
parseSummary <- function(summary,odsig=3,glm=FALSE,z=1.96){
  if(!is.data.frame(summary)){
    # get matrix
    summary = do.call(rbind,summary);			
    # order and data.frame
    if(glm){
      summary = data.frame(summary[order(rownames(summary),summary[,4]),])
    }else{
      summary = data.frame(summary[order(summary[,4]),])
    }
  }
  # odds and SE
  summary$'OR' = signif( exp(summary$Estimate) , odsig )
#  summary$'95 CI' = paste(   '(',   signif( exp(summary$Estimate - (z)*summary[,2]) ,odsig),' - ', signif( exp(summary$Estimate + (z)*summary[,2]) ,odsig), ')',sep='')
# exp(Estimate)  + 1.96 * exp(Estimate) * Std.err
  summary$'95 CI' = paste(   '(',   signif( exp(summary$Estimate)  - 1.96 * exp(summary$Estimate) * summary[,2] ,odsig),' - ', signif( exp(summary$Estimate)  + 1.96 * exp(summary$Estimate) * summary[,2] ,odsig), ')',sep='')
#  print(paste('95% CI =',exp((z)*summary[,2])))
  print(paste('95% CI = ',1.96 * exp(summary$Estimate) * summary[,2]))
  #pvalue adjust
  summary$'Pr(W)' = NEJM_pval(summary[,4])
  #rownames 
  rownames(summary) = gsub('_ug\\.mL', '' ,rownames(summary))
  # gen table 
  print(xtab <- xtable(summary[,5:7]),floating=FALSE)
  return(xtab)
}

#########
## aggregation functions

cumprob_f <- function(Prob,n){
  1-rollapply(1-Prob,2,function(x) prod(na.omit(x)),fill=NA)
}

geomeanprob_f <- function(Prob,n){
  rollapply(Prob,n,function(x) exp(mean(log(na.omit(x)))),fill=NA)
}

meanprob_f <- function(Prob,n){
  rollapply(Prob,n,function(x) mean(na.omit(x)),fill=NA)
}
meanprob_f_debug <- function(Prob,n){
  rollapply(Prob,n,function(x) print_n_func(mean,x),fill=NA)
}

print_n_func <- function(f,x){
  print(x)
  f_out <- f(x,na.rm = TRUE)
  print(f_out)
  f_out
}

library(gridExtra)
g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
}

#' cross_valid_resampling_sensitive
# cross validation method for sampling around subject resampling bias
cross_valid_resampling_sensitive <- function(data,formula,resp=all.vars(formula)[1],family,K=10,model=glm,resampled_var='ParticipantID',j=40,valid_prop=.2){
  K=K+1
  #ind = sample(1:nrow(data))
  #interval = floor(nrow(data)/K)
  #ranges = seq(1,nrow(data),interval)
  pred = list(); train = list(); labels = list()
  for(k in 1:(K-1)){ # k-th fold
      # iterate over j subjects without replacement
    ind = sapply( sample(unique(data[,resampled_var]),j), function(id){
      sample(which(data[,resampled_var]==id),1) }) # randomly choose a sample from each chosen subject

    # low = ranges[k] ; high = ranges[k+1]
    #valid_ind = ind[low:(high-1)]
    #train_ind = ind[!ind %in% valid_ind]
    
    #valid_max = round(length(ind)*valid_prop)
    #valid_ind = ind[1:valid_max]
    #train_ind = ind[(valid_max+1):length(ind)]
    #glm_i = model( formula , data=data[train_ind,] , family=family)
    #pred[[as.character(k)]] = predict(glm_i,newdata=data[,])
    pred[[as.character(k)]] = data[ind,all.vars(formula)[-1]]
    train[[as.character(k)]] = ind
    labels[[as.character(k)]] = data[ind,resp]
  }
  return(list(
    pred_out = do.call(cbind,pred),
    train_out = do.call(cbind,train),
    labels_out = do.call(cbind,labels)
  ))
}

####
# Get groupings
#' Returns case-control groupings 
#' 
#' @param ParticipantID A vector of participant IDs. Must alternate between cases and corresponding controls (e.g. case1,control1a,control1b...control1f,case2,control2a,control2b...)
#' @param cases A vector of cases.
#' @return A numeric vector denoting the group of each participant.
#' @examples
#' get_groupings(c('a','b','c','d'), c('a','c')) # -> c(1,1,2,2)
#' get_groupings(c('a1','a2','b','c','d1','d2'), c('a','c')) # -> c(1,1,1,2,2,2)
get_groupings <- function(ParticipantID,cases){
  #checks
  if(!all(cases %in% ParticipantID)) {warning(" 'ParticipantID' should contain all cases")}
  group=c()
  for(i in 1:length(ParticipantID)){ # for each participant ID
    if(ParticipantID[i] %in% cases){ # if ID is a case, begin a new group 
      group[i]=which(cases==ParticipantID[i])
    }else{
      group[i] = group[i-1]
    }
  }
  return(group)
}


clogit_screen<-function(){
  library(survival)
  summaryHMO <- list()
  variables = unique(c(grep('_ug',HMOs,value=T)))
  for(var in variables ){
    #### clogit Model for temporal HMOs and covariates
    form = as.formula( paste(c('NEC.binary ~ I(sqrt(DPP))+', 
                               ifelse(var %in% covariates,var,paste('scale(',var,')') ), '+ strata(matchedgroups)') ,collapse=' ') )
    print(form)
    NEC_data_tmp = unique(na.omit(NEC_data[,c('NEC.binary','DPP',var,'matchedgroups')])) # remove na and unique rows for var-specific data frame
    mod = clogit(form , data = NEC_data_tmp )
    summaryHMO[[var]] =  coef(summary(mod))[-1,-2] # HMOs
  }
  summaryHMO = data.frame(do.call(rbind,summaryHMO))
  colnames(summaryHMO)[1] = 'Estimate'
  summaryHMO = summaryHMO[order(summaryHMO[,4]),]
  summaryHMO.tab = parseSummary(summaryHMO);
}


########
# visualization

####
# Load Visualization Libraries
####

library(ggplot2)
library(directlabels)
library(wordcloud)
library(tm)

####
# Regression Volcano Plot
#' Volcano plot using wald probabilies and the corresponding beta statistics
#' 
#' @param pval A numeric vector containing p-values
#' @param beta A numeric vector containing corresponding beta values
#' @param name A character vector indicating the filename
#' @param thresh A number indicating the theshold for including lables
#' @param plot.directory, A directory, indicating where the plot should be saved.  If NULL, no plot is saved.
#' @examples
#' volcano.reg(c(.4,.04),c(.1,5))
volcano.reg<-function(pval,beta,nam,rnames,thresh=5,plot.directory=file.path(getwd(),'volcano.reg/'),exclude.outlier=TRUE){
  gees.df = data.frame(Beta=beta,Pr_wald=pval)  
  print(gees.df)
  rownames(gees.df) = rnames  # adding anova and 
  gees.df$q_wald     = p.adjust(gees.df$Pr_wald,'fdr')
  gees.df$negLog_10_wald = -log(gees.df$q_wald+1e-20,10)
# save(beta,pval,rnames,gees.df,file='tmp.rda')
  gees.df$HMO = rownames(gees.df)
  #qplot(x=Beta,y=negLog_10_wald,data=gees.df)

  label.dat = gees.df[gees.df$negLog_10_wald>=thresh,]

  gees.df = gees.df[gees.df$Beta<10,] # remove extreme values


  if(exclude.outlier){
    gees.df = gees.df[abs(gees.df$Beta)<(2*sd(gees.df$Beta)),]
  }

  #
# with(gees.df, textplot(Beta,negLog_10_wald,ifelse(negLog_10_wald>thresh,HMO,'')))
  #
# gees.df$HMO_lab <- ifelse(gees.df$negLog_10_wald>thresh,gees.df$HMO,'')
# direct.label(xyplot(Beta~negLog_10_wald,gees.df,groups=HMO_lab, col="black"))
  #
  g<-ggplot(gees.df, aes(Beta,negLog_10_wald)) + geom_point(color=gees.df$negLog_10_wald) 
    #geom_text(data = label.dat, aes(Beta,negLog_10_wald, label = HMO,color=negLog_10_wald), size = 2)
   # geom_dl(aes(label=HMO, colour=negLog_10_wald), list('last.bumpup', cex = 1.3, hjust = 1))
  #direct.label(g, list("last.points", cex=.7, hjust=1))
   # theme(axis.title.x=element_text('-log(Pr(wald))'))
    if(!is.null(plot.directory)){
    dir.create(plot.directory)
    ggsave(file.path(plot.directory,nam),plot=g)
    dev.off()
    write.csv(gees.df[order(gees.df$Pr_wald,-abs(gees.df$Beta)),],quote=FALSE,file=paste(file.path(plot.directory,nam),'.csv',sep=''))
  }
  return(gees.df)
}

####
# Spagetti Plot
#' Line plot for longitudinal data.
#' 
#' @param var.string A string refering to the column name of \code{var} in data.frame \code{dat}.
#' @param dat A data.frame containing a time variable (DPP), reponse variable (NEC).
#' @param scale A boolean indicating if the variable indicated by var.string should be scaled.
#' @param plot.directory, A directory, indicating where the plot should be saved.  If NULL, no plot is saved.
#' @examples
#' spagetify( 'var1',dat)
spagetify<-function(var.str,dat,scale=TRUE,plot.directory=file.path(getwd(),'spagetti/'),ymin=ifelse(scale,-3,0),ymax=ifelse(scale,3,1000)){
  if(!is.null(plot.directory)){
    dir.create(plot.directory)
    nam = paste(c('HMO',ifelse(scale,'scaled','raw'),var.str,'pdf') , collapse='.')
    pdf(file.path(plot.directory,nam))
  }
  id <- unique(dat$ParticipantID)
  plot(x=c(0, 30), y= c(ymin, ymax ) )
  for (ii in 1:length(id)){
      tmp <- subset(dat, ParticipantID==id[ii], c("DPP",var.str, "NEC.binary"))
      if(scale){
        tmp[,2] = scale(tmp[,2])
      }
      lines(tmp[, "DPP"], tmp[,var.str], col=tmp[,"NEC.binary"] + 1, type="l", pch=16)
  }
  if(!is.null(plot.directory)){
    dev.off()
  }
}

####
# GGPlot2 Lasagna Plot
#' Heatmap for longitudinal data.
#' 
#' @param var A numeric vector from the data.frame \code{dat}.
#' @param var.string A string refering to the column name of \code{var} in data.frame \code{dat}.
#' @param dat A data.frame containing a time variable (DPP), reponse variable (NEC).
#' @param smooth An integer, number of adjacent days to consider in generated mean and standard deviation.
#' @param maxCNT, An integer, the maximum repsonse variable to be considered Control.
#' @param plot.pval, A boolean, inclusion of z-test pvalue describing the change relative to the mean and standard deviation produced by the \code{smooth} window.
#' @param plot.fill, A string, indicating the fill of the heatmap: 'logfoldchange','foldChange'
#' @param plot.directory, A directory, indicating where the plot should be saved.  If NULL, no plot is saved.
#' @return A ggplot2 plot object.
#' @examples
#' lasagnify(dat$var1, 'var1',dat)
lasagnify<-function(var,var.str,dat,smooth=1,maxCNT=0,plot.pval=TRUE,plot.fill='logfoldChange',plot.directory=file.path(getwd(),'lasagnia/'),param=NULL){
  if(length(var) != nrow(dat)){stop("v")}
  # find the groupwise aggregated median
  datprep = data.frame(var=var,DPP=dat$DPP,MatchGroup=dat$matchedgroups,NEC=dat$NEC)
  d = aggregate(datprep$var,by=list(DPP=datprep$DPP,MatchGroup=datprep$MatchGroup,NEC=datprep$NEC),FUN=median,na.rm=TRUE)   # aggregates sampes, messes up standard deviation
  caseNEC=c( 2.2  , 2.1  , 1.2  , 3.2  , 3.1  , 2.4  , 2.3  , 3.3  , 1.1 , 2.5   )
  cases = c('A003','A029','A066','B032','C005','C024','C027','E002','D015','D024')
  d$NECgroups = paste(cases[d$MatchGroup],'- Bell: ',d$NEC)
  d$NECgroups2 = paste(paste('Stage:',caseNEC[d$MatchGroup]), cases[d$MatchGroup] )####
# d$NECgroups = paste("Group:",letters[d$MatchGroup],'- NEC:',d$NEC)

  d_m = c()
  d_s = c()
  d_md= c()
  Necdays = c()
  for(i in 1:nrow(d)){
    DPP_cnt = which( d$MatchGroup[i]  == dat$matchedgroups
                   & d$DPP[i]-smooth  <= dat$DPP
                   & d$DPP[i]+smooth  >= dat$DPP
                   & maxCNT     >= as.numeric(dat$NEC)    ) # find matching data in orgininal (unaggregated dat matrix)
    #DPP_cnt = which.min(  dat$DPP[DPP_cnt_pre])
    d_md[i] = median(na.omit(var[DPP_cnt]))
    d_m[i]= mean(na.omit(var[DPP_cnt]))
    d_s[i] = sd(na.omit(var[DPP_cnt]))
#   print(i)
#   print(DPP_cnt)
#   print(var[DPP_cnt])
  }

  #a<<-data.frame(d_s,d_m,d_md,d$x,d$NECgroups2,d$NECgroups,d$DPP)
  #print(a[a[,5]=='D015'])
# d_s = ifelse(  is.na(d_s) , mean(na.omit(d_s)) , d_s)
# d_m = ifelse(  is.na(d_m) , mean(na.omit(d_m)) , d_m)
# d_md = ifelse( is.na(d_md), mean(na.omit(d_md)), d_md)

  d$median = d_m  
  d$stdev  = d_s
  d$mean   = d_md
  #d$foldChange = d$x /  d$median
  d$foldChange = d$x /  d$mean
  d$logfoldChange = log(d$foldChange)
  d$z = (d$x - d$mean)/ d$stdev
  d$pValue_z =  2*pnorm(-abs(d$z))

  d = na.omit(d)

  NECdays=data.frame(DPP=1:4,NECgroups=1:8,bin=1)
  d$DPP[d$NEC==0] = -d$DPP[d$NEC==0]

  if(plot.pval){
    if(plot.fill=='logfoldChange'){
      plt <- ggplot(d, aes(x=DPP,y=NECgroups2,fill=logfoldChange,alpha=1-(pValue_z/2)))
      maxnum = max(abs(d$logfoldChange))
    }else if(plot.fill=='foldChange'){
      plt <- ggplot(d, aes(x=DPP,y=NECgroups2,fill=foldChange,alpha=1-(pValue_z/2)))
      maxnum = max(abs(d$foldChange))
    }
  }else{
    if(plot.fill=='logfoldChange'){
      plt <- ggplot(d, aes(x=DPP,y=NECgroups2,fill=logfoldChange))
      maxnum = max(abs(d$logfoldChange))
    }else if(plot.fill=='foldChange'){
      plt <- ggplot(d, aes(x=DPP,y=NECgroups2,fill=foldChange))
      maxnum = max(abs(d$foldChange))
    }
  }
# plt +  theme(panel.background = element_rect(fill = "grey"),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+geom_tile() + theme(axis.text.y=element_text(colour=c('blue','red')) ) + scale_fill_gradient2(low="red", mid = "white",high="blue" , breaks=seq(-4,2,.5),guide = guide_legend(title = paste('Log(fold change)\n',var.str) ))
# plt +  theme(panel.background = element_rect(fill = "white"),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+geom_tile() + theme(axis.text.y=element_text(colour=c('darkgreen','darkred')) ) + scale_fill_gradientn(colours=rev(c("darkgreen", "green", "yellow", "red", "darkred")) , breaks=seq(-4,2,.5),guide = guide_legend(title = paste('Log(fold change)\n',var.str) ))

  if(maxnum<5){maxnum = 2.5}

# plt <- ggplot(d, aes(x=DPP,y=NECgroups2,fill=logfoldChange))
  plt <- plt +  theme(panel.background = element_rect(fill = "grey"),panel.grid.major = element_blank(),panel.grid.minor = element_blank() )+
    geom_tile() + #theme(axis.text.y=element_text(colour=c('darkblue','white','red', 'darkred')[d$NEC]) ) +
    theme(axis.text.y= element_text(size = rel(.7)) ) + 
    scale_fill_gradientn(colours=rev(c("darkblue", "blue", "white", "red", "darkred")) , 
      breaks=seq(round(-maxnum-.5), round(maxnum+.5),1), 
      limits=c(-maxnum-.1, maxnum+.1), ##### fix scale
      #values=seq(-2,2,length=10), rescaler = function(x, ...) x, oob = identity,
      guide = 'colourbar' )+ #guide_legend(title = paste('Log(fold change)\n',var.str) )) +
    geom_tile() #+
#   geom_tile(onset_long, aes(x=DPP,y=ParticipantID,fill=NEC))
# plt

  ####
  # NEC onset
  NEC_onset = read.csv('~/Documents/Grad/Labs/Lewis/NEC_HMO/data/NEC_onset.csv')
  onset_long = as.data.frame(do.call(rbind,apply(NEC_onset,1,function(x) cbind(x[1],cbind(x[2],x[3]:x[4])))))
  colnames(onset_long)=c('ParticipantID','NEC', 'DPP')
  onset_long$DPP = as.numeric(onset_long$DPP)
# plt2 <- ggplot(onset_long, aes(x=DPP,y=ParticipantID,fill=NEC))

  if(!is.null(plot.directory)){
    dir.create(plot.directory)
    nam = paste(c('HMO',var.str,'smooth',smooth,'maxCNT',maxCNT,ifelse(plot.pval,'pval','withoutPval'),param,'eps') , collapse='.')
    ggsave(file.path(plot.directory,nam),plot=plt,height=2,width=7)
  }
  return(plt)
}








