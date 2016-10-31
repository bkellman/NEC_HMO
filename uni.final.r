######################################################
# 
#  LOAD DATA
#
######################################################

####
# Load Regression Libraries
library(glmmLasso)  # feature selection, linear
library(lme4)   # glmer, linear
library(geepack)  # geeglm, non-parametric
####
# Load Evaluation Libraries
library(pROC)
library(ROCR)
####

wd = '~/Documents/Grad/Labs/Lewis/NEC_HMO'
setwd(wd)
####
# Load Support Functions
source('uni.final.func.r')
####
# Load saved data
load('data/NEC_data_update.rda') # dat
NEC_data = dat
onset = read.csv('data/NEC_onset.csv')

setwd('odds')

##########################
## Set parameters for run
NECclass=2:3
CNTclass=0:1
#removeSecretorColumn=TRUE
#all = FALSE
#makePlots = FALSE

covarience = 'exchangeable'
##########################

####
# Binarize NEC classification
NEC_data$NEC.binary<- ifelse(NEC_data$NEC <= max(CNTclass), 0, ifelse(NEC_data$NEC >= min(NECclass) , 1 , NA )) 

####
# Initialize Time transforms
NEC_data$DPP.sqrt <-   sqrt(NEC_data$DPP)           # Square-root of the time variable ______
NEC_data$DPP.std  <-  scale(NEC_data$DPP)           # Scaled time variable _____
NEC_data$offset   <- ifelse(NEC_data$NEC.binary==0, 1, 5)   # Offset response helps account for the difference in the number of case and control subjects [1]

####
# Denote case-control matchings and covariates
cases = c('A003','A029' ,'A066' ,'B032' ,'C005' ,'C024' ,'C027' ,'E002' ,'D015' ,'D024' )
caseNEC=c( 2    , 2     , 1     , 3     , 3     , 2     , 2     , 3     , 1     , 2     )
case_df = data.frame(cases,caseNEC)
controls=unique(NEC_data$ParticipantID[which(!(NEC_data$ParticipantID%in%cases))])
NEC_data$matchedgroups = get_groupings(NEC_data$ParticipantID,cases)
covariates <- c("GA_calc","Birthweight","Gender","Delivery","Race_Ethnic","Location","HMO_Secretor","Diversity","Evenness", "Sum_ug.mL","SUM_nmol.mL")
external <- c("ParticipantID","DOB","GA_weeks","GA_days","NEC","SampleID" ,"DPP","NEC.binary","DPP.sqrt","DPP.std","offset","matchedgroups" )
HMOs <- colnames(NEC_data)[!(colnames(NEC_data) %in% c(covariates,external))]
temporal <- c(HMOs,"Diversity","Evenness", "Sum_ug.mL","SUM_nmol.mL")


gee0     <- geeglm( NEC.binary ~ Birthweight, id = ParticipantID, 
               data = NEC_data, family=binomial(link="logit"), 
                corstr = covarience, std.err="san.se") 

#######################
### Cohort Description
library(xtable)
library(Hmisc)
summary(NEC.binary ~ GA_calc + Birthweight + Gender + Delivery + Race_Ethnic + Diversity + Evenness + Sum_ug.mL, data=NEC_data)

######################################################
# 
#  INITIAL ASSESSMENT OF DSLNT AND SUM
#  FIGURE 1
#
####################################################### 
attach(NEC_data)

NEC_fact1 <- factor(ifelse(NEC==0,'CNT',ifelse(NEC==1,'Partial',ifelse(NEC==2|NEC==3,'NEC',NA))))

# t-test Normal 2-mean compare Normal
#t.test( DSLNT_ug.mL[NEC.binary==0] ,DSLNT_ug.mL[NEC.binary==1] ) # not valid because DSLNT is not normal, see shapiro.test
#t.test( SUM_nmol.mL[NEC.binary==0] ,SUM_nmol.mL[NEC.binary==1] ) # not valid because SUM is not normal, see shapiro.test
# mann-whitney-wilcoxon 2-mean compare Nonparametric 
wilcox.test( DSLNT_ug.mL[NEC==0] , DSLNT_ug.mL[NEC==1]  )           # reject equality, W = 6000, p-value = 0.004
wilcox.test( DSLNT_ug.mL[NEC==0] , DSLNT_ug.mL[NEC==2 | NEC==3]  )  # reject equality, W = 33000, p-value <2e-16
wilcox.test( DSLNT_ug.mL[NEC==1] , DSLNT_ug.mL[NEC==2 | NEC==3]  )  # reject equality, W = 920, p-value = 1e-08
wilcox.test( SUM_nmol.mL[NEC==0] , SUM_nmol.mL[NEC==1]  )           # accept equality, W = 4300, p-value = 0.8
wilcox.test( SUM_nmol.mL[NEC==0] , SUM_nmol.mL[NEC==2 | NEC==3]  )  # accept equality, W = 18000, p-value = 1
wilcox.test( SUM_nmol.mL[NEC==1] , SUM_nmol.mL[NEC==2 | NEC==3]  )  # accept equality, W = 400, p-value = 0.4
# kruskal-wallis n-mean compare Normal
kruskal.test( DSLNT_ug.mL ~ NEC_fact1  )  # reject equality, chi-squared = 130, df = 2, p-value <2e-16
kruskal.test(  SUM_nmol.mL ~ NEC_fact1 )  # accept equality, chi-squared = 0.082, df = 2, p-value = 1
# shapiro-wilk Normality
shapiro.test( DSLNT_ug.mL )   # reject normality, W = 0.92, p-value <2e-16
shapiro.test( SUM_nmol.mL )   # reject normality, W = 0.98, p-value = 1e-06

detach(NEC_data)

######################################################
# 
#  UNIVARIATE TEMPORAL GEES AND
#  UNIVARIATE GLMS
#  FIGURE 3A, 3B
#
####################################################### 
summaryHMO <- list()
summaryCov <- list()
summaryCov.tmp <- list()
variables = unique(c(grep('_ug',HMOs,value=T),covariates))
model_save <- list()

for(var in variables ){
	#### GEE Model for temporal HMOs and covariates
	if(var %in% temporal){  
		form = as.formula( paste(c('NEC.binary ~ I(sqrt(DPP))+', 
			ifelse(var %in% covariates,var,paste('scale(',var,')') ) ) , collapse=' ') )
		print(form)
		NEC_data_tmp = unique(na.omit(NEC_data[,c('NEC.binary','DPP',var,'ParticipantID')])) # remove na and unique rows for var-specific data frame
		gee     <- geeglm( form , id = ParticipantID, 
		               data = NEC_data_tmp, family=binomial(link="logit"), 
		               corstr = covarience, std.err="san.se")
    model_save[[var]] <- gee
		if(var %in% covariates){
			summaryCov.tmp[[var]] =  coef(summary(gee))[-c(1,2),] # diversity, evenness, sum
		}else{
			summaryHMO[[var]] =  coef(summary(gee))[-c(1,2),] # HMOs
		}
	}else{
	#### GLM Model for static covariates
		form = as.formula( paste(c('NEC.binary ~', var ) , collapse=' ') )
		print(form)
		NEC_data_tmp = unique(na.omit(NEC_data[,c('NEC.binary',var,'ParticipantID')])) # remove na and unique rows for var-specific data frame
		glm_i  <- glm( form ,  data = NEC_data_tmp, family=binomial(link="logit"))
		summaryCov[[var]] =  coef(summary(glm_i))[-1,]
		model_save[[var]] <- glm_i
		print (anova(glm_i,test='Chisq')$'Pr(>Chi)'[2]) # Chisq == LRT == Likelihood-Ratio test
	}
}
#summaryHMO_matched = summaryHMO

# calculate odds, adjust 
library(xtable)

#summaryHMO.tab = parseSummary(summaryHMO);	
#summaryCov.tmp.tab = parseSummary(summaryCov.tmp);			
#summaryCov.tab = parseSummary(summaryCov,glm=TRUE); 


######################################################
# 
#  FINAL MODEL
#  FIGURE 3C, TABLE S1, FIGURE S2
#
####################################################### 
library(MuMIn)

  gee0 = geeglm( NEC.binary ~ 1   , 
                  id = ParticipantID, data = NEC_data, family=binomial(link="logit"), 
                  corstr = covarience, std.err="san.se")
  gee0t = geeglm( NEC.binary ~ I(sqrt(DPP))   , 
                   id = ParticipantID, data = NEC_data, family=binomial(link="logit"), 
                   corstr = covarience, std.err="san.se") 

geeS <- geeglm( NEC.binary ~ I(sqrt(DPP)) + scale(DSLNT_ug.mL)   , 
	id = ParticipantID, data = NEC_data, family=binomial(link="logit"), 
	corstr = covarience, std.err="san.se") 

geeSN <- geeglm( NEC.binary ~ I(sqrt(DPP)) + scale(DSLNT_ug.mL) + scale(LNFP_I_ug.mL)   , 
	id = ParticipantID, data = NEC_data, family=binomial(link="logit"), 
	corstr = covarience, std.err="san.se")

geeSF <- geeglm( NEC.binary ~ I(sqrt(DPP)) + scale(DSLNT_ug.mL) + scale(DFLNT_ug.mL)   , 
	id = ParticipantID, data = NEC_data, family=binomial(link="logit"), 
	corstr = covarience, std.err="san.se") 
# Warning message:
# glm.fit: fitted probabilities numerically 0 or 1 occurred

geeSNF <- geeglm( NEC.binary ~ I(sqrt(DPP)) + scale(DSLNT_ug.mL) + scale(LNFP_I_ug.mL) + scale(DFLNT_ug.mL)   , 
	id = ParticipantID, data = NEC_data, family=binomial(link="logit"), 
	corstr = covarience, std.err="san.se")
# Warning message:
# glm.fit: fitted probabilities numerically 0 or 1 occurred

geeSN1FN3 <- geeglm( NEC.binary ~ I(sqrt(DPP)) + scale(DSLNT_ug.mL) + scale(LNFP_I_ug.mL) + scale(DFLNT_ug.mL) + scale(LNFP_III_ug.mL)   , 
                  id = ParticipantID, data = NEC_data, family=binomial(link="logit"), 
                  corstr = covarience, std.err="san.se")
# Warning message:
# glm.fit: fitted probabilities numerically 0 or 1 occurred

#model.sel(geeSN,geeSF,geeFN,geeS,geeSNF,geeSN1FN3,gee0,gee0t,rank=QIC)
model.sel(geeSN,geeSF,geeS,geeSNF,geeSN1FN3,gee0,gee0t,rank=QIC)

######################################################
# 
#  FINAL MODEL Selection: Backward Elimination
#  Table S3
#
####################################################### 

var4to3 <- list(
geeSN1FN3 = geeglm( NEC.binary ~ I(sqrt(DPP)) + scale(DSLNT_ug.mL) + scale(LNFP_I_ug.mL) + scale(DFLNT_ug.mL) + scale(LNFP_III_ug.mL)   , 
                     id = ParticipantID, data = NEC_data, family=binomial(link="logit"), 
                     corstr = covarience, std.err="san.se"),
geeN1FN3 = geeglm( NEC.binary ~ I(sqrt(DPP))  + scale(LNFP_I_ug.mL) + scale(DFLNT_ug.mL) + scale(LNFP_III_ug.mL)   , 
                     id = ParticipantID, data = NEC_data, family=binomial(link="logit"), 
                     corstr = covarience, std.err="san.se"),
geeSFN3 = geeglm( NEC.binary ~ I(sqrt(DPP)) + scale(DSLNT_ug.mL) + scale(DFLNT_ug.mL) + scale(LNFP_III_ug.mL)   , 
                     id = ParticipantID, data = NEC_data, family=binomial(link="logit"), 
                     corstr = covarience, std.err="san.se"),
geeSN1N3 = geeglm( NEC.binary ~ I(sqrt(DPP)) + scale(DSLNT_ug.mL) + scale(LNFP_I_ug.mL)  + scale(LNFP_III_ug.mL)   , 
                     id = ParticipantID, data = NEC_data, family=binomial(link="logit"), 
                     corstr = covarience, std.err="san.se"),
geeSN1F = geeglm( NEC.binary ~ I(sqrt(DPP)) + scale(DSLNT_ug.mL) + scale(LNFP_I_ug.mL) + scale(DFLNT_ug.mL)    , 
                     id = ParticipantID, data = NEC_data, family=binomial(link="logit"), 
                     corstr = covarience, std.err="san.se")
)
ms <- model.sel(var4to3 ,rank=QIC)
xtable(ms[,8:11])
#plot(ms,mar=c(5,5,9,9))

var3to2 <- list(
  geeSN1F = geeglm( NEC.binary ~ I(sqrt(DPP)) + scale(DSLNT_ug.mL) + scale(LNFP_I_ug.mL) + scale(DFLNT_ug.mL)    , 
                    id = ParticipantID, data = NEC_data, family=binomial(link="logit"), 
                    corstr = covarience, std.err="san.se"),
  geeN1F = geeglm( NEC.binary ~ I(sqrt(DPP)) + scale(LNFP_I_ug.mL) + scale(DFLNT_ug.mL)    , 
                    id = ParticipantID, data = NEC_data, family=binomial(link="logit"), 
                    corstr = covarience, std.err="san.se"),
  geeSF = geeglm( NEC.binary ~ I(sqrt(DPP)) + scale(DSLNT_ug.mL)  + scale(DFLNT_ug.mL)    , 
                    id = ParticipantID, data = NEC_data, family=binomial(link="logit"), 
                    corstr = covarience, std.err="san.se"),
  geeSN1 = geeglm( NEC.binary ~ I(sqrt(DPP)) + scale(DSLNT_ug.mL) + scale(LNFP_I_ug.mL)     , 
                    id = ParticipantID, data = NEC_data, family=binomial(link="logit"), 
                    corstr = covarience, std.err="san.se")
  )
ms <- model.sel(var3to2 ,rank=QIC)
xtable(ms[,7:10])
#plot(ms,mar=c(5,5,9,9))

var2to1to0 <- list(
  geeSF = geeglm( NEC.binary ~ I(sqrt(DPP)) + scale(DSLNT_ug.mL)  + scale(DFLNT_ug.mL)    , 
                  id = ParticipantID, data = NEC_data, family=binomial(link="logit"), 
                  corstr = covarience, std.err="san.se"),
  geeS = geeglm( NEC.binary ~ I(sqrt(DPP)) + scale(DSLNT_ug.mL)      , 
                  id = ParticipantID, data = NEC_data, family=binomial(link="logit"), 
                  corstr = covarience, std.err="san.se"),
  geeF = geeglm( NEC.binary ~ I(sqrt(DPP)) +  scale(DFLNT_ug.mL)    , 
                  id = ParticipantID, data = NEC_data, family=binomial(link="logit"), 
                  corstr = covarience, std.err="san.se"),
  gee0 = geeglm( NEC.binary ~ 1   , 
                  id = ParticipantID, data = NEC_data, family=binomial(link="logit"), 
                  corstr = covarience, std.err="san.se"),
  gee0t = geeglm( NEC.binary ~ I(sqrt(DPP))   , 
                   id = ParticipantID, data = NEC_data, family=binomial(link="logit"), 
                   corstr = covarience, std.err="san.se") 
  )
ms<-model.sel(var2to1to0 ,rank=QIC)
xtable(ms[,6:9])
#plot(ms,mar=c(5,5,9,9))

##
ms<-model.sel(c(var4to3,var3to2[-1],var2to1to0[-1]),rank=QIC)
xtable(ms[,8:11])
#plot(ms,mar=c(5,5,9,9))

model.sel(geeSNF,geeS,gee0,rank=QIC)


#library(survival)
#clogit(NEC.binary ~ I(sqrt(DPP)) + scale(DSLNT_ug.mL)  + scale(DFLNT_ug.mL) + scale(LNFP_I_ug.mL) + strata(matchedgroups), data=NEC_data)
#clogit(NEC.binary ~ I(sqrt(DPP)) +  scale(LNFP_II_ug.mL) + strata(matchedgroups), data=NEC_data)

######################################################
#
# Final Model Odds Visual
# FIGURE S4A, S4B
#
######################################################

require(dplyr)
require(ggplot2)
library(zoo)

# calculate odds, probability and cumulative probability at each time point
NEC_data_tmp = na.omit(NEC_data[,c('NEC.binary','NEC','DPP','DSLNT_ug.mL','LNFP_I_ug.mL','DFLNT_ug.mL','ParticipantID','DSLNT_percent','LNFP_I_percent','DFLNT_percent')])

# visualize each subject seperately
ggplot(NEC_data_tmp[NEC_data_tmp$NEC.binary>0,], aes(x=DPP,y=DFLNT_ug.mL,color=NEC)) + geom_point() + facet_grid(. ~ ParticipantID ) + geom_smooth(method = 'glm')
ggplot(NEC_data_tmp[NEC_data_tmp$NEC.binary>0,], aes(x=DPP,y=DSLNT_ug.mL,color=NEC)) + geom_point() + facet_grid(. ~ ParticipantID ) + geom_smooth(method = 'glm')
ggplot(NEC_data_tmp[NEC_data_tmp$NEC.binary>0,], aes(x=DPP,y=LNFP_I_ug.mL,color=NEC)) + geom_point() + facet_grid(. ~ ParticipantID ) + geom_smooth(method = 'glm')

use_final=FALSE
if(use_final){
  NEC_data_tmp$Odds1 = with(NEC_data_tmp, 
                          exp(  apply( cbind(1,sqrt(DPP) ,scale(DSLNT_ug.mL) , scale(LNFP_I_ug.mL) , scale(DFLNT_ug.mL)) %*% diag(coef(geeSNF)),1,sum ) ) )
}else{
  NEC_data_tmp$Odds1 = with(NEC_data_tmp, 
                          exp(  apply( cbind(1,sqrt(DPP) ,scale(DSLNT_ug.mL) ) %*% diag(coef(geeS)),1,sum ) ) )
}
NEC_data_tmp$prob1 = with(NEC_data_tmp, Odds1/(1 + Odds1) )

# local cumulative probabiliity
#NEC_data_tmp = NEC_data_tmp %>% group_by(ParticipantID) %>% arrange(DPP) %>% mutate(cumProb2 = 1-rollapply(1-prob1,2,function(x) prod(na.omit(x)),fill=NA))
#NEC_data_tmp = NEC_data_tmp %>% group_by(ParticipantID) %>% arrange(DPP) %>% mutate(cumProb3 = 1-rollapply(1-prob1,3,function(x) prod(na.omit(x)),fill=NA))
#NEC_data_tmp = NEC_data_tmp %>% group_by(ParticipantID) %>% arrange(DPP) %>% mutate(cumProb4 = 1-rollapply(1-prob1,4,function(x) prod(na.omit(x)),fill=NA))

# geometric mean of probability
NEC_data_tmp = NEC_data_tmp %>% group_by(ParticipantID) %>% arrange(DPP) %>% mutate(cumProb2 = rollapply(prob1,2,function(x) exp(mean(log(na.omit(x)))),fill=NA))
NEC_data_tmp = NEC_data_tmp %>% group_by(ParticipantID) %>% arrange(DPP) %>% mutate(cumProb3 = rollapply(prob1,3,function(x) exp(mean(log(na.omit(x)))),fill=NA))
NEC_data_tmp = NEC_data_tmp %>% group_by(ParticipantID) %>% arrange(DPP) %>% mutate(cumProb4 = rollapply(prob1,4,function(x) exp(mean(log(na.omit(x)))),fill=NA))

# joint probability
#NEC_data_tmp = NEC_data_tmp %>% group_by(ParticipantID) %>% arrange(DPP) %>% mutate(cumProb2 = rollapply(prob1,2,function(x) prod(na.omit(x)),fill=NA))
#NEC_data_tmp = NEC_data_tmp %>% group_by(ParticipantID) %>% arrange(DPP) %>% mutate(cumProb3 = rollapply(prob1,3,function(x) prod(na.omit(x)),fill=NA))
#NEC_data_tmp = NEC_data_tmp %>% group_by(ParticipantID) %>% arrange(DPP) %>% mutate(cumProb4 = rollapply(prob1,4,function(x) prod(na.omit(x)),fill=NA))

view.col=TRUE
if(view.col){
	col_fn <- scale_colour_hue( guide=guide_legend(title="Diagnosis"))
}else{
	col_fn <- scale_colour_manual(values=c('grey', 'black'), guide=guide_legend(title="Diagnosis"))
}
#GEE1
# spagetti plot odds
p1.1 <- ggplot(data = NEC_data_tmp, aes(x=DPP, y=Odds1, group = ParticipantID, color=ifelse(NEC.binary,'NEC','Control' ) )) +
	geom_line()+ xlab('')+ ylab('Odds') + 
	theme_bw() + col_fn
# NEC boxplots
p1.2 <- ggplot(data = NEC_data_tmp, aes(x=NEC, y=Odds1, group=NEC, color=ifelse(NEC.binary,'NEC','Control' )  )) +
	geom_boxplot() + xlab('Bell Stage') + ggtitle("Final Model 1 (equation 1)\n.") + ylab('Odds') + 
	theme_bw() + col_fn



# Cumulative probability of NEC onset
NEC_data_tmp = NEC_data_tmp %>% group_by(ParticipantID) %>% mutate(subjMean1 = median(na.omit(DSLNT_ug.mL)))
NEC_data_tmp = NEC_data_tmp %>% group_by(ParticipantID) %>% mutate(subjsd1 = sd(na.omit(DSLNT_ug.mL)))
order.df = unique(data.frame(NEC_data_tmp$ParticipantID,NEC_data_tmp$subjMean1))
NEC_data_tmp$ParticipantID = factor(NEC_data_tmp$ParticipantID , levels= order.df[order(order.df[,2]),1] )
# FIGURE S4A
conc1 <- ggplot(data = NEC_data_tmp, aes( y=DSLNT_ug.mL, x = ParticipantID, color=ifelse(NEC.binary,'NEC','Control' ) ))+
	geom_boxplot(outlier.size=1) + col_fn + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
	ggtitle("DSLNT of NEC given 1 Milk Sample") 
#	stat_summary(fun.y=sd, geom="line", aes(group=1),color='black')

NEC_data_tmp = NEC_data_tmp %>% group_by(ParticipantID) %>% mutate(subjMean1 = median(na.omit(prob1)))
order.df = unique(data.frame(NEC_data_tmp$ParticipantID,NEC_data_tmp$subjMean1))
NEC_data_tmp$ParticipantID = factor(NEC_data_tmp$ParticipantID , levels= order.df[order(order.df[,2]),1] )
cum1 <- ggplot(data = NEC_data_tmp, aes( y=prob1, x = ParticipantID, color=ifelse(NEC.binary,'NEC','Control' ) ))+
	geom_boxplot(outlier.size=1) + col_fn + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
	ggtitle("Probability of NEC given 1 Milk Sample")  
#	geom_hline(yintercept = 0.125) + 
#	stat_summary(fun.y=sd, geom="line", aes(group=1),color='black')

NEC_data_tmp = NEC_data_tmp %>% group_by(ParticipantID) %>% mutate(subjMean2 = median(na.omit(cumProb2)))
order.df = unique(data.frame(NEC_data_tmp$ParticipantID,NEC_data_tmp$subjMean2))
NEC_data_tmp$ParticipantID = factor(NEC_data_tmp$ParticipantID , levels= order.df[order(order.df[,2]),1] )
cum2 <- ggplot(data = NEC_data_tmp, aes( y=cumProb2, x = ParticipantID, color=ifelse(NEC.binary,'NEC','Control' ) ))+
	geom_boxplot(outlier.size=1) + col_fn + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
	ggtitle("Cumulative Probability of NEC given 2 Milk Samples") 
#	geom_hline(yintercept = 0.125) +
#	stat_summary(fun.y=sd, geom="line", aes(group=1),color='black')


NEC_data_tmp = NEC_data_tmp %>% group_by(ParticipantID) %>% mutate(subjMean4 = median(na.omit(cumProb4)))
order.df = unique(data.frame(NEC_data_tmp$ParticipantID,NEC_data_tmp$subjMean4))
NEC_data_tmp$ParticipantID = factor(NEC_data_tmp$ParticipantID , levels= order.df[order(order.df[,2]),1] )
cum4 <- ggplot(data = NEC_data_tmp, aes( y=cumProb4, x = ParticipantID, color=ifelse(NEC.binary,'NEC','Control' ) ))+
	geom_boxplot(outlier.size=1) + col_fn + theme(axis.text.x = element_text(angle = 90, hjust = 1))+
	ggtitle("Cumulative Probability of NEC given 4 Milk Samples") 
#	geom_hline(yintercept = 0.125) +
#	stat_summary(fun.y=sd, geom="line", aes(group=1),color='black')

# FIGURE S4B
model_used <- ifelse(use_final,'final_model','minimal_model' )
ggsave(plot=cum1,filename=paste('cumProb1.',model_used,'.eps',sep=''),height=15,width=20)
ggsave(plot=cum2,filename=paste('cumProb2.',model_used,'.eps',sep=''),height=15,width=20)
#ggsave(plot=cum3,filename=paste('cumProb3.',model_used,'.eps',sep=''),height=15,width=20)
ggsave(plot=cum4,filename=paste('cumProb4.',model_used,'.eps',sep=''),height=15,width=20)
       
library(gridExtra)
g<-arrangeGrob(conc1,cum1, cum2,cum4, ncol=2)
if(view.col){
	ggsave(plot=g,filename='jointProb.col.eps',height=15,width=20)
}else{
	ggsave(plot=g,filename='odds.bw.eps',height=15,width=15)
}

#write.csv(NEC_data_tmp,file='~/Documents/Grad/Labs/Lewis/NEC_HMO/odds/cumOdds/NEC_data_geometricMean.csv',quote=F,row.names=F)


######################################################
# 
#  Visualize OR
#  FIGURE 3B, 3C
#
####################################################### 
## visualize univariate OR
# subject grouped 
sumr = do.call(rbind,summaryHMO)

gee1 <- geeglm( NEC.binary ~ I(sqrt(DPP)) + scale(DSLNT_ug.mL) + scale(LNFP_I_ug.mL) + scale(DFLNT_ug.mL)    , 
	id = ParticipantID, data = NEC_data, family=binomial(link="logit"), 
	corstr = covarience, std.err="san.se") 
sumr = rbind(sumr,coef(summary(gee1))[-(1:2),])

sumr$Model = factor(c(rep("Prescreening",length(summaryHMO)),rep("Final Model",3) ),levels=rev(c('Final Model','Prescreening')))
sumr$HMO = gsub('_ug[.]mL|1|X|scale|[(]|[)]','',rownames(sumr))
sumr$HMO = factor(sumr$HMO,levels=sumr$HMO[order(-(sumr$'Pr(>|W|)'[1:length(summaryHMO)]))])

#limits <- aes(xmax = exp(sumr$Estimate + 1.96*sumr$Std.err), xmin=exp(sumr$Estimate - 1.96*sumr$Std.err))
p<- ggplot(sumr,aes(x=exp(Estimate),y=HMO))
#delta method: Std.Err(exp(\beta)) = exp(\beta) * Std.Err(\beta)
# CI = exp(0.1344)  +/- 1.96 * Std.Err(exp(\beta)
#p<-p + geom_point() + geom_errorbarh(aes(xmax = exp(Estimate + 1.96*Std.err), xmin=exp(Estimate - 1.96*Std.err)), height=0.0) + 
p<-p + geom_point() + geom_errorbarh(aes(xmax = exp(Estimate)  + 1.96 * exp(Estimate) * Std.err, xmin=exp(Estimate)  - 1.96 * exp(Estimate) * Std.err), height=0.0) + 
  facet_wrap(   ~ Model , nrow=1) +geom_vline(xintercept = 1) + xlab('Odds Ratio with 95% Confidence interval')
p

# FIGURE 3B, 3C
ggsave('ORvis.eps',p)
# DSLNT: 95% CI =  0.0449676823526655, LNFP1: 95% CI =  0.0638400157853984, DFLNT: 95% CI =  0.131122424069538 

######################################################
# 
#  Visualize HMO Conc
#
####################################################### 

# HMO Spaghetti
HMO_LNFP1 <- ggplot(data = NEC_data, aes(x=DPP, y=LNFP_I_ug.mL, group = ParticipantID, color=ifelse(NEC.binary,'NEC','Control' ) )) +
	geom_line()+ theme_bw() + scale_colour_hue( guide=guide_legend(title="Diagnosis")) + ylab('LNFP1') + xlab('')
HMO_DSLNT <- ggplot(data = NEC_data, aes(x=DPP, y=DSLNT_ug.mL, group = ParticipantID, color=ifelse(NEC.binary,'NEC','Control' ) )) +
	geom_line()+ theme_bw() + scale_colour_hue( guide=guide_legend(title="Diagnosis")) + ylab('DSLNT') + xlab('') + ggtitle("Selected and Suggested Contributors")
HMO_DFLNT <- ggplot(data = NEC_data, aes(x=DPP, y=DFLNT_ug.mL, group = ParticipantID, color=ifelse(NEC.binary,'NEC','Control' ) )) +
	geom_line()+ theme_bw() + scale_colour_hue( guide=guide_legend(title="Diagnosis")) + ylab('DFLNT')
HMO_LSTb <- ggplot(data = NEC_data, aes(x=DPP, y=LSTb_ug.mL, group = ParticipantID, color=ifelse(NEC.binary,'NEC','Control' ) )) +
	geom_line()+ theme_bw() + scale_colour_hue( guide=guide_legend(title="Diagnosis")) + ylab('LSTb') + xlab('') + ggtitle("Previously Rejected Contributors")
HMO_LNT <- ggplot(data = NEC_data, aes(x=DPP, y=LNT_ug.mL, group = ParticipantID, color=ifelse(NEC.binary,'NEC','Control' ) )) +
	geom_line()+ theme_bw() + scale_colour_hue( guide=guide_legend(title="Diagnosis")) + ylab('LNT') + xlab('')
HMO_Sum <- ggplot(data = NEC_data, aes(x=DPP, y=Sum_ug.mL, group = ParticipantID, color=ifelse(NEC.binary,'NEC','Control' ) )) +
	geom_line()+ theme_bw() + scale_colour_hue( guide=guide_legend(title="Diagnosis")) + ylab('Sum') 
g<-arrangeGrob(HMO_DSLNT, HMO_LSTb, HMO_LNFP1,HMO_LNT,HMO_DFLNT,HMO_Sum, ncol=2)
ggsave(plot=g,filename='HMO.spagetti.eps',height=15,width=15)

########################
#### Decision Tree
#library(RWeka)
#NEC_data_tmp = unique(na.omit(NEC_data[,c('NEC.binary','DSLNT_ug.mL','DFLNT_ug.mL','LNFP_I_ug.mL')])) 
#mod=J48(as.factor(NEC.binary) ~ . , data=NEC_data_tmp, control = Weka_control(M = 5))
#e <- evaluate_Weka_classifier(mod,cost = matrix(c(0,1,1,0), ncol = 2),numFolds = 10, complexity = TRUE, class = TRUE)

######################################################
# 
#  Decision Tree
#  FIGURE S3, 3, Table S2 
#
####################################################### 
library(partykit)

geeS <- geeglm( NEC.binary ~ I(sqrt(DPP)) + scale(DSLNT_ug.mL)   , 
	id = ParticipantID, data = NEC_data, family=binomial(link="logit"), 
	corstr = covarience, std.err="san.se") 
geeS.tab = parseSummary(coef(summary(geeS)))

geeSN <- geeglm( NEC.binary ~ I(sqrt(DPP)) + scale(DSLNT_ug.mL) + scale(LNFP_I_ug.mL)   , 
	id = ParticipantID, data = NEC_data, family=binomial(link="logit"), 
	corstr = covarience, std.err="san.se")
geeSN.tab = parseSummary(coef(summary(geeSN)))

geeSF <- geeglm( NEC.binary ~ I(sqrt(DPP)) + scale(DSLNT_ug.mL) + scale(DFLNT_ug.mL)   , 
	id = ParticipantID, data = NEC_data, family=binomial(link="logit"), 
	corstr = covarience, std.err="san.se") 
geeSF.tab = parseSummary(coef(summary(geeSF)))

geeSNF <- geeglm( NEC.binary ~ I(sqrt(DPP)) + scale(DSLNT_ug.mL) + scale(LNFP_I_ug.mL) + scale(DFLNT_ug.mL)   , 
	id = ParticipantID, data = NEC_data, family=binomial(link="logit"), 
	corstr = covarience, std.err="san.se") 
geeSNF.tab = parseSummary(coef(summary(geeSNF)))

geeSNF2 <- geeglm( NEC.binary ~ I(sqrt(DPP)) + scale(DSLNT_ug.mL) + scale(LNFP_I_ug.mL) + scale(DFLNT_ug.mL) + scale(LNFP_III_ug.mL)   , 
	id = ParticipantID, data = NEC_data, family=binomial(link="logit"), 
	corstr = covarience, std.err="san.se") 
geeSNF2.tab = parseSummary(coef(summary(geeSNF2)))

NEC_data_tmp = na.omit(NEC_data[,c('NEC.binary','NEC','DPP','DSLNT_ug.mL','LNFP_I_ug.mL','DFLNT_ug.mL','ParticipantID')])


combineby = geomeanprob_f # see uni.final.func.r

### remove minimum replicate samples for rolling average to work properly
NEC_data_tmp = NEC_data_tmp[! NEC_data_tmp$ParticipantID %in% c('E009','A064','C027'),]

# Models + cummulative odds and probability
NEC_data_tmp$Odds_Min_Model  =with(NEC_data_tmp, exp(  apply( cbind(1,sqrt(DPP) ,scale(DSLNT_ug.mL) ) %*% diag(coef(geeS)),1,sum ) ) )
NEC_data_tmp$Prob  = with(NEC_data_tmp, Odds_Min_Model/(1+Odds_Min_Model) )
#NEC_data_tmp = NEC_data_tmp %>% group_by(ParticipantID) %>% arrange(DPP) %>% mutate(cumProb2_Min_Model = 1-rollapply(1-Prob,2,function(x) prod(na.omit(x)),fill=NA))
NEC_data_tmp = NEC_data_tmp %>% group_by(ParticipantID) %>% arrange(DPP) %>% mutate(cumProb2_Min_Model = combineby(Prob,2))
NEC_data_tmp$cumOdds2_Min_Model = with(NEC_data_tmp, cumProb2_Min_Model/(1-cumProb2_Min_Model) )
NEC_data_tmp = NEC_data_tmp %>% group_by(ParticipantID) %>% arrange(DPP) %>% mutate(cumProb4_Min_Model = combineby(Prob,4))
NEC_data_tmp$cumOdds4_Min_Model = with(NEC_data_tmp, cumProb4_Min_Model/(1-cumProb4_Min_Model) )
NEC_data_tmp = NEC_data_tmp %>% group_by(ParticipantID) %>% arrange(DPP) %>% mutate(cumProb7_Min_Model = combineby(Prob,6))
NEC_data_tmp$cumOdds7_Min_Model = with(NEC_data_tmp, cumProb7_Min_Model/(1-cumProb7_Min_Model) )

#NEC_data_tmp2 = NEC_data_tmp2 %>% group_by(ParticipantID) %>% arrange(DPP) %>% mutate(tmp = combineby(Odds_Min_Model,6))

#stop()


NEC_data_tmp$Odds_Final_Model=with(NEC_data_tmp, exp(  apply( cbind(1,sqrt(DPP) ,scale(DSLNT_ug.mL)  , scale(LNFP_I_ug.mL), scale(DFLNT_ug.mL)) %*% diag(coef(geeSNF)),1,sum ) ) )
NEC_data_tmp$Prob  = with(NEC_data_tmp, Odds_Final_Model/(1+Odds_Final_Model) )
NEC_data_tmp = NEC_data_tmp %>% group_by(ParticipantID) %>% arrange(DPP) %>% mutate( cumProb2_Final_Model = combineby(Prob,2))
NEC_data_tmp$cumOdds2_Final_Model = with(NEC_data_tmp, cumProb2_Final_Model/(1-cumProb2_Final_Model) )
NEC_data_tmp = NEC_data_tmp %>% group_by(ParticipantID) %>% arrange(DPP) %>% mutate(cumProb4_Final_Model = combineby(Prob,4))
NEC_data_tmp$cumOdds4_Final_Model = with(NEC_data_tmp, cumProb4_Final_Model/(1-cumProb4_Final_Model) )
NEC_data_tmp = NEC_data_tmp %>% group_by(ParticipantID) %>% arrange(DPP) %>% mutate(cumProb7_Final_Model = combineby(Prob,6))
NEC_data_tmp$cumOdds7_Final_Model = with(NEC_data_tmp, cumProb7_Final_Model/(1-cumProb7_Final_Model) )

NEC_data_tmp$OddsSN =with(NEC_data_tmp, exp(  apply( cbind(1,sqrt(DPP) ,scale(DSLNT_ug.mL)  , scale(LNFP_I_ug.mL)) %*% diag(coef(geeSN)),1,sum ) ) )
NEC_data_tmp$OddsSF =with(NEC_data_tmp, exp(  apply( cbind(1,sqrt(DPP) ,scale(DSLNT_ug.mL)  , scale(DFLNT_ug.mL)) %*% diag(coef(geeSF)),1,sum ) ) )

#NEC_data_tmp$NEC.fact = factor( ifelse(NEC_data_tmp$NEC.binary,"NEC","Control") , levels=c('NEC','Control','NA'))
NEC_data_tmp$NEC.fact = as.factor(NEC_data_tmp$NEC)

library(RWeka)
mod=J48(as.factor(ifelse(NEC.binary,'NEC','Control')) ~ . ,
	data=NEC_data_tmp[,c('Odds_Min_Model','Odds_Final_Model','NEC.binary')], control = Weka_control(M = 20))
e <- evaluate_Weka_classifier(mod,cost = matrix(c(0,1,1,0), ncol = 2),numFolds = 10, complexity = TRUE, class = TRUE)
mod=J48(as.factor(ifelse(NEC.binary,'NEC','Control')) ~ . ,
	data=NEC_data_tmp[,c('DPP','DSLNT_ug.mL','LNFP_I_ug.mL','DFLNT_ug.mL','NEC.binary')], control = Weka_control(M = 10))
e <- evaluate_Weka_classifier(mod,cost = matrix(c(0,1,1,0), ncol = 2),numFolds = 10, complexity = TRUE, class = TRUE)
mod=J48(as.factor(ifelse(NEC.binary,'NEC','Control')) ~ . ,
	data=NEC_data_tmp[,c('DPP','DSLNT_ug.mL','LNFP_I_ug.mL','DFLNT_ug.mL','NEC.binary')], control = Weka_control(M = 5))
e <- evaluate_Weka_classifier(mod,cost = matrix(c(0,1,1,0), ncol = 2),numFolds = 10, complexity = TRUE, class = TRUE)

mod=J48(as.factor(ifelse(NEC.binary,'NEC','Control')) ~ . ,
	data=NEC_data_tmp[,c('cumOdds2_Min_Model','cumOdds2_Final_Model','NEC.binary')], control = Weka_control(M = 20))
e <- evaluate_Weka_classifier(mod,cost = matrix(c(0,1,1,0), ncol = 2),numFolds = 10, complexity = TRUE, class = TRUE)
mod=J48(as.factor(ifelse(NEC.binary,'NEC','Control')) ~ . ,
	data=NEC_data_tmp[,c('cumOdds4_Min_Model','cumOdds4_Final_Model','NEC.binary')], control = Weka_control(M = 20))
e <- evaluate_Weka_classifier(mod,cost = matrix(c(0,1,1,0), ncol = 2),numFolds = 10, complexity = TRUE, class = TRUE)
mod=J48(as.factor(ifelse(NEC.binary,'NEC','Control')) ~ . ,
	data=NEC_data_tmp[,c('cumOdds7_Min_Model','cumOdds7_Final_Model','NEC.binary')], control = Weka_control(M = 20))
e <- evaluate_Weka_classifier(mod,cost = matrix(c(0,1,1,0), ncol = 2),numFolds = 10, complexity = TRUE, class = TRUE)

# FIGURE S4
setwd(wd)
write.csv(NEC_data_tmp,file='odds/cumOdds/NEC_data_geomean.2.4.6.csv',quote=F,row.names=F)

# temporal decision trees (nlme) [not used]
if(F){
library(REEMtree)
NEC_data_tmp$NEC.bin.fact = as.factor(ifelse(NEC_data_tmp$NEC.binary,'nec','cnt' ))
reem.tree = REEMtree(NEC.bin.fact~I(sqrt(DPP))+scale(DSLNT_ug.mL)+scale(LNFP_I_ug.mL)+scale(DFLNT_ug.mL), data=NEC_data_tmp,
                     random=~ParticipantID,correlation=corSpher()) # spherical == exchangeable
reem.tree = REEMtree(NEC.fact~I(sqrt(DPP))+Odds_Min_Model+Odds_Final_Model, data=NEC_data_tmp,
                     random=~1|ParticipantID,correlation=corSpher()) # spherical == exchangeable

reem.tree = REEMtree(NEC.binary~scale(DSLNT_ug.mL)+scale(LNFP_I_ug.mL)+scale(DFLNT_ug.mL), data=NEC_data_tmp,
                     random=~1|ParticipantID,correlation=corSpher()) # spherical == exchangeable
reem.tree = REEMtree(NEC.binary~Odds_Min_Model+Odds_Final_Model, data=NEC_data_tmp,
                     random=~1|ParticipantID,correlation=corSpher()) # spherical == exchangeable

plot(reem.tree)
ranef(reem.tree) #random effects
}


######################################################
#
# Final Model Odds Visual
# FIGURE S4C
#
######################################################

# boxplot cumulative days x odds
df = data.frame( days=c(rep(1,nrow(NEC_data_tmp)),rep(2,nrow(NEC_data_tmp)),rep(4,nrow(NEC_data_tmp)),rep(6,nrow(NEC_data_tmp))),
                 NEC = rep(NEC_data_tmp$NEC.binary,4), ParticipantID = rep(NEC_data_tmp$ParticipantID,4),
                 Odds_Min_Model = with(NEC_data_tmp, c(Odds_Min_Model,cumOdds2_Min_Model,cumOdds4_Min_Model,cumOdds7_Min_Model)),
                 Odds_Final_Model = with(NEC_data_tmp, c(Odds_Final_Model,cumOdds2_Final_Model,cumOdds4_Final_Model,cumOdds7_Final_Model))
)
g1= ggplot(df,aes(x=as.factor(days),y=Odds_Min_Model,color=as.factor(NEC))) + geom_boxplot() 
g2= ggplot(df,aes(x=as.factor(days),y=Odds_Final_Model,color=as.factor(NEC))) + geom_boxplot()
# FIGURE S4C
ggsave(plot=g1,filename='odds/cumOdds/cumulativeOdds_boxplots_minModel.eps')
ggsave(plot=g2,filename='odds/cumOdds/cumulativeOdds_boxplots_FinalModel.eps')


### performance assessment
library(RegressionModelPipeline)
df$NEC=factor(df$NEC,levels=c('0','1'))
k=100
outL=list()
pdf('odds/perf.day1_min.pdf')
  outL[['day1_min']] =cross_assess_wrapper(data=na.omit(df[df$days==1,]),formula=formula('NEC~Odds_Min_Model'),resp='NEC',family='binomial',K=k,model=glm, cv_function=cross_valid_resampling_sensitive)
dev.off()
pdf('odds/perf.day1_full.pdf')
  outL[['day1_full']]=cross_assess_wrapper(data=na.omit(df[df$days==1,]),formula=formula('NEC~Odds_Final_Model'),resp='NEC',family='binomial',K=k,model=glm, cv_function=cross_valid_resampling_sensitive)
dev.off()
pdf('odds/perf.day2_min.pdf')
  outL[['day2_min']] =cross_assess_wrapper(data=na.omit(df[df$days==2,]),formula=formula('NEC~Odds_Min_Model'),resp='NEC',family='binomial',K=k,model=glm, cv_function=cross_valid_resampling_sensitive)
dev.off()
pdf('odds/perf.day2_full.pdf')
  outL[['day2_full']]=cross_assess_wrapper(data=na.omit(df[df$days==2,]),formula=formula('NEC~Odds_Final_Model'),resp='NEC',family='binomial',K=k,model=glm, cv_function=cross_valid_resampling_sensitive)
dev.off()
pdf('odds/perf.day4_min.pdf')
  try( outL[['day4_min']] <-cross_assess_wrapper(data=na.omit(df[df$days==4,]),formula=formula('NEC~Odds_Min_Model'),resp='NEC',family='binomial',K=k,model=glm, cv_function=cross_valid_resampling_sensitive)) 
dev.off()
pdf('odds/perf.day4_full.pdf')
  try(outL[['day4_full']]<-cross_assess_wrapper(data=na.omit(df[df$days==4,]),formula=formula('NEC~Odds_Final_Model'),resp='NEC',family='binomial',K=k,model=glm, cv_function=cross_valid_resampling_sensitive))
dev.off()
pdf('odds/perf.day6_min.pdf')
  outL[['day6_min']] =cross_assess_wrapper(data=na.omit(df[df$days==6,]),formula=formula('NEC~Odds_Min_Model'),resp='NEC',family='binomial',K=k,model=glm, cv_function=cross_valid_resampling_sensitive)
dev.off()
pdf('odds/perf.day6_full.pdf')
  outL[['day6_full']]=cross_assess_wrapper(data=na.omit(df[df$days==6,]),formula=formula('NEC~Odds_Final_Model'),resp='NEC',family='binomial',K=k,model=glm, cv_function=cross_valid_resampling_sensitive) 
dev.off()

df_pred=do.call(rbind, lapply(names(outL),function(x){
  data.frame(AUC=outL[[x]]$auc,
    MCC= unlist(lapply(outL[[x]]$mcc@y.values,function(y) max(na.omit(y[is.finite(y)]))) ),
    Consecutive_Samples=rep(strsplit(x,'_|day')[[1]][2],length(outL[[x]]$auc)) ,
    Model=toupper(rep(strsplit(x,'_')[[1]][2],length(outL[[x]]$auc) ))
      )
} ))
g1=ggplot(data=df_pred,aes(x=Consecutive_Samples,y=AUC,color=Model)) + geom_boxplot() + xlab('Consecutive Samples (days)') + theme(axis.text = element_text(size=12))
ggsave(plot=g1,filename='odds/perf.auc.eps')
g2=ggplot(data=df_pred,aes(x=Consecutive_Samples,y=MCC,color=Model)) + geom_boxplot() + xlab('Consecutive Samples (days)') + theme(axis.text = element_text(size=12))
ggsave(plot=g2,filename='odds/perf.mcc.eps')

g3=ggplot(data=df_pred,aes(x=AUC,y=MCC,color=Consecutive_Samples))+facet_grid(Model~Consecutive_Samples) + geom_density2d() + theme(axis.text = element_text(size=8))
ggsave(plot=g3,filename='odds/perf.aucVmcc.eps',width=12)

legend1 <- g_legend(g1)
legend3 <- g_legend(g3)

gs=list(g1+ theme(legend.position = 'none'),
        g2+ theme(legend.position = 'none'),
        g3+ theme(legend.position = 'none'),
        legend1,legend3)
lay <- rbind(c(1,2,4),
             c(3,3,5))
g=grid.arrange(grobs = gs, layout_matrix = lay)
ggsave(plot=g,filename='odds/perf.comb.eps',width=12,height=12)


ggplot(data=df,aes(x=Odds_Min_Model,color=NEC)) + geom_density() + facet_wrap(~days)

######################################################
#
# Fold change
# FIGURE 2
#
######################################################

# setwd(wd)
# source('dataVis/foldChange.r')