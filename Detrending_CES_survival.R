####  Example R code to produce detrended annual estimates of survival for synchrony models
####  Fits GLMM or GLM to produce annual estimates for each 100km square 
####  Fits GLM to detrend annual estimates

###Example species and country 

library(boot)
library(rjags)
library(jagsUI)
library(MCMCvis)

###Bring in capture history

CH<-read.csv('./Data/Capture_history_CES.csv')

### loop to estimate annual variation in count for each 100km square  
sqqs<-unique(CH$sqq)

### Object to save annual predictions in 
newdata_detrend<-NULL

for (i in 1:length(sqqs)){
  
  CHx<-CH[CH$sqq==sqqs[i],]
  
  ####Extract capture history data
  CHreal<-CHx[,3:I(CHx$nyears[1]+4)]   
  
  #### Extract variables for the model
  visit<-12-CHx$firstcapture_visit #### Visit on first capture
  last.year<-CHx$last.year ####Last year of capture for each indiviudal
  ff<-CHx$f #### First capture occasion for each individual
  CHsites<-unique(CHx$siteID) ####Number of unique sites

  #### Read in look up that indicates if sites were active in that year
  admat<-read.csv('./Data/activeLU_siteID.csv')
  ###Select sites and years in dataset
  admat<-admat[admat$X%in%CHsites,]
  admat<-admat[order(admat$X,CHsites),]
  admat<-admat[,colnames(admat)%in%colnames(CHreal)]
  
  ###Create a unique number ID for each site
  ssites<-data.frame(siteindex=1:length(unique(CHx$siteID)),siteID=unique(CHx$siteID))
  CHx_1<-merge(CHx,ssites,by='siteID')
  
  site<-CHx_1$siteindex  #### Unique ID for each site to be used in survival model
  
  #### Number of unique sites
  n.sites<-length(unique(CHx_1$siteindex))
  
  #### Number of unique ringing schemes
  nscheme<-1
  
  #### Function to create a matrix with information about known latent state z
  known.state.cjs <- function(ch){
    state <- ch 
    for (i in 1:dim(ch)[1]){
      n1 <- min(which(ch[i,]==1))  ###min 1
      n2 <- max(which(ch[i,]==1))  ###max 1
      state[i,n1:n2] <- 1  ###all inbetween min and max = 1
      state[i,n1] <- NA #####first cap =1
    }
    state[state==0] <- NA ###if not 1 NA
    return(state)
  }
  
  #### Bundle up JAGS data
  jags.data <- list(
    yy=CHreal, 
    ff=ff,
    n.individuals=nrow(CHreal), 
    n.years=ncol(CHreal)-1, 
    z=known.state.cjs(CHreal),
    visit=visit,
    site=site,
    n.sites=n.sites,
    last.year=last.year,
    admat=admat[,1:dim(admat)[2]],
    nscheme=nscheme
  )
  

  #### Initial values
  cjs.init.z <- function(ch,ff){ 
    for (i in 1:dim(ch)[1]){
      if (sum(ch[i,])==1) next
      n2 <- max(which(ch[i,]==1))
      ch[i,ff[i]:n2] <- NA
    }
    for (i in 1:dim(ch)[1]){
      ch[i,1:ff[i]] <- NA
    }
    return(ch)
  }
  
  #### initial values
  jags.inits <- function(){list(
    # z=cjs.init.z(CH,ff),
    beta.p=rep(0.4,n.sites),
    alpha.M=0.4,
    res=rep(0.6,nscheme),
    mean.sa=runif(ncol(CHreal)-2,0,1))
  }
  
  
  n.samples=5000000
  n.burnin=2000000
  n.thin=5000
  
  #### parameters to save 
  parameters=c('res','mu.sa','ep.sa','alpha.p','siteep','beta.p','tau.recapture','sigma.recapture','sigma2.recapture','sigma.survival.site','tau.survival.site','sigma2.survival.site')
  
  parameters_fixed=c('res','mu.sa','ep.sa','alpha.p','siteep','beta.p','tau.recapture','sigma.recapture','sigma2.recapture')
  
  ### Splits between GLMM and GLM - if more than 4 sites in 100km square then GLMM otherwise GLM
  if(length(CHsites)>4){
    
    out <- jags(data = jags.data,
                inits = jags.inits,
                parameters.to.save = parameters,
                model.file ='annual_survival_model.txt',
                n.chains = 2,
                n.adapt = 100,
                n.iter = n.samples,
                n.burnin = n.burnin,
                n.thin = n.thin,parallel=T)
  }
  
  
  if(length(CHsites)<5){
    
    out <- jags(data = jags.data,
                inits = jags.inits,
                parameters.to.save = parameters_fixed,
                  model.file ='annual_survival_model_fixed.txt',
                n.chains = 2,
                n.adapt = 100,
                n.iter = n.samples,
                n.burnin = n.burnin,
                n.thin = n.thin,parallel=T)
  }
  
  
  xx<-out$samples
  x <- mcmc.list(xx)
  
  ####Make annual estimates and convergence statistics
  ssummary<-MCMCsummary(out, round = 2)
  ssummary$parameter<-row.names(ssummary)
  
  predicted<-MCMCchains(out, params = 'mu.sa')
  predicted<-inv.logit(predicted)
  
  predicted2<-colMeans(predicted)
  predictedCRlo<-apply(predicted,2,quantile,0.025)
  predictedCRup<-apply(predicted,2,quantile,0.975)
  
  ###Detrend annual variation estimates###########
  ################################################
  
  yyear<-names(CHreal)[1:I(length(names(CHreal))-2)]
  yyear<-substr(yyear,2,5)
  yyear<-as.numeric(yyear)
  
  newdata<-data.frame(predicted=predicted2,predictedCRlo=predictedCRlo,predictedCRup=predictedCRup,year=yyear)
 
  ####Just predict and graph years in dataset
  newdata<-newdata[newdata$year%in%yyear,]
  
  newdata$year2<-newdata$year-(min(newdata$year)-1)
  newdata$year.c<-newdata$year2-mean(unique(newdata$year2)) ###Centered year
  newdata$year.c2<-newdata$year.c^2
  
  mod.detrend<-glm(predicted~year.c,data=newdata)
 
  ###Predictions
  newdata$predlinear<-predict(mod.detrend,newdata,type='response')
  newdata$residuals<-newdata$predlinear-newdata$predicted
  
  ###Graph residuals
  mmin<-min(c(newdata$predlinear,newdata$predicted+newdata$residual,newdata$predicted))
  mmax<-max(c(newdata$predlinear,newdata$predicted+newdata$residual,newdata$predicted))
  
  plot(predicted~year,data=newdata,col=1,pch=20,cex=2,ylab='Predicted survival',xlab='',main=sqqs[i],ylim=c(mmin,mmax))
  lines(predlinear~year,data=newdata,col=2)
  segments(newdata$year,newdata$predicted+newdata$residual,newdata$year,newdata$predicted)
  
  newdata$sqq<-sqqs[i]
  
  newdata_detrend<-rbind(newdata_detrend,newdata)
  
  print(i)
  
} ###end of i


####newdata_detrend contains detrended annual estimates for pairs analysis - see Pairs.R

write.csv(newdata_detrend,'./Data/detrended_survival_CES.csv',row.names=F)  



