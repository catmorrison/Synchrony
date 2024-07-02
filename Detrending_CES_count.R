####  Example R code to produce detrended annual estimates for synchrony models
####  Fits GLMM or GLM to produce annual estimates for each 100km square 
####  Fits GLM to detrend annual estimates

###Example species and country - for 39 randomly selected 100km squares###
######################################################################
library(lme4)
library(MuMIn)

count_dat<-read.csv('./Data/count_data_CES.csv')
count_dat$year.f<-as.factor(count_dat$year.f)

### loop to estimate annual variation in count for each 100km square       
sqqs<-unique(count_dat$sqq)

### Object to save annual predictions in 
newdata_detrend<-NULL
        
for (i in 1:length(sqqs)){
          
sqq.df<-count_dat[count_dat$sqq==sqqs[i],]

### Splits between GLMM and GLM - if more than 4 sites in 100km square then GLMM otherwise GLM
if (length(unique(sqq.df$siteID))>4){ 
mod_GLMM<-glmer(adults~year.f+(1|siteID),data=sqq.df,family='poisson')
newdata<-data.frame(year.f=sort(unique(sqq.df$year.f)))
newdata$pred<-predict(mod_GLMM,newdata,re.form=NA,type='response')
} ###end of if
else{####do a GLM
mod_GLM<-glm(adults~year.f,data=sqq.df,family='poisson')
newdata<-data.frame(year.f=sort(unique(sqq.df$year.f)))
newdata$pred<-predict(mod_GLM,newdata,type='response')
} ###end of else

###Detrend annual variation estimates###########
################################################

###Clean up variables
newdata$year<-as.numeric(as.character(newdata$year.f))
newdata$year2<-newdata$year-(min(newdata$year)-1)
newdata$year.c<-newdata$year2-mean(unique(newdata$year2)) ###Centered year

mod.detrend<-glm(pred~year.c,data=newdata)

###Predictions
  
newdata$predlinear<-predict(mod.detrend,newdata)
newdata$residuals<-newdata$predlinear-newdata$pred

###Graph residuals
          
minpreds<-min(c(newdata$pred,I(newdata$predlinear-newdata$residuals),I(newdata$pred-newdata$residuals),newdata$pred,newdata$predlinear))
maxpreds<-max(c(newdata$pred,I(newdata$predlinear+newdata$residuals),I(newdata$pred+newdata$residuals),newdata$pred,newdata$predlinear))

plot(pred~year,data=newdata,col=1,pch=20,cex=2,ylab='Predicted count',xlab='',main=sqqs[i],ylim=c(minpreds,maxpreds))
lines(predlinear~year,data=newdata,col=2)
segments(newdata$year,newdata$pred+newdata$residual,newdata$year,newdata$pred)
newdata$sqq<-sqqs[i]

###Save output for all squares          
newdata_detrend<-rbind(newdata_detrend,newdata)
          
} ###end of square loop i


####newdata_detrend contains detrended annual estimates for pairs analysis - see Pairs.R

write.csv(newdata_detrend,'./Data/detrended_count_CES.csv',row.names=F)        
   