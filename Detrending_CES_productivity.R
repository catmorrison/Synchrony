####  Example R code to produce detrended annual estimates of productivity for synchrony models
####  Fits GLMM or GLM to produce annual estimates for each 100km square 
####  Fits GLM to detrend annual estimates

###Example species and country 
######################################################################

library(lme4)
library(MuMIn)

###Bring in productivity data 

prod_dat<-read.csv('./Data/productivity_data_CES.csv')

prod_dat$response<-cbind(prod_dat$juveniles,prod_dat$adults)
prod_dat$year.f<-as.factor(prod_dat$year)

### loop to estimate annual variation in count for each 100km square  
sqqs<-unique(prod_dat$sqq)

### Object to save annual predictions in 
newdata_detrend<-NULL

for (i in 1:length(sqqs)){
  
sqq.df<-prod_dat[prod_dat$sqq==sqqs[i],]

### Splits between GLMM and GLM - if more than 4 sites in 100km square then GLMM otherwise GLM
if (length(unique(sqq.df$siteID))>4){  
mod_GLMM<-glmer(response~year.f+(1|siteID),data=sqq.df,family='binomial')
newdata<-data.frame(year.f=sort(unique(sqq.df$year.f)))
newdata$pred<-predict(mod_GLMM,newdata,re.form=NA)
newdata$predtrans<-exp(newdata$pred)
} ###end of if
else{####do a GLM
mod_GLM<-glm(response~year.f,data=sqq.df,family='binomial')
newdata<-data.frame(year.f=sort(unique(sqq.df$year.f)))
newdata$pred<-predict(mod_GLM,newdata)
newdata$predtrans<-exp(newdata$pred)
} ###end of else
  
###Detrend annual variation estimates###########
################################################

###Clean up variables
newdata$year<-as.numeric(as.character(newdata$year.f))
newdata$year2<-newdata$year-(min(newdata$year)-1)
newdata$year.c<-newdata$year2-mean(unique(newdata$year2)) ###Centered year
  
mod.detrend<-glm(predtrans~year.c,data=newdata)

###Predictions
newdata$predlinear<-predict(mod.detrend,newdata,type='response')
newdata$residuals<-newdata$predlinear-newdata$predtrans

###Graph residuals
plot(predtrans~year,data=newdata,col=1,pch=20,cex=2,ylab='Predicted productivity',xlab='',main=sqqs[i])
lines(predlinear~year,data=newdata,col=2)
segments(newdata$year,newdata$predtrans+newdata$residuals,newdata$year,newdata$predtrans)
  
newdata$sqq<-sqqs[i]
newdata_detrend<-rbind(newdata_detrend,newdata)
  
print(i)
} ###end of square loop i

####newdata_detrend contains detrended annual estimates for pairs analysis - see Pairs.R

write.csv(newdata_detrend,'./Data/detrended_productivity_CES.csv',row.names=F)   
