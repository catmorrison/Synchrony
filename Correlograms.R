####  Example R code to calculate correlograms
####  Pearson correlations produced in Pairs.R

library("mgcv")
library("ggplot2")
library("geosphere")
library(dplyr)
library(sf)

#dat<-read.csv('./Data/pairs_survival_CES.csv')
#dat<-read.csv('./Data/pairs_count_CES.csv')
dat<-read.csv('./Data/pairs_productivity_CES.csv')


###Filter out those pair-wise correlations based on less than 5 years
dat<-dat[dat$nyears>4,]
###Calculate distances between the squares

###Add on the centroids of the focal and secondary square
lats<-read.csv(paste0('./Data/sqq100centroids.csv'))

dat12<-merge(dat,lats,by.x='focal',by.y='sqq')
names(dat12)[6]<-'focal_long'
names(dat12)[5]<-'focal_lat'

dat13<-merge(dat12,lats,by.x='sqq',by.y='sqq')
names(dat13)[8]<-'sqq_long'
names(dat13)[7]<-'sqq_lat'

####Calculate distances between the 100km squares
dat13$dist<-distGeo(matrix(c(dat13$focal_long,dat13$focal_lat),ncol=2), matrix(c(dat13$sqq_long,dat13$sqq_lat),ncol=2))
dat13$dist2<-dat13$dist/1000

### loop to calculate max distances between each 100km square
sqqs<-unique(c(as.character(dat13$focal),as.character(dat13$sqq)))

maxdistance2<-NULL
for (ss in 1:length(sqqs)){
  ffocal<-dat13[dat13$focal%in%sqqs[ss],]
  ssqqs<-dat13[dat13$sqq%in%sqqs[ss],]
  dat1<-rbind(ffocal,ssqqs)
  maxdistance<-max(dat1$dist2)
  maxdistance2<-c(maxdistance2,maxdistance)
  print(ss)
} ###end of square loop

####Calculate the minimum maximum distance between the squares to use as a cut off for the scale
maxdistance3<-min(maxdistance2)  

### fit correlograms for each square

### Object to store data in
Ests2<-NULL

for (ss in 1:length(sqqs)){
  
  ffocal<-dat13[dat13$focal%in%sqqs[ss],]
  ssqqs<-dat13[dat13$sqq%in%sqqs[ss],]
  
  dat1<-rbind(ffocal,ssqqs)

 ###Fit the gam and make predictions and estimate confidence intervals
  m <- gam(cors ~ s(dist2), data = dat1)
  Vb <- vcov(m)  
  newd <-data.frame(dist2 = c(seq(0, max(dat1$dist2), length = 199),maxdistance3))
  newd<-sort(newd$dist2)
  newd<-data.frame(dist2=newd)
  
  pred <- predict(m, newd, se.fit = TRUE)
  se.fit <- pred$se.fit
  
  set.seed(42)
  N <- 10000
  
  rmvn <- function(n, mu, sig) {
    L <- mroot(sig)
    m <- ncol(L)
    t(mu + L %*% matrix(rnorm(m*n), m, n))
  }
  
  BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)), sig = Vb) 
  
  Cg <- predict(m, newd, type = "lpmatrix")  
  
  simDev <- Cg %*% t(BUdiff)
  
  absDev <- abs(sweep(simDev, 1, se.fit, FUN = "/"))
  
  masd <- apply(absDev, 2L, max)
  
  crit <- quantile(masd, prob = 0.95, type = 8)
  
  pred <- transform(cbind(data.frame(pred), newd),
                    uprP = fit + (2 * se.fit),
                    lwrP = fit - (2 * se.fit),
                    uprS = fit + (crit * se.fit),
                    lwrS = fit - (crit * se.fit))
  
  sims <- rmvn(N, mu = coef(m), sig = Vb)
  fits <- Cg %*% t(sims) 
  
  nrnd <-N
  rnd <- N
  stackFits <- stack(as.data.frame(fits[, rnd]))
  stackFits <- transform(stackFits, dist = rep(newd$dist, length(rnd)))

  #### Calculate the scale - the point at which predictions become < 0 
  predd<-sign(pred$fit)
  preddlo<-sign(pred$lwrS)
  preddup<-sign(pred$uprS)
  
  scale<-pred$dist2[which(predd==-1)[1]]
  scalelo<-pred$dist2[which(preddlo==-1)[1]]
  scaleup<-pred$dist2[which(preddup==-1)[1]]
  
  ###Cap everything at minimum maximum distance
  if(is.na(scale)==TRUE){scale<-maxdistance3}
  if(scale>maxdistance3){scale<-maxdistance3}
  
  ####Calculate the mean standard error of the scale (to use as a weight in GLMM)
  scale_se<-mean(I(pred$se.fit[1]:pred$se.fit[which(pred$dist2==scale)]))
  
  ####Calculate the strength
  newd3 <- with(dat1, data.frame(dist2 = 0))
  Cg <- predict(m, newd3, type = "lpmatrix") 
  fits2 <- Cg %*% t(sims)
  intercept<-apply(fits2,2,mean,na.rm=T)
  intercept<-mean(intercept,na.rm=T)
  int_CIs<-quantile(intercept,c(0.025,0.975),na.rm=T)
    
  ###Store results
  Ests<-data.frame(sqq=sqqs[ss],intercept=intercept,interceptlo=int_CIs[1],interceptup=int_CIs[2],scale=scale,scalelo=scalelo,scaleup=scaleup,scale_se=scale_se,minmaxdistance=maxdistance3)
  
  Ests2<-rbind(Ests2,Ests)
  
  ### Graph correlogram
  
  pp<-ggplot(data=pred[1:200,], aes(x = dist2, y = fit))+
    geom_ribbon(data=pred[1:200,],aes(ymin = lwrS, ymax = uprS), alpha = 0.2, fill = "red")+
    geom_path(lwd = 2) +
    labs(y = "Correlation coefficient",
         x = "Distance (km)")+
    geom_rug(data=dat1, aes(x =dist2), inherit.aes = F)+
    geom_hline(yintercept=0)+ 
    geom_vline(xintercept=scale,col=2,lwd=2,lty=2) ### add on scale 
  
  print(ss)
} ###End of square loop


###Ests2 contains the estimates of relative scale and strength for each square. 

#write.csv(Ests2,'./Data/correlogram_survival_CES.csv',row.names=F)  
#write.csv(Ests2,'./Data/correlogram_count_CES.csv',row.names=F)  
write.csv(Ests2,'./Data/correlogram_productivity_CES.csv',row.names=F)  

 