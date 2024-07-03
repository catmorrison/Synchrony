####  Example R code to calculate pearsons correlations between all pairs of 100km squares
####  Annual detrended estimates produced in Detrending.R

#detrended<-read.csv('./Data/detrended_survival_CES.csv')
#detrended<-read.csv('./Data/detrended_count_CES.csv')
detrended<-read.csv('./Data/detrended_productivity_CES.csv')

sqqs<-unique(detrended$sqq)

####Object to store pearsons correlations in for all squares
cors3<-NULL

###Start loop for each square
for (i in 1:length(sqqs)){ 
  
  ####Select focal square
  focal<-detrended[detrended$sqq%in%sqqs[i],]
  
  ####Find all secondary squares
  sqqs2<-unique(detrended$sqq)
  
  ####Object to store pearsons correlations in for each square
  cors2<-NULL
  
  ####Select a secondary square
  for (j in 1:length(sqqs2)){
    sqq<-detrended[detrended$sqq%in%sqqs2[j],]
    ####remove focal square from the main dataset
    restt<-detrended[!detrended$sqq%in%focal$sqq[1],]
    restt$sqq<-factor(restt$sqq)
    bboth<-merge(sqq,focal,by='year')
    ####Pearson correlation
    cors<-cor(bboth$residuals.y,bboth$residuals.x)
    cors<-data.frame(cors=cors,focal=focal$sqq[1],sqq=sqq$sqq[1],nyears=length(na.omit(bboth$year)))
    cors2<-rbind(cors2,cors)
  } ### End of j
  ####Remove focal square from the main dataset
  cors3<-rbind(cors3,cors2)
  detrended<-detrended[!detrended$sqq%in%focal$sqq,]
  detrended$sqq<-factor(detrended$sqq)
  if(length(unique(detrended$sqq))<2) break  ###End when only one site left in dataset
    print(i)
} ###End of i

####cors3 contains pearsons correlations for correlogram analysis- see Correlograms.R

#write.csv(cors3,'./Data/pairs_survival_CES.csv',row.names=F) 
#write.csv(cors3,'./Data/pairs_count_CES.csv',row.names=F)  
write.csv(cors3,'./Data/pairs_productivity_CES.csv',row.names=F)  







