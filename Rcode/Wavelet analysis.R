####Example R code to carry out wavelet analysis
###Example species and country 
######################################################################

library(wsyn)

detrend<-read.csv('/Users/catmorrison/Dropbox/UEA/Synchrony/Paper 1/Manuscript/R code for MS/Data/detrended_count_CES.csv')

###Create N (number of 100km squares * T (years) matrix) 
times<-max(detrend$year2)
mat<-matrix(NA,length(unique(detrend$sqq)),times) 

###Loop over all squares and years to fill up matrix
sqqs<-unique(detrend$sqq)

###Fill matrix with annual predictions
for (i in 1:length(sqqs)){ ###rows
for (j in 1:times){
matt<-detrend$pred[detrend$year2==j & detrend$sqq==sqqs[i]]
mat[i,j]<-ifelse(is.numeric(matt) && length(matt) == 0 ,NA,matt)
} ###end of j
} ###end of i

###Replace NAs with the median value (across all years) of the square
meds<-apply(mat,1,median,na.rm=T)
for (i in 1:dim(mat)[1]){ ###rows
for (j in 1:times){
mat[i,j]<-ifelse(is.na(mat[i,j]) ,meds[i],mat[i,j])
}
}

###Wavelet analysis
ttimes<-1:times ###generate time steps
 
### clean up the time series - clev=2, time series are linearly detrended and de-meaned
mat<-cleandat(mat,ttimes,2)$cdat
  
###Compute wavelet phasor mean field
wavelet_1<-wpmf(mat,ttimes,sigmethod="quick")
times <- get_times(wavelet_1)
timescales <- get_timescales(wavelet_1)
wavelet_1_result <- Mod(get_values(wavelet_1))
 
plotmag(wavelet_1)



