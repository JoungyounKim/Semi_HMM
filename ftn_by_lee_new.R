
#############################################
#### Original code is in "code_from_lee"  ###
#############################################
# For Weighted ECDF, use ewcdf() under library(spatstat)

library(spatstat)
library(Iso)
library(fdrtool)

# This function calculates the distance $\gamma  \
#d_n(\hat{F}_{s,n}^{\gamma},\check{F}_{s,n}^\gamma)$
# for grid of gamma values in [0,1].
# Input data is a numeric vector containing observations from the mixture model.
# Input w is the vector of weights.
# Bigger gridsize  gives more  accurate estimates of alpha.

EstMixMdl <- function(data, w, gridsize)
{
n <- length(data)
w<-w[order(data)]
data <- sort(data)
#plot(data, w)
data.1 <- unique(data)
Fn <- ewcdf(data, w)
Fn.1 <- Fn(data.1)
## Calculate the known F_b at the data points
## Use Standard Normal CDF
Fb <- pnorm(data.1)
## Compute the weights (= frequency/n) of the unique data values, i.e., dF_n
# Note that Fn.1 is already weighted
Freq <- diff(c(0,Fn.1))
distance <- rep(0,gridsize)
##distance[0]<- sqrt(t((Fn.1-Fb)^2)%*%Freq)
distance[0]<- sqrt(sum((Fn.1-Fb)^2*Freq, na.rm=T))
for(i in 1:gridsize)
{
  a <- i/gridsize               ## Assumes a value of the mixing proportion
  F.hat <- (Fn.1-(1-a)*Fb)/a     ## Computes the naive estimator of F_s
  F.is=pava(F.hat,Freq,decreasing=FALSE) ## Computes the Isotonic Estimator of F_s
  F.is[which(F.is<=0)]=0
  F.is[which(F.is>=1)]=1
  F.is[is.na(F.is)] <- 0
  ##distance[i] <- a*sqrt(t((F.hat-F.is)^2)%*%Freq);
  distance[i] <- a*sqrt(sum((F.hat-F.is)^2*Freq, na.rm=T))
}
  return(list("dist"=distance, "Fn.1"=Fn.1, "Fb"=Fb))
}

# The following function evaluates the numerical second derivative of any function

Comp_2ndDer <- function(dist.alpha, gridsize)
  {
  dder <- diff(dist.alpha)    ## Computes the 1st order differences
  dder <- diff(dder)      ## Computes the 2nd order differences
  dder <- c(0,0,dder)       ## The numerical double derivative vector

  return(dder)
}

# Compute CDF

CDFEst <- function(Fn.1, Fb, Est)
{
## Compute the weights (= frequency/n) of the unique data values, i.e., dF_n
Freq <- diff(c(0,Fn.1))
## Computes the naive estimator of F_s
Est.CDF.naive <- (Fn.1-(1-Est)*Fb)/Est
## Computes the Isotonic Estimator of F_s
Est.CDF=pava(Est.CDF.naive,Freq,decreasing=FALSE)
Est.CDF[which(Est.CDF<=0)]=0
Est.CDF[which(Est.CDF>=1)]=1
Est.CDF[is.na(Est.CDF)]=0
#plot(Est.CDF.naive, main="naive")
#plot(Est.CDF, main="final")
return(cbind(Est.CDF.naive,Est.CDF))
}

# Compute Density
# This was also written by Patra and Sen
# Not sure if we need this.

DensEst <- function(data, Fn.1, Fb, Est)
{
F.hat <- (Fn.1-(1-Est)*Fb)/Est
Freq <- diff(c(0,Fn.1))
F.is <- pava(F.hat,Freq,decreasing=FALSE)
F.is[which(F.is<=0)] <- 0
F.is[which(F.is>=1)] <- 1
F.check <- F.is
data.1 <- unique(sort(data))
x <- data.1
y <- F.check
ll <- gcmlcm(x,y, type="lcm")
xtemp=rep(ll$x.knots,each=2)
ytemp=c(0,rep(ll$slope.knots,each=2),0)
ans<-rbind(t(xtemp),t(ytemp))
return(ans)
}

cdf_finall<-function(data, w, gridsize){

	estMix<-EstMixMdl(data, w, gridsize)
	
	dist.alpha <- estMix$dist
	Fn.1 <- estMix$Fn.1
	Fb <- estMix$Fb

	dder <- Comp_2ndDer(dist.alpha, gridsize)
	Est <- which.max(dder)/gridsize

	tmp<-CDFEst (Fn.1, Fb, Est)

	##############
	### Added on 20170407
	############

	my_unique_data<-unique(sort(data))
	
	tmp<-tmp[,2]
	#plot(tmp)
	return(list("alpha"=Est, "tmp"=tmp))
}


pdf_finall<-function(data, w, gridsize){

	estMix<-EstMixMdl(data, w, gridsize)
	
	dist.alpha <- estMix$dist
	Fn.1 <- estMix$Fn.1
	Fb <- estMix$Fb

	dder <- Comp_2ndDer(dist.alpha, gridsize)
	Est <- which.max(dder)/gridsize

	#tmp<-CDFEst (Fn.1, Fb, Est)
	tmp<-DensEst(data, Fn.1, Fb, Est)

	##############
	### Added on 20170407
	############

	#my_unique_data<-unique(sort(data))
	

	return(tmp)
}







