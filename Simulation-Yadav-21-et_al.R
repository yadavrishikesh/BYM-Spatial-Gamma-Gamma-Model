##########################################
### LOAD LIBRARIES #######################
##########################################
rm(list=ls()) ### Deletes all current variables
library(mvtnorm)
library(MASS)
library(fields)
###########################
### SIMULATION SETTINGS ###
###########################

#DATA DIMENSION
nsites<-20#no of sites
ntime<-50# no of data at each site
n<-nsites*ntime #no of data points in space and time

#PARAMETERS for Gamma-Gamma model: Y(s)=X1(s)/X2(s), where:
# X1(s): iid Gamma(1,beta1), with rate 1 and shape beta1>0
# X2(s): marginally Gamma(alpha,beta2), with rate alpha>0 and shape beta2>0, dependence structure governed by Gaussian copula with exponential correlation function exp(-h/rho), h>=0, with range rho>0
alpha0<-1  # log scale
alpha1 <- 1# log scale
alpha2<-1 # log scale
alpha3<-1 # log scale
beta1 <- 5 #gamma SHAPE for X1 (numerator)
beta2 <- 5 #gamma SHAPE for X2 (denominator), related to the tail index as xi=1/beta2
rho <- 1 #RANGE parameter in exponential correlation function exp(-h/rho), h>=0, for the Gaussian copula in X2 (denominator)

#### tilde parametrization
alpha0.tilde<-alpha0+log(beta1)-log(beta2)
alpha1.tilde<-alpha1
alpha2.tilde<-alpha2
alpha3.tilde<-alpha3
beta1.tilde<-alpha0.tilde+log(beta1)
beta2.tilde<-log(1/beta2)
rho.tilde<-log(rho)

param <- c(alpha0,alpha1,alpha2,alpha3,beta1,beta2,rho) #Combines all hyperparameters in a single vector
log.tilde.param<-c(alpha0.tilde, alpha1.tilde, alpha2.tilde, alpha3.tilde,  beta1.tilde, beta2.tilde, rho.tilde) #Combines all transformed hyperparameters in a single vector

#SPATIAL LOCATIONS
set.seed(1)
loc<-cbind(runif(nsites),runif(nsites)) #locations generated uniformly in unit square [0,1]^2
#spatial covariates
Z1<-loc[,1]
Z2<-loc[,2]
set.seed(1)
### Simulating the third covarites by gaussian random filed with exponential covance function 
#DISTANCE MATRIX BETWEEN OUR SITES
dist.mat <- as.matrix(dist(loc)) #distance matrix of all locations (dimension: nsites x nsites)
### simulating the third covariates
range.param<-2
cov.grf<-exp(-dist.mat/range.param)
Z3<-mvtnorm::rmvnorm(n=1,mean = rep(0, nrow(cov.grf)), sigma =cov.grf)

###gamma RATE for X2 (denominator)
alpha<-exp(alpha0+alpha1*Z1+alpha2*Z2+alpha3*Z3)


#######################
### DATA SIMULATION ###
#######################

#sim.data: function to simulated data
# INPUTS:
#'   @ntime: integer, number of independent timen replicates
#'  @loc: matrix (dimension: nsites x 2), the i-th row gives the location of the i-th site
#'  @log.tilde.param: vector of transformed hyperparameters c(alpha0.tilde, alpha1.tilde, alpha2.tilde, alpha3.tilde,  beta1.tilde, beta2.tilde, rho.tilde)
# OUTPUT: list with the following elements
#'   @Y: simulated data matrix (dimension: ntime x nsites), the i-th row gives the i-th temporal replicate at all the sites
sim.data <- function(ntime,loc,log.tilde.param){
  
  nsites <- nrow(loc) #number of sites
  n <- nsites*ntime #total number of observations in space and time
  
  ## reverse transformation
  param <- c((2*log.tilde.param[1]-log.tilde.param[5]-log.tilde.param[6]),
             log.tilde.param[2],log.tilde.param[3],log.tilde.param[4],
            exp(log.tilde.param[5]-log.tilde.param[1]),
             exp(-log.tilde.param[6]),exp(log.tilde.param[7]))
  alpha.cov<-exp(param[1]+param[2]*Z1+param[3]*Z2+param[4]*Z3)
  beta1<-param[5]
  beta2 <- param[6]
  rho<-param[7]
  
  alpha <- matrix(rep(alpha.cov,ntime),nrow = ntime,ncol = nsites,byrow = TRUE) 
  dist.mat <- as.matrix(dist(loc)) #distance matrix of all locations (dimension: nsites x nsites)
  Sigma <- exp(-dist.mat/rho) #correlation matrix for Gaussian copula
  
  X1 <- matrix(rgamma(n,shape=beta1,rate=1),nrow=ntime,ncol=nsites) #process X1(s) in numerator as independent gamma random variables (dimension: ntime x nsites)
  Gauss <- mvtnorm::rmvnorm(n=ntime,mean=rep(0,nsites),Sigma) #simulation of ntime independent multivariate Gaussian vectors (dimension: ntime x nsites)
  X2 <- qgamma(pnorm(Gauss),shape=beta2,rate=alpha) #process X2(s) in denominator as independent vectors with gamma margins and Gaussian copula (dimension: ntime x nsites)
  
  Y <- X1/X2 #data simulation as gamma ratio
  return(list(Y=Y,X1=X1,X2=X2))
}
set.seed(1)
data.simulated <- sim.data(ntime,loc,log.tilde.param)
Y <- data.simulated$Y #data matrix (dimension: ntime x nsites)
X1 <- data.simulated$X1 #process X1(s) in numerator as independent gamma random variables (dimension: ntime x nsites)
X2 <- data.simulated$X2 #process X2(s) in denominator as independent vectors with gamma margins and Gaussian copula (dimension: ntime x nsites)

# #CHECK MARGINS OF SIMULATED DATA
# i <- 5#takes i-th location
# qqplot(Y[,i],alpha[i]*(beta1/beta2)*qf(c(1:ntime)/(ntime+1),df1=2*beta1,df2=2*beta2),log="xy")
# abline(0,1,col="red")


###########################################
### USEFUL FUNCTIONS FOR MCMC ALGORITHM ###
###########################################

#CORRELATION MATRIX
#Sigma.function: function to calculate the correlation matrix of the latent process X2(s) at some fixed sites
# INPUTS:
#'   @rho: positive real number, range parameter of exponential correlation function exp(-h/rho), h>=0
#'  @dist.mat: distance matrix between locations (dimension: nsites x nsites)
# OUTPUT:
#   correlation matrix Sigma (dimension: nsites x nsites)
Sigma.function <- function(rho,dist.mat){ # correlation matrix for some fixed sites
  return(exp(-dist.mat/rho))
}

# PC priors for beta2
PC.prior.beta2 <- function(beta,lambda){
  if(beta>1){
    l.d<-sqrt(2)*lambda*exp(-sqrt(2)*lambda*((beta*(beta-1)))^(-1/2))*(beta-1/2)*(beta*(beta-1))^(-3/2) 
  }else {l.d=10^(-50)}
  return(l.d)
}

# PC priors for beta1
KLD <- function(beta){
  return( (beta-1)*digamma(beta)-log(gamma(beta)) )	
}
dKLD <- function(beta){
  return( (beta-1)*trigamma(beta) )	
}
d.beta <- function(beta){
  return( sqrt(2*KLD(beta)) )
}
dd.beta <- function(beta){
  return( dKLD(beta)/sqrt(2*KLD(beta)) )
}
PC.prior.beta1 <- function(beta,lambda){
  res <- c()
  res[beta==0] <- 0
  res[beta==1] <- 0.5*lambda*mean(abs(c(dd.beta(1-10^(-6)),dd.beta(1+10^(-6)))))
  res[beta!=0 & beta!=1] <- 0.5*lambda*exp(-lambda*d.beta(beta[beta!=0 & beta!=1]))*abs(dd.beta(beta[beta!=0 & beta!=1])) ## The factor 0.5 is because we combine together the cases beta<1 and beta>1... (the PC prior is in fact a mixture between two priors defined over the intervals (0,1) and (1,Inf).)
  return( res )
}

# scale parameter in PC priors
lambda1=3  #for beta1
lambda2=3 #for beta2



#####  Creating artificial missingness in data for predictions
Y1<-Y
n.predict<-nsites/5
n.fit<-nsites-n.predict

prob<-0.75
fun.quan<-function(x){
  return(quantile(x[x>0],na.rm = TRUE,probs=prob))
}
quantile.matrix.fit<-matrix(rep(unname(apply(Y1[,-(1:n.predict)], 2, fun.quan)),ntime),nrow=ntime,ncol=n.fit,byrow = TRUE)

max.threhold<-10^3
## Creating the data full data Y (predict data + fit data) and threshold matrix u.thr (predict threshold + fit threshold)
## threshold for the stations where we fit th model (183x89) 
u.thr.fit<-matrix(NA,nrow=ntime,ncol=n.fit)
## data for fitting
Y.fit<-matrix(NA,nrow=ntime,ncol=n.fit)
for(i in 1: ntime){
  for (j in 1: n.fit) {
    if(is.na(Y1[,-(1:n.predict)][i,j])==TRUE){
      u.thr.fit[i,j]<-max(Y1,na.rm=TRUE)+max.threhold
      Y.fit[i,j]<- NA#-10^6 
    } else{
      u.thr.fit[i,j]<-quantile.matrix.fit[i,j]
      Y.fit[i,j]<-Y1[,-(1:n.predict)][i,j]
    }
  }
}

## threshold for missing sites
## threshold for missing sites assumed to be +10^3 than the maximum observed precipitation
u.thr.miss<-matrix(max(Y1,na.rm=TRUE)+max.threhold,nrow = ntime,ncol=n.predict)
## Threshold for all sites (predict+fit)
u.thr<-cbind(u.thr.miss,u.thr.fit)
data.to.predict<-matrix(NA,nrow = ntime,ncol=n.predict)
## data for all sites  (predict+fit)
Y<-cbind(data.to.predict,Y.fit)

#index for the censored observations
ind1<-Y<u.thr
ind1[is.na(ind1)] <- TRUE
## index for non-censored observations
ind2<-Y>=u.thr
ind2[is.na(ind2)] <- FALSE

#' @LOG-POSTERIOR-DENSITY 
# INPUTS:
#'   @Y: data matrix (dimension: ntime x nsites), the i-th row gives the i-th temporal replicate at all the sites
#'   @X2: latent parameters X2(s) in the denominator (dimension: ntime x nsites)
#'   @param: vector of hyperparameters 
#'   @ind1: censored index (dimension: ntime x nsites)
#'  @ind2: non-censored index (dimension: ntime x nsites)
#'  @u.thr: threshold matrix (dimension: ntime x nsites)
#'   @dist.mat: distance matrix between locations (dimension: nsites x nsites)
# OUTPUT:
#'  @log.post: value of the log-posterior density
log.posterior<-function(Y, X2, param, dist.mat, ind1, ind2,u.thr){
  ntime <- nrow(Y)
  nsites <- ncol(Y)
  n <- ntime*nsites
  alpha0 <- param[1]
  alpha1 <- param[2]
  alpha2 <- param[3]
  alpha3 <- param[4]
  beta1 <- param[5]
  beta2 <- param[6]
  rho <- param[7]
  alpha<-exp(matrix(rep((alpha0+alpha1*Z1+alpha2*Z2+alpha3*Z3),ntime),nrow = ntime,ncol = nsites,byrow = TRUE))
  
  log.prior.alpha0 <- dnorm(alpha0,0,10,log=TRUE) # log-prior for alpha0
  log.prior.alpha1<- dnorm(alpha1,0,10,log=TRUE) # log-prior for alpha1
  log.prior.alpha2 <- dnorm(alpha2,0,10,log=TRUE) # log-prior for alpha2
  log.prior.alpha3 <- dnorm(alpha3,0,10,log=TRUE)
  # log.prior.beta1 <- dgamma(beta1,shape = 1/100,rate=1/100,log=TRUE) # log-prior for beta1
  # log.prior.beta2 <- dgamma(beta2,shape = 1/100,rate=1/100,log=TRUE) # log-prior for beta2
  log.prior.beta1 <- log(PC.prior.beta1(beta=beta1,lambda=lambda1)) # log-prior for beta1
  log.prior.beta2<-log(PC.prior.beta2(beta=beta2,lambda=lambda2)) # log-prior for beta2
  log.prior.rho <- dgamma(rho,shape = 1/100,rate=1/100,log=TRUE) # log-prior for rho
  
  log.density.X1 <- sum(dgamma(Y[ind2], shape = beta1, rate = X2[ind2], log = TRUE))+sum(pgamma(u.thr[ind1]*X2[ind1],shape = beta1,rate = 1,log.p=TRUE)) # log-density for X1(s) (numerator) given X2(s) (denominator)
  x.arg<- qnorm(matrix(pgamma(X2, shape = beta2, rate = alpha),nrow = ntime, ncol = nsites))
  log.density.X2 <- sum(mvtnorm::dmvnorm(x.arg, mean=rep(0,nsites), sigma=Sigma.function(rho=rho,dist.mat=dist.mat),log = TRUE)) + sum(dgamma(X2, shape = beta2, rate = alpha, log = TRUE))-sum(dnorm(x.arg, log=TRUE)) # log-density for X2(s) (denominator)
  
  
  X1 <- matrix(rgamma(n,shape=beta1,rate=1),nrow=ntime,ncol=nsites) #process X1(s) in numerator as independent gamma random variables (dimension: ntime x nsites)
  Gauss <- mvtnorm::rmvnorm(n=ntime,mean=rep(0,nsites),Sigma) #simulation of ntime independent multivariate Gaussian vectors (dimension: ntime x nsites)
  X2 <- qgamma(pnorm(Gauss),shape=beta2,rate=alpha) #process X2(s) in denominator as independent vectors with gamma margins and Gaussian copula (dimension: ntime x nsites)
  

  log.post<-log.prior.alpha0+log.prior.alpha1+log.prior.alpha2+log.prior.alpha3+log.prior.beta1+log.prior.beta2+log.prior.rho+log.density.X1+log.density.X2
  return(log.post)
}

#'@Tilde.LOG-POSTERIOR DENSITY based on transformed hyperparameters
#'  @Y: data matrix (dimension: ntime x nsites), the i-th row gives the i-th temporal replicate at all the sites
#'   @log.X2: latent parameters X2(s) in the denominator  at log scale (dimension: ntime x nsites)
#'   @log.tilde.param: vector of transformed hyperparameters 
#'  @ind1: censored  index (dimension: ntime x nsites)
#'  @ind2: censored  index (dimension: ntime x nsites)
#'  @u.thr: threshold matrix (dimension: ntime x nsites)
#'  @dist.mat: distance matrix between locations (dimension: nsites x nsites)
# OUTPUT:
#'  @log.post.tilde: value of the log-posterior density based on transformed hyperparameters
tilde.log.posterior<-function(Y, log.X2, log.tilde.param, dist.mat,ind1,ind2,u.thr){
  ntime <- nrow(Y)
  nsites <- ncol(Y)
  n <- ntime*nsites
  X2<-exp(log.X2)
  param <- c((2*log.tilde.param[1]-log.tilde.param[5]-log.tilde.param[6]),log.tilde.param[2],log.tilde.param[3],log.tilde.param[4],exp(log.tilde.param[5]-log.tilde.param[1]),exp(-log.tilde.param[6]),exp(log.tilde.param[7]))
  log.post.tilde<-log.posterior(Y, X2, param, dist.mat,ind1,ind2,u.thr) + sum(log.X2) + log.tilde.param[5] - log.tilde.param[1] - log.tilde.param[6]+log.tilde.param[7]
  return(log.post.tilde)
}

#To check the tilde.function is correctly parametrized
# tilde.log.posterior(Y, log.X2=log(X2), log.tilde.param, dist.mat,ind1,ind2,u.thr)
# log.posterior(Y, X2, param, dist.mat,ind1,ind2,u.thr)

#' @GRADIENT OF LOG-POSTERIOR DENSITY WITH RESPECT TO LATENT PARAMETRS
#' @gradient.log.posterior.numer: function to NUMERICALLY calculate the gradient of the log-posterior density
# INPUTS:
#' @Y: data matrix (dimension: ntime x nsites), the i-th row gives the i-th temporal replicate at all the sites
#' @X2: latent parameters X2(s) in the denominator (dimension: ntime x nsites)
#' @param: vector of hyperparameters (alpha,beta1,beta2,rho)
#' @dist.mat: distance matrix between locations (dimension: nsites x nsites)
#' @delta: positive real number, the precision for the numerical derivative
#' @ind1: censored  index (dimension: ntime x nsites)
#' @ind2: non-censored  index (dimension: ntime x nsites)
#' @u.thr: threshold matrix (dimension: ntime x nsites)
# OUTPUT:
#   gradient of the log-posterior density with respect to latent parameters (dimension: ntime x nsites)
gradient.log.posterior.numer<-function(Y, X2, param, dist.mat, delta,ind1,ind2,u.thr){
  ntime <- nrow(Y)
  nsites <- ncol(Y)
  n <- ntime*nsites
  
  gradient.log.post <- matrix(nrow=ntime,ncol=nsites)
  
  for(i in 1:ntime){
    for(j in 1:nsites){
      X2.plus <- X2.minus <- X2
      X2.plus[i,j] <- X2[i,j]+delta
      X2.minus[i,j] <- X2[i,j]-delta
      gradient.log.post[i,j] <- (log.posterior(Y=Y,X2=X2.plus,param=param,dist.mat=dist.mat,ind1,ind2,u.thr)-log.posterior(Y=Y,X2=X2.minus,param=param,dist.mat=dist.mat,ind1,ind2,u.thr))/(2*delta)
    }
  }
  
  return(gradient.log.post)
}

#' @gradient.log.posterior.theor: function to THEORETICALLY calculate the gradient of the log-posterior density
# INPUTS:
#'  @Y: data matrix (dimension: ntime x nsites), the i-th row gives the i-th temporal replicate at all the sites
#'  @X2: latent parameters X2(s) in the denominator (dimension: ntime x nsites)
#'  @param: vector of hyperparameters (alpha,beta1,beta2,rho)
#'  @dist.mat: distance matrix between locations (dimension: nsites x nsites)
#'  @ind1: censored  index (dimension: ntime x nsites)
#'  @u.thr: threshold matrix (dimension: ntime x nsites)
# OUTPUT:
#' @log.gradient: gradient of the log-posterior density with respect to latent parameters (dimension: ntime x nsites)
gradient.log.posterior.theor<-function(Y, X2, param, dist.mat,ind1,u.thr){ 
  ntime <- nrow(Y)
  nsites <- ncol(Y)
  n <- ntime*nsites
  
  alpha0 <- param[1]
  alpha1 <- param[2]
  alpha2 <- param[3]
  alpha3 <- param[4]
  beta1 <- param[5]
  beta2 <- param[6]
  rho <- param[7]
  alpha<-exp(matrix(rep((alpha0+alpha1*Z1+alpha2*Z2+alpha3*Z3),ntime),nrow = ntime,ncol = nsites,byrow = TRUE))
  
  log.gradient.X1<-ifelse(ind1, (u.thr*dgamma(u.thr*X2,shape = beta1, rate=1))/(pgamma(u.thr*X2,shape=beta1,rate=1)),(beta1/X2)-Y)
  
  latent.1<-((beta2-1)/X2)-alpha
  quantile.norm<-qnorm(pgamma(X2, shape = beta2, rate = alpha))
  grad.quantile.norm<-dgamma(X2,shape=beta2, rate=alpha)/dnorm(quantile.norm)
  latent.2<--(grad.quantile.norm)*t(solve(Sigma.function(rho,dist.mat))%*%t(quantile.norm))+quantile.norm*grad.quantile.norm
  log.gradient <- log.gradient.X1+latent.1+latent.2  
  return(log.gradient)
}
#To check of the gradinet is coded correctly; compare with numerical gradients 
# delta <- 10^(-4)
# gradient.log.posterior.numer(Y=Y,X2=X2,param=param,dist.mat=dist.mat,delta=delta,ind1=ind1,ind2=ind2,u.thr=u.thr)[1:10]
# gradient.log.posterior.theor(Y,X2,param,dist.mat,ind1,u.thr)[1:10]


#' @tilde.gradient.log.posterior.theor: function to THEORETICALLY calculate the gradient of the log-posterior density in terms of log scale of the latent varibles
# INPUTS:
#'   @Y: data matrix (dimension: ntime x nsites), the i-th row gives the i-th temporal replicate at all the sites
#'   @log.X2: latent parameters X2(s) in the denominator  at log scale (dimension: ntime x nsites)
#'  @log.tilde.param: vector of transformed hyperparameters 
#'  @ind1: censored  index (dimension: ntime x nsites)
#'  @u.thr: threshold matrix (dimension: ntime x nsites)
#'  @dist.mat: distance matrix between locations (dimension: nsites x nsites)
# OUPUT:
#  theoretical gradients of log posterior wrt to log.X2
tilde.gradient.log.posterior.theor<-function(Y, log.X2, log.tilde.param, dist.mat,ind1,u.thr){ 
  ntime <- nrow(Y)
  nsites <- ncol(Y)
  n <- ntime*nsites
  X2<-exp(log.X2)
  param <- c((2*log.tilde.param[1]-log.tilde.param[5]-log.tilde.param[6]),log.tilde.param[2],log.tilde.param[3],log.tilde.param[4],exp(log.tilde.param[5]-log.tilde.param[1]),exp(-log.tilde.param[6]),exp(log.tilde.param[7]))
  
  return(gradient.log.posterior.theor(Y, X2, param, dist.mat,ind1,u.thr)*X2+1)
}


#' @tilde.gradient.log.posterior.theor: function to NUMERICALLY to calculate the gradient of the log-posterior density in terms of log scale of the latent varibles i.e., log.X2=log(X2)
# INPUTS:
#'   @Y: data matrix (dimension: ntime x nsites), the i-th row gives the i-th temporal replicate at all the sites
#'  @log.X2: latent parameters X2(s) in the denominator  at log scale (dimension: ntime x nsites)
#'  @log.tilde.param: vector of transformed hyperparameters 
#'  @ind1: censored  index (dimension: ntime x nsites)
#'  @u.thr: threshold matrix (dimension: ntime x nsites)
#'  @dist.mat: distance matrix between locations (dimension: nsites x nsites)
# OUPUT:
# numerical gradients of log posterior wrt to log.X2
tilde.gradient.log.posterior.numer<-function(Y, log.X2, log.tilde.param, dist.mat, delta,ind1,ind2,u.thr){
  ntime <- nrow(Y)
  nsites <- ncol(Y)
  n <- ntime*nsites
  X2<-exp(log.X2)
  param <- c((2*log.tilde.param[1]-log.tilde.param[5]-log.tilde.param[6]),log.tilde.param[2],log.tilde.param[3],log.tilde.param[4],exp(log.tilde.param[5]-log.tilde.param[1]),exp(-log.tilde.param[6]),exp(log.tilde.param[7]))
  
  return(gradient.log.posterior.numer(Y, X2, param, dist.mat,delta,ind1,ind2,u.thr)*X2+1)
}

#Other tilde function...
tilde.gradient.log.posterior.numer2<-function(Y, log.X2, log.tilde.param, dist.mat, delta,ind1,ind2,u.thr){
  
  ntime <- nrow(Y)
  nsites <- ncol(Y)
  n <- ntime*nsites
  
  gradient.log.post <- matrix(nrow=ntime,ncol=nsites)
  
  for(i in 1:ntime){
    for(j in 1:nsites){
      log.X2.plus <- log.X2.minus <- log.X2
      log.X2.plus[i,j] <- log.X2[i,j]+delta
      log.X2.minus[i,j] <- log.X2[i,j]-delta
      gradient.log.post[i,j] <- (tilde.log.posterior(Y=Y,log.X2=log.X2.plus,log.tilde.param=log.tilde.param,dist.mat=dist.mat, ind1,ind2,u.thr)-tilde.log.posterior(Y=Y,log.X2=log.X2.minus,log.tilde.param=log.tilde.param,dist.mat=dist.mat,ind1,ind2,u.thr))/(2*delta)
    }
  }
  return(gradient.log.post)
}

#To check above function (should give the same result)
# delta <- 10^(-4)
# tilde.gradient.log.posterior.numer(Y,log.X2=log(X2),log.tilde.param,dist.mat,delta,ind1,ind2,u.thr)[1:10]
# tilde.gradient.log.posterior.numer2(Y,log.X2=log(X2),log.tilde.param,dist.mat,delta,ind1,ind2,u.thr)[1:10]
# tilde.gradient.log.posterior.theor(Y,log.X2=log(X2),log.tilde.param,dist.mat,ind1,u.thr)[1:10]


#' @proposal.mala.hyper: RANDOM-WALK-PROPOSALS FOR HYPERPARAMETERS BLOCKS
# INPUTS:
#' @log.tilde.param: vector of transformed hyperparameters 
#' @sigma.hyper2: tuning parameter for random walk updates of hyperparameters
# OUTPUTS:
#' @proposals: proposed values for hyperparameters 
proposal.mala.hyper<-function(log.tilde.param, sigma.hyper2){
  means.mala <- log.tilde.param
  variances.mala <- rep(sigma.hyper2,length(log.tilde.param)) 
  proposals <- rnorm(length(log.tilde.param),mean=means.mala,sd=sqrt(variances.mala))
  return(proposals)
}

#' @proposal.mala: MALA-PROPOSAL FOR  HYPER-PARAMETERS: N(log hyper-parametrs + (sigma.latent2^2/2) *Gradient of log(pi), sigma^2 D), D=I
#proposal.mala: function to sample latent parameters log(X2) from the MALA proposal distribution
# INPUTS:
#'   @Y: data matrix (dimension: ntime x nsites), the i-th row gives the i-th temporal replicate at all the sites
#'   @log.X2: current value of the logarithm of latent parameters X2(s) in the denominator (dimension: ntime x nsites)
#' @log.tilde.param: vector of  transfomred hyperparameters 
#'   @ist.mat: distance matrix between locations (dimension: nsites x nsites)
#'   @sigma.latent2: tuning parameter for MALA proposal
#'   @delta: the precision for the numerical derivative (if using the numerical gradient and/or hessian)
#'  @ind1: censored  index (dimension: ntime x nsites)
#'   @u.thr: threshold matrix (dimension: ntime x nsites)
# OUTPUT:
#'   @proposals: proposed values for each latent parameter X2 (dimension: ntime x nsites)

proposal.mala<-function(Y, log.X2, log.tilde.param, dist.mat, sigma.latent2,ind1,u.thr){
  ntime <- nrow(Y)
  nsites <- ncol(Y)
  n <- ntime*nsites
  means.mala <- log.X2+(sigma.latent2/2)*tilde.gradient.log.posterior.theor(Y, log.X2, log.tilde.param, dist.mat,ind1,u.thr)
  variances.mala <- rep(sigma.latent2,n) 
  proposals <- matrix(rnorm(n,mean=as.numeric(means.mala),sd=as.numeric(sqrt(variances.mala))),nrow=ntime,ncol=nsites)
  return(proposals)
}

#' @log.proposal.density.mala: function to calculate the log-density for the MALA proposal
# INPUTS:
#'   @Y: data matrix (dimension: ntime x nsites), the i-th row gives the i-th temporal replicate at all the sites
#'   @log.X2: current value of the logarithm of latent parameters X2(s) in the denominator (dimension: ntime x nsites)
#'  @log.X2.star: candidate value of the logarithm of latent parameters X2(s) in the denominator (dimension: ntime x nsites)
#'  @param: vector of hyperparameters (alpha,beta1,beta2,rho)
#'  @dist.mat: distance matrix between locations (dimension: nsites x nsites)
#'  @sigma.latent2: tuning parameter for MALA proposal
#'  @delta: the precision for the numerical derivative (if using the numerical gradient and/or hessian)
# OUTPUT:
#' @densities.mala:  value of the MALA log-density
log.proposal.density.mala<-function(Y, log.X2, log.X2.star, log.tilde.param, dist.mat, sigma.latent2,ind1,u.thr){
  ntime <- nrow(Y)
  nsites <- ncol(Y)
  n <- ntime*nsites
  means.mala <- log.X2+(sigma.latent2/2)*tilde.gradient.log.posterior.theor(Y, log.X2, log.tilde.param, dist.mat,ind1,u.thr)
  variances.mala <- rep(sigma.latent2,n) 
  densities.mala <- sum(matrix(dnorm(x=log.X2.star,mean=as.numeric(means.mala),sd=as.numeric(sqrt(variances.mala)),log=TRUE),nrow=ntime,ncol=nsites))
  return(densities.mala)
}

### To check above functions
# sigma.latent2 <- 0.001
# delta <- 10^(-4)
# log.X2 <- log(X2)
# thr.prob <- 0.75 ### NOT USED IN THE FUNCTION...
# log.X2.star <- proposal.mala(Y, log.X2, param, dist.mat, sigma.latent2,ind1,u.thr)
# log.proposal.density.mala(Y=Y, log.X2=log.X2, log.X2.star=log.X2.star, log.tilde.param, dist.mat, sigma.latent2,ind1,u.thr)
# log.proposal.density.mala(Y=Y, log.X2=log.X2.star, log.X2.star=log.X2, log.tilde.param, dist.mat, sigma.latent2,ind1,u.thr)

######################
### MCMC ALGORITHM ###
######################
# INITIAL VALUES FOR FIRST CHAIN
alpha0.init1 <- 0.1 #gamma RATE for X2 (denominator) (i.e., controls the overall scale of the data)
alpha1.init1 <- 0
alpha2.init1 <- 0
alpha3.init1 <- 0
beta1.init1 <- 2 
beta2.init1 <- 2
rho.init1<-0.5 
# tilde parametrization
alpha0.tilde.init1<-alpha0.init1+log(beta1.init1)-log(beta2.init1)
alpha1.tilde.init1<-alpha1.init1
alpha2.tilde.init1<-alpha2.init1
alpha3.tilde.init1<-alpha3.init1
beta1.tilde.init1<-alpha0.tilde.init1+log(beta1.init1)
beta2.tilde.init1<-log(1/beta2.init1)
rho.tilde.init1<-log(rho.init1)
tilde.param.init1 <- c(alpha0.tilde.init1,alpha1.tilde.init1,alpha2.tilde.init1,alpha3.tilde.init1,beta1.tilde.init1,beta2.tilde.init1,rho.tilde.init1) # "random" vector of all initial hyperparameters

sim.data.init<- function(ntime,loc,log.tilde.param){
  nsites <- nrow(loc) #number of sites
  n <- nsites*ntime #total number of observations in space and time
  param <- c((2*log.tilde.param[1]-log.tilde.param[5]-log.tilde.param[6]),log.tilde.param[2],log.tilde.param[3],log.tilde.param[4],exp(log.tilde.param[5]-log.tilde.param[1]),exp(-log.tilde.param[6]),exp(log.tilde.param[7]))
  alpha.cov<-exp(param[1]+param[2]*Z1+param[3]*Z2+param[4]*Z3)
  beta1<-param[5]
  beta2 <- param[6]
  rho<-param[7]
  
  alpha <- matrix(rep(alpha.cov,ntime),nrow = ntime,ncol = nsites,byrow = TRUE) 
  dist.mat <- as.matrix(dist(loc)) #distance matrix of all locations (dimension: nsites x nsites)
  Sigma <- exp(-dist.mat/rho) #correlation matrix for Gaussian copula
  
  X1 <- matrix(rgamma(n,shape=beta1,rate=1),nrow=ntime,ncol=nsites) #process X1(s) in numerator as independent gamma random variables (dimension: ntime x nsites)
  Gauss <- mvtnorm::rmvnorm(n=ntime,mean=rep(0,nsites),Sigma) #simulation of ntime independent multivariate Gaussian vectors (dimension: ntime x nsites)
  X2 <- qgamma(pnorm(Gauss),shape=beta2,rate=alpha) #process X2(s) in denominator as independent vectors with gamma margins and Gaussian copula (dimension: ntime x nsites)
  
  Y <- X1/X2 #data simulation as gamma ratio
  return(list(Y=Y,X1=X1,X2=X2))
}
set.seed(123)
data.simulated.init1 <- sim.data.init(ntime=ntime,loc=loc,log.tilde.param=tilde.param.init1)
X2.init1 <- data.simulated.init1$X2 #process X2(s) in denominator as independent vectors with gamma margins and Gaussian copula (dimension: ntime x nsites)
log.X2.init1 <- log(X2.init1)
init1 <- c(log.X2.init1,tilde.param.init1) ## Initial values for MCMC: log(X2) and log(hyperparameters)

# INITIAL VALUES FOR SECOND CHAIN
alpha0.init2 <- 1.5 #gamma RATE for X2 (denominator) (i.e., controls the overall scale of the data)
alpha1.init2 <- 2
alpha2.init2 <- 2
alpha3.init2 <- 2
beta1.init2<- 8 #gamma SHAPE for X1 (numerator) (i.e., controls departure from GP)
beta2.init2 <- 8#gamma SHAPE for X2 (denominator), related to the tail index as xi=1/beta2
rho.init2<-3 #RANGE parameter in exponential correlation function exp(-h/rho), h>=0, for the Gaussian copula in X2 (denominator)

# tilde parametrization
alpha0.tilde.init2<-alpha0.init2+log(beta1.init2)-log(beta2.init2)
alpha1.tilde.init2<-alpha1.init2
alpha2.tilde.init2<-alpha2.init2
alpha3.tilde.init2<-alpha3.init2
beta1.tilde.init2<-alpha0.tilde.init2+log(beta1.init2)
beta2.tilde.init2<-log(1/beta2.init2)
rho.tilde.init2<-log(rho.init2)
tilde.param.init2 <- c(alpha0.tilde.init2,alpha1.tilde.init2,alpha2.tilde.init2,alpha3.tilde.init2,beta1.tilde.init2,beta2.tilde.init2,rho.tilde.init2) # "random" vector of all initial hyperparameters

set.seed(123)
data.simulated.init2 <- sim.data.init(ntime=ntime,loc=loc,log.tilde.param=tilde.param.init2)
X2.init2 <- data.simulated.init2$X2 
log.X2.init2 <- log(X2.init2)

init2 <- c(log.X2.init2,tilde.param.init2) ## Initial values for MCMC: log(X2) and log(hyperparameters)

#MAIN FUNCTION
#' @MCMC.gamma.gamma: function to run the MCMC algorithm for the Gamma ratio model
# INPUTS:
#'   @N.MCMC: integer, number of MCMC iterations
#'   @Y: data matrix (dimension: ntime x nsites), the i-th row gives the i-th temporal replicate at all the sites
#'  @init: vector of length (ntime x nsites + length of hypeparameters) containing the initial values for log(X2) (latent parameters) and log(hyperparameters)
#'  @dist.mat: distance matrix between locations (dimension: nsites x nsites)
#'  @sigma.latent2: tuning parameter for MALA proposal
#'  @by: store (and plot) samples every 'by' iterations
#'  @adapt: the number of run after which the adaptation in tuning parameytrs takes place
#'  @theta.hyper : the fixed parameter for hyperparameters in the adaptation strategy
#'   @theta.latent : the fixed parameter for latent parameters in the adaptation strategy
#'   @burn_in1   : the first burning period 
#'   @burn_in2  : the second burning period 
#'  @ind1: censored  index (dimension: ntime x nsites)
#'  @ind2: censored  index (dimension: ntime x nsites)
#'  @u.thr: threshold matrix (dimension: ntime x nsites)
# OUTPUT:
#'  @samples:  matrix (dimension: N.MCMC x (ntime x nsites + 4)), containing the Markov chains for each latent variable X2 and hyperparameters
#'  @tuning_param_x_hyper: matrix  with two coloumns, where the first row is the tuning parameters values for latent parameter and seond row is for hyperparameter vector
#'  @Acc.rate-hyper: acc rate for the hyperparameters vector blocks
#'  @Acc.rate.latent: acc rate for the latent vector blocks

start_time<-Sys.time()
MCMC.gamma.gamma=function(N.MCMC, Y, init, dist.mat, sigma.hyper2, sigma.latent2, by, adapt, burn_in1, burn_in2, theta.hyper, theta.latent, ind1, ind2, u.thr)
{
  samples <- matrix(nrow=floor(N.MCMC/by),ncol=length(init),byrow=TRUE)
  sigma.matrix<-matrix(nrow=floor(N.MCMC/adapt),ncol=2,byrow=TRUE)
  sigma.matrix[1,]<-c(sigma.latent2,sigma.hyper2)
  samples[1,]<-init
  cur.samples.hyper <- init[(n+1):(n+7)]
  cur.samples.latent<-init[1:n]
  
  ntime <- nrow(Y)
  nsites <- ncol(Y)
  n <- ntime*nsites
  j<-1
  l<-1
  m<-1
  k<-1
  for (i in 1:(N.MCMC-1)){
    if((i%%by)-1==0){
      par(mfrow=c(2,5),oma=c(0,0,2,0),mar=c(4,5,1,1))
      if((i%%(adapt))-1==0 & i< (burn_in1+burn_in2+2)){ #to calculate the accpetance rate based on only current samples, burning+2 to calculate the accpertance rate after the burning samples
        rate.latent<-0
        rate.hyper<-0
      }
      if(i< (burn_in1+burn_in2+2)){
        print(paste("Iteration: ",i,"; Acc rate hyp=",rate.hyper/((i%%(adapt))+1),"; Acc rate lat=",rate.latent/((i%%(adapt))+1),"; sigma hyper=",sigma.hyper2,"; sigma latent=",sigma.latent2,sep=""))
      } else{
        print(paste("Iteration: ",i,"; Acceptance rate hyper=",rate.hyper/(i-(burn_in1+burn_in2+2)),"; Acceptance rate latent=",rate.latent/(i-(burn_in1+burn_in2+2)),";sigma hyper=",sigma.hyper2,"; sigma latent=",sigma.latent2,sep=""))
      }
      
      #ploteps <- 1
      #ploteps <- c(-ploteps,ploteps)
      plot(by*c(0:(l-1))+1,samples[1:l,n+1],type = "l",xlab="MCMC iteration",ylab=expression(tilde(alpha)[0])) # Plot for alpha.tilde
      abline(h=alpha0.tilde,col="4")
      
      #ploteps <- 1
      #ploteps <- c(-ploteps,ploteps)
      plot(by*c(0:(l-1))+1,samples[1:l,n+2],type = "l",xlab="MCMC iteration",ylab=expression(tilde(alpha)[1])) # Plot for beta1.tilde
      abline(h=alpha1.tilde,col="4")
      
      #ploteps <- 2
      #ploteps <- c(-ploteps,ploteps)
      plot(by*c(0:(l-1))+1,samples[1:l,n+3],type = "l",xlab="MCMC iteration",ylab=expression(tilde(alpha)[2])) # Plot for beta2
      abline(h=alpha2.tilde,col="4")
      
      plot(by*c(0:(l-1))+1,samples[1:l,n+4],type = "l",xlab="MCMC iteration",ylab=expression(tilde(alpha)[3])) # Plot for beta2
      abline(h=alpha3.tilde,col="4")
      
      #ploteps <- 0.3
      #ploteps <- c(-ploteps,ploteps)
      plot(by*c(0:(l-1))+1,samples[1:l,n+5],type = "l",xlab="MCMC iteration",ylab=expression(tilde(beta)[1])) # Plot for rho
      abline(h=beta1.tilde,col="4")
      
      #ploteps <- 0.3
      #ploteps <- c(-ploteps,ploteps)
      plot(by*c(0:(l-1))+1,samples[1:l,n+6],type = "l",xlab="MCMC iteration",ylab=expression(tilde(beta2))) # Plot for rho
      abline(h=beta2.tilde,col="4")
      
      plot(by*c(0:(l-1))+1,samples[1:l,n+7],type = "l",xlab="MCMC iteration",ylab=expression(tilde(rho))) # Plot for rho
      abline(h=rho.tilde,col="4")
      
      
      plot(by*c(0:(l-1))+1,samples[1:l,1],type = "l",xlab="MCMC iteration",ylab=expression(log(X[1]))) # Plot for log(X2[1])
      abline(h=as.vector(log(X2))[1],col="4")
      
      plot(by*c(0:(l-1))+1,samples[1:l,2],type = "l",xlab="MCMC iteration",ylab=expression(log(X[2]))) # Plot for log(X2[1])
      abline(h=as.vector(log(X2))[2],col="4")
      
      plot(by*c(0:(l-1))+1,samples[1:l,3],type = "l",xlab="MCMC iteration",ylab=expression(log(X[3]))) # Plot for log(X2[1])
      abline(h=as.vector(log(X2))[3],col="4")
      
      
      
      # plot(by*c(0:(l-1))+1,samples[1:l,500],type = "l",xlab="MCMC iteration",ylab=expression(log(X2[500]))) # Plot for log(X2[1])
      # abline(h=as.vector(log(X2))[500],col="4")
      # 
      # plot(by*c(0:(l-1))+1,samples[1:l,1000],type = "l",xlab="MCMC iteration",ylab=expression(log(X2[1000]))) # Plot for log(X2[1])
      # abline(h=as.vector(log(X2))[1000],col="4")
      if((i%%adapt)-1==0){
        #plot(sigma.matrix[1:k,1],sigma.matrix[1:k,2],xlab="sigma.latent",ylab="sigma.hyper2") #plot for the scale parameters chosen adaptively
        k<-k+1
      }
      
      l<-l+1
    }
    
    ### Extracting current parameters (log of latent parameters X2(s) + hyperparameters)
    # cur <- cur.samples # Current sample of the Markov chain
    cur.log.tilde.param <- cur.samples.hyper # Current hypeparameters (on log scale)
    cur.log.X2<-matrix(cur.samples.latent, nrow = ntime, ncol = nsites) # Current latent parameters X2(s) (on log scale)
    
    ### Proposing new parameters
    prop.log.tilde.param <- proposal.mala.hyper(log.tilde.param=cur.log.tilde.param,sigma.hyper2 = sigma.hyper2)  # Proposed hyperparameters using uniform random walks
    
    ### Terms in acceptance ratio (using tilde parametrization)
    cur.log.posterior.hyper <- tilde.log.posterior(Y=Y,log.X2=cur.log.X2,log.tilde.param=cur.log.tilde.param,dist.mat=dist.mat,ind1=ind1,ind2=ind2,u.thr=u.thr)
    prop.log.posterior.hyper <- tilde.log.posterior(Y=Y,log.X2=cur.log.X2,log.tilde.param=prop.log.tilde.param,dist.mat=dist.mat,ind1=ind1,ind2=ind2,u.thr=u.thr)
    
    log.num.hyper <- prop.log.posterior.hyper
    log.denom.hyper<- cur.log.posterior.hyper
    log.ratio.hyper <- log.num.hyper-log.denom.hyper
    
    
    if (log(runif(1)) < log.ratio.hyper & is.na(log.ratio.hyper)==FALSE){
      cur.samples.hyper <-prop.log.tilde.param
      rate.hyper<-rate.hyper+1
    } else{
      cur.samples.hyper<- cur.log.tilde.param
    }
    
    
    if (i<(burn_in1+2)){
      if(i%%adapt==0 & i%%(2*adapt)!=0){
        sigma.hyper2<-exp(((rate.hyper)/(adapt)-0.23)/theta.hyper)*sigma.hyper2
      }
    }
    if ((burn_in1+1) < i  & i< (burn_in1+burn_in2+2)){
      if ((rate.hyper/(adapt))>0.30 | (rate.hyper/(adapt))<0.15){ #burning+2 to calculate the accpertance rate after the burning samples
        if(i%%adapt==0 & i%%(2*adapt)!=0){
          sigma.hyper2<-exp(((rate.hyper)/(adapt)-0.23)/theta.hyper)*sigma.hyper2
        }
      } else {
        sigma.hyper2<-sigma.hyper2
      }
    }
    
    ### Terms in acceptance ratio (using tilde parametrization)
    
    prop.log.X2 <- proposal.mala(Y=Y,log.X2=cur.log.X2,log.tilde.param=cur.samples.hyper,dist.mat=dist.mat,sigma.latent2=sigma.latent2,ind1=ind1,u.thr=u.thr) # Proposed latent parameters X2(s) using the MALA
    
    cur.log.posterior.latent<- tilde.log.posterior(Y=Y,log.X2=cur.log.X2,log.tilde.param=cur.samples.hyper,dist.mat=dist.mat,ind1=ind1,ind2=ind2,u.thr=u.thr)
    prop.log.posterior.latent <- tilde.log.posterior(Y=Y,log.X2=prop.log.X2,log.tilde.param=cur.samples.hyper,dist.mat=dist.mat,ind1=ind1,ind2=ind2,u.thr=u.thr)
    cur.log.proposal.mala.latent <- log.proposal.density.mala(Y=Y,log.X2=prop.log.X2,log.X2.star=cur.log.X2,log.tilde.param=cur.samples.hyper,dist.mat=dist.mat,sigma.latent2=sigma.latent2,ind1=ind1,u.thr=u.thr)
    prop.log.proposal.mala.latent <- log.proposal.density.mala(Y=Y,log.X2=cur.log.X2,log.X2.star=prop.log.X2,log.tilde.param=cur.samples.hyper,dist.mat=dist.mat,sigma.latent2=sigma.latent2,ind1=ind1,u.thr=u.thr)
    
    log.num.latent <- prop.log.posterior.latent+ cur.log.proposal.mala.latent
    log.denom.latent<- cur.log.posterior.latent + prop.log.proposal.mala.latent
    log.ratio.latent <- log.num.latent-log.denom.latent
    
    if (log(runif(1)) < log.ratio.latent & is.na(log.ratio.latent)==FALSE){
      cur.samples.latent <- c(as.vector(prop.log.X2))
      rate.latent<-rate.latent+1
    } else{
      cur.samples.latent<- cur.log.X2
    }
    
    if((i%%by)-1==0){
      samples[j,]<- c(cur.samples.latent,cur.samples.hyper)
      j=j+1
    }
    # #ADAPTIVE ALGORITHM 3 (not conditioning in the range (0.50,0.65))
    if (i<(burn_in1+2)){
      if(i%%(2*adapt)==0){
        sigma.latent2<-exp(((rate.latent)/(adapt)-0.57)/theta.latent)*sigma.latent2
      }
    }
    if ((burn_in1+1) < i  & i< (burn_in1+burn_in2+2)){
      if ((rate.latent/(adapt))>0.65 | (rate.latent/(adapt))<0.50){ #burning+2 to calculate the accpertance rate after the burning samples
        if(i%%(2*adapt)==0){
          sigma.latent2<-exp(((rate.latent)/(adapt)-0.57)/theta.latent)*sigma.latent2
        }
      } else {
        sigma.latent2<-sigma.latent2
      }
    }
    
    if((i%%adapt)-1==0){ # to save allexp the scale parameter of the MALA
      sigma.matrix[m,]<- c(sigma.latent2,sigma.hyper2)
      m=m+1
    }
  }
  return(list("samples"=samples,"tuning_param_x_hyper"=sigma.matrix,"Acc.rate hyper"=rate.hyper/(N.MCMC-(burn_in1+burn_in2)),"Acc.rate.latent"=rate.latent/(N.MCMC-(burn_in1+burn_in2))))
}

#TUNING PARAMETERS
sigma.hyper2<-0.001 # tuning parameter for mean/variance of the MALA proposals(hyperparameter)
sigma.latent2<-0.0001 # tuning parameter for mean/variance of the MALA proposals
#MCMC ITERATIONS
N.MCMC<-1.5*10^6
by<-250  #  the thinning samples
adapt<-1000 ## number of samples aftyer which we upadte the tuning parameters
theta.hyper<-0.4 
theta.latent<-0.4 

burn_in1<-2.5*10^5
burn_in2<-5*10^5


# RUNNING FIRST CHAIN
set.seed(2)
start_time1<-Sys.time()
MCMC1<- MCMC.gamma.gamma(N.MCMC=N.MCMC,Y=Y,init=init1,dist.mat=dist.mat,sigma.hyper2=sigma.hyper2,sigma.latent2=sigma.latent2,by=by,adapt=adapt,burn_in1 = burn_in1,burn_in2=burn_in2,theta.hyper=theta.hyper,theta.latent = theta.latent,ind1=ind1,ind2=ind2,u.thr=u.thr)
end_time1<-Sys.time()
MCMC.out1<-MCMC1$samples
MCMC_tuning_param_x_chain1<-MCMC1$tuning_param_x_hyper
MCMC_Acc.rate_chain1<-MCMC1$Acc.rate
#repoerting runtime for chain2
time1<-print(end_time1-start_time1)
#save.image("Mydata_chain1.RData")

# RUNNING SECOND CHAIN
set.seed(3)
start_time2<-Sys.time()
MCMC2<- MCMC.gamma.gamma(N.MCMC=N.MCMC,Y=Y,init=init2,dist.mat=dist.mat,sigma.hyper2=sigma.hyper2,sigma.latent2=sigma.latent2,by=by,adapt=adapt,burn_in1 = burn_in1,burn_in2=burn_in2,theta.hyper=theta.hyper,theta.latent = theta.latent,ind1=ind1,ind2=ind2,u.thr=u.thr)
end_time2<-Sys.time()
MCMC.out2<-MCMC2$samples
MCMC_tuning_param_x_chain2<-MCMC2$tuning_param_x_hyper
MCMC_Acc.rate_chain2<-MCMC2$Acc.rate
#repoerting runtime for chain2
time2<-print(end_time2-start_time2)
#save.image("Mydata_chain2.RData")
#save.image("BothChainMoreIterations.RData")

