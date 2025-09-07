library("boot", include.only = "inv.logit")
library(mixtools)
library(coda)
library(dplyr) 
library(MCMCpack,include.only = "riwish")
library(DirichletReg, include.only =  "rdirichlet")
library(mixtools)
library(mvtnorm, include.only =  "dmvnorm")
library(MASS, include.only = "mvrnorm")
library(scales, include.only = "alpha")
library(invgamma)
library(truncdist)

############################################################################
## Because this method requires extensive Gibbs sampling , these simulations are 
## very computationally Expensive and were originally run on a HPC
#############################################################################


setwd(".../Sim_Files")

iterations = 250 #number of times model is run on new data

truth_runs = 10000 # number of simulations used to obtain true effect values

n.imps = 1000 # number of used samples from gibbs sampler per iteration
thin = 5 # amount of thinning
burn=200 # burn in period
Sims <- burn + thin*n.imps # Set number of MCMC draws taken per simulation
samps_use <- burn + 1:n.imps*thin #Sample index to collect in MCMC

# Read in data and obtain standardized covariates for simulation
Covariates<- read.csv("CC_Sim_Covariates.csv")


### Make Data ###
J=60 # number of providers in simulation
Nj = 20 # number of people per providers
K = 2 # number of latent compliance strata
pi_k <- c(0.5,0.5) # probability of observing a facility in each strata
P = 1 # number of facility covariates used in mixture
L=1 # number of facility compliance measures using in mixture
YCN = ncol(Covariates) # number of response covariates in covariate matrix
DCN = ncol(Covariates)# number of individual compliance covariates 
N = Nj*J # total number of individuals generated for simulation study




# Coverage, Bias and Interval width for effect estimates
ATECover <- rep(NA,iterations)
ATEBias <- rep(NA,iterations)
ATEIW <- rep(NA,iterations)

CACECover <- rep(NA,iterations)
CACEBias <- rep(NA,iterations)
CACEIW <- rep(NA,iterations)

ATE1Cover <- rep(NA,iterations)
ATE1Bias <- rep(NA,iterations)
ATE1IW <- rep(NA,iterations)

ATE2Cover <- rep(NA,iterations)
ATE2Bias <- rep(NA,iterations)
ATE2IW <- rep(NA,iterations)

CACE1Cover <- rep(NA,iterations)
CACE1Bias <- rep(NA,iterations)
CACE1IW <- rep(NA,iterations)

CACE2Cover <- rep(NA,iterations)
CACE2Bias <- rep(NA,iterations)
CACE2IW <- rep(NA,iterations)


# Additional metric tracking matricies
Delta0Bias <- matrix(NA, nrow=2,ncol = iterations)
Sigma2YBias <- matrix(NA, nrow=2,ncol = iterations)
Tau2PhiYBias <- matrix(NA, nrow=2,ncol = iterations)
Tau2PhiDBias <- matrix(NA, nrow=2,ncol = iterations)
Beta1Bias <- matrix(NA, nrow=7,ncol = iterations)
Beta2Bias <- matrix(NA, nrow=7,ncol = iterations)
Alpha1Bias <- matrix(NA, nrow=2,ncol = iterations)
Alpha2Bias <- matrix(NA, nrow=2,ncol = iterations)
muPhiDBias <- matrix(NA, nrow=2,ncol = iterations)
SampCover <- rep(NA,iterations)
PhiY_cover <-  matrix(NA, nrow=J,ncol = iterations)
PhiYBias <- matrix(NA,nrow=J,ncol=iterations)
PhiD_cover <-  matrix(NA, nrow=J,ncol = iterations)
PhiDBias <- matrix(NA,nrow=J,ncol=iterations)




# Define skew and degrees of freedom for multivariate skewed t that
# generates facility covariates and compliance 
# (skew=0 and df = Inf implies correctly specified MVN)
skew<-rep(0,K)
mix.df <- Inf


# Decide whether to use Burr link or probit to generate individual 
# compliance probability
# (useBurr = 1 implies Burr link, useBurr =0 implies correctly specified probit)
useBurr=0
burr.link <- function(x,c=0.5){
  return(1-(1+exp(x))^(-1*c))
}

# Decide on interaction effects that are present when generating the outcome data
# 0's imply correctly specified
interaction_coeff_treat <- c(0,0)
interaction_coeff_control <- c(0,0)

# Decide on the effects of compliance on control outcomes 
# D0 <- rep(0,K)
# useD0 <- ifelse(D0 == 0, 0, 1)



# Randomly sample individuals from "population" to be included in study
set.seed(1)
X_full <- Covariates[sample(1:nrow(Covariates),N,replace = T),]
rownames(X_full) <- 1:nrow(X_full)
set.seed(NULL)



# Define vectors an matricies to hold data on providers (facilities) and outcomes
# for simulation of true effects
Facility <- sort(rep(1:J,Nj))
Person <- 1:(N)
Y1 <- rep(NA,length(Person))  
Y0 <- rep(NA,length(Person))
D <- rep(NA,length(Person))
treat_Fac_bin <- ifelse(1:J<=J/2,0,1)
Treat <- rep(treat_Fac_bin,each=Nj) # assignment vector
OutcomeData <- as.data.frame(cbind(Person, Facility, Treat, Y0,Y1,D,Strata=NA)) 
treat_Fac_bool <- treat_Fac_bin ==1
treat_Fac_ind <- which(treat_Fac_bin ==1)
treat_Person_bin <- OutcomeData$Treat
treat_Person_bool <- treat_Person_bin==1
treat_person_ind <- which(treat_Person_bin==1)

#Vectors to hold differences in genreated outcomes that correspond to effects of interest
Diff_ATE2 <- rep(NA,truth_runs)
Diff_CACE1 <- rep(NA,truth_runs)
Diff_ATE1 <- rep(NA,truth_runs)
Diff_CACE2 <- rep(NA,truth_runs)
Diff_CACE <- rep(NA,truth_runs)
Diff_ATE <- rep(NA,truth_runs)

# Simulate outcome data for true effects
for (t in 1:truth_runs) {

  X <- as.matrix(X_full[sample(1:nrow(X_full),N,replace = F),])
  rownames(X) <- 1:nrow(X)
  
  # Creating person by Person-level strata membership
  Facility_Strata <- sample(1:K,J,prob = pi_k,replace = T) #Strata of Facilities (providers) in Simulation
  OutcomeData$Strata <- rep(Facility_Strata ,each=Nj)
  Facility_df <- data.frame(Facility= Facility)
  Facility_Strata <- Facility_Strata
  
  
  
  # Generate true compliance mean and true characteristic mean vectors for facilites
  mu.C <- matrix(NA, nrow = L, ncol = K)
  mu.Z <- matrix(NA, nrow = P, ncol = K)
  
  
  for (k in 1:K) {
    for(l in 1:L){
      mu.C[l,k] <- (1+k*4)
    }
    for(p in 1:P){
      mu.Z[p,k] <- (1+k*4)
    }
  }
  
  mu.C <- mu.C -mean(mu.C)
  mu.Z <- mu.Z -mean(mu.Z)
  
  
  mu <- as.matrix(rbind(mu.C,mu.Z))
  
  #Generate True covariance matrix for compliance and facility covariates
  rho_Sigma <- runif(1,-0.8,0.8)
  Sigma2_C <- runif(1,0.5,2)
  Sigma2_Z <- runif(1,0.5,2)
  
  Sigma <- matrix(c(Sigma2_C, rho_Sigma*sqrt(Sigma2_C*Sigma2_Z),
                    rho_Sigma*sqrt(Sigma2_C*Sigma2_Z),Sigma2_Z), nrow = 2, byrow = T)

  
  # Sample a Matrix of Compliance and Characterisitcs for facilities using true parameters
  True_CZ <- matrix(NA, nrow = J, ncol = L+P)
  
  for (k in 1:K) {
    fac_in_k <- Facility_Strata==k
    True_CZ[fac_in_k,] <-sn::rmst(n=sum(fac_in_k), xi=mu[,k],
                                  Omega=Sigma, alpha=skew, nu=mix.df)
  }
  
  True_C <- as.matrix(True_CZ[,1:L])
  True_Z <- as.matrix(True_CZ[,(L+1):(L+P)])
  # pairs(True_CZ)
  
  
  #############################
  ### Individual Compliance ###
  #############################
  
  # Create covariate matrix for individual compliance 
  X_D <- as.matrix(cbind(1,X))
  
  # Create true coefficients for individual compliance measure
  Coeff_Alpha <- matrix(NA, nrow = ncol(X_D), ncol=K)
  # mu_PhiD <- rep(NA,K)
  PhiD_Fac <- rep(NA,J)
  
  #Create variance for facility effect (used for compliance data)
  Tau2_PhiD <- 0.25
  
  for (k in 1:K) {
    Coeff_Alpha[,k] <- c((k-1)/2,rep(-0.25*k,ncol(X_D)-1))
    # mu_PhiD[k] <- k/2
    PhiD_Fac[Facility_Strata==k] <- rnorm(sum(Facility_Strata==k),0 ,sqrt(Tau2_PhiD))
  }
  
  PhiD_Person <- PhiD_Fac[OutcomeData$Facility]
  OutcomeData$PhiD <- PhiD_Person 
  
  #Generate true sample individual compliance probababilities and observed compliance
  if(useBurr==1){
    D_Probs <- burr.link(diag(X_D%*%Coeff_Alpha[,OutcomeData$Strata])+PhiD_Person) 
  }else{ D_Probs <-pnorm(diag(X_D%*%Coeff_Alpha[,OutcomeData$Strata])+PhiD_Person)}
  OutcomeData$D <- rbinom(length(D_Probs),1, D_Probs)
  
  ###########################
  ### Outcome of interest ###
  ###########################
  
  # Create response matricies for outcome variable
  
  X_Main_obs <- cbind(1,as.matrix(X))
  X_Main_obs <- cbind(X_Main_obs,(1- OutcomeData$Treat)*OutcomeData$D*1)
  X_Main_obs <- cbind(X_Main_obs,(OutcomeData$Treat)*OutcomeData$D* cbind(as.matrix(X),1))
  
  X_Main_treat <- cbind(1,as.matrix(X))
  X_Main_treat <- cbind(X_Main_treat,(1- 1)*OutcomeData$D*1)
  X_Main_treat <- cbind(X_Main_treat,(1)*OutcomeData$D* cbind(as.matrix(X),1))
  
  X_Main_control <- cbind(1,as.matrix(X))
  X_Main_control <- cbind(X_Main_control,(1- 0)*OutcomeData$D*1)
  X_Main_control <- cbind(X_Main_control,(0)*OutcomeData$D* cbind(as.matrix(X),1))
  
  X_Main_control_comp <- cbind(1,as.matrix(X))
  X_Main_control_comp <- cbind(X_Main_control_comp,(1- 0)*1*1)
  X_Main_control_comp <- cbind(X_Main_control_comp,(0)*1* cbind(as.matrix(X),1))
  
  X_Main_control_nocomp <- cbind(1,as.matrix(X))
  X_Main_control_nocomp <- cbind(X_Main_control_nocomp,(1- 0)*0*1)
  X_Main_control_nocomp <- cbind(X_Main_control_nocomp,(0)*0* cbind(as.matrix(X),1))
  
  
  #Create coefficients for outcome, and mean for facility effect
  Coeff_BetaD0 <- matrix(NA, nrow = YCN,ncol = K)
  Coeff_BetaD1 <- matrix(NA, nrow = YCN,ncol = K)
  Coeff_D0 <- rep(NA,K)
  Coeff_D1 <- rep(NA,K)
  Coeff_muY <- rep(NA,K)
  PhiY_Fac <- rep(NA,J)
  
  # Create variance for facility effect (used for outcome data)
  Tau2_PhiY <- 9
  
  for (k in 1:K) {
    Coeff_BetaD0[,k] <-rep(1*k,ncol(X))
    Coeff_BetaD1[,k] <-rep(1*k,ncol(X))
    Coeff_D0[k] <- k
    Coeff_D1[k] <- k + 4.5 + (k-1)
    Coeff_muY[k] <- 2*k
    PhiY_Fac[Facility_Strata==k] <- rnorm(sum(Facility_Strata==k), 0,sqrt(Tau2_PhiY))
  }
  
  # Coeff_FullY  <- rbind(Coeff_muY,Coeff_BetaD0, Coeff_BetaD1,  Coeff_D1)
  Coeff_FullY  <- rbind(Coeff_muY,Coeff_BetaD0,Coeff_D0, Coeff_BetaD1,  Coeff_D1)
  
  PhiY_Person <- PhiY_Fac[OutcomeData$Facility]
  OutcomeData$PhiY <- PhiY_Person
  
  #Create variance for outcome data
  Sigma2Y <- c(16,16)
  
  
  
  #Sample true outcome data
  Y0Means <- rep(NA,N)
  Y1Means <- rep(NA,N)
  
  Y0Means <- diag(X_Main_control%*%Coeff_FullY[,OutcomeData$Strata])+PhiY_Person +
    diag(X_Main_control[,2]*X_Main_control[,3])%*% interaction_coeff_control[OutcomeData$Strata]
  Y1Means <- diag(X_Main_treat%*%Coeff_FullY[,OutcomeData$Strata])+PhiY_Person+
    diag(X_Main_control[,2]*X_Main_control[,3])%*% interaction_coeff_control[OutcomeData$Strata]+
    OutcomeData$D*diag(X_Main_treat[,2]*X_Main_treat[,3])%*% interaction_coeff_treat[OutcomeData$Strata]
  
  OutcomeData$Y0 <- rnorm(length(Y0Means), Y0Means,sqrt(Sigma2Y[OutcomeData$Strata]))
  OutcomeData$Y1 <- rnorm(length(Y1Means), Y1Means,sqrt(Sigma2Y[OutcomeData$Strata]))
  
  Obs_CACE <- mean(OutcomeData$Y1[OutcomeData$D==1]-OutcomeData$Y0[OutcomeData$D==1])
  Obs_ATE <- mean(OutcomeData$Y1-OutcomeData$Y0)
  Obs_CACE1 <- mean(OutcomeData$Y1[OutcomeData$D==1&OutcomeData$Strata==1]-OutcomeData$Y0[OutcomeData$D==1&OutcomeData$Strata==1])
  Obs_CACE2 <- mean(OutcomeData$Y1[OutcomeData$D==1&OutcomeData$Strata==2]-OutcomeData$Y0[OutcomeData$D==1&OutcomeData$Strata==2])
  Obs_ATE1 <- mean(OutcomeData$Y1[OutcomeData$Strata==1]-OutcomeData$Y0[OutcomeData$Strata==1])
  Obs_ATE2 <- mean(OutcomeData$Y1[OutcomeData$Strata==2]-OutcomeData$Y0[OutcomeData$Strata==2])
  
  Diff_CACE[t] <- Obs_CACE
  Diff_CACE1[t] <- Obs_CACE1
  Diff_CACE2[t] <- Obs_CACE2
  Diff_ATE[t] <- Obs_ATE
  Diff_ATE1[t] <- Obs_ATE1
  Diff_ATE2[t] <- Obs_ATE2
  
  if(t%%1000 ==0){print(t)}
  
}

# Set True effect values
TrueCACE <- mean(Diff_CACE)
TrueCACE1 <- mean(Diff_CACE1)
TrueCACE2 <- mean(Diff_CACE2)
TrueATE <- mean(Diff_ATE)
TrueATE1 <- mean(Diff_ATE1)
TrueATE2 <- mean(Diff_ATE2)


# Re-Define vectors an matricies to hold data on providers (facilities) and outcomes
# for simulation study of method
Facility <- sort(rep(1:J,Nj))
Person <- 1:(N)
Y1 <- rep(NA,length(Person))  
Y0 <- rep(NA,length(Person))
D <- rep(NA,length(Person))
treat_Fac_bin <- ifelse(1:J<=J/2,0,1)
Treat <- rep(treat_Fac_bin,each=Nj) # assignment vector
OutcomeData <- as.data.frame(cbind(Person, Facility, Treat, Y0,Y1,D,Strata=NA)) 



for (m in 1:iterations) {

  # Set random seed based on time and sample individuals for simulation
  set.seed(as.numeric(format(Sys.time(), "%OS3"))*1000)
  X <- as.matrix(X_full[sample(1:nrow(X_full),N,replace = F),])
  rownames(X) <- 1:nrow(X)
  
  # Creating person by Person-level strata membership
  Facility_Strata <- sample(1:K,J,prob = pi_k,replace = T) #Strata of Facilities (providers) in Simulation
  OutcomeData$Strata <- rep(Facility_Strata ,each=Nj)
  Facility_df <- data.frame(Facility= Facility)
  Facility_Strata <- Facility_Strata
  
  
  
  # Generate true compliance mean and true characteristic mean vectors for facilites
  mu.C <- matrix(NA, nrow = L, ncol = K)
  mu.Z <- matrix(NA, nrow = P, ncol = K)
  
  
  for (k in 1:K) {
    for(l in 1:L){
      mu.C[l,k] <- (1+k*4)
    }
    for(p in 1:P){
      mu.Z[p,k] <- (1+k*4)
    }
  }
  
  mu.C <- mu.C -mean(mu.C)
  mu.Z <- mu.Z -mean(mu.Z)
  
  
  mu <- as.matrix(rbind(mu.C,mu.Z))
  
  #Generate True covariance matrix for compliance and facility covariates
  rho_Sigma <- runif(1,-0.8,0.8)
  Sigma2_C <- runif(1,0.5,2)
  Sigma2_Z <- runif(1,0.5,2)
  
  Sigma <- matrix(c(Sigma2_C, rho_Sigma*sqrt(Sigma2_C*Sigma2_Z),
                    rho_Sigma*sqrt(Sigma2_C*Sigma2_Z),Sigma2_Z), nrow = 2, byrow = T)
  Sigma
  
  
  # Sample a Matrix of Compliance and Characterisitcs for facilities using true parameters
  True_CZ <- matrix(NA, nrow = J, ncol = L+P)
  
  for (k in 1:K) {
    fac_in_k <- Facility_Strata==k
    True_CZ[fac_in_k,] <-sn::rmst(n=sum(fac_in_k), xi=mu[,k],
                                  Omega=Sigma, alpha=skew, nu=mix.df)
  }
  
  True_C <- as.matrix(True_CZ[,1:L])
  True_Z <- as.matrix(True_CZ[,(L+1):(L+P)])
  # pairs(True_CZ)
  
  
  #############################
  ### Individual Compliance ###
  #############################
  
  # Create covariate matrix for individual compliance 
  X_D <- as.matrix(cbind(1,X))
  
  # Create true coefficients for individual compliance measure
  Coeff_Alpha <- matrix(NA, nrow = ncol(X_D), ncol=K)
  # mu_PhiD <- rep(NA,K)
  PhiD_Fac <- rep(NA,J)
  
  #Create variance for facility effect (used for compliance data)
  Tau2_PhiD <- 0.25
  
  for (k in 1:K) {
    Coeff_Alpha[,k] <- c((k-1)/2,rep(-0.25*k,ncol(X_D)-1))
    # mu_PhiD[k] <- k/2
    PhiD_Fac[Facility_Strata==k] <- rnorm(sum(Facility_Strata==k),0 ,sqrt(Tau2_PhiD))
  }
  
  PhiD_Person <- PhiD_Fac[OutcomeData$Facility]
  OutcomeData$PhiD <- PhiD_Person 
  
  #Generate true sample individual compliance probababilities and observed compliance
  if(useBurr==1){
    D_Probs <- burr.link(diag(X_D%*%Coeff_Alpha[,OutcomeData$Strata])+PhiD_Person) 
  }else{ D_Probs <-pnorm(diag(X_D%*%Coeff_Alpha[,OutcomeData$Strata])+PhiD_Person)}
  OutcomeData$D <- rbinom(length(D_Probs),1, D_Probs)
  
  ###########################
  ### Outcome of interest ###
  ###########################
  
  # Create response matricies for outcome variable
  
  X_Main_obs <- cbind(1,as.matrix(X))
  X_Main_obs <- cbind(X_Main_obs,(1- OutcomeData$Treat)*OutcomeData$D*1)
  X_Main_obs <- cbind(X_Main_obs,(OutcomeData$Treat)*OutcomeData$D* cbind(as.matrix(X),1))
  
  X_Main_treat <- cbind(1,as.matrix(X))
  X_Main_treat <- cbind(X_Main_treat,(1- 1)*OutcomeData$D*1)
  X_Main_treat <- cbind(X_Main_treat,(1)*OutcomeData$D* cbind(as.matrix(X),1))
  
  X_Main_control <- cbind(1,as.matrix(X))
  X_Main_control <- cbind(X_Main_control,(1- 0)*OutcomeData$D*1)
  X_Main_control <- cbind(X_Main_control,(0)*OutcomeData$D* cbind(as.matrix(X),1))
  
  X_Main_control_comp <- cbind(1,as.matrix(X))
  X_Main_control_comp <- cbind(X_Main_control_comp,(1- 0)*1*1)
  X_Main_control_comp <- cbind(X_Main_control_comp,(0)*1* cbind(as.matrix(X),1))
  
  X_Main_control_nocomp <- cbind(1,as.matrix(X))
  X_Main_control_nocomp <- cbind(X_Main_control_nocomp,(1- 0)*0*1)
  X_Main_control_nocomp <- cbind(X_Main_control_nocomp,(0)*0* cbind(as.matrix(X),1))
  
  
  #Create coefficients for outcome, and mean for facility effect
  Coeff_BetaD0 <- matrix(NA, nrow = YCN,ncol = K)
  Coeff_BetaD1 <- matrix(NA, nrow = YCN,ncol = K)
  Coeff_D0 <- rep(NA,K)
  Coeff_D1 <- rep(NA,K)
  Coeff_muY <- rep(NA,K)
  PhiY_Fac <- rep(NA,J)
  
  # Create variance for facility effect (used for outcome data)
  Tau2_PhiY <- 9
  
  for (k in 1:K) {
    Coeff_BetaD0[,k] <-rep(1*k,ncol(X))
    Coeff_BetaD1[,k] <-rep(1*k,ncol(X))
    Coeff_D0[k] <- k
    Coeff_D1[k] <- k + 4.5 + (k-1)
    Coeff_muY[k] <- 2*k
    PhiY_Fac[Facility_Strata==k] <- rnorm(sum(Facility_Strata==k), 0,sqrt(Tau2_PhiY))
  }
  
  # Coeff_FullY  <- rbind(Coeff_muY,Coeff_BetaD0, Coeff_BetaD1,  Coeff_D1)
  Coeff_FullY  <- rbind(Coeff_muY,Coeff_BetaD0,Coeff_D0, Coeff_BetaD1,  Coeff_D1)
  
  PhiY_Person <- PhiY_Fac[OutcomeData$Facility]
  OutcomeData$PhiY <- PhiY_Person
  
  #Create variance for outcome data
  Sigma2Y <- c(16,16)
  
  
  
  #Sample true outcome data
  Y0Means <- rep(NA,N)
  Y1Means <- rep(NA,N)
  
  Y0Means <- diag(X_Main_control%*%Coeff_FullY[,OutcomeData$Strata])+PhiY_Person +
    diag(X_Main_control[,2]*X_Main_control[,3])%*% interaction_coeff_control[OutcomeData$Strata]
  Y1Means <- diag(X_Main_treat%*%Coeff_FullY[,OutcomeData$Strata])+PhiY_Person+
    diag(X_Main_control[,2]*X_Main_control[,3])%*% interaction_coeff_control[OutcomeData$Strata]+
    OutcomeData$D*diag(X_Main_treat[,2]*X_Main_treat[,3])%*% interaction_coeff_treat[OutcomeData$Strata]
  
  OutcomeData$Y0 <- rnorm(length(Y0Means), Y0Means,sqrt(Sigma2Y[OutcomeData$Strata]))
  OutcomeData$Y1 <- rnorm(length(Y1Means), Y1Means,sqrt(Sigma2Y[OutcomeData$Strata]))
  
  CZ_Matrix <- True_CZ
  CZ_Matrix[treat_Fac_bin==0,1:L] <- NA
  C_Matrix <- as.matrix(CZ_Matrix[,1:L])
  Z_Matrix <- as.matrix(CZ_Matrix[,(L+1):(L+P)])
  
  
  
  #######################
  #### Model Fitting ####
  #######################
  
  
  C = C_Matrix; Z=Z_Matrix;
  K=2;
  D=OutcomeData$D; X.use=X;
  Y1 = OutcomeData$Y1; Y0 = OutcomeData$Y0;
  Facility= as.numeric(OutcomeData$Facility);
  L <- ncol(C_Matrix)
  P <- ncol(Z_Matrix)
  if(is.null(L)){L=1}
  if(is.null(P)){P=1}
  J <- nrow(CZ_Matrix)
  
  
  ### Get initial values for mixture using Maximum Liklihood for treat_Fac_bool providers
  CZ <- CZ_Matrix[treat_Fac_bin==1,]
  tot_S <- rep(0,K)
  while (min(tot_S)<3) {
    check_mix <- try(mvnormalmixEM(CZ,arbvar = F, k=K ))
    init_clust_treat <- apply(check_mix$posterior,1,which.max)
    init_clust<-rep(NA,J)
    init_clust[treat_Fac_bin==1]<-init_clust_treat
    
    for (k in 1:K) {
      tot_S[k] <- sum(init_clust_treat==k)
    }
  }
  
  # Set value for control providers
  cont_fac <-which(treat_Fac_bin==0)
  for (i in 1:sum(treat_Fac_bin==0)) {
    
    mean1 <- abs(CZ_Matrix[cont_fac[i],2] -check_mix$mu[[1]][2])
    mean2 <- abs(CZ_Matrix[cont_fac[i],2] -check_mix$mu[[2]][2])
    
    init_clust[cont_fac[i]] <- ifelse(mean1<mean2,1,2)
  }
  
  pairs(CZ, col = apply(check_mix$posterior,1,which.max))
  

  
  
  # Create matricies that will hold posterior draws for each parameter
  
  mu_Sim <- list()
  for (k in 1:K) {
    mu_Sim[[k]] <- matrix(NA,nrow = L+P, Sims)
  }
  
  Sigma_Sim <- list()
  
  C_Sim <- list()
  for (l in 1:L) {
    C_Sim[[l]] <- matrix(NA, nrow = J,ncol = Sims)
    C_Sim[[l]][treat_Fac_bin==1,1:Sims] <- C_Matrix[treat_Fac_bin==1,l]
  }
  
  S_Sim <- matrix(NA, J, Sims)
  pi_Sim <- matrix(NA, nrow = K, ncol= Sims)
  
  
  ## Initial Values ###
  
  S_Sim[,1] <- init_clust
  
  for (k in 1:K) {
    mu_Sim[[k]][,1] <- colMeans(CZ_Matrix[treat_Fac_bin==1 &S_Sim[,1]==k,])
  }
  
  Sigma_Sim[[1]] <- diag(0, nrow = L+P)
  for (k in 1:K) {
    Sigma_Sim[[1]] <- Sigma_Sim[[1]] + cov(CZ_Matrix[treat_Fac_bin==1 & S_Sim[,1]==k,])
  }
  
  for (l in 1:L) {
    for (k in 1:K) {
      C_Sim[[l]][treat_Fac_bin==0 & init_clust ==k ,1] <- mean(C_Matrix[treat_Fac_bin==1 & init_clust ==k,l])
    }
  }
  
  
  pi_Sim[,1] <- 1/K
  
  
  #Compliance and Outcome Simulation Matricies
  D_Sim <- matrix(NA, nrow = N, ncol = Sims)

  
  D_Sim[treat_Person_bool,1]<- D[treat_Person_bool]
  for (s in 1:Sims) {
    D_Sim[treat_Person_bool,s] <- D[treat_Person_bool]
  }
  
  
  # Facility_df <- data.frame(Facility= Facility)
  provider_Fac <- unique(Facility)
  person_Fac <- as.numeric(Facility_df$Facility)
  treat_Fac_bool <- 1:J %in% unique(OutcomeData$Facility[OutcomeData$Treat==1])
  treat_Fac_ind <- as.numeric(unique(OutcomeData$Facility[OutcomeData$Treat==1]))
  
  Y0_Sim <- matrix(NA, nrow = N, ncol = Sims)
  Y1_Sim <- matrix(NA, nrow = N, ncol = Sims)
  
  
  muY_Sim <- matrix(NA,nrow = K, ncol = Sims)
  PhiY_Sim <- matrix(NA, nrow = J, ncol = Sims)
  Sigma2_Sim <- matrix(NA, nrow = K, ncol = Sims)
  # Tau2_PhiY_Sim <- matrix(NA,nrow = K,ncol = Sims)
  Tau2_PhiY_Sim <- rep(NA, Sims)
  PhiY_Person <- rep(0,nrow(X_Main_obs))
  
  Y1_Sim[treat_Person_bool,1] <-Y1[treat_Person_bool]
  Y0_Sim[!treat_Person_bool,1] <- Y0[!treat_Person_bool]
  for (s in 1:Sims) {
    
    Y1_Sim[treat_Person_bool,s]<-Y1[treat_Person_bool]
    Y0_Sim[!treat_Person_bool,s] <- Y0[!treat_Person_bool]
  }
  
  PhiY_Sim[,1] <- 0
  
  Tau2_PhiY_Sim[1] <- 1
  
  Sigma2_Sim[,1] <- 1
  
  
  PhiD_Sim <- matrix(0,J,Sims)
  # mu_PhiD_Sim <- matrix(0,nrow = 2, ncol= Sims)
  Tau2_PhiD_Sim <- rep(0.1,Sims)
  PhiD_Person <- rep(0,nrow(X_D))
  
  Strata_Fac <- S_Sim[,1]
  Strata_Person <- Strata_Fac[Facility]
  
  Alpha_Sim <- list()
  for (k in 1:K) {
    Alpha_Sim[[k]] <- matrix(NA,nrow = ncol(X_D), ncol = Sims)
  }
  
  #Get Initial values for compliance coefficients
  CheckData_comp <- as.data.frame(cbind(D=D_Sim[,1], X_D))
  CheckData_comp$D[!treat_Person_bool]<-NA
  prob_D <- rep(NA,length(D))
  for (k in 1:K) {
    reg_compk <- glm(D~.-1, data = CheckData_comp[Strata_Person==k,], family = binomial("probit"))
    Alpha_Sim[[k]][,1]<- summary(reg_compk)$coefficients[,1]
    prob_D[Strata_Person==k] <- predict(reg_compk,newdata = CheckData_comp[Strata_Person==k,], type = "response")
  }
  
  D_Sim[!treat_Person_bool,1] <- rbinom(sum(!treat_Person_bool), 1, prob_D[!treat_Person_bool])
  
  # Initial values for probit regression
  U_Sim <- matrix(nrow = nrow(X_D), ncol = Sims)
  U_Sim[D_Sim[,1]==1,1] <- truncnorm::rtruncnorm(sum(D_Sim[,1]), a=0, b=Inf)
  U_Sim[D_Sim[,1]==0,1] <- truncnorm::rtruncnorm(sum(1-D_Sim[,1]), a=-Inf, b=0)
  X_D<- as.matrix(X_D)
  
  X_Main_obs <- cbind(1,as.matrix(X))
  X_Main_obs <- cbind(X_Main_obs,(1- OutcomeData$Treat)*D_Sim[,1]*1)
  X_Main_obs <- cbind(X_Main_obs,(OutcomeData$Treat)*D_Sim[,1]* cbind(as.matrix(X),1))
  colnames(X_Main_obs) <- paste0("V",1:ncol(X_Main_obs))
  
  X_Main_treat <- cbind(1,as.matrix(X))
  X_Main_treat <- cbind(X_Main_treat,(0)*D_Sim[,1]*1)
  X_Main_treat <- cbind(X_Main_treat,(1)*D_Sim[,1]* cbind(as.matrix(X),1))
  colnames(X_Main_treat) <- paste0("V",1:ncol(X_Main_treat))
  
  X_Main_control <- cbind(1,as.matrix(X))
  X_Main_control <- cbind(X_Main_control,(1)*D_Sim[,1]*1)
  X_Main_control <- cbind(X_Main_control,(0)*D_Sim[,1]* cbind(as.matrix(X),1))
  colnames(X_Main_control) <- paste0("V",1:ncol(X_Main_control))
  
  
  YCN <- ncol(X)
  BetaD0_Sim <- list()
  BetaD1_Sim <- list()
  Delta0_Sim <- matrix(NA,nrow = K,ncol=Sims)
  Delta1_Sim <- matrix(NA,nrow = K,ncol=Sims)
  Beta_full_Sim <- list()
  
  for (k in 1:K) {
    BetaD0_Sim[[k]] <-  matrix(NA,nrow = YCN,ncol = Sims)
    BetaD1_Sim[[k]] <-  matrix(NA,nrow = YCN,ncol = Sims)
    Beta_full_Sim[[k]] <- matrix(NA,nrow=ncol(X_Main_obs),ncol=Sims)
  }
  
  Y_obs <- rep(NA,length(OutcomeData))
  Y_obs[treat_Person_bool]<-Y1[treat_Person_bool]
  Y_obs[!treat_Person_bool]<-Y0[!treat_Person_bool]
  
  CheckData <- as.data.frame(cbind(Y_obs, X_Main_obs))
  colnames(CheckData)[1]<-"Y"
  
  
  #Get initial values for outcome model
  
  for (k in 1:K) {
    reg1k <- lm(Y~.-1,data = CheckData[Strata_Person==k,])
    muY_Sim[k,1]<-summary(reg1k)$coefficients[1,1]
    BetaD0_Sim[[k]][,1] <-summary(reg1k)$coefficients[(1+(1:nrow(BetaD0_Sim[[1]]))),1]
    Delta0_Sim[k,1] <- summary(reg1k)$coefficients[(2+nrow(BetaD0_Sim[[1]])),1]
    BetaD1_Sim[[k]][,1] <-summary(reg1k)$coefficients[(2+nrow(BetaD0_Sim[[1]])+(1:nrow(BetaD1_Sim[[1]]))),1]
    Delta1_Sim[k,1]<-summary(reg1k)$coefficients[ncol(X_Main_obs),1]
    Y1_Sim[Strata_Person==k &treat_Person_bool==0,1] <- predict(reg1k, newdata= as.data.frame(X_Main_treat[Strata_Person==k&treat_Person_bool==0,]))
    Y0_Sim[Strata_Person==k &treat_Person_bool==1,1] <- predict(reg1k, newdata= as.data.frame(X_Main_control[Strata_Person==k&treat_Person_bool==1,]))
    Beta_full_Sim[[k]][,1] <- c( muY_Sim[k,1],BetaD0_Sim[[k]][,1],Delta0_Sim[k,1] ,
                                 BetaD1_Sim[[k]][,1],Delta1_Sim[k,1])
  }
  
  
  ######################################################################
  ######################################################################
  ###################### Gibbs Sampling ################################
  ######################################################################
  ######################################################################
  library(invgamma)
  
  Y1_use <- matrix(NA, nrow = N, ncol= n.imps)
  Y0_use <- matrix(NA, nrow = N, ncol= n.imps)
  D_Samps_use <- matrix(NA, nrow = N/2, ncol= n.imps)
  Beta_full_use1 <- matrix(NA,nrow=nrow(Beta_full_Sim[[1]]),ncol = n.imps)
  Beta_full_use2 <- matrix(NA,nrow=nrow(Beta_full_Sim[[2]]),ncol = n.imps)
  ATE_use <- rep(NA,n.imps)
  CACE_use <- rep(NA,n.imps)
  ATE1_use <- rep(NA,n.imps)
  ATE2_use <- rep(NA,n.imps)
  CACE1_use <- rep(NA,n.imps)
  CACE2_use <- rep(NA,n.imps)
  PhiY_use <- matrix(NA,nrow = J,ncol = n.imps)
  PhiD_use <- matrix(NA,nrow = J,ncol = n.imps)
  D_use <- matrix(NA,nrow = length(D), ncol = n.imps)
  
  
  V_0_inv=NULL;
  m_0 =NULL;
  nu_0 = 7;
  S_0 = NULL; 
  lambda_pk =3;
  if(is.null(V_0_inv)){
    V_0_inv <- diag(1/100,nrow = L+P)}
  if(is.null(m_0)){
    m_0 <- matrix(0,nrow = L+P, ncol = 1)}
  if(is.null(S_0)){
    S_0 <- diag(1/100,nrow = L+P) }
  m_alpha <- rep(0,3)
  v_alpha_inv <- diag(1/100, 3)
  m_beta <- rep(0,ncol(X_Main_obs))
  v_beta_inv <- diag(1/100, ncol(X_Main_obs) )
  v_tauD_mu <- 25
  m_tauD_mu <- 0
  a_sigmas <- 1
  b_sigmas <-1
  max_tauPhiD <- 10
  max_tauPhiY <- 225
  alt_phiD <- matrix(0,nrow = J,ncol = K)
  alt_phiY <- matrix(0,nrow = J,ncol = K)
  
  
  CZ_obs <- CZ_Matrix[treat_Fac_bool,]
  C_obs <- as.matrix(CZ_Matrix[treat_Fac_bool,1:L])
  Z_obs <- as.matrix(CZ_Matrix[treat_Fac_bool,(L+1):(L+P)])
  treat_index <- which(treat_Fac_bool)
  J_treat <- sum(treat_Fac_bool)
  control_index <- which(!treat_Fac_bool)
  J_control <- sum(!treat_Fac_bool)
  
  
  
  ######################################################################
  ######################################################################
  ###################### Gibbs Sampling ################################
  ######################################################################
  ######################################################################
  samps_use <- burn+(1:n.imps*thin)
  S_use <- matrix(NA,nrow=J,ncol = n.imps)
  C_use <- list()
  
  
  C_use <- list()
  for (l in 1:L) {
    C_use[[l]] <- matrix(nrow = J, ncol = n.imps)
  }
  
  
  Coeff_Alpha_est <- matrix(NA,nrow=nrow(Alpha_Sim[[1]]),ncol=K)
  # mu_PhiD_est <- rep(NA,K)
  Tau2_PhiD_est <- NA
  pi_est <- rep(NA,K)
  
  Coeff_BetaD1_est <-matrix(NA,nrow=nrow(BetaD1_Sim[[1]]),ncol=K)
  Coeff_D1_est<- rep(NA,K)
  Coeff_D0_est <- rep(NA,K)
  
  
  
  Strata_Fac <- S_Sim[,1]
  Strata_Person <- S_Sim[person_Fac,1]
  PhiY_Person <- PhiY_Sim[person_Fac,1]
  PhiD_Person <- PhiD_Sim[person_Fac,1]
  
  CZ_i <- CZ_Matrix
  C_i <-  C_Matrix
  Z_i <-Z_Matrix
  
  for (i in 1:(Sims-1)) {

    
    new = i+1
    
    Sigma_inv <- solve(Sigma_Sim[[i]])
    
    for (l in 1:L) {
      CZ_i[treat_Fac_bin==0,l] <- C_Sim[[l]][treat_Fac_bin==0,i] 
      C_i[treat_Fac_bin==0,l] <- C_Sim[[l]][treat_Fac_bin==0,i] 
    }
    
    ### Sampling mu_ks ###
    
    # for (k in 1:K) {
    #   mu_Sim[[k]][,new] <-  mu[,k]
    # }
    
    new_mu <- matrix(NA, nrow = L+P, ncol = K)
    
    for (k in 1:K) {

      fac_in_k <- S_Sim[,i] == k
      J_k <- sum(fac_in_k)

      VJSig_inv <- solve(V_0_inv + J_k*Sigma_inv)

      if(J_k>0){
        CZ_k_bar <- as.matrix(colMeans(matrix(CZ_i[fac_in_k,],ncol=2)))
        mu_hat_k <- VJSig_inv%*%(V_0_inv%*%m_0 + J_k*Sigma_inv%*%CZ_k_bar)
      }else{mu_hat_k <- VJSig_inv%*%(V_0_inv%*%m_0)
      }
      new_mu[,k] <- MASS::mvrnorm(1,mu_hat_k,VJSig_inv)
    }

    # new_mu <- new_mu[,order(new_mu[2,])]
    for (k in 1:K) {
      mu_Sim[[k]][,new] <-  new_mu[,k]
    }
    
  
    
    
    ### Sampling Sigma ###
    
    Sigma_Sim[[new]]<-Sigma
    
    S_J <- S_0
    for (k in 1:K) {

      fac_in_k <- S_Sim[,i] ==k
      J_k <- sum(fac_in_k)

      if(J_k >0){
        CZ_k <- rbind(t(C_i[fac_in_k,]),t(Z_i[fac_in_k,]))
        CZ_k_minus_mu_k <- sweep(CZ_k,1,mu_Sim[[k]][,new])
        S_J = S_J +
          Reduce('+',lapply(seq_len(J_k), function(j) (CZ_k_minus_mu_k[,j])%*%t(CZ_k_minus_mu_k[,j])))}


    }



    Sigma_Sim[[new]]<-riwish(nu_0+J, S_J)

    
    
    ### Sampling p's ###
    
    # pi_Sim[,new] <- pi_k
    
    num_fac_each_k <-c()

    for (k in 1:K) {
      J_k<-sum(Strata_Fac==k)
      num_fac_each_k <- c(num_fac_each_k,J_k)
    }

    pk_post_dir_param <- lambda_pk+num_fac_each_k

    pi_Sim[,new] <- rdirichlet(1,pk_post_dir_param)
  
    
    
    #Determining compliers and noncompliers for probit regression
    
    D1s <- list()
    D0s <- list()
    for (k in 1:K) {

      if(sum(Strata_Fac==k)>0){

        D1s[[k]] <- which(D_Sim[,i]==1 & Strata_Person==k)
        D0s[[k]] <-  which(D_Sim[,i]==0 & Strata_Person==k)


        ## Sampling latent U's


        if(length(D1s[[k]])>0){
          U_Sim[D1s[[k]],new] <- truncnorm::rtruncnorm(length(D1s[[k]]), a=0, b= Inf,
                                                       mean = X_D[D1s[[k]],]%*% Alpha_Sim[[k]][,i]+PhiD_Person[D1s[[k]]])}
        if(length(D0s[[k]])>0){
          U_Sim[D0s[[k]],new] <- truncnorm::rtruncnorm(length(D0s[[k]]), a=-Inf, b= 0,
                                                       mean = X_D[D0s[[k]],]%*% Alpha_Sim[[k]][,i]+PhiD_Person[D0s[[k]]])}
      }
    }

    
    
    
    ### Sample PhiD
    # PhiD_Sim[,new] <- PhiD_Fac
    
    for(k in 1:K){

      people_in_k <- Strata_Person==k
      fac_in_k <- provider_Fac[Strata_Fac==k]
      num_fac_in_k <- length(fac_in_k)
      peop_fac <- person_Fac[people_in_k]

      if(num_fac_in_k>0){
        U_phi <- U_Sim[people_in_k,new] - X_D[people_in_k,]%*%Alpha_Sim[[k]][,i]
        var_phi <- (Tau2_PhiD_Sim[i])/(Nj*Tau2_PhiD_Sim[i]+1)
        mean_den_phi <- Nj + (1/Tau2_PhiD_Sim[i])
        mean_num_phi <- tapply( U_phi ,peop_fac ,sum)
        mu_phi <- matrix(mean_num_phi/mean_den_phi,nrow=num_fac_in_k,ncol=1)

        PhiD_Sim[fac_in_k,new] <- MASS::mvrnorm(1, mu=mu_phi, Sigma = diag(var_phi ,length(fac_in_k)))}


    }
    

    PhiD_Person <- PhiD_Sim[person_Fac,new]
    
    
    
    
    ### Sample Variance for Compliance Facility Effects
    # Tau2_PhiD_Sim[new] <- Tau2_PhiD
 
        Rk_shape <-  1/2 * J - 1/2
        Rk_rate <- 1/2 * sum((PhiD_Sim[,new])^2)
        t <- try(Tau2_PhiD_Sim[new] <- rtrunc(1, "invgamma",a=0,b=  max_tauPhiD,
                                                shape=Rk_shape ,rate=Rk_rate), silent = T)
        if("try-error" %in% class(t) | t==0) Tau2_PhiD_Sim[new] <- max_tauPhiD



    
    
    
    ###  Sample Mean of Random Effects
     
    # for (k in 1:K) {
    # 
    #   people_in_k <- Strata_Person==k
    #   fac_in_k <- Strata_Fac == k
    #   J_k <- sum(fac_in_k)
    # 
    #   if(sum(fac_in_k)>0){
    #     P_tauD_mu <- 1/(J_k/Tau2_PhiD_Sim[k,new] +1/v_tauD_mu)
    #     p_tauD_mu <- sum(PhiD_Sim[fac_in_k,new])/Tau2_PhiD_Sim[k,new] +m_tauD_mu/v_tauD_mu
    #   }else{
    #     P_tauD_mu <- v_tauD_mu
    #     p_tauD_mu <- m_tauD_mu/v_tauD_mu
    #   }
    #   mu_PhiD_Sim[k,new] <- rnorm(1,P_tauD_mu*p_tauD_mu,sqrt(P_tauD_mu))
    # }
  
    
    # mu_PhiD_Sim[,new] <- mu_PhiD
    
    
    
    ## Sampling Compliance Coefficient Alphas
    
    for(k in 1:K){
      
      people_in_k <- Strata_Person==k

      if(sum(people_in_k)>0){
        X_in_k<-X_D[people_in_k,]
        X_in_k_t <-t(X_in_k)

        XtX_inv_k<- chol2inv(chol(X_in_k_t%*%X_in_k +v_alpha_inv))
        a_mu_k<- XtX_inv_k%*%(X_in_k_t%*%(U_Sim[people_in_k,new]-PhiD_Person[people_in_k])+
                                v_alpha_inv%*%m_alpha)
      }else{
        XtX_inv_k <- solve(v_alpha_inv)
        a_mu_k <- m_alpha
      }


      Alpha_Sim[[k]][,new] <- MASS::mvrnorm(1, mu =a_mu_k,
                                            Sigma = XtX_inv_k)

      # Alpha_Sim[[k]][,new] <- Coeff_Alpha[,k]
    }
    
    
    
    ### Sample PhiY
    
        # PhiY_Sim[,new] <- PhiY_Fac
        
        
    for(k in 1:K){

      people_in_k <- Strata_Person==k
      fac_in_k <- provider_Fac[Strata_Fac==k]
      num_fac_in_k <- length(fac_in_k)
      peop_fac <- person_Fac[people_in_k]

      if(num_fac_in_k>0){
        Y_phi <- Y_obs[people_in_k] - X_Main_obs[people_in_k,]%*%Beta_full_Sim[[k]][,i]

        var_phi <- (Sigma2_Sim[k,i]*Tau2_PhiY_Sim[i])/(Nj*Tau2_PhiY_Sim[i]+Sigma2_Sim[k,i])
        mean_den_phi <- Nj + (Sigma2_Sim[k,i]/Tau2_PhiY_Sim[i])
        mean_num_phi <- tapply( Y_phi ,peop_fac ,sum)
        mu_phi <- matrix(mean_num_phi/mean_den_phi,nrow=num_fac_in_k,ncol=1)

        PhiY_Sim[fac_in_k,new] <- MASS::mvrnorm(1, mu=mu_phi, Sigma = diag(var_phi ,length(fac_in_k)))}

    }
    

    PhiY_Person <- PhiY_Sim[person_Fac,new]
    
    
    
    
    #### Sample Variance for Outcome Facility Effects
    # Tau2_PhiY_Sim[new] <- Tau2_PhiY
        
        Rk_shape <-  1/2 * length(fac_in_k) -1/2
        Rk_rate <- 1/2 * sum((PhiY_Sim[fac_in_k,new])^2)
        t <- try(Tau2_PhiY_Sim[new] <- rtrunc(1, "invgamma",a=0,b=  max_tauPhiY,
                                                shape=Rk_shape ,rate=Rk_rate), silent = T)
        if("try-error" %in% class(t) | t==0) Tau2_PhiY_Sim[new] <- max_tauPhiY
        


    ### Sampling Betas
    
    
    for (k in 1:K) {
      
      people_in_k <- Strata_Person==k
      num_in_k <- sum(people_in_k)

      if(num_in_k>0){
        X_in_k <-X_Main_obs[people_in_k,]
        X_in_k_t <- t(X_in_k)
        Sig_inv <- 1/Sigma2_Sim[k,i]

        vXkXk_inv <- chol2inv(chol(Sig_inv*X_in_k_t%*%X_in_k + v_beta_inv))

        Betak_hat <- vXkXk_inv%*%(X_in_k_t%*%(Sig_inv*(Y_obs[people_in_k]-PhiY_Person[people_in_k]))
                                  +v_beta_inv%*%m_beta)
        Sigmak_hat <- vXkXk_inv
      }else{
        Betak_hat <- m_beta
        Sigmak_hat <- solve(v_beta_inv)
      }

      Beta_full_Sim[[k]][,new] <- MASS::mvrnorm(1, mu=Betak_hat, Sigma = Sigmak_hat)
      
      # Beta_full_Sim[[k]][,new] <- c(Coeff_muY[k],Coeff_BetaD0[,k], Coeff_D0[k],
      #                               Coeff_BetaD1[,k],Coeff_D1[k])
      
      muY_Sim[k,new]<-Beta_full_Sim[[k]][1,new]
      BetaD0_Sim[[k]][,new] <- Beta_full_Sim[[k]][1+1:YCN,new]
      Delta0_Sim[k,new] <- Beta_full_Sim[[k]][2+YCN,new]
      BetaD1_Sim[[k]][,new] <- Beta_full_Sim[[k]][((2+YCN)+1:YCN),new]
      Delta1_Sim[k,new] <- Beta_full_Sim[[k]][ncol(X_Main_obs),new]
    }
    
    
    ### Sampling Sigmas
    
        # Sigma2_Sim[,new] <- Sigma2Y
    
    for (k in 1:K) {

      people_in_k <- Strata_Person==k
      Sigma2_shape <- a_sigmas+sum(people_in_k)/2

      if(sum(people_in_k) >0){
        X_in_k <-X_Main_obs[people_in_k,]
        Sigma2_SS = sum((Y_obs[people_in_k]-
                           X_in_k%*%as.matrix(Beta_full_Sim[[k]][,new])-PhiY_Person[people_in_k])^2)
        Sigma2_rate <- b_sigmas+1/2*Sigma2_SS
      }else{
        Sigma2_rate <- b_sigmas
      }

      Sigma2_Sim[k,new] <- rinvgamma(1, Sigma2_shape, Sigma2_rate)


    }
    

    
    
     ## Sampling S ###
    
    CZ_S <- CZ_i
    Sk_probs <- matrix(NA, nrow =J , ncol=K)
    
    for (k in 1:K) {

      people_in_k =Strata_Person==k
      fac_in_k = Strata_Fac==k
      
      # Strata_Person =S_Sim[person_Fac,i-1]
      # Strata_Fac= S_Sim[,i-1]
      
      CZ_dens <- apply(CZ_i, 1, dmvnorm, mean= mu_Sim[[k]][,new],
                       sigma = Sigma_Sim[[new]], log=T)
      
      Y_mean_k <- X_Main_obs%*%Beta_full_Sim[[k]][,new]+
        PhiY_Person
      
      Y_dens <- dnorm(Y_obs, Y_mean_k, sqrt(Sigma2_Sim[k,new]), log = T)
      Y_dens_fac <- aggregate(Y_dens, by=list(person_Fac), FUN=sum)$x
      
      
      D_dens <- dbinom(D_Sim[,i],1, pnorm(X_D%*%Alpha_Sim[[k]][,new]+
                                                   PhiD_Person), log = T)
      
      D_dens_fac <- aggregate(D_dens, by=list(person_Fac), FUN=sum)$x
      
      
      Sk_probs[,k] <- CZ_dens+Y_dens_fac+D_dens_fac+log(pi_Sim[k,new])
      
      # Sk_probs[,k] <- CZ_dens+log(pi_Sim[k,new])
    }
    
    Sk_probs <- exp(Sk_probs-apply(Sk_probs,1,max))
    Sk_probs <- Sk_probs/rowSums(Sk_probs)
    # round(Sk_probs,3)
    
    for (j in 1:J) {
      S_Sim[j,new] <- which(rmultinom(1,1,Sk_probs[j,])==1)
    }
    
    

    Strata_Fac <- S_Sim[,new]
    Strata_Person <- Strata_Fac[Facility]
    
    
    #Imputing Control C
    
    SigmaCC <- as.matrix(Sigma_Sim[[new]][1:L,1:L])
    SigmaCZ <- matrix(Sigma_Sim[[new]][1:L,(L+1):(L+P)], nrow=L, ncol = P)
    SigmaZC <- t(SigmaCZ)
    SigmaZZ <- as.matrix(Sigma_Sim[[new]][(L+1):(L+P),(L+1):(L+P)])
    SigmaZZ_inv <- solve(SigmaZZ)

    S_control <- S_Sim[control_index,new]
    C_control_new <- matrix(NA,nrow = J_control, ncol = L)

    for (k in 1:K) {
      control_and_k <- S_Sim[,new]==k & !treat_Fac_bool
      J_k_cont <- sum(control_and_k )
      if(J_k_cont >0){
      mu_Ck <- mu_Sim[[k]][1:L,new]
      mu_Zk <- mu_Sim[[k]][(L+1):(L+P),new]

      Ck_hat <- mu_Ck +
        SigmaCZ%*%SigmaZZ_inv%*%t(as.matrix(sweep(as.matrix(Z_Matrix[control_and_k,]),2,mu_Zk)))

      Ck_hat_list <- lapply(seq_len(ncol(Ck_hat)), function(i) Ck_hat[,i])
      Sigma_mu_C_hat <- SigmaCC - SigmaCZ%*%SigmaZZ_inv%*%SigmaZC
      sampled_Cs <- matrix(unlist(lapply(Ck_hat_list , mvrnorm, n=1, Sigma=Sigma_mu_C_hat)),
                           nrow=J_k_cont, ncol= L, byrow=T)

      for (l in 1:L) {
        C_Sim[[l]][control_and_k,new] <- sampled_Cs[,l]
      }
      }
    }
    
    
    
    
    
    
    
    
    ## Sample Unobserved D
    D_samp_probs <- rep(NA,N)
    
    
    
    for (k in 1:K) {
      
      if(sum(Strata_Fac==k)>0){
        UnObs_D_k <- Strata_Person==k & !treat_Person_bool
        Y_UnObs_D_k <- Y_obs[UnObs_D_k]
        Y_Mean_UnObs_D_k_comp <- X_Main_control_comp[UnObs_D_k,]%*%Beta_full_Sim[[k]][,new] + PhiY_Person[UnObs_D_k]
        Y_Mean_UnObs_D_k_nocomp <- X_Main_control_nocomp[UnObs_D_k,]%*%Beta_full_Sim[[k]][,new] + PhiY_Person[UnObs_D_k]
        
        comp_Dens <- dnorm(Y_UnObs_D_k, Y_Mean_UnObs_D_k_comp,sqrt(Sigma2_Sim[k,new]))*
          pnorm(X_D[UnObs_D_k,]%*%Alpha_Sim[[k]][,new] + PhiD_Person[UnObs_D_k])
        
        noncomp_Dens <- dnorm(Y_UnObs_D_k, Y_Mean_UnObs_D_k_nocomp,sqrt(Sigma2_Sim[k,new]))*
          (1- pnorm(X_D[UnObs_D_k,]%*%Alpha_Sim[[k]][,new] + PhiD_Person[UnObs_D_k]))
        
        D_samp_probs[UnObs_D_k]<- comp_Dens/(comp_Dens+noncomp_Dens)
        
        D_Sim[UnObs_D_k, new] <- rbinom(sum(UnObs_D_k),1,D_samp_probs[UnObs_D_k])
      }
      
    }
    
    # D_Sim[,new] <-OutcomeData$D
    
    #Creating new response matrix
    X_Main_obs <- cbind(1,as.matrix(X))
    X_Main_obs <- cbind(X_Main_obs,(1- OutcomeData$Treat)*D_Sim[,new]*1)
    X_Main_obs <- cbind(X_Main_obs,(OutcomeData$Treat)*D_Sim[,new]* cbind(as.matrix(X),1))
    colnames(X_Main_obs) <- paste0("V",1:ncol(X_Main_obs))
    
    
    
    
    
    
    
    
    if(i %% 500 ==0){print(i)
      print(table(S_Sim[,i],treat_Fac_bin))}
    
    if(new %in% samps_use){
      thin_ind <- which(samps_use==new)
      S_use[,thin_ind]<-S_Sim[,new]
      for (l in 1:L) {
        C_use[[l]][,thin_ind]<-C_Sim[[l]][,new]
      }
      
      D_use[,thin_ind]<-D_Sim[,new]
      Y1_use[,thin_ind]<-Y1_Sim[,new]
      Y0_use[,thin_ind]<-Y0_Sim[,new]
      # Beta11_use[,thin_ind] <- BetaD1_Sim[[1]][,new]
      # Beta12_use[,thin_ind] <- BetaD1_Sim[[2]][,new]
      # Beta01_use[,thin_ind] <- BetaD0_Sim[[1]][,new]
      # Beta02_use[,thin_ind] <- BetaD0_Sim[[2]][,new]
      Beta_full_use1[,thin_ind] <- Beta_full_Sim[[1]][,new]
      Beta_full_use2[,thin_ind] <- Beta_full_Sim[[2]][,new]
      # D_Samps_use[,thin_ind] <- D_samp_probs[!treat_Person_bool]
      PhiY_use[,thin_ind] <- PhiY_Sim[,new]
      PhiD_use[,thin_ind] <- PhiD_Sim[,new]
      
      
      for (k in 1:K) {
        Coeff_Alpha_est[,k] <- Alpha_Sim[[k]][,new]
        Tau2_PhiD_est<- Tau2_PhiD_Sim[new]
        pi_est[k] <- pi_Sim[k,new]
        
        Coeff_BetaD1_est[,k] <- BetaD1_Sim[[k]][,new]
        Coeff_D1_est[k] <- Delta1_Sim[k,new]
        Coeff_D0_est[k] <- Delta0_Sim[k,new]
        
      }
      
      
      X_D_Alpha_Phi<-as.matrix(X_D)%*%Coeff_Alpha_est
      # X_D_Alpha_Phi <-sweep(X_D_Alpha,2,mu_PhiD_est,"+")
      # X_D_Alpha_Phi_Mean <- sweep(X_D_Alpha_Phi,2,sqrt(Tau2_PhiD_est+1),"/")
      X_D_Alpha_Phi_Mean <- X_D_Alpha_Phi/sqrt(Tau2_PhiD_est+1)
      Strata_D_Probs <- pnorm(X_D_Alpha_Phi_Mean)
      pi_k_Strata_D_Probs <- sweep(Strata_D_Probs,2,pi_est,"*")
      sum_pi_k_Strata_D_Probs <- rowSums(pi_k_Strata_D_Probs)
      X_Main_times_sum_pi_k_Strata_D_Probs <- sweep(cbind(as.matrix(X),1),1,sum_pi_k_Strata_D_Probs ,"*")
      mean_X_Main_times_sum_pi_k_Strata_D_Probs <- colMeans(X_Main_times_sum_pi_k_Strata_D_Probs)
      
      means_Strata_D_Probs <- colMeans(Strata_D_Probs)
      Beta_diff <- Coeff_BetaD1_est
      Beta_diff_phi_D <- rbind(Beta_diff,(Coeff_D1_est-Coeff_D0_est))
      pi_k_Beta_diff <- sweep(Beta_diff_phi_D,2,pi_est,"*")
      pi_k_Beta_diff_times_means_Strata_D_Probs <- rowSums(sweep(pi_k_Beta_diff,2,means_Strata_D_Probs ,"*"))
      
      
      CACE_Num <- mean_X_Main_times_sum_pi_k_Strata_D_Probs %*%pi_k_Beta_diff_times_means_Strata_D_Probs
      CACE_Den <- mean(sum_pi_k_Strata_D_Probs)^2
      CACE_use[thin_ind] <- c(CACE_Num /CACE_Den)
      
      X_Main_beta_diff <-  as.matrix(cbind(X,1)) %*% Beta_diff_phi_D 
      X_Main_beta_diff_D_diff <-  X_Main_beta_diff *Strata_D_Probs
      pi_k_X_Main_beta_diff_D_diff <- sweep(X_Main_beta_diff_D_diff ,2,pi_est,"*")
      
      ATE_use[thin_ind] <- mean(rowSums(pi_k_X_Main_beta_diff_D_diff))
      
      S_1<- which.min(c(mu_Sim[[1]][1,new] , mu_Sim[[2]][1,new]))
      S_2<- which.max(c(mu_Sim[[1]][1,new] , mu_Sim[[2]][1,new]))
      
      CACE1_use[thin_ind] <- sum(X_Main_beta_diff_D_diff[,S_1]) / (sum(Strata_D_Probs[,S_1]))
      CACE2_use[thin_ind] <- sum(X_Main_beta_diff_D_diff[,S_2]) / (sum(Strata_D_Probs[,S_2]))
      
      ATE1_use[thin_ind] <- mean(X_Main_beta_diff_D_diff[,S_1])
      ATE2_use[thin_ind] <- mean(X_Main_beta_diff_D_diff[,S_2])
    }
    
  }

  ATECover[m] <- as.numeric(TrueATE >= quantile(ATE_use,0.025, na.rm = T)&TrueATE <= quantile(ATE_use,0.975, na.rm = T))
  ATEBias[m] <- mean(ATE_use, na.rm = T)-TrueATE
  ATEIW[m] <- quantile(ATE_use,0.975, na.rm = T)-quantile(ATE_use,0.025, na.rm = T)
  
  CACECover[m] <- as.numeric(TrueCACE >= quantile(CACE_use,0.025, na.rm = T)&TrueCACE <= quantile(CACE_use,0.975, na.rm = T))
  CACEBias[m] <- mean(CACE_use, na.rm = T)-TrueCACE
  CACEIW[m] <- quantile(CACE_use,0.975, na.rm = T)-quantile(CACE_use,0.025, na.rm = T)

  ATE1Cover[m] <- as.numeric(TrueATE1 >= quantile(ATE1_use,0.025, na.rm = T)&TrueATE1 <= quantile(ATE1_use,0.975, na.rm = T))
  ATE1Bias[m] <- mean(ATE1_use, na.rm = T)-TrueATE1
  ATE1IW[m] <- quantile(ATE1_use,0.975, na.rm = T)-quantile(ATE1_use,0.025, na.rm = T)
    
  ATE2Cover[m] <- as.numeric(TrueATE2 >= quantile(ATE2_use,0.025, na.rm = T)&TrueATE2 <= quantile(ATE2_use,0.975, na.rm = T))
  ATE2Bias[m] <- mean(ATE2_use, na.rm = T)-TrueATE2
  ATE2IW[m] <- quantile(ATE2_use,0.975, na.rm = T)-quantile(ATE2_use,0.025, na.rm = T)
  
  CACE1Cover[m] <- as.numeric(TrueCACE1 >= quantile(CACE1_use,0.025, na.rm = T)&TrueCACE1 <= quantile(CACE1_use,0.975, na.rm = T))
  CACE1Bias[m] <- mean(CACE1_use, na.rm = T)-TrueCACE1
  CACE1IW[m] <- quantile(CACE1_use,0.975, na.rm = T)-quantile(CACE1_use,0.025, na.rm = T)
  
  CACE2Cover[m] <- as.numeric(TrueCACE2 >= quantile(CACE2_use,0.025, na.rm = T)&TrueCACE2 <= quantile(CACE2_use,0.975, na.rm = T))
  CACE2Bias[m] <- mean(CACE2_use, na.rm = T)-TrueCACE2
  CACE2IW[m] <- quantile(CACE2_use,0.975, na.rm = T)-quantile(CACE2_use,0.025, na.rm = T)
  

  
  
  print(m)
  
  print(mean(ATECover,na.rm=T))
  print(mean(CACECover,na.rm=T))
  print(mean(ATE1Cover,na.rm=T))
  print(mean(ATE2Cover,na.rm=T))
  print(mean(CACE1Cover,na.rm=T))
  print(mean(CACE2Cover,na.rm=T))
  print(mean(ATEBias,na.rm=T))
  print(mean(CACEBias,na.rm=T))
  print(mean(ATE1Bias,na.rm=T))
  print(mean(ATE2Bias,na.rm=T))
  print(mean(CACE1Bias,na.rm=T))
  print(mean(CACE2Bias,na.rm=T))
  
  
  if(m==1 | m%%10 ==0){ 

  
  
  EstimandData <- data.frame(ATEBias,
                             ATECover,
                             ATEIW,
                             ATE1Bias,
                             ATE1Cover,
                             ATE1IW,
                             ATE2Bias,
                             ATE2Cover,
                             ATE2IW,
                             CACEBias,
                             CACECover,
                             CACEIW,
                             CACE1Bias,
                             CACE1Cover,
                             CACE1IW,
                             CACE2Bias,
                             CACE2Cover,
                             CACE2IW)


  }
  
  

  
}

Estimand_Data %>% round(3)



