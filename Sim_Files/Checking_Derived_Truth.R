runs <- 2000

Diff_ATE1 <- rep(NA,runs)
Diff_ATE2 <- rep(NA,runs)
Diff_CACE1 <- rep(NA,runs)
Diff_CACE2 <- rep(NA,runs)
Diff_CACE <- rep(NA,runs)
Tot_Comp <- rep(NA,runs) 
Diff_ATE <- rep(NA,runs)
Diff_CACE_Y1<- rep(NA,runs)
Diff_CACE_Y0<- rep(NA,runs)
Diff_CACE1_Y1<- rep(NA,runs)
Diff_CACE1_Y0<- rep(NA,runs)
Diff_CACE2_Y1<- rep(NA,runs)
Diff_CACE2_Y0<- rep(NA,runs)
Diff_ATE_Y1<- rep(NA,runs)
Diff_ATE_Y0<- rep(NA,runs)
Diff_ATE1_Y1<- rep(NA,runs)
Diff_ATE1_Y0<- rep(NA,runs)
Diff_ATE2_Y1<- rep(NA,runs)
Diff_ATE2_Y0<- rep(NA,runs)

setwd(".../Sim_Files")
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



Diff_ATE_Var <- rep(NA,runs)
Diff_ATE_CI <- matrix(NA,nrow =runs, ncol=2)
Complier_Prop <- rep(NA,runs)
X.D.Comps <- matrix(NA,nrow =runs, ncol=ncol(Covariates))
Coeff_Diffs <- matrix(NA,nrow =runs, ncol=ncol(Covariates))
Phi1.mean.check<- rep(NA,runs)
Phi0.mean.check<- rep(NA,runs)
PhiD.mean.check<- rep(NA,runs)
S_check <- matrix(NA,nrow = J,ncol=runs)
pi_k <- c(0.5,0.5)
Y1.D.mean.check<- rep(NA,runs)
Y0.D.mean.check<- rep(NA,runs)
Samp_ATE <- rep(NA,runs)
Samp_Y1_Mean<- rep(NA,runs)
Samp_Y0_Mean<- rep(NA,runs)
check_var_Y11 <- rep(NA,runs)
check_cov_Y12 <- rep(NA,runs)
check_var_S1 <- rep(NA,runs)
check_var_S2 <- rep(NA,runs)
Total_Var_D <- rep(NA,runs)
Total_Var_Y <- rep(NA,runs)
Treat_CACE <- rep(NA,runs)
CACE_D_Rep <- rep(NA,runs)







set.seed(1)
X_full <- as.matrix(Covariates[sample(1:nrow(Covariates),N,replace = T),])
# X.full <- as.matrix(Covariates[rep(1,1200),])
rownames(X_full) <- 1:nrow(X_full)
set.seed(NULL)

# Define skew and degrees of freedom for multivariate skewed t that
# generates facility covariates and compliance 
# (skew=0 and df = Inf imply coorectly specified MVN)
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
# 0 's imply correctly specified
interaction_coeff_treat <- c(0,0)
interaction_coeff_control <- c(0,0)

for (t in 1:runs) {

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
  
  Obs_CACE <- mean(OutcomeData$Y1[OutcomeData$D==1]-OutcomeData$Y0[OutcomeData$D==1])
  Obs_ATE <- mean(OutcomeData$Y1-OutcomeData$Y0)
  Obs_CACE1 <- mean(OutcomeData$Y1[OutcomeData$D==1&OutcomeData$Strata==1]-OutcomeData$Y0[OutcomeData$D==1&OutcomeData$Strata==1])
  Obs_CACE2 <- mean(OutcomeData$Y1[OutcomeData$D==1&OutcomeData$Strata==2]-OutcomeData$Y0[OutcomeData$D==1&OutcomeData$Strata==2])
  Obs_ATE1 <- mean(OutcomeData$Y1[OutcomeData$Strata==1]-OutcomeData$Y0[OutcomeData$Strata==1])
  Obs_ATE2 <- mean(OutcomeData$Y1[OutcomeData$Strata==2]-OutcomeData$Y0[OutcomeData$Strata==2])
  Obs_CACE_Y1 <-  mean(OutcomeData$Y1[OutcomeData$D==1])
  Obs_CACE_Y0 <-  mean(OutcomeData$Y0[OutcomeData$D==1])
  Obs_ATE_Y1 <- mean(OutcomeData$Y1)
  Obs_ATE_Y0 <- mean(OutcomeData$Y0)
  Obs_CACE1_Y1 <- mean(OutcomeData$Y1[OutcomeData$D==1&OutcomeData$Strata==1])
  Obs_CACE1_Y0 <- mean(OutcomeData$Y0[OutcomeData$D==1&OutcomeData$Strata==1])
  Obs_CACE2_Y1 <- mean(OutcomeData$Y1[OutcomeData$D==1&OutcomeData$Strata==2])
  Obs_CACE2_Y0 <- mean(OutcomeData$Y0[OutcomeData$D==1&OutcomeData$Strata==2])
  Obs_ATE1_Y1 <- mean(OutcomeData$Y1[OutcomeData$Strata==1])
  Obs_ATE1_Y0 <- mean(OutcomeData$Y0[OutcomeData$Strata==1])
  Obs_ATE2_Y1 <- mean(OutcomeData$Y1[OutcomeData$Strata==2])
  Obs_ATE2_Y0 <- mean(OutcomeData$Y0[OutcomeData$Strata==2])
  
  
  
  Diff_CACE[t] <- Obs_CACE
  Diff_CACE1[t] <- Obs_CACE1
  Diff_CACE2[t] <- Obs_CACE2
  Diff_ATE[t] <- Obs_ATE
  Diff_ATE1[t] <- Obs_ATE1
  Diff_ATE2[t] <- Obs_ATE2
  Diff_ATE_Y1[t] <- Obs_ATE_Y1
  Diff_ATE_Y0[t] <- Obs_ATE_Y0
  Diff_CACE_Y1[t] <- Obs_CACE_Y1
  Diff_CACE_Y0[t]<- Obs_CACE_Y0
  Diff_CACE1_Y1[t] <- Obs_CACE1_Y1
  Diff_CACE1_Y0[t]<- Obs_CACE1_Y0
  Diff_CACE2_Y1[t] <- Obs_CACE2_Y1
  Diff_CACE2_Y0[t]<- Obs_CACE2_Y0
  Diff_ATE1_Y1[t] <- Obs_ATE1_Y1
  Diff_ATE1_Y0[t]<- Obs_ATE1_Y0
  Diff_ATE2_Y1[t] <- Obs_ATE2_Y1
  Diff_ATE2_Y0[t]<- Obs_ATE2_Y0
  
  if(t%%1000 ==0){print(t)}
  
}

## Full Empirical ###

### ATE_Y1

XD_Alpha<- cbind(1,as.matrix(X_full))%*%Coeff_Alpha
p_ijk <- pnorm(XD_Alpha/sqrt(Tau2_PhiD+1))
mu_k_plus_X_ij_Beta0_k <- sweep(X_full%*%Coeff_BetaD0,2,Coeff_muY,"+")
delta1_k_plus_X_ij_Beta1_k <- sweep(X_full%*%Coeff_BetaD1,2,Coeff_D1,"+")
p_ijk_times_delta1_k_plus_X_ij_Beta1_k <- p_ijk*delta1_k_plus_X_ij_Beta1_k
ATE_Y1_P1_plus_ATE_Y1_P2<- mu_k_plus_X_ij_Beta0_k+p_ijk_times_delta1_k_plus_X_ij_Beta1_k 
ATE_Y1k <- sweep(ATE_Y1_P1_plus_ATE_Y1_P2,2,pi_k,"*")
True_ATE_Y1 <- mean(rowSums(ATE_Y1k))

mean(Diff_ATE_Y1)
True_ATE_Y1


### ATE_Y0

p_ijk_times_delta0_k <- sweep(p_ijk,2,Coeff_D0,"*")
ATE_Y0_P1_plus_ATE_Y0_P2<- mu_k_plus_X_ij_Beta0_k+p_ijk_times_delta0_k
ATE_Y0k <- sweep(ATE_Y0_P1_plus_ATE_Y0_P2,2,pi_k,"*")
True_ATE_Y0 <- mean(rowSums(ATE_Y0k))

mean(Diff_ATE_Y0)
True_ATE_Y0


### CACE_Y1

den_CACE <- sum(pi_k*colMeans(p_ijk))
p_ijk_weighted_X_Beta_Both <- colMeans(p_ijk*(X_full%*%(Coeff_BetaD0+Coeff_BetaD1)))/ colMeans(p_ijk)
mu_k_delta1_k_plus_weighted_XBeta <-Coeff_muY+ p_ijk_weighted_X_Beta_Both  +Coeff_D1                                       
pijk_times_weighted_out1 <-   sweep(p_ijk,2,mu_k_delta1_k_plus_weighted_XBeta,"*")
num_CACE_Y1 <- sum(pi_k*colMeans(pijk_times_weighted_out1))
True_CACE_Y1 <- num_CACE_Y1/den_CACE  
True_CACE_Y1
mean(Diff_CACE_Y1)

### CACE_Y0


p_ijk_weighted_X_Beta0 <- colMeans(p_ijk*(X_full%*%(Coeff_BetaD0)))/ colMeans(p_ijk)
mu_k_delta0_k_plus_weighted_XBeta0 <-Coeff_muY+ p_ijk_weighted_X_Beta0 +Coeff_D0                                       
pijk_times_weighted_out0 <-   sweep(p_ijk,2,mu_k_delta0_k_plus_weighted_XBeta0,"*")
num_CACE_Y0 <- sum(pi_k*colMeans(pijk_times_weighted_out0))
True_CACE_Y0 <- num_CACE_Y0/den_CACE  
True_CACE_Y0
mean(Diff_CACE_Y0)
  


### ATE1_Y1

mu1_plus_Mean_XBeta01 <- Coeff_muY[1]+mean(X_full%*%Coeff_BetaD1[,1])
mean_pij1_times_XBeta11_plus_Delta11 <- mean(p_ijk[,1]*(X_full%*%Coeff_BetaD1[,1]+Coeff_D1[1]))
True_ATE1_Y1 <- mu1_plus_Mean_XBeta01 +mean_pij1_times_XBeta11_plus_Delta11
True_ATE1_Y1
mean(Diff_ATE1_Y1)


### ATE1_Y0

mu1_plus_Mean_XBeta01 <- Coeff_muY[1]+mean(X_full%*%Coeff_BetaD1[,1])
mean_pij1_times_Delta01 <- mean(p_ijk[,1]*(Coeff_D0[1]))
True_ATE1_Y0 <- mu1_plus_Mean_XBeta01 +mean_pij1_times_Delta01
True_ATE1_Y0 
mean(Diff_ATE1_Y0)


### CACE1_Y1
True_CACE1_Y1 <- mu_k_delta1_k_plus_weighted_XBeta[1]
True_CACE1_Y1
mean(Diff_CACE1_Y1)

### CACE1_Y0
True_CACE1_Y0 <- mu_k_delta0_k_plus_weighted_XBeta0[1]
True_CACE1_Y0
mean(Diff_CACE1_Y0)


#### Strata Specific X


#  ATE_k  Y1_k and Y0_k
Y1_k <- rep(NA,K)
Y0_k <- rep(NA,K)
ATE_k <-rep(NA,K)
for (k in 1:K) {
  
  XD_Alpha<- cbind(1,as.matrix(X_full))%*%Coeff_Alpha[,k]
  p_ijk <- pnorm(XD_Alpha/sqrt(Tau2_PhiD+1))
  muY_weighted_XBeta0_k <- Coeff_muY[k] +mean(X_full%*%Coeff_BetaD0[,k])
  weighted_XBeta1_k_delta1_k <- mean(p_ijk*(X_full%*%Coeff_BetaD1[,k]+ Coeff_D1[k])) 
  weighted_delta0_k <- mean(p_ijk*(Coeff_D0[k]))
  
  Y1_k[k] <-  muY_weighted_XBeta0_k + weighted_XBeta1_k_delta1_k
  Y0_k[k] <-  muY_weighted_XBeta0_k +weighted_delta0_k
  ATE_k[k] <-Y1_k[k]-Y0_k[k]
}

ATE1_True_S <- ATE_k[1]
mean(Diff_ATE1)
ATE1_True_S

ATE1_Y1_True_S <- Y1_k[1]
mean(Diff_ATE1_Y1)
ATE1_Y1_True_S

ATE1_Y0_True_S <- Y0_k[1]
mean(Diff_ATE1_Y0)
ATE1_Y0_True_S

ATE2_True_S <- ATE_k[2]
mean(Diff_ATE2)
ATE2_True_S

ATE2_Y1_True_S <- Y1_k[2]
mean(Diff_ATE2_Y1)
ATE2_Y1_True_S

ATE2_Y0_True_S <- Y0_k[2]
mean(Diff_ATE2_Y0)
ATE2_Y0_True_S


#  CACE_k  Y1_D_k and Y0_D_k
p_ijk_mean <- rep(NA,K)
Y1_D_k <- rep(NA,K)
Y0_D_k <- rep(NA,K)
CACE_k <-rep(NA,K)
for (k in 1:K) {

  XD_Alpha<- cbind(1,as.matrix(X_full))%*%Coeff_Alpha[,k]
  p_ijk <- pnorm(XD_Alpha/sqrt(Tau2_PhiD+1))
  p_ijk_mean[k] <- mean(p_ijk)
  weighted_XBeta0_kBeta1_k <- mean( p_ijk *X_full%*%(Coeff_BetaD0[,k]+Coeff_BetaD1[,k]))
  weighted_XBeta0_k <- mean(p_ijk *X_full%*%(Coeff_BetaD0[,k]))
  Y1_D_k[k] <-  Coeff_muY[k] + weighted_XBeta0_kBeta1_k/p_ijk_mean[k] + Coeff_D1[k]
  Y0_D_k[k] <- Coeff_muY[k] + weighted_XBeta0_k/p_ijk_mean[k] + Coeff_D0[k]
  CACE_k[k] <- Y1_D_k[k]-Y0_D_k[k]
}

CACE1_True_S <- CACE_k[1]
mean(Diff_CACE1)
CACE1_True_S 

CACE1_Y1_True_S <- Y1_D_k[1]
mean(Diff_CACE1_Y1)
CACE1_Y1_True_S 

CACE1_Y0_True_S <- Y0_D_k[1]
mean(Diff_CACE1_Y0)
CACE1_Y0_True_S 

CACE2_True_S  <- CACE_k[2]
mean(Diff_CACE2)
CACE2_True_S 

CACE2_Y1_True_S <- Y1_D_k[2]
mean(Diff_CACE2_Y1)
CACE2_Y1_True_S 

CACE2_Y0_True_S <- Y0_D_k[2]
mean(Diff_CACE2_Y0)
CACE2_Y0_True_S 

#ATE and CACE

ATE <- 0
ATE_Y1 <-0
ATE_Y0 <-0
CACE_num <- 0
CACE_den <- 0
CACE_Y1_num <- 0
CACE_Y1_den <- 0
CACE_Y0_num <- 0
CACE_Y0_den <- 0

for (k in 1:K) {
  ATE <- ATE+ pi_k[k]*ATE_k[k]
  ATE_Y1 <- ATE_Y1+pi_k[k]*Y1_k[k]
  ATE_Y0 <- ATE_Y0+ pi_k[k]*Y0_k[k]
  
  CACE_num <-  CACE_num+ pi_k[k]*CACE_k[k]*p_ijk_mean[k]
  CACE_den <- CACE_den+pi_k[k]*p_ijk_mean[k]
  CACE_Y1_num <-  CACE_Y1_num+ pi_k[k]*Y1_D_k[k]*p_ijk_mean[k]
  CACE_Y0_num <-  CACE_Y0_num+ pi_k[k]*Y0_D_k[k]*p_ijk_mean[k]
} 
CACE_Y1_True_S <- CACE_Y1_num/CACE_den
mean(Diff_CACE_Y1)
CACE_Y1_True_S

CACE_Y0_True_S <- CACE_Y0_num/CACE_den
mean(Diff_CACE_Y0)
CACE_Y0_True_S

CACE_True_S <- CACE_num/CACE_den
mean(Diff_CACE)
CACE_True_S

ATE_Y1_True_S  <- ATE_Y1
mean(Diff_ATE_Y1)
ATE_Y1_True_S

ATE_Y0_True_S  <- ATE_Y0
mean(Diff_ATE_Y0)
ATE_Y0_True_S

ATE_True_S  <- ATE 
mean(Diff_ATE)
ATE_True_S


