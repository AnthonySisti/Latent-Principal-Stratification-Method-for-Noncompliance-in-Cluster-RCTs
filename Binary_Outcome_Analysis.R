library("boot", include.only = "inv.logit")
library("mixtools")
library(coda)
library(MCMCpack,include.only = "riwish")
library(DirichletReg, include.only =  "rdirichlet")
library(mixtools)
library(mvtnorm, include.only =  "dmvnorm")
library(MASS, include.only = "mvrnorm")
library(scales, include.only = "alpha")
library(invgamma)
library(truncdist)
library(dplyr)

# Set Gibbs sampler parameters 
burn=500
thin = 3
n.imps=1000
Sims <- burn + thin*n.imps
samps_use <- burn+(1:n.imps*thin)

#Load in Data
CZ_Matrix <-  read.csv(".../Data/CZ_Matrix_Synth.csv")
X <- read.csv(".../Data/X_Synth.csv")
OutcomeData <- read.csv( ".../Data/Outcome_Synth_Binary.csv")

#Define outcome, individual and cluster level covariate and compliance matricies
X <- data.matrix(X)
X_D <- data.matrix(cbind(Int=1,X))
C_Matrix <- CZ_Matrix %>% select(C)
Z_Matrix <- CZ_Matrix %>% select(Z)

# Treatment indicators
treat_Fac_bool <- !is.na(CZ_Matrix$C)
treat_Fac_bin <- as.numeric(treat_Fac_bool)
treat_Person_bool <- OutcomeData$Treat==1
treat_Person_bin <- OutcomeData$Treat
provider_Fac <- unique(OutcomeData$Facility)
person_Fac <- as.numeric(OutcomeData$Facility)



#K = 2 strata
K=2;N =nrow(X)
#Set outcomes and covariates
D=OutcomeData$D; X.use=X;
Y1 = OutcomeData$Y1; Y0 = OutcomeData$Y0;
Facility= as.numeric(OutcomeData$Facility);
#Set number of compliance metrics and covariates
L <- ncol(C_Matrix)
P <- ncol(Z_Matrix)
if(is.null(L)){L=1}
if(is.null(P)){P=1}
CZ_obs <- as.matrix(CZ_Matrix[treat_Fac_bool,])
C_obs <- as.matrix(CZ_Matrix[treat_Fac_bool,1:L])
Z_obs <- as.matrix(CZ_Matrix[treat_Fac_bool,(L+1):(L+P)])

# Number of treatment and control facilities
J_treat <- sum(treat_Fac_bool)
treat_index <- which(treat_Fac_bool)
J_control <- sum(!treat_Fac_bool)
control_index <- which(!treat_Fac_bool)



# Get initial values for cluster compliance strata
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

#Plot initial cluster compliance strata
pairs(CZ, col = apply(check_mix$posterior,1,which.max))


# Cluster Strata and probability Gibbs sample holders
S_Sim <- matrix(NA, J, Sims)
S_Sim[,1] <- init_clust
Strata_Fac <- S_Sim[,1]
Strata_Person <- Strata_Fac[person_Fac]
pi_Sim <- matrix(NA, nrow = K, ncol= Sims)
pi_Sim[,1] <- 1/K


# Variance Matrix of cluster compliance/characteristic
Sigma_Sim <- list()
#initial value
Sigma_Sim[[1]] <- diag(0, nrow = L+P)
for (k in 1:K) {
  Sigma_Sim[[1]] <- Sigma_Sim[[1]] + cov(CZ_Matrix[treat_Fac_bin==1 & S_Sim[,1]==k,])
}


# Means of cluster compliance/characteristic mixture
mu_Sim <- list()
for (k in 1:K) {
  mu_Sim[[k]] <- matrix(NA,nrow = L+P, Sims)
  #initial value
  mu_Sim[[k]][,1] <- colMeans(CZ_Matrix[treat_Fac_bin==1 &S_Sim[,1]==k,])
}


# Cluster compliance values, observed and imputed
C_Sim <- list()
for (l in 1:L) {
  C_Sim[[l]] <- matrix(NA, nrow = J,ncol = Sims)
  C_Sim[[l]][treat_Fac_bin==1,1:Sims] <- C_Matrix[treat_Fac_bin==1,l]
  for (k in 1:K) {
    C_Sim[[l]][treat_Fac_bin==0 & init_clust ==k ,1] <- mean(C_Matrix[treat_Fac_bin==1 & init_clust ==k,l])
  }
}

# Individual Compliance and Outcome Simulation Matricies
D_Sim <- matrix(NA, nrow = N, ncol = Sims)
D_Sim[treat_Person_bool,1]<- D[treat_Person_bool]
for (s in 1:Sims) {
  D_Sim[treat_Person_bool,s] <- D[treat_Person_bool]
}

Y0_Sim <- matrix(NA, nrow = N, ncol = Sims)
Y1_Sim <- matrix(NA, nrow = N, ncol = Sims)
Y1_Sim[treat_Person_bool,1] <-Y1[treat_Person_bool]
Y0_Sim[!treat_Person_bool,1] <- Y0[!treat_Person_bool]

for (s in 1:Sims) {
  Y1_Sim[treat_Person_bool,s]<-Y1[treat_Person_bool]
  Y0_Sim[!treat_Person_bool,s] <- Y0[!treat_Person_bool]
}



# Mean and facility effects simulation matricies
muY_Sim <- matrix(NA,nrow = K, ncol = Sims)

PhiY_Sim <- matrix(NA, nrow = J, ncol = Sims)
PhiY_Sim[,1] <- 0

Tau2_PhiY_Sim <- rep(NA, Sims)
Tau2_PhiY_Sim[1] <- 1
PhiY_Person <- rep(0,N)

Tau2_PhiD_Sim <- rep(NA,Sims)
Tau2_PhiD_Sim[1] <- 0.1 
PhiD_Person <- rep(0,N)

PhiD_Sim <- matrix(0,J,Sims)



# Coefficients for compliance regression model
DCN <- ncol(X_D)
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

D_Sim[is.na(OutcomeData$D),1] <- rbinom(sum(is.na(OutcomeData$D)), 1, prob_D[is.na(OutcomeData$D)])

U_Sim <- matrix(nrow = nrow(X_D), ncol = Sims)
U_Sim[D_Sim[,1]==1,1] <- truncnorm::rtruncnorm(sum(D_Sim[,1]), a=0, b=Inf)
U_Sim[D_Sim[,1]==0,1] <- truncnorm::rtruncnorm(sum(1-D_Sim[,1]), a=-Inf, b=0)



# Coefficients for outcome regression models
YCN <- ncol(X)
BetaD0_Sim <- list()
BetaD1_Sim <- list()
Delta0_Sim <- matrix(NA,nrow = K,ncol=Sims)
Delta1_Sim <- matrix(NA,nrow = K,ncol=Sims)
Beta_full_Sim <- list()

# Matrices used for imputations
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


for (k in 1:K) {
  BetaD0_Sim[[k]] <-  matrix(NA,nrow = YCN,ncol = Sims)
  BetaD1_Sim[[k]] <-  matrix(NA,nrow = YCN,ncol = Sims)
  Beta_full_Sim[[k]] <- matrix(NA,nrow=ncol(X_Main_obs),ncol=Sims)
}


#Get initial values for outcome model

Y_obs <- matrix(NA,nrow=nrow(OutcomeData),ncol = Sims)

for (s in 1:Sims) {
  Y_obs[treat_Person_bool,s]<-Y1[treat_Person_bool]
  Y_obs[!treat_Person_bool,s]<-Y0[!treat_Person_bool]
}

CheckData <- as.data.frame(cbind(Y_obs[,1], X_Main_obs))
colnames(CheckData)[1]<-"Y"

for (k in 1:K) {
  
  reg1k <- glm(Y~.-1,data = CheckData[Strata_Person==k,], family = binomial("probit"))
  muY_Sim[k,1]<-summary(reg1k)$coefficients[1,1]
  BetaD0_Sim[[k]][,1] <-summary(reg1k)$coefficients[(1+(1:nrow(BetaD0_Sim[[1]]))),1]
  Delta0_Sim[k,1] <- summary(reg1k)$coefficients[(2+nrow(BetaD0_Sim[[1]])),1]
  BetaD1_Sim[[k]][,1] <-summary(reg1k)$coefficients[(2+nrow(BetaD0_Sim[[1]])+(1:nrow(BetaD1_Sim[[1]]))),1]
  Delta1_Sim[k,1]<-summary(reg1k)$coefficients[ncol(X_Main_obs),1]
  
  Y_obs[treat_Person_bool &is.na(Y1),1] <- as.numeric(predict(reg1k, newdata= as.data.frame(X_Main_treat[treat_Person_bool &is.na(Y1),]), type ="response")>0.5)
  Y_obs[!treat_Person_bool &is.na(Y0),1] <- as.numeric(predict(reg1k, newdata= as.data.frame(X_Main_control[!treat_Person_bool &is.na(Y0),]), type ="response")>0.5) 
  
  Y1_Sim[Strata_Person==k &!treat_Person_bool,1] <- as.numeric(predict(reg1k, newdata= as.data.frame(X_Main_treat[Strata_Person==k&treat_Person_bool==0,]), type="response")>0.5)
  Y0_Sim[Strata_Person==k &treat_Person_bool,1] <- as.numeric(predict(reg1k, newdata= as.data.frame(X_Main_control[Strata_Person==k&treat_Person_bool==1,]), type="response")>0.5)
  Beta_full_Sim[[k]][,1] <- c( muY_Sim[k,1],BetaD0_Sim[[k]][,1],Delta0_Sim[k,1] ,
                               BetaD1_Sim[[k]][,1],Delta1_Sim[k,1])
}

Y1_Sim[treat_Person_bool,1] <- Y_obs[treat_Person_bool,1]
Y0_Sim[!treat_Person_bool,1]<- Y_obs[!treat_Person_bool,1]

Q_Sim <- matrix(nrow = nrow(X_Main_obs), ncol = Sims)
Q_Sim[Y_obs[,1]==1,1] <- truncnorm::rtruncnorm(sum(Y_obs[,1]), a=0, b=Inf)
Q_Sim[Y_obs[,1]==0,1] <- truncnorm::rtruncnorm(sum(1-Y_obs[,1]), a=-Inf, b=0)





### Prior Parameters ###

covsX <- c(2:(YCN+1), (YCN+3):(ncol(X_Main_obs)-1))
noncovsX <- c(1,YCN+2,ncol(X_Main_obs))
V_0_inv=NULL;
m_0 =NULL;
nu_0 = 7;
S_0 = NULL; 
lambda_pk =3;
V_0_inv <- diag(c(1/1000,1/100000))
if(is.null(m_0)){
  m_0 <- matrix(0,nrow = L+P, ncol = 1)}
if(is.null(S_0)){
  S_0 <- diag(1/10000,nrow = L+P) }
m_alpha <- rep(0,DCN)
v_alpha_inv <- diag(1/25, DCN)
m_beta <- rep(0,length(covsX))
m_beta_k <- rep(0,length(noncovsX))
v_beta_inv <- diag(1/100,length(covsX))
v_beta_k_inv<-diag(c(1/10000, 1/16, 1/100))
v_tauD_mu <- 9
m_tauD_mu <- 0
a_sigmas <- 1
b_sigmas <-1
max_tauPhiD <- 10
max_tauPhiY <- 225
alt_phiD <- matrix(0,nrow = J,ncol = K)
alt_phiY <- matrix(0,nrow = J,ncol = K)
if(is.null(V_0_inv)){
  V_0_inv <- diag(1/1000000,nrow = L+P)}
if(is.null(m_0)){
  m_0 <- matrix(0,nrow = L+P, ncol = 1)}
if(is.null(S_0)){
  S_0 <- diag(1/10000,nrow = L+P) }



######################################################################
######################################################################
###################### Gibbs Sampling ################################
######################################################################
######################################################################

Y1_use <- matrix(NA, nrow = N, ncol= n.imps)
Y0_use <- matrix(NA, nrow = N, ncol= n.imps)
D_use <- matrix(NA,nrow = length(D), ncol = n.imps)

Beta_full_use1 <- matrix(NA,nrow=nrow(Beta_full_Sim[[1]]),ncol = n.imps)
Beta_full_use2 <- matrix(NA,nrow=nrow(Beta_full_Sim[[2]]),ncol = n.imps)
ATE_use <- rep(NA,n.imps)
ATE_Y1_use <- rep(NA,n.imps)
ATE_Y0_use <-rep(NA,n.imps)
CACE_use <- rep(NA,n.imps)
CACE_Y1_use <- rep(NA,n.imps)
CACE_Y0_use  <- rep(NA,n.imps)
ATE1_use <- rep(NA,n.imps)
ATE1_Y1_use <- rep(NA,n.imps)
ATE1_Y0_use <-rep(NA,n.imps)
ATE2_use <- rep(NA,n.imps)
ATE2_Y1_use <- rep(NA,n.imps)
ATE2_Y0_use <-rep(NA,n.imps)
CACE1_use <- rep(NA,n.imps)
CACE1_Y1_use <- rep(NA,n.imps)
CACE1_Y0_use <- rep(NA,n.imps)
CACE2_use <- rep(NA,n.imps)
CACE2_Y1_use <- rep(NA,n.imps)
CACE2_Y0_use <- rep(NA,n.imps)
PhiY_use <- matrix(NA,nrow = J,ncol = n.imps)
PhiD_use <- matrix(NA,nrow = J,ncol = n.imps)
D_use <- matrix(NA,nrow = length(D), ncol = n.imps)

S_use <- matrix(NA,nrow=J,ncol = n.imps)
C_use <- list()
C_use <- list()
for (l in 1:L) {
  C_use[[l]] <- matrix(nrow = J, ncol = n.imps)
}

Coeff_Alpha_est <- matrix(NA,nrow=nrow(Alpha_Sim[[1]]),ncol=K)
mu_PhiD_k_est <- rep(NA,K)
Tau2_PhiD_est <- rep(NA,K)
pi_est <- rep(NA,K)

Coeff_BetaD1_est <-matrix(NA,nrow=nrow(BetaD1_Sim[[1]]),ncol=K)
Coeff_BetaD0_est <-matrix(NA,nrow=nrow(BetaD0_Sim[[1]]),ncol=K)
Coeff_D1_est<- rep(NA,K)
Coeff_D0_est <- rep(NA,K)
muY_est <- rep(NA,K)
X_k <- list()
XD_k <- list()

missingYs <- (treat_Person_bool&is.na(OutcomeData$Y1))|(!treat_Person_bool&is.na(OutcomeData$Y0))
Nj <-table(person_Fac)

CZ_i <- CZ_Matrix
C_i <-  C_Matrix
Z_i <-Z_Matrix

PPcheckD <- matrix(NA, nrow=N,ncol=n.imps)
PPcheckY <- matrix(NA, nrow=N,ncol=n.imps)
PPcheckS <- matrix(NA, nrow=J,ncol=n.imps)
PPcheckS_person <- matrix(NA, nrow=N,ncol=n.imps)

X_Main_control_comp <- cbind(1,as.matrix(X))
X_Main_control_comp <- cbind(X_Main_control_comp,(1- 0)*1*1)
X_Main_control_comp <- cbind(X_Main_control_comp,(0)*1* cbind(as.matrix(X),1))

X_Main_control_nocomp <- cbind(1,as.matrix(X))
X_Main_control_nocomp <- cbind(X_Main_control_nocomp,(1- 0)*0*1)
X_Main_control_nocomp <- cbind(X_Main_control_nocomp,(0)*0* cbind(as.matrix(X),1))

for (i in 1:(Sims-1)) {
  
  new = i+1
  
  Sigma_inv <- solve(Sigma_Sim[[i]])
  for (l in 1:L) {
    CZ_i[treat_Fac_bin==0,l] <- C_Sim[[l]][treat_Fac_bin==0,i] 
    C_i[treat_Fac_bin==0,l] <- C_Sim[[l]][treat_Fac_bin==0,i] 
  }
  
  
  ### Sampling mu_ks ###
  new_mu <- matrix(NA, nrow = L+P, ncol = K)
  
  for (k in 1:K) {
    
    fac_in_k <- Strata_Fac ==k
    J_k <- sum(fac_in_k)
    
    VJSig_inv <- solve(V_0_inv + J_k*Sigma_inv)
    
    if(J_k>0){
      CZ_k_bar <- as.matrix(colMeans(CZ_i[fac_in_k,]))
      mu_hat_k <- VJSig_inv%*%(V_0_inv%*%m_0 + J_k*Sigma_inv%*%CZ_k_bar)
    }else{mu_hat_k <- VJSig_inv%*%(V_0_inv%*%m_0)
    }
    new_mu[,k] <- MASS::mvrnorm(1,mu_hat_k,VJSig_inv)
  }
  
  for (k in 1:K) {
    mu_Sim[[k]][,new] <-  new_mu[,k]
  }
  
  
  ### Sampling Sigma ###
  S_J <- S_0
  for (k in 1:K) {
    fac_in_k <- Strata_Fac ==k
    J_k <- sum(fac_in_k)
    
    if(J_k >0){
      CZ_k <- rbind(t(C_i[fac_in_k,]),t(Z_i[fac_in_k,]))
      CZ_k_minus_mu_k <- sweep(CZ_k,1,mu_Sim[[k]][,new])
      S_J = S_J +
        Reduce('+',lapply(seq_len(J_k), function(j) (CZ_k_minus_mu_k[,j])%*%t(CZ_k_minus_mu_k[,j])))}
  }
  Sigma_Sim[[new]]<-riwish(nu_0+J, S_J)
  
  
  ### Sampling p's ###
  num_fac_each_k <-c()
  
  for (k in 1:K) {
    J_k<-sum(Strata_Fac==k)
    num_fac_each_k <- c(num_fac_each_k,J_k)
  }
  
  pk_post_dir_param <- lambda_pk+num_fac_each_k
  
  pi_Sim[,new] <- rdirichlet(1,pk_post_dir_param)
  
  
  ### Determining compliers and noncompliers for probit regression ###
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
  
  
  ### Sample PhiD ###
  for(k in 1:K){
    
    people_in_k <- Strata_Person==k
    fac_in_k <- provider_Fac[Strata_Fac==k]
    num_fac_in_k <- length(fac_in_k)
    peop_fac <- person_Fac[people_in_k]
    
    if(num_fac_in_k>0){
      U_phi <- U_Sim[people_in_k,new] - X_D[people_in_k,]%*%Alpha_Sim[[k]][,i]
      Njk <- Nj[fac_in_k]
      
      var_phi <- (Tau2_PhiD_Sim[i])/(Njk*Tau2_PhiD_Sim[i]+1)
      mean_den_phi <- Njk + (1/Tau2_PhiD_Sim[i])
      mean_num_phi <- tapply( U_phi ,peop_fac ,sum)
      mu_phi <- matrix(mean_num_phi/ mean_den_phi,nrow=num_fac_in_k,ncol=1)
      
      PhiD_Sim[fac_in_k,new] <- MASS::mvrnorm(1, mu=mu_phi, Sigma = diag(var_phi,num_fac_in_k))}
    
  }
  
  PhiD_Person <- PhiD_Sim[person_Fac,new]
  
  
  ### Sample Variance for Compliance Facility Effects ###
  Rk_shape <-  1/2 * J - 1/2
  Rk_rate <- sum(PhiD_Sim[,new]^2)
  Rk_rate <- 1/2* Rk_rate
  t <- try(Tau2_PhiD_Sim[new] <- rtrunc(1, "invgamma",a=0,b=  max_tauPhiD, 
                                        shape=Rk_shape ,rate=Rk_rate), silent = T)
  if("try-error" %in% class(t) | t==0) Tau2_PhiD_Sim[new] <- max_tauPhiD
  
  
  ### Sampling Compliance Coefficient Alphas ###
  X_in_k<-X_D
  X_in_k_t <-t(X_in_k)
  XtX_inv_k<- chol2inv(chol(X_in_k_t%*%X_in_k +v_alpha_inv))
  a_mu_k<- XtX_inv_k%*%(X_in_k_t%*%(U_Sim[,new]-PhiD_Person)+
                          v_alpha_inv%*%m_alpha)
  
  Alpha_new<- MASS::mvrnorm(1, mu =a_mu_k,
                            Sigma = XtX_inv_k)
  for(k in 1:K){
    Alpha_Sim[[k]][,new] <-Alpha_new
  }
  
  
  ### Determining 1s and 0s for probit regression ###
  Y1s <- list()
  Y0s <- list()
  for (k in 1:K) {
    if(sum(Strata_Fac==k)>0){
      
      Y1s[[k]] <- which(Y_obs[,i]==1 & Strata_Person==k)
      Y0s[[k]] <-  which(Y_obs[,i]==0 & Strata_Person==k)
      
      ## Sampling latent U's
    
      if(length(Y1s[[k]])>0){
        Q_Sim[Y1s[[k]],new] <- truncnorm::rtruncnorm(length(Y1s[[k]]), a=0, b= Inf,
                                                     mean = X_Main_obs[Y1s[[k]],]%*% Beta_full_Sim[[k]][,i]+PhiY_Person[Y1s[[k]]])}
      if(length(Y0s[[k]])>0){
        Q_Sim[Y0s[[k]],new] <- truncnorm::rtruncnorm(length(Y0s[[k]]), a=-Inf, b= 0,
                                                     mean = X_Main_obs[Y0s[[k]],]%*% Beta_full_Sim[[k]][,i]+PhiY_Person[Y0s[[k]]])}
    }
  }     
  
  
  ### Sample PhiY
  for(k in 1:K){
    
    people_in_k <- Strata_Person==k
    fac_in_k <- provider_Fac[Strata_Fac==k]
    num_fac_in_k <- length(fac_in_k)
    peop_fac <- person_Fac[people_in_k]
    
    if(num_fac_in_k>0){
      Njk <- Nj[fac_in_k]
      Q_phi <- Q_Sim[people_in_k,new] -  X_Main_obs[people_in_k,]%*%Beta_full_Sim[[k]][,i]
      
      var_phi <- (Tau2_PhiY_Sim[i])/(Njk*Tau2_PhiY_Sim[i]+1)
      mean_den_phi <- Njk + (1/Tau2_PhiY_Sim[i])
      mean_num_phi <- tapply( Q_phi ,peop_fac ,sum)
      mu_phi <- matrix(mean_num_phi/ mean_den_phi,nrow=num_fac_in_k,ncol=1)
  
      PhiY_Sim[fac_in_k,new] <- MASS::mvrnorm(1, mu=mu_phi, Sigma = diag(var_phi,num_fac_in_k))}
   }
  
  PhiY_Person <- PhiY_Sim[person_Fac,new]
  

  #### Sample Variance for Outcome Facility Effects
  Rk_shape <-  1/2 * J -1/2
  Rk_rate <- 1/2 * sum((PhiY_Sim[,new])^2)
  t <- try(Tau2_PhiY_Sim[new] <- rtrunc(1, "invgamma",a=0,b=  max_tauPhiY, 
                                        shape=Rk_shape ,rate=Rk_rate), silent = T)
  if("try-error" %in% class(t) | t==0) Tau2_PhiY_Sim[new] <- max_tauPhiY
  
  
  ### Sampling muY, delta0 and delta1 (complier specific) ###
  for (k in 1:K) {
    
    people_in_k <- Strata_Person==k
    num_in_k <- sum(people_in_k)
    
    if(num_in_k>0){
      X_in_k <-X_Main_obs[people_in_k,noncovsX]
      X_in_k_t <- t(X_in_k)
      
      vXkXk_inv <- chol2inv(chol(X_in_k_t%*%X_in_k + v_beta_k_inv))
      
      Xc_in_k <- X_Main_obs[people_in_k,covsX]
      Betac <- Beta_full_Sim[[k]][covsX,i]
      Betak_hat <- vXkXk_inv%*%(X_in_k_t%*%(Q_Sim[people_in_k,new]- Xc_in_k%*% Betac-PhiY_Person[people_in_k])
                                +v_beta_k_inv%*%m_beta_k)
      
      Sigmak_hat <- vXkXk_inv
    }else{
      Betak_hat <- m_beta_k
      Sigmak_hat <- solve(v_beta_k_inv)
    }
    
    muy_deltak_new <- MASS::mvrnorm(1, mu=Betak_hat, Sigma = Sigmak_hat)
    
    muY_Sim[k,new] <- muy_deltak_new[1]
    Delta0_Sim[k,new] <- muy_deltak_new[2]
    Delta1_Sim[k,new] <- muy_deltak_new[3]
  
  }
  
  
  #Sampling regression coefficient beta vector, constant across strata
  X_in_k <-X_Main_obs[,covsX]
  X_in_k_t <- t(X_in_k)
  vXkXk_inv <- chol2inv(chol(X_in_k_t%*%X_in_k + v_beta_inv))
  Xc_in_k <- X_Main_obs[,noncovsX]
  Betac <- rbind(muY_Sim[Strata_Person,new],Delta0_Sim[Strata_Person,new],Delta1_Sim[Strata_Person,new])
  Betak_hat <- vXkXk_inv%*%(X_in_k_t%*%(Q_Sim[,new]- diag(Xc_in_k%*% Betac)-  PhiY_Person)
                            +v_beta_inv%*%m_beta)
  Sigmak_hat <- vXkXk_inv
  
  Beta_covs_new <- MASS::mvrnorm(1, mu=Betak_hat, Sigma = Sigmak_hat)
  
  for (k in 1:K) {
    Beta_full_Sim[[k]][,new] <- c(muY_Sim[k,new],
                                  Beta_covs_new[1:YCN],
                                  Delta0_Sim[k,new],
                                  Beta_covs_new[(YCN+1):(2*YCN)],
                                  Delta1_Sim[k,new])
  }
  
  
  ### Sampling S ###
  Sk_probs <- matrix(NA, nrow =J , ncol=K)
  CZ_S <- CZ_i
  
  for (k in 1:K) {
    
    people_in_k =Strata_Person==k
    fac_in_k = Strata_Fac==k
    
    CZ_dens <- apply(CZ_S, 1, dmvnorm, mean= mu_Sim[[k]][,new],
                     sigma = Sigma_Sim[[new]], log=T)
    
    Y_mean_k <- pnorm(X_Main_obs%*%Beta_full_Sim[[k]][,new]+ PhiY_Person)
    Y_dens <- dbinom(Y_obs[,i], 1, Y_mean_k, log = T)
    
    Y_dens_fac <- aggregate(Y_dens, by=list(person_Fac), FUN=sum)$x
    
    D_dens <- dbinom(D_Sim[,i],1, pnorm(X_D%*%Alpha_Sim[[k]][,new]+ PhiD_Person),log = T)
    
    D_dens_fac <- aggregate(D_dens, by=list(person_Fac), FUN=sum)$x
    
    Sk_probs[,k] <-  CZ_dens+Y_dens_fac+D_dens_fac+log(pi_Sim[k,new])
    
  }
  
  Sk_probs <- exp(Sk_probs-apply(Sk_probs,1,max))
  Sk_probs <- Sk_probs/rowSums(Sk_probs)
  round(Sk_probs,3)
  
  for (j in 1:J) {
    S_Sim[j,new] <- which(rmultinom(1,1,Sk_probs[j,])==1)
  }
  
  Strata_Fac <- S_Sim[,new]
  Strata_Person <- Strata_Fac[Facility]
  
  
  ### Imputing Control C ###
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
  
  
  ### Sample Unobserved D ###
  D_samp_probs <- rep(NA,N)
  for (k in 1:K) {
    
    if(sum(Strata_Fac==k)>0){
      UnObs_D_k <- Strata_Person==k & is.na(OutcomeData$D)
      Y_UnObs_D_k <- Y_obs[UnObs_D_k, i]
      Y_Mean_UnObs_D_k_comp <- pnorm(X_Main_control_comp[UnObs_D_k,]%*%Beta_full_Sim[[k]][,new] + PhiY_Person[UnObs_D_k])
      Y_Mean_UnObs_D_k_nocomp <- pnorm(X_Main_control_nocomp[UnObs_D_k,]%*%Beta_full_Sim[[k]][,new] + PhiY_Person[UnObs_D_k])
      
      comp_Dens <- dbinom(Y_UnObs_D_k,1, Y_Mean_UnObs_D_k_comp)*
        pnorm(X_D[UnObs_D_k,]%*%Alpha_Sim[[k]][,new] + PhiD_Person[UnObs_D_k])
      
      noncomp_Dens <- dbinom(Y_UnObs_D_k, 1, Y_Mean_UnObs_D_k_nocomp)*
        (1- pnorm(X_D[UnObs_D_k,]%*%Alpha_Sim[[k]][,new] + PhiD_Person[UnObs_D_k]))
      
      D_samp_probs[UnObs_D_k]<- comp_Dens/(comp_Dens+noncomp_Dens)
      
      D_Sim[UnObs_D_k, new] <- rbinom(sum(UnObs_D_k),1,D_samp_probs[UnObs_D_k])
    }
    
  }
  
  
  ### Creating new response matrix ###
  X_Main_obs <- cbind(1,as.matrix(X))
  X_Main_obs <- cbind(X_Main_obs,(1- OutcomeData$Treat)*D_Sim[,new]*1)
  X_Main_obs <- cbind(X_Main_obs,(OutcomeData$Treat)*D_Sim[,new]* cbind(as.matrix(X),1))
  colnames(X_Main_obs) <- paste0("V",1:ncol(X_Main_obs))
  
  ### Imputing Missing Y's ###
  for (k in 1:K) {
    
    missingYks <- missingYs&Strata_Person==k
    Ymiss_mean <- pnorm(X_Main_obs[missingYks,]%*%Beta_full_Sim[[k]][,new] +
                          PhiY_Person[missingYks])
    
    Y_obs[missingYks,new] <- rbinom(sum(missingYks),1,Ymiss_mean)
  }
  
  Y1_Sim[treat_Person_bool,new] <- Y_obs[treat_Person_bool,new]
  Y0_Sim[!treat_Person_bool,new] <- Y_obs[!treat_Person_bool,new]
  
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
    Beta_full_use1[,thin_ind] <- Beta_full_Sim[[1]][,new]
    Beta_full_use2[,thin_ind] <- Beta_full_Sim[[2]][,new]
    PhiY_use[,thin_ind] <- PhiY_Sim[,new]
    PhiD_use[,thin_ind] <- PhiD_Sim[,new]
    
    for (k in 1:K) {
      Coeff_Alpha_est[,k] <- Alpha_Sim[[k]][,new]
      Tau2_PhiD_est[k] <- Tau2_PhiD_Sim[new]
      pi_est[k] <- pi_Sim[k,new]
      Coeff_BetaD1_est[,k] <- Beta_full_Sim[[k]][((2+YCN)+1:YCN),new]
      Coeff_D1_est[k] <- Delta1_Sim[k,new]
      Coeff_D0_est[k] <- Delta0_Sim[k,new]
    }
    

    ### Using Procedure 3 to approximate posterior distribution of the estimands ###
    R=1000
    Y1s<- list()
    Y0s<-list()
    Ss <- list()
    Ds <- list()
    XD_k <- list()
    X_k <- list()
    
    for (r in 1:R) {
      Y1s_imp <- rep(NA,N)
      Y0s_imp <- rep(NA,N)
      Ds_imp <- rep(NA,N)
      S_imp <- rbinom(N, 1, pi_Sim[2,new]) +1
      
      for(k in 1:K){
        N_k <- sum(S_imp==k)
        S_ind <- sample(which(Strata_Person==k),N_k,replace = T)
        
        imp_PersonPhiD <- rnorm(N_k, mean = 0, sd = sqrt(Tau2_PhiD_Sim[new]))
        imp_PersonPhiY <- rnorm(N_k, mean = 0, sd = sqrt(Tau2_PhiY_Sim[new]))
        
        D_imps <- rbinom(N_k, 1, pnorm(X_D[S_ind,]%*%Alpha_Sim[[k]][,new]+imp_PersonPhiD))
        
        X_Main_treat_imp <- cbind(1,as.matrix(X[S_ind,]))
        X_Main_treat_imp <- cbind(X_Main_treat_imp,(0)*D_imps*1)
        X_Main_treat_imp <- cbind(X_Main_treat_imp,(1)*D_imps* cbind(as.matrix(X[S_ind,]),1))
        colnames(X_Main_treat_imp) <- paste0("V",1:ncol(X_Main_treat_imp))
        
        X_Main_control_imp <- cbind(1,as.matrix(X[S_ind,]))
        X_Main_control_imp <- cbind(X_Main_control_imp,(1)*D_imps*1)
        X_Main_control_imp <- cbind(X_Main_control_imp,(0)*D_imps* cbind(as.matrix(X[S_ind,]),1))
        colnames(X_Main_control_imp) <- paste0("V",1:ncol(X_Main_control_imp))
        
        Y1s_imp[S_imp==k] <- rbinom(N_k, 1, pnorm(X_Main_treat_imp%*%Beta_full_Sim[[k]][,new]+imp_PersonPhiY))
        Y0s_imp[S_imp==k] <- rbinom(N_k, 1, pnorm(X_Main_control_imp%*%Beta_full_Sim[[k]][,new]+imp_PersonPhiY))
        Ds_imp[S_imp==k] <- D_imps
        
      }
      
      Y1s[[r]] <- Y1s_imp
      Y0s[[r]] <- Y0s_imp
      Ss[[r]] <- S_imp
      Ds[[r]] <- Ds_imp
    }
    
    Y1s_calc <- unlist(Y1s)
    Y0s_calc <- unlist(Y0s)
    Ss_calc <-  unlist(Ss)
    Ds_calc <- unlist(Ds)
    
    ATE_Y1_use[thin_ind] <- mean(Y1s_calc)
    ATE_Y0_use[thin_ind] <- mean(Y0s_calc)
    ATE_use[thin_ind] <- mean(Y1s_calc - Y0s_calc)
    
    ATE1_Y1_use[thin_ind] <- mean(Y1s_calc[Ss_calc==1])
    ATE1_Y0_use[thin_ind] <- mean(Y0s_calc[Ss_calc==1])
    ATE1_use[thin_ind] <- mean(Y1s_calc[Ss_calc==1] - Y0s_calc[Ss_calc==1])
    
    ATE2_Y1_use[thin_ind] <- mean(Y1s_calc[Ss_calc==2])
    ATE2_Y0_use[thin_ind] <- mean(Y0s_calc[Ss_calc==2])
    ATE2_use[thin_ind] <- mean(Y1s_calc[Ss_calc==2] - Y0s_calc[Ss_calc==2])
    
    CACE_Y1_use[thin_ind] <- mean(Y1s_calc[Ds_calc==1])
    CACE_Y0_use[thin_ind] <- mean(Y0s_calc[Ds_calc==1])
    CACE_use[thin_ind] <- mean(Y1s_calc[Ds_calc==1] - Y0s_calc[Ds_calc==1])
    
    CACE1_Y1_use[thin_ind] <- mean(Y1s_calc[Ds_calc==1& Ss_calc==1])
    CACE1_Y0_use[thin_ind] <- mean(Y0s_calc[Ds_calc==1& Ss_calc==1])
    CACE1_use[thin_ind] <- mean(Y1s_calc[Ds_calc==1 & Ss_calc==1] - Y0s_calc[Ds_calc==1& Ss_calc==1])
    
    CACE2_Y1_use[thin_ind] <- mean(Y1s_calc[Ds_calc==1& Ss_calc==2])
    CACE2_Y0_use[thin_ind] <- mean(Y0s_calc[Ds_calc==1& Ss_calc==2])
    CACE2_use[thin_ind] <- mean(Y1s_calc[Ds_calc==1 & Ss_calc==2] - Y0s_calc[Ds_calc==1& Ss_calc==2])
    
    
    ##PPchecks
    PPcheckS[,thin_ind] <- rbinom(J,1, pi_Sim[2,i])+1
    PPcheckS_person[,thin_ind] <-  PPcheckS[,thin_ind][Facility]
    
    for (k in 1:K) {
      
      people_in_k <- PPcheckS_person[,thin_ind]==k
      PP_phiD <- rnorm(J, mean = 0, sd = sqrt(Tau2_PhiD_Sim[new]))
      PP_phiY <- rnorm(J, mean = 0, sd = sqrt(Tau2_PhiY_Sim[new]))
      PP_phiD_person <- PP_phiD[Facility]
      PP_phiY_person <- PP_phiY[Facility]
      
      PP_Dmean <- pnorm(X_D[people_in_k,]%*%Alpha_Sim[[k]][,new] + PP_phiD_person[people_in_k])
      
      PPcheckD[people_in_k,thin_ind] <- rbinom(sum(people_in_k),1,PP_Dmean)
      
      PP_X_Main_obsk <- cbind(1,as.matrix(X[people_in_k,]))
      PP_X_Main_obsk <- cbind(PP_X_Main_obsk,(1- treat_Person_bin[people_in_k])*PPcheckD[people_in_k,thin_ind]*1)
      PP_X_Main_obsk  <- cbind(PP_X_Main_obsk ,treat_Person_bin[people_in_k]*PPcheckD[people_in_k,thin_ind]* cbind(as.matrix(X[people_in_k,]),1))
      colnames(PP_X_Main_obsk ) <- paste0("V",1:ncol(PP_X_Main_obsk ))
      
      PP_Ymean <- pnorm(PP_X_Main_obsk%*%Beta_full_Sim[[k]][,new] +
                          PP_phiY_person[people_in_k])
      
      PPcheckY[people_in_k,thin_ind] <- rbinom(sum(people_in_k),1,PP_Ymean)
    }
  }
}






out_df <- data.frame(ATE = ATE_use,
                     ATE1 = ATE1_use,
                     ATE2 = ATE2_use,
                     ATE1_ATE2 = ATE1_use-ATE2_use,
                     CACE = CACE_use,
                     CACE1 = CACE1_use,
                     CACE2 = CACE2_use,
                     CACE1_CACE2 = CACE1_use-CACE2_use)



out_df %>% summarise(estimand = colnames(.),
                     mean = colMeans(.),
                     sd = apply(.,2,sd),
                     L_95_CrI =apply(.,2, quantile,0.025),
                     U_95_CrI = apply(.,2, quantile,0.975)) %>% mutate(across(where(is.numeric), ~round(.x,3)))



plot(1:n.imps, ATE_use,type="l")
plot(1:n.imps, ATE1_use,type="l")
plot(1:n.imps, ATE2_use,type="l")
plot(1:n.imps, CACE_use,type="l")
plot(1:n.imps, CACE1_use,type="l")
plot(1:n.imps, CACE2_use,type="l")
plot(1:n.imps, Beta_full_use2[5,],type="l")
plot(1:n.imps, Alpha_Sim[[1]][4,samps_use],type="l")



hist(colMeans(PPcheckY[treat_Person_bool & !is.na(OutcomeData$Y1),]),breaks = 30)
abline(v=mean(OutcomeData$Y1,na.rm=T),col="red")

hist(colMeans(PPcheckY[!treat_Person_bool& !is.na(OutcomeData$Y0),]),breaks = 30)
abline(v=mean(OutcomeData$Y0,na.rm=T),col="red")

hist(colMeans(PPcheckD[!is.na(OutcomeData$D),]),breaks=20)
abline(v=mean(OutcomeData$D,na.rm=T),col="red")


## PPcheck 

PP_Y_treat_comp_S1 <- list()
PP_Y_treat_comp_S2 <- list()
PP_Y_control_comp_S1 <- list()
PP_Y_control_comp_S2 <- list()

N_Y_treat_comp_S1 <- rep(NA, ncol(PPcheckY))
N_Y_treat_comp_S2 <- rep(NA, ncol(PPcheckY))
N_Y_control_comp_S1 <- rep(NA, ncol(PPcheckY))
N_Y_control_comp_S2 <- rep(NA, ncol(PPcheckY))

PP_signal_S1 <- rep(NA, ncol(PPcheckY))
PP_signal_S2 <- rep(NA, ncol(PPcheckY))

PP_noise_S1 <- rep(NA, ncol(PPcheckY))
PP_noise_S2 <- rep(NA, ncol(PPcheckY))


for (i in 1:ncol(PPcheckY)) {
  
  
  PP_Y_treat_comp_S1[[i]] <- PPcheckY[(treat_Person_bool & PPcheckD[,i] ==1 &
                                         PPcheckS_person[,i]==1),i]
  N_Y_treat_comp_S1[i] <- sum((treat_Person_bool &PPcheckD[,i] ==1 &
                                 PPcheckS_person[,i]==1))
  
  PP_Y_treat_comp_S2[[i]] <- PPcheckY[(treat_Person_bool & PPcheckD[,i] ==1 &
                                         PPcheckS_person[,i]==2),i]
  
  N_Y_treat_comp_S2[i] <- sum((treat_Person_bool &PPcheckD[,i] ==1 &
                                 PPcheckS_person[,i]==2))
  
  PP_Y_control_comp_S1[[i]] <- PPcheckY[(!treat_Person_bool &PPcheckD[,i] ==1 &
                                           PPcheckS_person[,i]==1),i]
  
  N_Y_control_comp_S1[i] <- sum((!treat_Person_bool &PPcheckD[,i] ==1 &
                                   PPcheckS_person[,i]==1))
  
  PP_Y_control_comp_S2[[i]] <- PPcheckY[( !treat_Person_bool &PPcheckD[,i] ==1 &
                                            PPcheckS_person[,i]==2),i]
  
  N_Y_control_comp_S2[i] <- sum(!treat_Person_bool &PPcheckD[,i] ==1 &
                                  PPcheckS_person[,i]==2)
  
  
  PP_signal_S1[i] <- abs(mean(PP_Y_treat_comp_S1[[i]]) -
                           mean(PP_Y_control_comp_S1[[i]]   ))  
  PP_signal_S2[i] <- abs(mean(PP_Y_treat_comp_S2[[i]]) -
                           mean(PP_Y_control_comp_S2[[i]]   ))  
  
  PP_noise_S1[i] <- sqrt(var(PP_Y_treat_comp_S1[[i]])/N_Y_treat_comp_S1[i] +
                           var(PP_Y_control_comp_S1[[i]])/N_Y_control_comp_S1[i])
  
  PP_noise_S2[i] <- sqrt(var(PP_Y_treat_comp_S2[[i]])/N_Y_treat_comp_S2[i] +
                           var(PP_Y_control_comp_S2[[i]])/N_Y_control_comp_S2[i])
}





Obs_treat_comp_S1 <- list()
Obs_control_comp_S1 <- list()

Obs_treat_comp_S2 <- list()
Obs_control_comp_S2 <- list()

Obs_treat_N_S1 <- rep(NA,ncol(PPcheckY))
Obs_control_N_S1 <- rep(NA,ncol(PPcheckY))

Obs_treat_N_S2 <- rep(NA,ncol(PPcheckY))
Obs_control_N_S2 <- rep(NA,ncol(PPcheckY))

Obs_signal_S1 <- rep(NA,ncol(PPcheckY))
Obs_signal_S2 <- rep(NA,ncol(PPcheckY))

Obs_noise_S1 <- rep(NA,ncol(PPcheckY))
Obs_noise_S2 <- rep(NA,ncol(PPcheckY))




for(i in 1:ncol(PPcheckY)){
  
  Obs_treat_comp_S1[[i]] <-  Y_obs[treat_Person_bool & 
                                     D_use[,i]==1 &
                                     S_use[,i][Facility]==1,i]
  
  Obs_control_comp_S1[[i]] <-   Y_obs[! treat_Person_bool & 
                                        D_use[,i]==1 &
                                        S_use[,i][Facility]==1,i]
  
  Obs_treat_comp_S2[[i]] <-  Y_obs[treat_Person_bool & 
                                     D_use[,i]==1 &
                                     S_use[,i][Facility]==2,i]
  
  Obs_control_comp_S2[[i]] <-   Y_obs[! treat_Person_bool & 
                                        D_use[,i]==1 &
                                        S_use[,i][Facility]==2,i]
  
  Obs_treat_N_S1[i] <-  sum(treat_Person_bool & 
                              D_use[,i]==1 &
                              S_use[,i][Facility]==1)
  
  Obs_control_N_S1[i] <-   sum(! treat_Person_bool & 
                                 D_use[,i]==1 &
                                 S_use[,i][Facility]==1)
  
  Obs_treat_N_S2[i] <-  sum(treat_Person_bool & 
                              D_use[,i]==1 &
                              S_use[,i][Facility]==2)
  
  Obs_control_N_S2[i] <-   sum(! treat_Person_bool & 
                                 D_use[,i]==1 &
                                 S_use[,i][Facility]==2)
  
  
  
  Obs_signal_S1[i] <- abs(mean(Obs_treat_comp_S1[[i]]) -
                            mean(Obs_control_comp_S1[[i]]   ))  
  Obs_signal_S2[i] <- abs(mean(Obs_treat_comp_S2[[i]]) -
                            mean(Obs_control_comp_S2[[i]]   ))  
  
  Obs_noise_S1[i] <- sqrt(var(Obs_treat_comp_S1[[i]])/Obs_treat_N_S1[i] +
                            var(Obs_control_comp_S1[[i]])/Obs_control_N_S1[i])
  
  Obs_noise_S2[i] <- sqrt(var(Obs_treat_comp_S2[[i]])/Obs_treat_N_S2[i] +
                            var(Obs_control_comp_S2[[i]])/Obs_control_N_S2[i])
  
}


Obs_signal_S1_use <- mean(Obs_signal_S1, na.rm=T)
Obs_signal_S2_use <- mean(Obs_signal_S2, na.rm=T)

Obs_noise_S1_use <- mean(Obs_noise_S1, na.rm=T)
Obs_noise_S2_use <- mean(Obs_noise_S2, na.rm=T)

hist(PP_signal_S1)
abline(v= Obs_signal_S1_use)
mean(PP_signal_S1 >Obs_signal_S1, na.rm=T)

hist(PP_signal_S2)
abline(v= Obs_signal_S2_use)
mean(PP_signal_S2 >Obs_signal_S2, na.rm=T)

hist(PP_noise_S1)
abline(v= Obs_noise_S1_use)
mean(PP_noise_S1 >Obs_noise_S1, na.rm=T)

hist(PP_noise_S2)
abline(v= Obs_noise_S2_use)
mean(PP_noise_S2 >Obs_noise_S2, na.rm=T)

hist(PP_signal_S1/PP_noise_S1)
abline(v= Obs_signal_S1_use/Obs_noise_S1_use)
mean(PP_signal_S1/PP_noise_S1 >Obs_signal_S1/Obs_noise_S1, na.rm=T)

hist(PP_signal_S2/PP_noise_S2)
abline(v= Obs_signal_S2_use/Obs_noise_S2_use)
mean(PP_signal_S2/PP_noise_S2 >Obs_signal_S2/Obs_noise_S2, na.rm=T)


PP_check_results <- data.frame(Outcome= "Continuous",
                               `Compliance Level` ="30 Minutes",
                               `Signal S=1` =   mean(PP_signal_S1 >Obs_signal_S1, na.rm=T),
                               `Signal S=2` =   mean(PP_signal_S2 >Obs_signal_S2, na.rm=T),
                               `Noise S=1` =    mean(PP_noise_S1 >Obs_noise_S1, na.rm=T),
                               `Noise S=2` =    mean(PP_noise_S2 >Obs_noise_S2, na.rm=T),
                               `Signal/Noise S=1` = mean(PP_signal_S1/PP_noise_S1 >Obs_signal_S1/Obs_noise_S1, na.rm=T),
                               `Signal/Noise S=2` = mean(PP_signal_S2/PP_noise_S2 >Obs_signal_S2/Obs_noise_S2, na.rm=T))




post_data_C <- matrix(NA, nrow =J, ncol = n.imps)
post_data_Z <- matrix(NA, nrow =J, ncol = n.imps)

mean_1C_samps <- matrix(NA, nrow =5000, ncol =  n.imps)
mean_2C_samps <- matrix(NA, nrow =5000, ncol =  n.imps)
mean_1Z_samps <- matrix(NA, nrow =5000, ncol = n.imps)
mean_2Z_samps <- matrix(NA, nrow =5000, ncol =  n.imps)


PPCs <-c()
PPZs <-c()

for (i in 1:n.imps) {
  
  ind <-samps_use[i]
  
  S1 <- S_Sim[,ind]==1
  # 
  mv_samp1 <-rmvnorm(sum(S1), mu = mu_Sim[[1]][,ind], sigma = Sigma_Sim[[ind]])
  mv_samp2 <-rmvnorm(sum(!S1), mu = mu_Sim[[2]][,ind], sigma = Sigma_Sim[[ind]])
  
  post_data_C[S1,i] <- inv.logit(mv_samp1[,1])
  post_data_C[!S1,i] <- inv.logit(mv_samp2[,1])
  post_data_Z[S1,i] <- mv_samp1[,2] 
  post_data_Z[!S1,i] <-mv_samp2[,2]
  
  
  samp1 <-rmvnorm(5000, mu = mu_Sim[[1]][,ind], sigma = Sigma_Sim[[ind]])
  samp2 <-rmvnorm(5000, mu = mu_Sim[[2]][,ind], sigma = Sigma_Sim[[ind]])
  mean_1C_samps[,i] <-inv.logit(samp1[,1])
  mean_2C_samps[,i] <- inv.logit(samp2[,1])
  mean_1Z_samps[,i] <- samp1[,2]
  mean_2Z_samps[,i] <- samp2[,2]
  
  
  
  
}


par(mfrow=c(1,2))
plot(inv.logit(C_Matrix$C),Z_Matrix$Z, pch=16, main="Obseved Data",
     ylab="Z", xlab= "C",
     ylim =c(-70,100), xlim = c(0,1))
plot(post_data_C[,1:200], post_data_Z[,1:200], pch = 16,cex=0.5, 
     col= alpha(S_Sim[,samps_use[1:200]]*2,0.2), main="Posterior Predictive Distribution",
     ylab="Z", xlab= "C",
     ylim =c(-70,100), xlim = c(0,1))



C_name <- "C"
Z_name <- "Z"

check_dose_mean <- data.frame(Means = c(colMeans(mean_1C_samps),colMeans(mean_1Z_samps),
                                        colMeans(mean_2C_samps), colMeans(mean_2Z_samps)),
                              Strata= c(rep(1,n.imps*2),rep(2,n.imps*2)),
                              Var = rep(c(rep(C_name,n.imps),rep(Z_name,n.imps)),2)
)

check_dose_mean$Var <- as.factor(check_dose_mean$Var)
check_dose_mean$Strata <- as.factor(check_dose_mean$Strata)

Var_name <-levels(check_dose_mean$Var)

library(ggplot2)

ggplot(check_dose_mean, aes(x=Means, fill = Strata)) +   
  geom_density(alpha =0.6, adjust=1.6)+ 
  facet_wrap(~Var, scales = "free") +theme_bw() +
  ggtitle("Posterior Means within each Latent Compliance Strata")+
  scale_fill_manual(breaks = c(1,2), values = c("red","deepskyblue3"))


S1_CMix_Mean <- colMeans(mean_1C_samps)
S2_CMix_Mean <- colMeans(mean_2C_samps)

S1_ZMix_Mean <- colMeans(mean_1Z_samps)
S2_ZMix_Mean <- colMeans(mean_2Z_samps)

Diff_CMix_Mean <- colMeans(mean_1C_samps-mean_2C_samps)
Diff_ZMix_Mean <- colMeans(mean_1Z_samps-mean_2Z_samps)


Mean_Mix_results <- data.frame(Mixture_Parameter = c("mu_C", "mu_Z"), 
                               Mean_S1 = c(round(mean(S1_CMix_Mean),2),
                                           round(mean(S1_ZMix_Mean),1)),
                               SD_S1 = c(round(sd(S1_CMix_Mean),2),
                                         round(sd(S1_ZMix_Mean) ,1)),
                               CRI_95_S1 = c(paste("[",round(quantile(S1_CMix_Mean,0.025),2),", ",
                                                   round(quantile(S1_CMix_Mean,0.975),2),"]", sep=""),
                                             paste("[",round(quantile(S1_ZMix_Mean,0.025),1),", ",
                                                   round(quantile(S1_ZMix_Mean,0.975),1),"]", sep="")),
                               Mean_S2 = c(round(mean(S2_CMix_Mean),2),
                                           round(mean(S2_ZMix_Mean),1)),
                               SD_S2 = c(round(sd(S2_CMix_Mean),2),
                                         round(sd(S2_ZMix_Mean) ,1)),
                               CRI_95_S2 = c(paste("[",round(quantile(S2_CMix_Mean,0.025),2),", ",
                                                   round(quantile(S2_CMix_Mean,0.975),2),"]", sep = ""),
                                             paste("[",round(quantile(S2_ZMix_Mean,0.025),1),", ",
                                                   round(quantile(S2_ZMix_Mean,0.975),1),"]", sep = "")),
                               Mean_Diff = c(round(mean(S1_CMix_Mean-S2_CMix_Mean),2),
                                             round(mean(S1_ZMix_Mean-S2_ZMix_Mean) ,1)),
                               SD_Diff = c(round(sd(S1_CMix_Mean-S2_CMix_Mean),2),
                                           round(sd(S1_ZMix_Mean-S2_ZMix_Mean) ,1)),
                               CRI_95_Diff = c(paste("[",round(quantile(S1_CMix_Mean-S2_CMix_Mean,0.025),2),", ",
                                                     round(quantile(S1_CMix_Mean-S2_CMix_Mean,0.975),2),"]", sep=""),
                                               paste("[",round(quantile(S1_ZMix_Mean-S2_ZMix_Mean,0.025),1),", ",
                                                     round(quantile(S1_ZMix_Mean-S2_ZMix_Mean,0.975),1),"]", sep="")))


Mean_Mix_results






