# Simulation Study

## Data Generation

```
trainMA <- function(n.study, n.AE, n.cancer, n.drug, ss, spar, b, mu_u){

  study=rep(1:n.study, n.AE)
  cancer=rep(1:n.cancer, n.study)
  drug=rep(rep(1:n.drug, each=n.study/n.drug),n.AE)
  AE = rep(1:n.AE, each=n.study)
  
  # construct design matrices
  study.m=model.matrix(~-1+factor(study))    # 
  cancer.m=model.matrix(~-1+factor(cancer))  # 
  drug.m=model.matrix(~-1+factor(drug))      # 
  AE.m = model.matrix(~-1+factor(AE))        # 
  
  # Design matrix for random effects
  Z = cbind(study.m,cancer.m, drug.m, AE.m) # 
  colnames(Z) <- c(paste("study", 1:n.study, sep=""), # 
                   paste("cancer", 1:n.cancer, sep=""),
                   paste("drug", 1:n.drug, sep=""),
                   paste("AE", 1:n.AE, sep=""))
  
  AEcancer = as.numeric(factor(cancer*100+AE))      # 
  AEdrug = as.numeric(factor(drug*100+AE))          # 
  cancerDrug <- as.numeric(factor(cancer*10+drug))  #
  
  AEcancer.m <- model.matrix(~-1+factor(AEcancer))
  AEdrug.m <- model.matrix(~-1+factor(AEdrug))
  cancerDrug.m <- model.matrix(~-1+factor(cancerDrug))
  
  # Design matrix for interaction effects  
  X <- cbind(AEcancer.m, AEdrug.m, cancerDrug.m)  #      
  colnames(X) <- c(paste("AEcancer", 1:(n.AE*n.cancer), sep=""),  
                   paste("AEdrug", 1:(n.AE*n.drug), sep=""),
                   paste("cancerDrug", 1:(n.cancer*n.drug), sep=""))   # 
  
  nobs = nrow(Z)
  ss = ss # sample size per study: e.g. 100
  nME = ncol(Z)
  nIA = ncol(X)
  
  # coefficients of interaction terms: beta
  b1 = b2 = b3 = b # 1
  beta = rep(0, nIA)
  n.b1 = length(table(AEcancer))
  n.b2 = length(table(AEdrug))
  n.b3 = length(table(cancerDrug))
  
  id.b1_s = sort(sample(1:n.b1, n.b1*spar))
  id.b2_s = sort(sample(1:n.b2, n.b2*spar)) + n.b1
  id.b3_s = sort(sample(1:n.b3, n.b3*spar)) + n.b1 + n.b2

  beta[id.b1_s] <- rep(b1, n.b1*spar)
  beta[id.b2_s] <- rep(b2, n.b2*spar)
  beta[id.b3_s] <- rep(b3, n.b3*spar)

  flag1 = as.numeric(ifelse(AEcancer.m%*%beta[1:n.b1]==1, 1, 0))
  flag2 = as.numeric(ifelse(AEdrug.m%*%beta[n.b1+1:n.b2]==1, 1, 0))
  flag3 = as.numeric(ifelse(cancerDrug.m%*%beta[n.b1+n.b2+1:n.b3]==1, 1, 0))
  
  # coefficients of marginal effects: u
  u <- rnorm(n=nME, mean=mu_u, sd = 0.4)  # -0.5, -1
  
  library(gtools)
  incidence = inv.logit(Z%*%u+X%*%beta)
  Y <- rbinom(n=nobs, size=ss, prob=incidence)
  
  N = rep(ss, length(Y))
  id=1:nobs
  
 
  df=data.frame(study, N, drug, AE, cancer, AEcancer,
                AEdrug, cancerDrug, Y, id, flag1, flag2, flag3);
  # by descending order 
  df <- df[order(-Y),]
}
```

## model fitting 
### scenario: 0% censoring
```
horseshoeJAGS0 <- function(df){

  attach(df)
  J1 <- length(Y)
  n.drug <- length(table(drug))
  n.AE <- length(table(AE))
  n.study <- length(table(study))   
  n.cancer <-length(table(cancer))
  n.AEcancer <-length(table(AEcancer))
  n.AEdrug <- length(table(AEdrug))
  n.cancerDrug <- length(table(cancerDrug))
  
  library(rjags)
  a <- 0.0016; #scale.prior: A = 25
  data <- list (Y=Y, N=N, J1=J1, AE=AE, n.AE=n.AE, study=study, 
                n.study=n.study, cancer=cancer, n.cancer=n.cancer, 
                n.AEcancer=n.AEcancer, AEcancer=AEcancer,
                AEdrug=AEdrug, cancerDrug=cancerDrug, 
                n.AEdrug=n.AEdrug, n.cancerDrug=n.cancerDrug,
                drug=drug, n.drug=n.drug, a=a)
  
  inits <- function() {list (theta.study=rnorm(n.study,-1,0.4), 
                             theta.AE=rnorm(n.AE,-1,0.4),
                             theta.cancer=rnorm(n.cancer,-1,0.4), 
                             theta.drug=rnorm(n.drug,-1,0.4)) }
  parameters <- c("p.beta1.adj", "p.beta2.adj","p.beta3.adj", 
                  "mu.adj", "p.study.adj", "p.drug.adj", 
                  "p.cancer.adj", "p.AE.adj",
                  "beta1", "tau1",
                  "beta2", "tau2", "beta3", "tau3",
                  "sigma.study","sigma.AE", "sigma.cancer", "sigma.drug")
                  
  #system.time(model.1 <- jags.model("model.file.txt", data, inits, n.chains=3, n.adapt=1000))
  #system.time(update(model.1,30000)) 
  #system.time(m1 <- coda.samples(model=model.1, variable.names=parameters, n.iter=30000, thin=3))
  
  #Parallel computing with JAGS
  library(dclone)
  n.cores=3
  cl <- makePSOCKcluster(n.cores)
  tmp <- clusterEvalQ(cl, library(dclone))
  m1 <- jags.parfit(cl = cl, data = data, params = parameters, inits=inits,
                    model = "model.file.txt", 
                    n.chains = 3, 
                    n.adapt = 1000, 
                    n.update = 30000,
                    n.iter = 30000, 
                    thin = 3)
  stopCluster(cl)
  m1
}
```

###  scenario: 80% censoring
```
horseshoeJAGS2 <- function(df){
  attach(df)
  cut = df$Y[nrow(df)*0.2+1]
  # indicator
  df$delta <- c(rep(1, nrow(df)*0.2),
                rep(0, nrow(df)*0.8)) 
  
  df.obs <- df[which(df$delta==1), ]
  df.mis <- df[which(df$delta==0), ]
  df.all <- rbind(df.obs, df.mis); 
  
  Y <- df.obs$Y
  W <- ifelse(df.mis$delta==0, 1, 0)  # left-censored AEs
  cut = rep(cut, nrow(df.mis)) 
  J1 <- length(Y)
  J2 <- length(W)
  
  drug <- df.all$drug
  AE <- df.all$AE
  study <- df.all$study
  cancer <- df.all$cancer
  AEcancer <- df.all$AEcancer
  AEdrug <- df.all$AEdrug
  cancerDrug <- df.all$cancerDrug
  
  n.drug <- length(table(drug))
  n.AE <- length(table(AE))
  n.study <- length(table(study))   
  n.cancer <-length(table(cancer))
  n.AEcancer <-length(table(AEcancer))
  n.AEdrug <- length(table(AEdrug))
  n.cancerDrug <- length(table(cancerDrug))
  
  library(rjags)
  a <- 0.0016; #scale.prior: A = 25
  data <- list (Y=Y, N=N, J1=J1, W=W, cut=cut, J2=J2, 
                AE=AE, n.AE=n.AE, study=study, 
                n.study=n.study, cancer=cancer, n.cancer=n.cancer, 
                n.AEcancer=n.AEcancer, AEcancer=AEcancer,
                AEdrug=AEdrug, cancerDrug=cancerDrug, 
                n.AEdrug=n.AEdrug, n.cancerDrug=n.cancerDrug,
                drug=drug, n.drug=n.drug, a=a)
  
  inits <- function() {list (theta.study=rnorm(n.study,-1,0.4), 
                             theta.AE=rnorm(n.AE,-1,0.4),
                             theta.cancer=rnorm(n.cancer,-1,0.4), 
                             theta.drug=rnorm(n.drug,-1,0.4)) }
  parameters <- c("p.beta1.adj", "p.beta2.adj","p.beta3.adj",
                  "mu.adj", "p.study.adj", "p.drug.adj", "p.cancer.adj", "p.AE.adj",
                  "beta1", "tau1", "beta2", "tau2", "beta3", "tau3",
                  "sigma.study","sigma.AE", "sigma.cancer", "sigma.drug")
              
  #system.time(model.2 <- jags.model("model.file.txt", data, inits, n.chains=3, n.adapt=1000))
  #system.time(update(model.2,30000)) 
  #system.time(m2 <- coda.samples(model=model.2, variable.names=parameters, n.iter=30000, thin=3))
  
  #Parallel computing with JAGS
  library(dclone)
  n.cores=3
  #timer <- proc.time()
  cl <- makePSOCKcluster(n.cores)
  tmp <- clusterEvalQ(cl, library(dclone))
  m2 <- jags.parfit(cl = cl, data = data, params = parameters, inits=inits,
                    model = "model.file.txt", 
                    n.chains = 3, 
                    n.adapt = 1000, 
                    n.update = 30000,
                    n.iter = 30000, 
                    thin = 3)
  stopCluster(cl)
  #time.taken <- proc.time() - timer
  #timings <- time.taken[3] 
  m2
}

```
