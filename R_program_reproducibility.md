### The following R code is to reproduce the case study

```
# load the spreadsheet
#setwd()

data <- read.csv("FinalAEs.csv", header=T)
data35 <- data[ ,c(1:12, seq(13,162, by=2)+1)]
AE.group <- read.csv("AElist.csv", h=T)[, "Group2"]
group.name = rep(AE.group, each=nrow(data35))

library(reshape2)
datal35 <- melt(data35, id.vars = c(1,4:10),
                measure.vars = 13:87,
                variable.name = "AE",       
                value.name = "Y")
datal35 = cbind(datal35, group.name)

datal35$cutoff.for.grades.3.or.higher[datal35$Trial==15 & datal35$group.name=="Endocrine"] <- 1 
datal35$cutoff.for.grades.3.or.higher[datal35$Trial==82 & datal35$group.name=="Endocrine"] <- 0 
datal35$cutoff.for.grades.3.or.higher[datal35$Trial==84 & datal35$group.name=="Endocrine"] <- 0 
datal35$cutoff.for.grades.3.or.higher[datal35$Trial==86 & datal35$group.name=="Endocrine"] <- 2 
datal35$cutoff.for.grades.3.or.higher[datal35$Trial==15 & datal35$group.name=="Others"] <- 1
datal35$cutoff.for.grades.3.or.higher[datal35$Trial==82 & datal35$group.name=="Others"] <- 0
datal35$cutoff.for.grades.3.or.higher[datal35$Trial==84 & datal35$group.name=="Others"] <- 0
datal35$cutoff.for.grades.3.or.higher[datal35$Trial==86 & datal35$group.name=="Others"] <- 2

datal35$Y[is.na(datal35$Y) & datal35$cutoff.for.grades.3.or.higher==0] <- 0
datal35$Y[is.na(datal35$Y) & datal35$cutoff.for.grades.3.or.higher>=1] <- NA

datal35$R <- ifelse(is.na(datal35$Y), 0, 1)
datal35$L <- ifelse(datal35$R==0, datal35$cutoff.for.grades.3.or.higher, -1)
datal35$U <-  datal35$No..of.treated.patients

data.obs <- datal35[!is.na(datal35$Y), ] # 
data.mis <- datal35[is.na(datal35$Y), ]  # 
data.all <- rbind(data.obs, data.mis)    #

W <- rep(1, dim(data.mis)[1]) # RC: W=0, LC: W=1
cut <- data.mis$L 

Y <- as.numeric(data.obs$Y) # true number of AEs (observed)
W <- as.numeric(W)          # "missing" AEs
cut <- as.numeric(cut)
J1 <- length(Y)
J2 <- length(W)
N <- as.numeric(data.all$No..of.treated.patients)  # number of pts

study <- as.numeric(as.factor(data.all$Trial))
drug <- as.numeric(factor(data.all$Drug, 
                          level = c("Nivolumab", "Pembrolizumab",
                                    "Atezolizumab", "Avelumab",
                                    "Durvalumab"))) # 5 levels
AE <- as.numeric(as.factor(as.numeric(data.all$AE))) #
cancer <- as.numeric(factor(data.all$Cancer.type, 
                            level = c("Lung cancer","GU cancer",
                                      "Melanoma", "Other cancer",
                                      "Mixed cancer types", "GI cancer", 
                                      "Hematologic malignancy"))) # 

AEcancer = as.numeric(factor(cancer*100+AE))   # 
AEdrug = as.numeric(factor(drug*100+AE))       # 
cancerDrug <- as.numeric(factor(cancer*10+drug))  #

n.study <- length(table(study))   
n.drug <- length(table(drug))    
n.AE <- length(table(AE))        
n.cancer <-length(table(cancer))  

n.AEcancer <-length(table(AEcancer))
n.AEdrug <- length(table(AEdrug))
n.cancerDrug <- length(table(cancerDrug))

set.seed(123)
library(rjags)
a <- 0.0016; #scale.prior: A = 25
data <- list (Y=Y, N=N, J1=J1, W=W, cut=cut, J2=J2, 
              AE=AE, n.AE=n.AE, study=study, 
              n.study=n.study, cancer=cancer, n.cancer=n.cancer, 
              n.AEcancer=n.AEcancer, AEcancer=AEcancer,
              AEdrug=AEdrug, cancerDrug=cancerDrug, 
              n.AEdrug=n.AEdrug, n.cancerDrug=n.cancerDrug,
              drug=drug, n.drug=n.drug, a=a)
inits <- function() {list (theta.study=rnorm(n.study,0,1), 
                           theta.AE=rnorm(n.AE,0,1),
                           theta.cancer=rnorm(n.cancer,0,1), 
                           theta.drug=rnorm(n.drug,0,1)) }
parameters <- c("p.beta1.adj", "p.beta2.adj","p.beta3.adj",
                "mu.adj", "p.study.adj", "p.drug.adj", 
                "p.cancer.adj", "p.AE.adj",
                "beta1", "tau1","beta2", "tau2", "beta3", "tau3",
                "sigma.study","sigma.AE", "sigma.cancer", "sigma.drug")

system.time(model.1 <- jags.model ("model.txt", data, inits, n.chains=3, n.adapt=1000))
system.time(update(model.1,30000)) # 30000
system.time(fit.1 <- coda.samples(model=model.1, variable.names=parameters,n.iter=30000, thin=3))

```
