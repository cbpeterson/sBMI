# Bayesian sparse modeling of censored safety data in meta-analysis


### Objective:
A sparse Bayesian model with interaciton selection using horsehsoe prior (sBMI) is proposed to simultaneously identify nonzero interactions and high-risk groups with an elevated probability of adverse events and to address a key challenge in the meta-analysis of safety data.

### Model specification in JAGS
```
model {
        # Observed data
        for (j in 1:J1) {
        # Likelihood
        Y[j] ~ dbin(theta[j], N[j])
        logit(theta[j]) <- theta.v1[v1[j]] + theta.v2[v2[j]] + theta.v3[v3[j]] +
                            theta.v4[v4[j]] + beta1[v5[j]] + beta2[v6[j]] + beta3[v7[j]]
        }
        
        # Censored data
        for (j in 1:J2) {
        # Likelihood
        W[j] ~ dbern(p[j])
        p[j] <- pbin(cut[j], theta[j+J1], N[j+J1]) # Pr(Y<=cut)
        logit(theta[j+J1]) <- theta.v1[v1[j+J1]] + theta.v2[v2[j+J1]] + theta.v3[v3[j+J1]] +
                              theta.v4[v4[j+J1]] + beta1[v5[j+J1]] + beta2[v6[j+J1]] + beta3[v7[j+J1]]
        }
        
        for (i1 in 1:n.v1) {
        theta.v1[i1] ~ dnorm(0, prec.v1)
        } # a half-Cauchy prior on standard deviation
        prec.v1 <- pow(sigma.v1, -2)
        # a=1/A^2, where scale parameter A=25
        sigma.v1 ~ dt(0, a, 1)T(0,)
        
        ... # same priors on other main coefficients: v2, v3, v4
        
        for (k1 in 1:n.v5) {
        beta1[k1] ~ dnorm(0, prec1[k1])
        # precision = 1/variance
        prec1[k1] <- pow(sigma1[k1], -2)
        sigma1[k1] <- lambda1[k1] * tau1
        # Local shrinkage parameters of horseshoe prior
        lambda1[k1] ~ dt(0, 1, 1)T(0,)
        }
        # Global shrinkage parameter of horseshoe prior
        tau1 ~ dt(0, 1, 1)T(0,)
        
        ... # same priors on other interaction coefficients: v6, v7
}

```
The R code to apply the proposed method, sBMI, to simulation and case study is available in the ``R_program.md``.
