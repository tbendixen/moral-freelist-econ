########################################################################
### "The role of ``morality'' in moralistic supernatural punishment"
###
### Pre-registration; Parameter recovery for RAG
###
### Script by Theiss Bendixen & Benjamin G. Purzycki

rm(list=ls())

### Load packages
library(rethinking) # prepare Stan code and data
library(rstanarm) # fitting Stan model
library(cmdstanr) # fitting Stan model
library(bayesplot) # convenience plotting
library(mice) # to ampute simulated data
library(loo) # model comparison

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

memory.limit(size=56000)

### unblock to load simulated data and fitted models
# base::load("binom_sim_testdat.RData")

### Simulate data set
set.seed(1992)

N <- 2600 # total number of observations: no. of ids*2, since each participant played two rounds (self and local) of the RAG
n <- N/2 # number of individuals

# simulate covariates
P <- as.integer(sample(c(0,1,2), n, replace = TRUE)) |> rep(2) # punishment
O <- as.integer(sample(c(0,1,2), n, replace = TRUE)) |> rep(2) # omniscience
C <- as.numeric(sample(0:8, n, replace = TRUE)) |> rep(2) # kids
S <- as.integer(rbinom(n, size = 1, inv_logit(0.01+0.1*C))) |> rep(2) # food security
M <- as.numeric(sample(seq(0,0.8, by=0.2), n, replace = TRUE)) + 0.2*S # free-listed morality
check <- as.integer(sample(0:1, n, replace = TRUE, prob = c(0.9, 0.1)))  |> rep(2)# no effect of check
order <- as.integer(sample(0:1, n, replace = TRUE))  |> rep(2) # no effect of order

id <- rep(seq(1:n), 2) # no effect of id
type <- c(rep(0, n), rep(1, n)) # no effect of type
group <- as.integer(c(rep(seq(1, 10), length.out=n), rep(seq(1, 10), length.out=n))) # no effect of group

simdat <- data.frame(id=id,
                     group=group,
                     P=P,
                     O=O,
                     C=C,
                     S=S,
                     M=M,
                     check=check,
                     type=type,
                     order=order)
# setting coefs
b0 <- 0.5 # intercept
bM <- 0.3
bO <- 0.3
bP <- 0.3
bC <- -0.3
bS <- -0.3
bMP <- 0.3
bMO <- 0.3

y <- with(simdat, 
          rbinom(N, size = 30, inv_logit(b0 + bM * M + bO * O + bP * P + bC * C + bS * S + bMP*M*P + bMO*M*O))
          )

hist(y)

### Simulate missing data
y_notmiss <- which( !is.na(y)) # we do not simulate missing outcome values, but in the real data set we'll exclude missing NAs in the outcome
dat_list_full <- with(simdat, list(
  M = M[y_notmiss],
  P = P[y_notmiss],
  O = O[y_notmiss],
  S = S[y_notmiss],
  C = C[y_notmiss],
  order = order[y_notmiss],
  check = check[y_notmiss]
))

dat_df <- as.data.frame(dat_list_full)

dat_df_amp <- ampute(dat_df, prop = 0.1)

dat_list <- as.list(dat_df_amp$amp)

dat_list$id <- id[y_notmiss]
dat_list$group <- group[y_notmiss]
dat_list$y <- y[y_notmiss]
dat_list$type <- as.integer(type[y_notmiss])
dat_list$N <- length(y[y_notmiss])

str(dat_list)

### quick'n'dirty parameter recovery
freq_dat <- with(dat_list, data.frame(y=y, M=M, O=O, P=P, C=C, S=S, order=order, check=check, type=type))

freq_dat$y_prop <- freq_dat$y/30

# full interaction model
quick_int_binom <- glm(y_prop ~ M*O*P+S+C+order+check+type, data = freq_dat, 
                   family = "quasibinomial") # quasibinomial analyzes the proportion of 'succeses'/'trials

summary(quick_int_binom)
confint(quick_int_binom)

# additive model
quick_noint_binom <- glm(y_prop ~ M+O+P+S+C+order+check+type, data = freq_dat, 
                       family = "quasibinomial") # quasibinomial analyzes the proportion of 'succeses'/'trials

summary(quick_noint_binom)
confint(quick_noint_binom)

### Prepare interaction model and model data with rethinking

m_binom_int_prep <- map2stan(
  alist(
    ## coin model
    y ~ dbinom(30,p),

    logit(p) <- a + zi[id]*sigma_id + aj[group] + # z standardizes adaptive prior for varying effects on individuals
                (bM+bMj[group])*M + (bO+bOj[group])*O + (bP+bPj[group])*P +
                (bMP+bMPj[group])*M*P + (bMO+bMOj[group])*M*O +
                (bPO+bPOj[group])*P*O + (bMPO+bMPOj[group])*M*P*O +
                bC*C + bS*S +
                border*order + btype*type + bcheck*check,

    ## morality model
    M ~ normal(M_mu,M_sd),
    M_mu <- Mavg[group],
    M_sd ~ exponential(1),
    Mavg[group] ~ normal(Mu_Mavg,sigmaMavg),
    Mu_Mavg ~ normal(0.5,2),
    sigmaMavg ~ exponential(1),
    
    ## pun model
    P ~ normal(P_mu,P_sd),
    P_mu <- Pavg[group],
    P_sd ~ exponential(1),
    Pavg[group] ~ normal(Mu_Pavg,sigmaPavg),
    Mu_Pavg ~ normal(1,2),
    sigmaPavg ~ exponential(1),
    
    ## omni model
    O ~ normal(O_mu,O_sd),
    O_mu <- Oavg[group],
    O_sd ~ exponential(1),
    Oavg[group] ~ normal(Mu_Oavg,sigmaOavg),
    Mu_Oavg ~ normal(1,2),
    sigmaOavg ~ exponential(1),
    
    ## children
    C ~ normal(C_mu,C_sd), # >=0 constraint imposed later
    C_mu <- C[group],
    C_sd ~ exponential(10),
    Cavg[group] ~ normal(Mu_Cavg,sigmaCavg),
    Mu_Cavg ~ normal(1,2),
    sigmaCavg ~ exponential(1),
    
    ## priors
    a ~ normal(0,1.5),
    c(bM, bP, bO, bMP, bMO, bPO, bMPO, bC, bS, border, bcheck, btype, bMavg, bPavg, bOavg) ~ normal(0,1),
    ## varying intercepts and slopes
    c(aj, bMj , bPj , bOj, bMPj, bMOj, bPOj, bMPOj)[group] ~ dmvnormNC( Sigmaj , Rhoj ),
    Sigmaj ~ dexp(1),
    Rhoj ~ dlkjcorr(4),
    
    # individual varying intercepts
    zi[id] ~ normal(0,1),
    sigma_id ~ exponential(1),
    
    ## imputation distributions below
    S ~ bernoulli(phi_S),
    phi_S ~ beta(1,1),
    order ~ bernoulli(phi_order),
    phi_order ~ beta(1,1),
    check ~ bernoulli(phi_check),
    phi_check ~ beta(1,1),
    type ~ bernoulli(phi_type),
    phi_type ~ beta(1,1)
  ),
  
  constraints=list(
    sigma_id="lower=1",
    phi_check="lower=0,upper=1",
    phi_order="lower=0,upper=1",
    phi_type="lower=0,upper=1",
    phi_S="lower=0,upper=1",
    C_impute="lower=0",
    order_impute="lower=0,upper=1",
    check_impute="lower=0,upper=1",
    type_impute="lower=0,upper=1",
    S_impute="lower=0,upper=1",
    M_impute="lower=0,upper=1",
    P_impute="lower=0,upper=2",
    O_impute="lower=0,upper=2"
  ),
  data=dat_list ,
  sample=FALSE ,
  do_discrete_imputation = TRUE)

### Add to Stan code:
# in generated quantities (top):
# vector[N] yrep;
# vector[N] log_lik;
# in generated quantities (pi needs to be defined/computed before it can be used here)
# for ( i in 1:N ) yrep[i] = binomial_rng(30, inv_logit(p[i]));
# for ( i in 1:N ) log_lik[i] = binomial_logit_lpmf( y[i] | 30 , p[i] );

# extract raw Stan code
m_binom_int_code <- rethinking::stancode(m_binom_int_prep)

# amend raw Stan code
m_binom_int_code_gq <- gsub("vector[N] S_impute;", 
                     "vector[N] S_impute;
                      vector[N] yrep;
                      vector[N] log_lik;", m_binom_int_code, fixed=TRUE)

m_binom_int_stan <- gsub("check_impute[i] = check[i];\n        }\n    }//i", 
                     "check_impute[i] = check[i];\n        }\n    }//i
                      for ( i in 1:N ) yrep[i] = binomial_rng(30, inv_logit(p[i]));
                      for ( i in 1:N ) log_lik[i] = binomial_logit_lpmf( y[i] | 30 , p[i] );", m_binom_int_code_gq, fixed=TRUE)

### fit model in CmdStan
# library(cmdstanr) # fitting Stan model # load cmdstanr
# set_cmdstan_path("C:\\cmdstan\\cmdstan-2.29.0") # set path to whereever cmdstan in installed
# m_binom_int_cmdstan <- write_stan_file(m_binom_int_stan)
# m_binom_int_cmdstan_model <- cmdstan_model(m_binom_int_cmdstan, compile=TRUE)
# m_binom_int_cmdstan_fit <- m_binom_int_cmdstan_model$sample(
#  data = m_binom_int_prep$data,
#  seed = 123,
#  chains = 4,
#  parallel_chains = 4)

### fit model in Rstan
iter <- 2000
chains <- 4

m_binom_int_fit <- rstan::stan( 
  model_code = m_binom_int_stan, 
  data = m_binom_int_prep$data,
  init = "0", # help sampler get started
  chains = chains, iter = iter, cores = 4,
  control = list(adapt_delta=0.99))

### quick MCMC diagnostic check
rstan::check_hmc_diagnostics(m_binom_int_fit)

### posterior predictive check
m_binom_int_yrep <- rstan::extract(m_binom_int_fit)$yrep

# global
int_ppc_rag <- ppc_bars(y = m_binom_int_prep$data[["y"]],
         yrep = m_binom_int_yrep,
         size = 0.2,
         freq = FALSE) + # easier to compare groups with proportions, instead of counts
  xlab("Number of coins allocated to DISTANT cup") +
  theme(legend.position = "none")# remove legend

# site-specific
ppc_bars_grouped(y = m_binom_int_prep$data[["y"]],
                 yrep = m_binom_int_yrep,
                 group = m_binom_int_prep$data[["group"]],
                 size = 0.2,
                 facet_args = list(ncol=3), # set number of columns
                 freq = FALSE) + # easier to compare groups with proportions, instead of counts
                 xlab("Number of coins allocated to DISTANT cup") +
                 theme(legend.position = "none")# remove legend

### coef plot
int_pars <- c("bM", "bO", "bP", "bMP", "bMO", "bMPO", "bS", "bC")
int_posterior <- as.data.frame(m_binom_int_fit)[int_pars]
int_post <- int_posterior
int_post <- as.matrix(int_post)

int_coef_rag <- mcmc_areas(int_post,
                           pars = int_pars,
                           prob = 0.8) + xlab("Raw coefficient (log odds)") + ylab("Focal parameters")

### model comparison metrics
m_binom_int_loo <- loo(m_binom_int_fit)

### Prepare non-interaction model and model data with rethinking

m_binom_noint_prep <- map2stan(
  alist(
    ## coin model
    y ~ dbinom(30,p),
    
    logit(p) <- a + zi[id]*sigma_id + aj[group] + # z standardizes adaptive prior for varying effects on individuals
      (bM+bMj[group])*M + (bO+bOj[group])*O + (bP+bPj[group])*P +
      bC*C + bS*S +
      border*order + btype*type + bcheck*check,
    
    ## morality model
    M ~ normal(M_mu,M_sd),
    M_mu <- Mavg[group],
    M_sd ~ exponential(1),
    Mavg[group] ~ normal(Mu_Mavg,sigmaMavg),
    Mu_Mavg ~ normal(0.5,2),
    sigmaMavg ~ exponential(1),
    
    ## pun model
    P ~ normal(P_mu,P_sd),
    P_mu <- Pavg[group],
    P_sd ~ exponential(1),
    Pavg[group] ~ normal(Mu_Pavg,sigmaPavg),
    Mu_Pavg ~ normal(1,2),
    sigmaPavg ~ exponential(1),
    
    ## omni model
    O ~ normal(O_mu,O_sd),
    O_mu <- Oavg[group],
    O_sd ~ exponential(1),
    Oavg[group] ~ normal(Mu_Oavg,sigmaOavg),
    Mu_Oavg ~ normal(1,2),
    sigmaOavg ~ exponential(1),
    
    ## children
    C ~ normal(C_mu,C_sd), # >=0 constraint imposed later
    C_mu <- C[group],
    C_sd ~ exponential(10),
    Cavg[group] ~ normal(Mu_Cavg,sigmaCavg),
    Mu_Cavg ~ normal(1,2),
    sigmaCavg ~ exponential(1),
    
    ## priors
    a ~ normal(0,1.5),
    c(bM, bP, bO, bC, bS, border, bcheck, btype, bMavg, bPavg, bOavg) ~ normal(0,1),
    ## varying intercepts and slopes
    c(aj, bMj , bPj , bOj)[group] ~ dmvnormNC( Sigmaj , Rhoj ),
    Sigmaj ~ dexp(1),
    Rhoj ~ dlkjcorr(4),
    
    # individual varying intercepts
    zi[id] ~ normal(0,1),
    sigma_id ~ exponential(1),
    
    ## imputation distributions below
    S ~ bernoulli(phi_S),
    phi_S ~ beta(1,1),
    order ~ bernoulli(phi_order),
    phi_order ~ beta(1,1),
    check ~ bernoulli(phi_check),
    phi_check ~ beta(1,1),
    type ~ bernoulli(phi_type),
    phi_type ~ beta(1,1)
  ),
  
  constraints=list(
    sigma_id="lower=1",
    phi_check="lower=0,upper=1",
    phi_order="lower=0,upper=1",
    phi_type="lower=0,upper=1",
    phi_S="lower=0,upper=1",
    C_impute="lower=0",
    order_impute="lower=0,upper=1",
    check_impute="lower=0,upper=1",
    type_impute="lower=0,upper=1",
    S_impute="lower=0,upper=1",
    M_impute="lower=0,upper=1",
    P_impute="lower=0,upper=2",
    O_impute="lower=0,upper=2"
  ),
  data=dat_list ,
  sample=FALSE ,
  do_discrete_imputation = TRUE)

### Add to Stan code:
# in generated quantities (top):
# vector[N] yrep;
# vector[N] log_lik;
# in generated quantities (pi needs to be defined/computed before it can be used here)
# for ( i in 1:N ) yrep[i] = binomial_rng(30, inv_logit(p[i]));
# for ( i in 1:N ) log_lik[i] = binomial_logit_lpmf( y[i] | 30 , p[i] );

# extract raw Stan code
m_binom_noint_code <- rethinking::stancode(m_binom_noint_prep)

# amend raw Stan code
m_binom_noint_code_gq <- gsub("vector[N] S_impute;", 
                              "vector[N] S_impute;
                              vector[N] yrep;
                              vector[N] log_lik;", m_binom_noint_code, fixed=TRUE)

m_binom_noint_stan <- gsub("check_impute[i] = check[i];\n        }\n    }//i", 
                         "check_impute[i] = check[i];\n        }\n    }//i
                         for ( i in 1:N ) yrep[i] = binomial_rng(30, inv_logit(p[i]));
                         for ( i in 1:N ) log_lik[i] = binomial_logit_lpmf( y[i] | 30 , p[i] );", m_binom_noint_code_gq, fixed=TRUE)

### fit model in CmdStan
# library(cmdstanr) # fitting Stan model # load cmdstanr
# set_cmdstan_path("C:\\cmdstan\\cmdstan-2.29.0") # set path to whereever cmdstan in installed
# m_binom_noint_cmdstan <- write_stan_file(m_binom_noint_stan)
# m_binom_noint_cmdstan_model <- cmdstan_model(m_binom_noint_cmdstan, compile=TRUE)
# m_binom_noint_cmdstan_fit <- m_binom_noint_cmdstan_model$sample(
#  data = m_binom_noint_prep$data,
#  seed = 123,
#  chains = 4,
#  parallel_chains = 4)

### fit model in Rstan
m_binom_noint_fit <- rstan::stan( 
  model_code = m_binom_noint_stan, 
  data = m_binom_noint_prep$data,
  init = "0", # help sampler get started
  chains = chains, iter = iter, cores = 4,
  control = list(adapt_delta=0.99))

### posterior predictive check
m_binom_noint_yrep <- extract(m_binom_noint_fit)$yrep

# global
noint_ppc_rag <- ppc_bars(y = m_binom_noint_prep$data[["y"]],
         yrep = m_binom_noint_yrep,
         size = 0.2,
         freq = FALSE) + # easier to compare groups with proportions, instead of counts
  xlab("Number of coins allocated to DISTANT cup") +
  theme(legend.position = "none")# remove legend

# site-specific
ppc_bars_grouped(y = m_binom_noint_prep$data[["y"]],
                 yrep = m_binom_noint_yrep,
                 group = m_binom_noint_prep$data[["group"]],
                 size = 0.2,
                 facet_args = list(ncol=3), # set number of columns
                 freq = FALSE) + # easier to compare groups with proportions, instead of counts
                 xlab("Number of coins allocated to DISTANT cup") +
                 theme(legend.position = "none")# remove legend

### coef plot
noint_pars <- c("bM", "bO", "bP", "bS", "bC")
noint_posterior <- as.data.frame(m_binom_noint_fit)[noint_pars]
noint_post <- noint_posterior
noint_post <- as.matrix(noint_post)

noint_coef_rag <- mcmc_areas(noint_post,
           pars = noint_pars,
           prob = 0.8) + xlab("Raw coefficient (log odds)") + ylab("Focal parameters")

### model comparison metrics
m_binom_noint_loo <- loo(m_binom_noint_fit)

# compare models
model_comparison <- loo_compare(m_binom_int_loo, m_binom_noint_loo)

### export environment
save.image("binom_sim_testdat.RData")

### END
