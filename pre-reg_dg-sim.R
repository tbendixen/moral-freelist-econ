########################################################################
### "The role of ``morality'' in moralistic supernatural punishment"
###
### Pre-registration; Parameter recovery for DG
###
### Script by Theiss Bendixen & Benjamin G. Purzycki

rm(list=ls())

### Load packages
library(rethinking) # prepare Stan code and data
library(rstanarm) # fitting Stan model
library(bayesplot) # convenience plotting
library(mice) # to ampute simulated data
library(ordinal) # for quick 'frequentist' parameter recovery
library(loo) # model comparison

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

memory.limit(size=56000)

### unblock to load simulated data and fitted models
base::load("ord_cat_sim_testdat.RData")

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
bM <- 0.3
bO <- 0.3
bP <- 0.3
bC <- -0.3
bS <- -0.3
bMP <- 0.3
bMO <- 0.3

# simulate outcome
# credit to: https://stats.stackexchange.com/questions/321770/simulating-data-for-an-ordered-logit-model
b01 <- 2
b02 <- 1.5
b03 <- -0.3
b04 <- -0.5
b05 <- -1
b06 <- -2.2
b07 <- -2.3
b08 <- -2.4
b09 <- -2.5
b10 <- -3

logodds1 <- b01 + bM * M + bO * O + bP * P + bC * C + bS * S + bMP*M*P + bMO*M*O
logodds2 <- b02 + bM * M + bO * O + bP * P + bC * C + bS * S + bMP*M*P + bMO*M*O
logodds3 <- b03 + bM * M + bO * O + bP * P + bC * C + bS * S + bMP*M*P + bMO*M*O
logodds4 <- b04 + bM * M + bO * O + bP * P + bC * C + bS * S + bMP*M*P + bMO*M*O
logodds5 <- b05 + bM * M + bO * O + bP * P + bC * C + bS * S + bMP*M*P + bMO*M*O
logodds6 <- b06 + bM * M + bO * O + bP * P + bC * C + bS * S + bMP*M*P + bMO*M*O
logodds7 <- b07 + bM * M + bO * O + bP * P + bC * C + bS * S + bMP*M*P + bMO*M*O
logodds8 <- b08 + bM * M + bO * O + bP * P + bC * C + bS * S + bMP*M*P + bMO*M*O
logodds9 <- b09 + bM * M + bO * O + bP * P + bC * C + bS * S + bMP*M*P + bMO*M*O
logodds10 <- b10 + bM * M + bO * O + bP * P + bC * C + bS * S + bMP*M*P + bMO*M*O

# cumulative probs
prob_2to11 <- inv_logit(logodds1)
prob_3to11 <- inv_logit(logodds2)
prob_4to11 <- inv_logit(logodds3)
prob_5to11 <- inv_logit(logodds4)
prob_6to11 <- inv_logit(logodds5)
prob_7to11 <- inv_logit(logodds6)
prob_8to11 <- inv_logit(logodds7)
prob_9to11 <- inv_logit(logodds8)
prob_10to11 <- inv_logit(logodds9)
prob_11 <- inv_logit(logodds10)

# individual response prob
prob_1 <- 1 - prob_2to11
prob_2 <- prob_2to11 - prob_3to11
prob_3 <- prob_3to11 - prob_4to11
prob_4 <- prob_4to11 - prob_5to11
prob_5 <- prob_5to11 - prob_6to11
prob_6 <- prob_6to11 - prob_7to11
prob_7 <- prob_7to11 - prob_8to11
prob_8 <- prob_8to11 - prob_9to11
prob_9 <- prob_9to11 - prob_10to11
prob_10 <- prob_9to11 - prob_11

y <- c()
for (i in 1:N) {
  y[i] <- sample(
    x = c(1:11), 
    size = 1, 
    prob = c(prob_1[i], prob_2[i], prob_3[i], prob_4[i], prob_5[i],
             prob_6[i], prob_7[i], prob_8[i], prob_9[i], prob_10[i], 
             prob_11[i])
  )
}

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

dat_list$id <- as.integer(id[y_notmiss])
dat_list$group <- group[y_notmiss]
dat_list$y <- y[y_notmiss]
dat_list$type <- as.integer(type[y_notmiss])
dat_list$N <- as.integer(length(y[y_notmiss]))
dat_list$K <- as.integer(max(y[y_notmiss]))

str(dat_list)

### quick'n'dirty parameter recovery
freq_dat <- with(dat_list, data.frame(y=factor(y, ordered = TRUE), M=M, O=O, P=P, C=C, S=S, order=order, check=check, type=type))

# full interaction model
quick_int_clm <- clm(y ~ M*O*P+S+C+order+check+type, data = freq_dat)
summary(quick_int_clm)
confint(quick_int_clm)

# additive model
quick_noint_clm <- clm(y ~ M+O+P+S+C+order+check+type, data = freq_dat)
summary(quick_noint_clm)
confint(quick_noint_clm)

### Prepare interaction model and model data with rethinking

m_cat_int_prep <- map2stan(
  alist(
    ## coin model
    y ~ dordlogit(eta, cutpoints),

    eta <- zi[id]*sigma_id + # z standardizes adaptive prior for varying effects on individuals
      (bM+bMj[group])*M + (bO+bOj[group])*O + (bP+bPj[group])*P +
      (bMP+bMPj[group])*M*P + (bMO+bMOj[group])*M*O +
      (bPO+bPOj[group])*P*O + (bMPO+bMPOj[group])*M*P*O +
      bC*C + bS*S +
      border*order + btype*type + bcheck*check,
    
    ## cutpoints
    cutpoints ~ normal(0,2),
    
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
    c(bM, bP, bO, bMP, bMO, bPO, bMPO, bC, bS, border, bcheck, btype, bMavg, bPavg, bOavg) ~ normal(0,1),
    ## varying intercepts and slopes
    c(bMj , bPj , bOj, bMPj, bMOj, bPOj, bMPOj)[group] ~ dmvnormNC( Sigmaj , Rhoj ),
    Sigmaj ~ dexp(1),
    Rhoj ~ dlkjcorr(4),
    
    ## individual varying effects
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
  sample=FALSE, 
  do_discrete_imputation = TRUE)

### additions/modifications to rethinking-generated Stan code:
# in data block:
# int<lower=2> K;
# in parameter block: 
# ordered[K-1] cutpoints;
# in generated quantities block:
# add to the other vectors:
#   vector[N] yrep;
#   vector[N] log_lik;
# add at the end of the gq block (eta needs to be defined/computed before it can be used here):
#   for (i in 1:N) {
#   yrep[i] = ordered_logistic_rng(eta[i], cutpoints);
#   }
#   for (i in 1:N) {
#   log_lik[i] = ordered_logistic_lpmf( y[i] | eta[i] , cutpoints );
#   }
# 
# This is automated in the following by extracting the model code string and then using gsub():

# extract raw Stan code
m_cat_int_code <- rethinking::stancode(m_cat_int_prep)

# amend raw Stan code
m_cat_int_code_db <- gsub("int<lower=1> N;", 
                     "int<lower=1> N;
                      int<lower=2> K;", m_cat_int_code, fixed=TRUE)

m_cat_int_code_pb <- gsub("ordered cutpoints;", 
                      "ordered[K-1] cutpoints;", m_cat_int_code_db, fixed=TRUE)

m_cat_int_code_gq <- gsub("generated quantities{", 
                      "generated quantities{ 
                      vector[N] yrep; 
                      vector[N] log_lik;", m_cat_int_code_pb, fixed=TRUE)

m_cat_int_stan <- gsub("check_impute[i] = check[i];\n        }\n    }//i", 
                   "check_impute[i] = check[i];\n        }\n    }//i 
                      for (i in 1:N) yrep[i] = ordered_logistic_rng(eta[i], cutpoints);
                      for (i in 1:N) log_lik[i] = ordered_logistic_lpmf( y[i] | eta[i] , cutpoints );", 
                   m_cat_int_code_gq, fixed=TRUE)

### fit Stan model in CmdStan
#library(cmdstanr) # fitting Stan model
#set_cmdstan_path("C:\\cmdstan\\cmdstan-2.29.0")
#m_cat_cmdstan_file <- write_stan_file(m_cat_stan)
#m_cat_cmdstan_model <- cmdstan_model(m_cat_cmdstan_file, compile=TRUE)
#cat_cmdstanfit <- m_cat_cmdstan_model$sample(
#  data = m_cat_sim$data,
#  seed = 123,
#  chains = 4,
#  parallel_chains = 4)

### fit Stan model in RStan
iter <- 2000
chains <- 4

m_cat_int_fit <- rstan::stan( 
  model_code = m_cat_int_stan, 
  data = m_cat_int_prep$data,
  init = "0", # help sampler get started
  chains = chains, iter = iter, cores = 4,
  seed = 2021,
  control = list(adapt_delta=0.99))

### quick MCMC diagnostic check
rstan::check_hmc_diagnostics(m_cat_int_fit)

### posterior predictive checks
m_cat_int_yrep <- extract(m_cat_int_fit)$yrep

# global
int_ppc_dg <- ppc_bars(y = m_cat_int_prep$data[["y"]],
         yrep = m_cat_int_yrep,
         size = 0.2,
         freq = FALSE) + # easier to compare groups with proportions, instead of counts
  xlab("Number of coins allocated to DISTANT cup") +
  theme(legend.position = "none")# remove legend

# site-specific
ppc_bars_grouped(y = m_cat_int_prep$data[["y"]],
                 yrep = m_cat_int_yrep,
                 group = m_cat_int_prep$data[["group"]],
                 size = 0.2,
                 facet_args = list(ncol=3), # set number of columns
                 freq = FALSE) + # easier to compare groups with proportions, instead of counts
                 xlab("Number of coins allocated to DISTANT cup") +
                 theme(legend.position = "none")# remove legend

### coef plot
int_pars <- c("bM", "bO", "bP", "bMP", "bMO", "bMPO", "bS", "bC")
int_posterior <- as.data.frame(m_cat_int_fit)[int_pars]
int_post <- as.matrix(int_posterior)

int_coef_dg <- mcmc_areas(int_post,
           pars = int_pars,
           prob = 0.8) + xlab("Raw coefficient (log odds)") + ylab("Focal parameters")

### plot MCMC diagnostics
posterior_chains <- as.array(m_cat_int_fit)[ndraws, int_pars]
fargs <- list(ncol = 2, labeller = label_parsed)
chains_trace <- mcmc_trace(posterior_chains, pars = int_pars,
                           n_warmup = 300, facet_args = fargs)

### model comparison metrics
m_cat_int_loo <- loo(m_cat_int_fit)

### Prepare non-interaction model and model data with rethinking

m_cat_noint_prep <- map2stan(
  alist(
    ## coin model
    y ~ dordlogit(eta, cutpoints),
    
    eta <- zi[id]*sigma_id + # z standardizes adaptive prior for varying effects on individuals
      (bM+bMj[group])*M + (bO+bOj[group])*O + (bP+bPj[group])*P +
      bC*C + bS*S+
      border*order + btype*type + bcheck*check,
    
    ## cutpoints
    cutpoints ~ normal(0,2),
    
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
    c(bM, bP, bO, bC, bS, border, bcheck, btype, bMavg, bPavg, bOavg) ~ normal(0,1),
    ## varying slopes for groups and h (individual response) on y
    c(bMj , bPj , bOj)[group] ~ dmvnormNC( Sigmaj , Rhoj ),
    Sigmaj ~ dexp(1),
    Rhoj ~ dlkjcorr(4),
    
    ## individual varying effects
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

### additions/modifications to rethinking-generated Stan code:

# in data block:
# int<lower=2> K;

# in parameter block: 
# ordered[K-1] cutpoints;

# in generated quantities block:
# add to the other vectors:
#   vector[N] yrep;
#   vector[N] log_lik;
# add at the end of the gq block (eta needs to be defined/computed before it can be used here):
#   for (i in 1:N) {
#   yrep[i] = ordered_logistic_rng(eta[i], cutpoints);
#   }
#   for (i in 1:N) {
#   log_lik[i] = ordered_logistic_lpmf( y[i] | eta[i] , cutpoints );
#   }

# This is automated in the following by extracting the model code string and then using gsub():

# extract raw Stan code
m_cat_noint_code <- rethinking::stancode(m_cat_noint_prep)

# amend raw Stan code
m_cat_noint_code_db <- gsub("int<lower=1> N;", 
                          "int<lower=1> N;
                      int<lower=2> K;", m_cat_noint_code, fixed=TRUE)

m_cat_noint_code_pb <- gsub("ordered cutpoints;", 
                          "ordered[K-1] cutpoints;", m_cat_noint_code_db, fixed=TRUE)

m_cat_noint_code_gq <- gsub("generated quantities{", 
                          "generated quantities{ 
                      vector[N] yrep; 
                      vector[N] log_lik;", m_cat_noint_code_pb, fixed=TRUE)

m_cat_noint_stan <- gsub("check_impute[i] = check[i];\n        }\n    }//i", 
                       "check_impute[i] = check[i];\n        }\n    }//i 
                      for (i in 1:N) yrep[i] = ordered_logistic_rng(eta[i], cutpoints);
                      for (i in 1:N) log_lik[i] = ordered_logistic_lpmf( y[i] | eta[i] , cutpoints );", 
                       m_cat_noint_code_gq, fixed=TRUE)

### fit Stan model in CmdStan
# library(cmdstanr) # fitting Stan model
# set_cmdstan_path("C:\\cmdstan\\cmdstan-2.29.0")
# m_cat_cmdstan_file <- write_stan_file(m_cat_stan)
# m_cat_cmdstan_model <- cmdstan_model(m_cat_cmdstan_file, compile=TRUE)
# cat_cmdstanfit <- m_cat_cmdstan_model$sample(
#  data = m_cat_sim$data,
#  seed = 123,
#  chains = 4,
#  parallel_chains = 4)

### fit Stan model in RStan
m_cat_noint_fit <- rstan::stan( 
  model_code = m_cat_noint_stan, 
  data = m_cat_noint_prep$data,
  init = "0", # help sampler get started
  chains = chains, iter = iter, cores = 4,
  seed = 2021,
  control = list(adapt_delta=0.99))

### quick MCMC diagnostic check
rstan::check_hmc_diagnostics(m_cat_noint_fit)

### posterior predictive checks
m_cat_noint_yrep <- extract(m_cat_noint_fit)$yrep

# global
noint_ppc_dg <- ppc_bars(y = m_cat_noint_prep$data[["y"]],
                       yrep = m_cat_noint_yrep,
                       size = 0.2,
                       freq = FALSE) + # easier to compare groups with proportions, instead of counts
  xlab("Number of coins allocated to DISTANT cup") +
  theme(legend.position = "none")# remove legend

# site-specific
ppc_bars_grouped(y = m_cat_noint_prep$data[["y"]],
                 yrep = m_cat_noint_yrep,
                 group = m_cat_noint_prep$data[["group"]],
                 freq = FALSE) # easier to compare groups with proportions, instead of counts

### coef plot
noint_pars <- c("bM", "bO", "bP", "bS", "bC")
noint_posterior <- as.data.frame(m_cat_noint_fit)[noint_pars]
noint_post <- as.matrix(noint_posterior)

noint_coef_dg <- mcmc_areas(noint_post,
           pars = noint_pars,
           prob = 0.8) + xlab("Raw coefficient (log odds)") + ylab("Focal parameters")

### model comparison metrics
m_cat_noint_loo <- loo(m_cat_noint_fit)

# compare models
model_comparison <- loo_compare(m_cat_int_loo, m_cat_noint_loo)

### save environment
save.image("ord_cat_sim_testdat.RData")

### END 
