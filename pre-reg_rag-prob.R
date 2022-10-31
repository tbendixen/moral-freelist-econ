library(tidyr)
library(dplyr)
library(ggplot2)
library(tidybayes)

rm(list=ls())

memory.limit(size=56000)

set.seed(1992)

load("binom_sim_testdat.RData")

fit <- m_binom_int_fit

ndraws <- 1:nrow(fit) # all draws from the posterior

a <- rstan::extract(fit,"a")$a[ndraws]
aj <- rstan::extract(fit,"aj")$aj[ndraws,]
bM <- rstan::extract(fit,"bM")$bM[ndraws]
bMj <- rstan::extract(fit,"bMj")$bMj[ndraws,]
bO <- rstan::extract(fit,"bO")$bO[ndraws]
bOj <- rstan::extract(fit,"bOj")$bOj[ndraws,]
bP <- rstan::extract(fit,"bP")$bP[ndraws]
bPj <- rstan::extract(fit,"bPj")$bPj[ndraws,]
bMO <- rstan::extract(fit,"bMO")$bMO[ndraws]
bMOj <- rstan::extract(fit,"bMOj")$bMOj[ndraws,]
bMP <- rstan::extract(fit,"bMP")$bMP[ndraws]
bMPj <- rstan::extract(fit,"bMPj")$bMPj[ndraws,]
bPO <- rstan::extract(fit,"bPO")$bPO[ndraws]
bPOj <- rstan::extract(fit,"bPO")$bPO[ndraws]
bMPO <- rstan::extract(fit,"bMPO")$bMPO[ndraws]
bMPOj <- rstan::extract(fit,"bMPOj")$bMPOj[ndraws,]
border <- rstan::extract(fit,"border")$border[ndraws]
btype <- rstan::extract(fit,"btype")$btype[ndraws]
bcheck <- rstan::extract(fit,"bcheck")$bcheck[ndraws]
bC <- rstan::extract(fit,"bC")$bC[ndraws]
bS <- rstan::extract(fit,"bS")$bS[ndraws]

### Marginalizing over predictors
# 1) (marginal) prediction matrix (for "average treatment effect"/g-computation-ish calculation)

M_pred <- with(as.data.frame(dat_list),
               c(rep(0, length(id)), 
                 rep(1, length(id))))

O_pred <- with(as.data.frame(dat_list),
               c(rep(0, length(id)), 
                 rep(2, length(id))))

P_pred <- with(as.data.frame(dat_list),
               c(rep(0, length(id)), 
                 rep(2, length(id))))

levels <- length(unique(M_pred)) # how many levels of the focal predictor?

pred_ate_group <-  with(as.data.frame(dat_list),
                        data.frame(
                          id = rep(id, levels),
                          M = M_pred, # levels of focal predictor
                          group = rep(group, levels),
                          O = rep(O, levels), # leave omniscience rating as observed
                          P = rep(P, levels), # leave punishment rating as observed
                          #O = O_pred, # include min. and max. omniscience in focal g-computation; will throw an error when imputing below
                          #P = P_pred,#  include min. and max. punishment in focal g-computation; will throw an error when imputing below
                          C = rep(C, levels),
                          S = rep(S, levels),
                          check = rep(check, levels),
                          order = rep(order, levels),
                          type = rep(type, levels))
)

# plug posterior mean of imputations for variables with missingness into prediction grid,
# so that we get predictions for these individuals even if they had missing values

O_impute <- rstan::extract(fit,"O_impute")$"O_impute"[ndraws,]
O_impute_mean <- apply(O_impute, 2, mean)
O_impute_idx <- which(is.na(pred_ate_group$O))
pred_ate_group[O_impute_idx,]$O <- O_impute_mean

P_impute <- rstan::extract(fit,"P_impute")$P_impute[ndraws,]
P_impute_mean <- apply(P_impute, 2, mean)
P_impute_idx <- which(is.na(pred_ate_group$P))
pred_ate_group[P_impute_idx,]$P <- P_impute_mean

C_impute <- rstan::extract(fit,"C_impute")$C_impute[ndraws,]
C_impute_mean <- apply(C_impute, 2, mean)
C_impute_idx <- which(is.na(pred_ate_group$C))
pred_ate_group[C_impute_idx,]$C <- C_impute_mean

S_impute <- rstan::extract(fit,"S_impute")$S_impute[ndraws,]
S_miss_idx <- which(is.na(dat_list$S))
S_impute_mean <- apply(S_impute[,S_miss_idx], 2, mean)
S_impute_idx <- which(is.na(pred_ate_group$S))
pred_ate_group[S_impute_idx,]$S <- S_impute_mean

order_impute <- rstan::extract(fit,"order_impute")$order_impute[ndraws,]
order_miss_idx <- which(is.na(dat_list$order))
order_impute_mean <- apply(order_impute[,order_miss_idx], 2, mean)
order_impute_idx <- which(is.na(pred_ate_group$order))
pred_ate_group[order_impute_idx,]$order <- order_impute_mean

check_impute <- rstan::extract(fit,"check_impute")$check_impute[ndraws,]
check_miss_idx <- which(is.na(dat_list$check))
check_impute_mean <- apply(check_impute[,check_miss_idx], 2, mean)
check_impute_idx <- which(is.na(pred_ate_group$check))
pred_ate_group[check_impute_idx,]$check <- check_impute_mean

# check number of missing values; should be zero
sum(is.na(pred_ate_group))

# 2) (marginal) prediction matrix

# the multilevel - multilevel linear model, group-level, ignoring id random effect

pred_ate_group <- pred_ate_group[rep(seq_len(nrow(pred_ate_group)), each = length(ndraws)), ] # expand prediction matrix

group_levels <- length(unique(dat_list$group))
times <- nrow(pred_ate_group)/length(ndraws)
timesj <- (nrow(pred_ate_group)/group_levels/length(ndraws))

# interaction model
pred_coef_group <- data.frame(a=rep(a, times), aj=rep(as.vector(aj),timesj),
                              bM=rep(bM, times), bO=rep(bO, times), bP=rep(bP, times),
                              bMj=rep(as.vector(bMj), timesj),bOj=rep(as.vector(bOj), timesj),bPj=rep(as.vector(bPj), timesj),
                              bMOj=rep(as.vector(bMOj), timesj),bMPj=rep(as.vector(bMPj), timesj), bPOj=rep(bPOj, times), bMPOj=rep(bMPOj, times),
                              bMO=rep(bMO, times), bMP=rep(bMP, times), bPO=rep(bPO, times), bMPO=rep(bMPO, times),
                              bC=rep(bC, times), bS=rep(bS, times), 
                              border=rep(border, times), btype=rep(btype, times), bcheck=rep(bcheck,times),
                              .draw=rep(1:length(ndraws), times))

predp_fun <- function(x) {
  p <- with(x, plogis(a + aj + (bM+bMj)*M + 
                (bO+bOj)*O + (bP + bPj)*P +
                (bMP+bMPj)*M*P + (bMO+bMOj)*M*O + (bPO+bPOj)*P*O +
                (bMPO+bMPOj)*M*P*O +
                bC*C + bS*S +
                border*order + btype*type + bcheck*check))
}

# additive model
# pred_coef_group <- data.frame(a=rep(a, times), aj=rep(as.vector(aj),timesj),
#                              bM=rep(bM, times), bO=rep(bO, times), bP=rep(bP, times),
#                              bMj=rep(as.vector(bMj), timesj),bOj=rep(as.vector(bOj), timesj),bPj=rep(as.vector(bPj), timesj),
#                              bC=rep(bC, times), bS=rep(bS, times),
#                              border=rep(border, times), btype=rep(btype, times), bcheck=rep(bcheck,times),
#                              .draw=rep(1:length(ndraws), times))
# predp_fun <- function(x) {
#  p <- with(x, plogis(a + aj + (bM+bMj)*M + 
#                (bO+bOj)*O + (bP + bPj)*P +
#                bC*C + bS*S +
#                border*order + btype*type + bcheck*check))
# }


# get predictions 
pred_mat_group <- cbind(pred_ate_group, pred_coef_group, row.names = NULL)

for (group in unique(pred_mat_group$group)){  # for each category that is unique in the culture variable
  sub <- subset(pred_mat_group, group == group) |> as.data.frame() # subset
  sub$p <- predp_fun(sub) # compute probabiltiy of coin allocation
  return(sub) # return data frame
}

# prep plot 

postpredcoin_select <- sub %>% select(c(M, id, group, .draw, p))

M0coins <- postpredcoin_select[postpredcoin_select$M==0,]
colnames(M0coins) <- paste0("M0", colnames(M0coins))
M1coins <- postpredcoin_select[postpredcoin_select$M==1,]
colnames(M1coins) <- paste0("M1", colnames(M1coins))

ridgeplotcoin <- cbind(M0coins, M1coins)
ridgeplotcoin <- ridgeplotcoin |> group_by(M0id, M0.draw) |> mutate(diff = M1p-M0p) |> ungroup()

cairo_pdf("binridge_global.pdf",
          width = 7, height = 5) # start print to pdf

# ridge plot (normalized): contrast of posterior predicted probabilities in percentage points difference
ggplot(ridgeplotcoin, aes(x = diff)) +
  
  stat_halfeye(aes(slab_alpha = stat(f),
               fill = stat(x > 0)), 
               normalize="all",
               orientation = "horizontal",
               .width = 0.95,
               size = 1,
               point_interval = "mean_hdi", # summarize the distribution with its mean and highest posterior density
               show.legend=FALSE) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("gray40", "#3182BD")) +
  theme_classic() + 
  #theme(legend.position = "none") +
  scale_x_continuous(breaks = seq(-0.25,0.50, by=0.25), labels = c("-25%","0%","25%","50%")) +
  xlab("Contrast of posterior predicted probabilities of coin to DISTANT cup:\nNot free-listing ''Morality'' vs. free-listing only ''Morality''") +
  ylab("Density (normalized)") +
  annotate("text", x = 0.4, y = .8, 
           label = "Higher probability\nwhen free-listing only ''Morality''", 
           color = "#3182BD",
           alpha = 1,
           fontface = 2, size = 3.5) + 
  annotate("text", x = 0.4, y = .8, # duplicate to increase readability
           label = "Higher probability\nwhen free-listing only ''Morality''", 
           color = "#3182BD",
           alpha = 1,
           fontface = 2, size = 3.5) +
  annotate("text", x = -0.175, y = .8, 
           label = "Higher probability\nwhen not free-listing ''Morality''", 
           color = "gray40",
           fontface = 2, size = 3.5) +
  coord_cartesian(xlim=c(-0.3, 0.55)) # reduced margin at top (ylim) and center (xlim)

dev.off() # end print to pdf

