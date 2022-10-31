library(tidyr)
library(dplyr)
library(ggplot2)
library(tidybayes)

rm(list=ls())

memory.limit(size=56000)

set.seed(1992)

load("ord_cat_sim_testdat.RData")

fit <- m_cat_int_fit

ndraws <- 1:nrow(fit) # all draws from the posterior

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
C <- rstan::extract(fit,"cutpoints")$cutpoints[ndraws,]

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

levels <- length(unique(M_pred)) # # how many levels of the focal predictor?

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
pred_coef_group <- data.frame(bM=rep(bM, times), bO=rep(bO, times), bP=rep(bP, times),
                              bMj=rep(as.vector(bMj), timesj),bOj=rep(as.vector(bOj), timesj),bPj=rep(as.vector(bPj), timesj),
                              bMOj=rep(as.vector(bMOj), timesj),bMPj=rep(as.vector(bMPj), timesj), bPOj=rep(bPOj, times), bMPOj=rep(bMPOj, times),
                              bMO=rep(bMO, times), bMP=rep(bMP, times), bPO=rep(bPO, times), bMPO=rep(bMPO, times),
                              bC=rep(bC, times), bS=rep(bS, times), 
                              border=rep(border, times), btype=rep(btype, times), bcheck=rep(bcheck,times),
                              .draw=rep(1:length(ndraws), times))

predeta_fun <- function(x) {
  eta <- with(x, (bM+bMj)*M + 
                (bO+bOj)*O + (bP + bPj)*P +
                (bMP+bMPj)*M*P + (bMO+bMOj)*M*O + (bPO+bPOj)*P*O +
                (bMPO+bMPOj)*M*P*O +
                bC*C + bS*S +
                border*order + btype*type + bcheck*check)
}

# additive model
# pred_coef_group <- data.frame(bM=rep(bM, times), bO=rep(bO, times), bP=rep(bP, times),
#                              bMj=rep(as.vector(bMj), timesj),bOj=rep(as.vector(bOj), timesj),bPj=rep(as.vector(bPj), timesj),
#                              bC=rep(bC, times), bS=rep(bS, times),
#                              border=rep(border, times), btype=rep(btype, times), bcheck=rep(bcheck,times),
#                              .draw=rep(1:length(ndraws), times))
# predeta_fun <- function(x) {
#  eta <- with(x, (bM+bMj)*M + 
#                          (bO+bOj)*O + (bP + bPj)*P +
#                         bC*C + bS*S +
#                          border*order + btype*type + bcheck*check)
# }

pred_mat_group <- cbind(pred_ate_group, pred_coef_group, row.names = NULL)

predcoin_fun <- function(eta, C){
  coin0 <- 1 - plogis((eta)-C[,1])
  coin1 <- plogis((eta)-C[,1]) - plogis((eta)-C[,2])
  coin2 <- plogis((eta)-C[,2]) - plogis((eta)-C[,3])
  coin3 <- plogis((eta)-C[,3]) - plogis((eta)-C[,4])
  coin4 <- plogis((eta)-C[,4]) - plogis((eta)-C[,5])
  coin5 <- plogis((eta)-C[,5]) - plogis((eta)-C[,6])
  coin6 <- plogis((eta)-C[,6]) - plogis((eta)-C[,7])
  coin7 <- plogis((eta)-C[,7]) - plogis((eta)-C[,8])
  coin8 <- plogis((eta)-C[,8]) - plogis((eta)-C[,9])
  coin9 <- plogis((eta)-C[,9]) - plogis((eta)-C[,10])
  coin10 <- plogis((eta)-C[,10])
  return(as.data.frame(cbind(coin0,coin1,coin2,coin3,coin4,coin5,coin6,coin7,coin8,coin9,coin10)))
}

# calculate eta
for (group in unique(pred_mat_group$group)){  # for each category that is unique in the culture variable
  sub <- subset(pred_mat_group, group == group) |> as.data.frame() # subset
  sub$eta <- predeta_fun(sub) 
  subeta <- sub
  return(subeta) # return data frame
}

# calculate probability of each coin allocation
for (group in unique(subeta$group)){  # for each category that is unique in the culture variable
  prepcoin <- subset(subeta, group == group)
  subcoin <- predcoin_fun(eta=prepcoin$eta, C=C) # compute probability of coin allocations
  return(postpredcoin <- cbind(prepcoin, subcoin)) # return data frame
}

# get predictions in format for plotting
postpredcoin_select <- postpredcoin %>% select(c(M, id, group, .draw, starts_with("coin")))

postpredcoin_select_summary_l <- postpredcoin_select %>% pivot_longer(5:15, names_to = "coins")

M0coins <- postpredcoin_select_summary_l[postpredcoin_select_summary_l$M==0,]
M0coins$coins <- sub("^coin", "", M0coins$coins) |> as.numeric()
colnames(M0coins) <- paste0("M0", colnames(M0coins))
M1coins <- postpredcoin_select_summary_l[postpredcoin_select_summary_l$M==1,]
M1coins$coins <- sub("^coin", "", M1coins$coins) |> as.numeric()
colnames(M1coins) <- paste0("M1", colnames(M1coins))

ridgeplotcoin <- cbind(M0coins, M1coins)
ridgeplotcoin <- ridgeplotcoin |> group_by(M0id, M0.draw, M0coins) |> mutate(diff = M1value-M0value) |> ungroup()

cairo_pdf("catridge.pdf",
          width = 7, height = 7) # start print to pdf

# ridge plot (normalized): contrast of posterior predicted probabilities in percentage points difference
ggplot(ridgeplotcoin, aes(y = M0coins, x = diff)) +
  
  stat_halfeye(aes(fill = stat(x > 0)), 
               normalize="all",
               orientation = "horizontal",
               .width = 0.95,
               size = 1,
               point_interval = "median_hdi", # summarize the distribution with its median and highest posterior density
               show.legend=FALSE) +

  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("gray40", "#3182BD")) +
  theme(panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major.y = element_line( size=.05, color="grey" ),
        legend.position = "none") +
  scale_y_continuous(breaks = seq(0, 10, by = 1)) + 
  scale_x_continuous(breaks = seq(-0.15,0.15, by=0.05), labels = c("-15%","-10%","-5%","0%","5%","10%","15%")) +
  xlab("Contrast of posterior predicted probabilities:\nNot free-listing ''Morality'' vs. free-listing only ''Morality''") +
  ylab("Number of coins allocated to the DISTANT cup") +
  annotate("text", x = 0.1, y = 7.5, 
           label = "Higher probability\nwhen free-listing only ''Morality''", 
           color = "#3182BD",
           alpha = 1,
           fontface = 2) + 
  annotate("text", x = 0.1, y = 7.5, # duplicate to increase readability
           label = "Higher probability\nwhen free-listing only ''Morality''", 
           color = "#3182BD",
           alpha = 1,
           fontface = 2) +
  annotate("text", x = -0.1, y = 7.5, 
           label = "Higher probability\nwhen not free-listing ''Morality''", 
           color = "gray40",
           fontface = 2) +
  coord_cartesian(ylim=c(0,10), xlim=c(-0.15, 0.15)) # reduced margin at top (ylim) and center (xlim)

dev.off() # end print to pdf

# "line-ribbon" plot 
postpredcoin_select_summary_l
postpredcoin_select_summary_l$coins <- factor(postpredcoin_select_summary_l$coins,
                               levels = c("coin0", "coin1", "coin2", "coin3", "coin4",
                                          "coin5", "coin6", "coin7", "coin8", "coin9", "coin10"),
                               ordered=TRUE)

cairo_pdf("catlinerib.pdf",
          width = 7, height = 7) # start print to pdf

result_plot <- 
  ggplot(postpredcoin_select_summary_l, 
         aes(x = M, y = value, color = coins, fill = coins)) +
  tidybayes::stat_lineribbon(point_interval = tidybayes::median_qi, .width = 0.1, alpha = 0.5) +
  theme_classic() + 
  ylab("Probability of coin allocation to DISTANT cup") + xlab("Proportion of ''Morality'' items in free-lists")

result_plot

dev.off() # end print to pdf

### END
