########################################################################
### "The role of ``morality'' in moralistic supernatural punishment"
###
### Pre-registration: Posterior predictive checks, coef plots and "effect plots"
###
### Script by Theiss Bendixen & Benjamin G. Purzycki

library(tidyverse)
library(ggplot2)
library(tidybayes)

memory.limit(size=56000)

### unblock to load simulated data and fitted models
base::load("ord_cat_sim_testdat.RData")
base::load("binom_sim_testdat.RData")

### Panel: Posterior predictive checks and coefficient plots
cairo_pdf("ppc_coef.pdf",
          width = 7, height = 11) # start print to pdf

(int_ppc_rag + ggtitle("(A)\nPosterior predictive check RAG")) + 
  (int_ppc_dg + ylab(NULL) + ggtitle("(B)\nPosterior predictive check DG") + 
     scale_x_continuous(breaks = seq(1,11, by=1), 
                      labels = c("0","1","2","3","4","5","6","7","8","9","10"))) +
  (int_coef_rag + xlab(NULL) + ggtitle("(C)\nCoefficients RAG interaction model") + xlim(-1.4,1.4)) + 
  (int_coef_dg + ylab(NULL) + xlab(NULL) + ggtitle("(D)\nCoefficients DG interaction model") + xlim(-1.4,1.4)) +
  (noint_coef_rag + ggtitle("(E)\nCoefficients RAG additive model") + xlim(-1.4,1.4)) + 
  (noint_coef_dg + ylab(NULL) + ggtitle("(F)\nCoefficients DG additive model") + xlim(-1.4,1.4)) +
  patchwork::plot_layout(ncol = 2)

dev.off() # end print to pdf
