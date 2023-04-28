library(ggtext)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidybayes)

### Load prepped data
full_df <- read.csv("data_prepped.csv")

### Free-listed Morality

cairo_pdf("m_raw.pdf",
          width = 6, height = 3.5) # start print to pdf

# exclude non-freelist sites
full_df[!full_df$GROUP %in% c("China-Buddhists", "China-Christians"),] |>
  
  ggplot(aes(BGD.MORPROP)) +
  stat_halfeye(aes(alpha = 0.5,
                   fill = "#3182BD"),
               n = 1001,
               orientation = "horizontal",
               point_interval = NULL, # remove interval
               show.legend=FALSE) + 
  scale_fill_manual(values = c("#3182BD")) +
  facet_wrap(~GROUP, nrow=3) +
  theme(panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 6),
        strip.text.x = element_text(size = 6.5),
        plot.title = element_text(size = 10),
        plot.subtitle = element_text(size = 8)) +
  xlab(NULL) + 
  ylab(NULL) + 
  scale_x_continuous(breaks = seq(0,1, by=0.5), labels = c(0,0.5,1)) +
  theme(plot.title = element_markdown(), plot.subtitle = element_markdown()) +
  labs(title = "Supernatural moral beliefs across cultures")

dev.off() # end print to pdf

### RAG
y_notmiss <- which(full_df$RAG==1) # Only RAG SELF and non-missing outcome

rag_dat_long <- with(full_df,
                          data.frame(
                            id = NEW.ID[y_notmiss],
                            group = GROUP[y_notmiss],
                            y = Y[y_notmiss],
                            SELF = SELF[y_notmiss]
                          )
)

rag_dat <- tidyr::pivot_wider(rag_dat_long, id_cols = c("id", "group"), names_from = "SELF", values_from = "y") |>
  rename("LOCAL" = "0", "SELF" = "1")

rag_raw <- rag_dat[!rag_dat$group %in% c("China-Buddhists", "China-Christians"),] |>
  ggplot() + 
  geom_line(aes(LOCAL), stat="density", color = "#d35400cc", alpha = 0.8, size = 1) +
  geom_line(aes(SELF), stat="density", color = "#0072B2", alpha = 0.8, size = 1) +
  geom_vline(xintercept = 15, linetype = "dashed", alpha = 0.5) +
  facet_wrap(~ group, nrow = 3) +
  theme(plot.title = element_text(size = 10), 
        plot.subtitle = element_text(size = 8),
        legend.position = "none",
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="white"),
        axis.text.x = element_text(size = 6),
        strip.text.x = element_text(size = 6.5)) +
  theme(plot.title = element_markdown(), plot.subtitle = element_markdown()) +
  xlab(NULL) +
  ylab(NULL) +
  labs(title = "Random Allocation Games",
       subtitle = "Raw data distributions of coins to DISTANT cup comparing <b style='color:#0072B2'>SELF</b> and <b style='color:#d35400;'>LOCAL</b> games")

### DG
y_notmiss <- which(full_df$RAG==0) # Only RAG SELF and non-missing outcome

dg_dat_long <- with(full_df,
                     data.frame(
                       id = NEW.ID[y_notmiss],
                       group = GROUP[y_notmiss],
                       y = Y[y_notmiss],
                       SELF = SELF[y_notmiss]
                     )
)

dg_dat <- tidyr::pivot_wider(dg_dat_long, id_cols = c("id", "group"), names_from = "SELF", values_from = "y") |>
  rename("LOCAL" = "0", "SELF" = "1")

dg_raw <- 
  # exclude wave 1 and non-freelist sites
  dg_dat[!rag_dat$group %in% c("China-Buddhists", "China-Christians", "Coastal Tanna", "Lovu Fiji", "Marajo", "Samburu", "Tyva Republic", "Yasawa Fiji"),] |> 
  ggplot() + 
  geom_line(aes(LOCAL), stat="density", color = "#d35400cc", alpha = 0.8, size = 1) +
  geom_line(aes(SELF), stat="density", color = "#0072B2", alpha = 0.8, size = 1) +
  geom_vline(xintercept = 5, linetype = "dashed", alpha = 0.5) +
  facet_wrap(~ group, ncol = 3) +
  theme(plot.title = element_text(size = 10), 
        plot.subtitle = element_text(size = 8),
        legend.position = "none",
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="white"),
        axis.text.x = element_text(size = 6),
        strip.text.x = element_text(size = 6.5)) +
  theme(plot.title = element_markdown(), plot.subtitle = element_markdown()) +
  xlab(NULL) +
  ylab(NULL) +
  scale_x_continuous(breaks = c(0,5,10)) +
  labs(title = "Dictator Games",
       subtitle = "Raw data distributions of coins to DISTANT cup comparing <b style='color:#0072B2'>SELF</b> and <b style='color:#d35400;'>LOCAL</b> games")

### Panel of raw game plots
cairo_pdf("game_raw.pdf",
          width = 6, height = 7) # start print to pdf

rag_raw + 
  dg_raw + 
  patchwork::plot_layout(ncol = 1, nrow = 2)

dev.off() # end print to pdf
