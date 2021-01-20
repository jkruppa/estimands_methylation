## ------------------------------------------------------------
pacman::p_load(simstudy, lumi, betareg, tidyverse, plyr, faraway)
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## ------------------------------------------------------------
## https://www.rdatagen.net/page/simstudy/

## ------------------------------------------------------------
## by J.Kruppa on Tuesday, January 21, 2020 (08:29)
## beta regression

b0_prob_lower_vec <- c(0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1,
                       0.15, 0.2, 0.3, 0.4, 0.5)
b0_logit_vec <- logit(b0_prob_lower_vec) 
b1_prob <- 0.1
b1_logit_vec <- logit(b0_prob_lower_vec + b1_prob) - logit(b0_prob_lower_vec)

## chek if we have everything right
ilogit(b0_logit_vec + b1_logit_vec)
(b1_logit_vec - b0_logit_vec)[7]

## start simulation
n_sim <- 100
n_sample <- 10
sim_beta_logit_df <- ldply(seq_along(b0_logit_vec), function(i) {
  options(show.error.messages = FALSE)
  b0_logit <- b0_logit_vec[i]
  b1_logit <- b1_logit_vec[i]
  b1_obs_vec <- laply(1:n_sim, function(...) {
    def <- defData(varname = "trt", dist = "binary", 
                   formula = 0.5 , id = "pat_id")
    def <- defData(def, varname = "b_values", dist = "beta",
                   formula = str_c(b0_logit, " + ", b1_logit, " * trt"), variance = 1, link = "logit")
    beta_data <- genData(n_sample, def) %>%
    as_tibble %>%
    mutate(trt = as.factor(trt),
           b_values = b_values)
    beta_fit <- suppressWarnings(try(betareg(b_values ~ trt, beta_data), silent = TRUE))
    if(class(beta_fit) == "try-error") {
      b1_obs <- NA
    } else {
      b1_obs <- beta_fit %>%
        coef %>%
        pluck("trt1")
    }
    return(b1_obs)
  })
  eff_out <- tibble(b1_obs_vec,
                    b0 = ilogit(b0_logit),
                    b1 = b1_logit,
                    bias = b1 - b1_obs_vec,
                    bias_perc = ((b1 - b1_obs_vec) / b1))
  options(show.error.messages = TRUE)
  return(eff_out)
}, .progress = "text") 

## if this runs ours, it is better to save it
## write_rds(sim_beta_logit_df, "sim_beta_logit_df.rds")

sim_beta_logit_tbl <- sim_beta_logit_df %>%
  as_tibble %>% 
  mutate(b0 = as.factor(round(b0, 4)))

conv_tbl <- sim_beta_logit_tbl %>%
  ddply("b0", summarise, conv_rate = 1-sum(is.na(bias))/n_sim)

pdf(file.path("sim_beta_reg.pdf"), width = 6, height = 4.5)
sim_beta_logit_tbl %>% 
  ggplot(., aes(b0, bias, fill = b0)) +
  geom_boxplot() +
  geom_line(data = conv_tbl, aes(b0, conv_rate, group = 1), size = 1,
            color = cbbPalette[1]) +
  geom_point(data = conv_tbl, aes(b0, conv_rate, group = 1), size = 2,
             color = cbbPalette[1], shape = 4) +
  theme_bw()  +
  labs(x = expression(paste(beta[0])),
       y = "Covergence rate") +
  scale_y_continuous(sec.axis = sec_axis(~., name = expression(paste("Bias (", hat(beta)[1], " - ", beta[1], ")"))),
                     limits = c(-.5, 1)) +
  theme(legend.position = "none") 
dev.off()

