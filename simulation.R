## ------------------------------------------------------------
## Further information
## https://www.rdatagen.net/page/simstudy/
## https://cran.r-project.org/web/packages/simstudy/vignettes/simstudy.html

## ------------------------------------------------------------
pacman::p_load(tidyverse, simstudy, lumi, plotly, mosaic, colorspace, limma,
               minfi, xtable, plyr)
extract <- function(...) magrittr::extract(...)
## ------------------------------------------------------------

## ------------------------------------------------------------
## R Code for the generation of table 4  
eff <- 0.1
eff_tbl <- tibble(b0 = seq(0.001, 1.001 - eff, 0.1),
                  m0 = round(beta2m(b0), 2),
                  b1 = ifelse(b0 + eff <= 1, b0 + eff, 0.999),
                  m1 = round(beta2m(b1), 2),
                  db = round(b1 - b0, 1),
                  mb = m1 - m0,
                  formula = str_c(m0, " + ", mb, " * trt"))
write_csv(eff_tbl, "beta_mvalues.csv")
print(xtable(eff_tbl), include.rownames = FALSE)
## ------------------------------------------------------------

## ------------------------------------------------------------
## R code for case 1: one treatment with two levels

## number of simulations
n_sim <- 10 
## effect size in difference of Beta-values between groups
eff_sim <- c(0.005, 0.01, 0.05, 0.1, 0.2, 0.3) 
## number of overall samples
n_sample <- 100 

sim_case_1_tbl <- ldply(eff_sim, function(eff) {
   eff_tbl <- tibble(b0 = seq(0.001, 1.001 - eff, 0.1),
                    b1 = ifelse(b0 + eff <= 1, b0 + eff, 0.999),
                    m0 = round(beta2m(b0), 4),
                    m1 = round(beta2m(b1), 4),
                    mb = m1 - m0,
                    db = round(b1 - b0, 1),
                    formula = str_c(m0, " + ", mb, " * trt"))
  ##
  eff_vec_tbl <- ldply(1:n_sim, function(...) {
    sim_formula <- sample(eff_tbl$formula, 1)
    def <- defData(varname = "trt", dist = "binary", 
                   formula = 0.5 , id = "pat_id")
    def <- defData(def, varname = "m_values", dist = "normal",
                   formula = sim_formula, variance = 1)
    case_1_data <- genData(n_sample, def) %>%
      as_tibble %>%
      mutate(trt = as.factor(trt),
             b_values = m2beta(m_values))
    mean_beta <- lm(b_values ~ trt, data = case_1_data) %>%
      coef %>%
      extract(2)
    mean_beta_tbl <- tibble(mean_beta,
                            b0 = str_replace(sim_formula, "(\\d+)\\s\\+.*", "\\1")) %>%
      mutate(b0 = round(m2beta(as.numeric(b0)), 4))
    return(mean_beta_tbl)
  })
   eff_out <- eff_vec_tbl %>%
     as_tibble %>% 
     mutate(bias = mean_beta - eff,
            eff = eff,
            bias_perc = ((mean_beta - eff) / eff) * 100) 
   return(eff_out)
}, .progress = "text")  

pdf(file.path("case_1.pdf"), width = 6, height = 4.5)
sim_case_1_tbl %>%
  as_tibble %>%
  mutate(eff = as.factor(eff),
         b0 = round(b0, 1),
         b0 = as.factor(b0)) %>%
  filter(b0 %in% c(0.1, 0.3, 0.5, 0.7, 0.9)) %>% 
  ggplot(., aes(eff, bias_perc, fill = b0)) + geom_boxplot() +
  theme_bw() +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  labs(
  x = expression(paste("Predefined Beta-value differences from ", Trt['Placebo'], " to ",
                       Trt['Verum'],
                       " (", Delta['Beta'], ")")),
  y = expression(paste("Percentage deviation from predifined to estimated ", Delta['Beta'], "")),
  fill = expression(paste("", Trt['Placebo']))
  ) +
  scale_y_continuous(breaks = c(-400, -300, -200, -100, 0, 100, 200, 300, 400),
                     limits = c(-400, 400),
                     labels = c("-400%", "-300%", "-200%", "-100%", 0, "100%", "200%", "300%", "400%")) +
  theme(legend.position = c(0.925, 0.75)) 
dev.off()

## ------------------------------------------------------------
## R code for case 2: one treatment with two levels with two confounder

sim_case_2_tbl <- ldply(eff_sim, function(eff) {
  confounder_mat <- tibble(trt_eff = seq(from = 100, to = 80, by = -10),
                           age_eff = seq(from = 0, to = 10, length.out = length(trt_eff)),
                           gender_eff = seq(from = 0, to = 10, length.out = length(trt_eff))
                           ) %>%
  as.matrix(.)/100 
  ## run over the eff_matrix
  confounder_df <- ldply(1:nrow(confounder_mat), function(i) {
    eff_tbl <- tibble(b0 = seq(0.01, 1.01 - eff, 0.1),
                      b1 = ifelse(b0 + eff <= 1, b0 + eff, 0.999),
                      m0 = round(beta2m(b0), 4),
                      m1 = round(beta2m(b1), 4),
                      dm = m1 - m0,
                      db = round(b1 - b0, 1))
    ## run the simulation
    eff_vec_tbl <- ldply(1:n_sim, function(...) {
      ##
      eff_tbl_row <- eff_tbl %>%
        sample_n(1)
      eff_formula <- confounder_mat[i,] * eff_tbl_row$dm
      sim_formula <- str_c(eff_tbl_row$m0, "+",
                           eff_formula["trt_eff"], "* trt +",
                           eff_formula["age_eff"], "* age +",
                           eff_formula["gender_eff"], "* gender", sep = " ")
      ## 
      def <- defData(varname = "trt", dist = "binary", 
                     formula = 0.5 , id = "pat_id")
      def <- defData(def, varname = "age", dist = "normal", 
                     formula = 60 , variance = 5)
      def <- defData(def, varname = "gender", dist = "binary", 
                     formula = 0.5)
      def <- defData(def, varname = "m_values", dist = "normal",
                     formula = sim_formula,
                     variance = 1)
      case_2_data <- genData(n_sample, def) %>%
        as_tibble %>%
        mutate(trt = as.factor(trt),
               b_values = m2beta(m_values))
      mean_beta <- lm(b_values ~ trt, data = case_2_data) %>%
        coef %>%
        extract(2)
      ##
      mean_beta_tbl <- tibble(mean_beta,
                              b0 = str_replace(sim_formula, "(\\d+)\\s\\+.*", "\\1")) %>%
        mutate(b0 = round(m2beta(as.numeric(b0)), 4))
      mean_beta_tbl <- cbind(mean_beta_tbl, t(as.data.frame(confounder_mat[i,]))) %>%
        as_tibble
      return(mean_beta_tbl)
    })
    eff_out <- eff_vec_tbl %>%
      mutate(bias = mean_beta - eff * trt_eff,
             eff = eff,
             bias_perc = ((mean_beta - (eff * trt_eff)) / (eff * trt_eff)) * 100) 
    return(eff_out)
  })
}, .progress = "text")  


pdf("case_2.pdf", width = 8, height = 4.5)
sim_case_2_tbl %>%
  filter(trt_eff >= 0.8) %>%
  filter(eff >= 0.01 & eff <= 0.3) %>%
  as_tibble %>% 
  mutate(eff = as.factor(eff),
         b0 = round(b0, 1),
         b0 = as.factor(b0),
         age_eff = as.factor(age_eff)) %>%
  filter(b0 %in% c(0.1, 0.3, 0.5, 0.7, 0.9)) %>% 
  ggplot(., aes(eff, bias_perc, fill = b0)) +
  facet_wrap(~ ordered(trt_eff, levels = c(1, 0.9, 0.8),
                       labels = c("0% confounder effect",
                                  "10% confounder effect",
                                  "20% confounder effect"))) +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  geom_boxplot() +
  labs(
    x = expression(paste("Predefined Beta-value differences from ", Trt['Placebo'], " to ",
                         Trt['Verum'],
                         " (", Delta['Beta'], ")")),
    y = expression(paste("Percentage deviation from predifined to estimated ", Delta['Beta'], "")),
    fill = expression(paste("", Trt['Placebo']))
  ) +
  scale_y_continuous(breaks = c(-400, -300, -200, -100, 0, 100, 200, 300, 400),
                     limits = c(-400, 400),
                     labels = c("-400%", "-300%", "-200%", "-100%", 0, "100%", "200%", "300%", "400%")) +
  theme_bw() +
  theme(legend.position = c(0.925, 0.75)) 
dev.off()

## ------------------------------------------------------------
## by J.Kruppa on Monday, February 10, 2020 (19:01)
## R code for case 3: one treatment with two levels with two confounder using limma

sim_case_3_tbl <- ldply(eff_sim, function(eff) {
  confounder_mat <- tibble(trt_eff = seq(from = 100, to = 80, by = -10),
                           age_eff = seq(from = 0, to = 10, length.out = length(trt_eff)),
                           gender_eff = seq(from = 0, to = 10, length.out = length(trt_eff))
                           ) %>%
  as.matrix(.)/100 
  ## run over the eff_matrix
  confounder_df <- ldply(1:nrow(confounder_mat), function(i) {
    eff_tbl <- tibble(b0 = seq(0.01, 1.01 - eff, 0.1),
                      b1 = ifelse(b0 + eff <= 1, b0 + eff, 0.999),
                      m0 = round(beta2m(b0), 4),
                      m1 = round(beta2m(b1), 4),
                      dm = m1 - m0,
                      db = round(b1 - b0, 1))
    ## run the simulation
    eff_vec_tbl <- ldply(1:n_sim, function(...) {
      ##
      eff_tbl_row <- eff_tbl %>%
        sample_n(1)
      eff_formula <- confounder_mat[i,] * eff_tbl_row$dm
      sim_formula <- str_c(eff_tbl_row$m0, "+",
                           eff_formula["trt_eff"], "* trt +",
                           eff_formula["age_eff"], "* age +",
                           eff_formula["gender_eff"], "* gender", sep = " ")
      ## 
      def <- defData(varname = "trt", dist = "binary", 
                     formula = 0.5 , id = "pat_id")
      def <- defData(def, varname = "age", dist = "normal", 
                     formula = 60 , variance = 5)
      def <- defData(def, varname = "gender", dist = "binary", 
                     formula = 0.5)
      def <- defData(def, varname = "m_values", dist = "normal",
                     formula = sim_formula,
                     variance = 1)
      case_2_data <- genData(n_sample, def) %>%
        as_tibble %>%
        mutate(trt = as.factor(trt),
               b_values = m2beta(m_values))
      ## limma with intercept
      lim_fit <- lmFit(t(case_2_data$m_values),
                       with(case_2_data, model.matrix(~ trt + age + gender))) %>%
        eBayes
      lim_fit_coef <- eBayes(lim_fit)$coefficients
      ##
      m0_eff <- lim_fit_coef[, "(Intercept)"] 
      m1_eff <- m0_eff + lim_fit_coef[, "trt1"] 
      mean_beta <- abs(m2beta(m0_eff) - m2beta(m1_eff))
      ## go on here
      mean_beta_tbl <- tibble(mean_beta,
                              b0 = str_replace(sim_formula, "(\\d+)\\s\\+.*", "\\1")) %>%
        mutate(b0 = round(m2beta(as.numeric(b0)), 4))
      mean_beta_tbl <- cbind(mean_beta_tbl, t(as.data.frame(confounder_mat[i,]))) %>%
        as_tibble
      return(mean_beta_tbl)
    })
    eff_out <- eff_vec_tbl %>%
      mutate(bias = mean_beta - (eff * trt_eff),
             eff = eff,
             bias_perc = ((mean_beta - (eff * trt_eff)) / (eff * trt_eff)) * 100) 
    return(eff_out)
  })
}, .progress = "text")  

pdf("case_3.pdf", width = 8, height = 4.5)
sim_case_3_tbl %>%
  filter(trt_eff >= 0.8) %>%
  filter(eff >= 0.01 & eff <= 0.3) %>%
  as_tibble %>% 
  mutate(eff = as.factor(eff),
         b0 = round(b0, 1),
         b0 = as.factor(b0),
         age_eff = as.factor(age_eff)) %>%
  filter(b0 %in% c(0.1, 0.3, 0.5, 0.7, 0.9)) %>% 
  ggplot(., aes(eff, bias_perc, fill = b0)) +
  facet_wrap(~ ordered(trt_eff, levels = c(1, 0.9, 0.8),
                       labels = c("0% confounder effect",
                                  "10% confounder effect",
                                  "20% confounder effect"))) +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  geom_boxplot() +
  labs(
    x = expression(paste("Predefined Beta-value differences from ", Trt['Placebo'], " to ",
                         Trt['Verum'],
                         " (", Delta['Beta'], ")")),
    y = expression(paste("Percentage deviation from predifined to estimated ", Delta['Beta'], "")),
    fill = expression(paste("", Trt['Placebo']))
  ) +
  scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400),
                     limits = c(-100, 400),
                     labels = c("-100%", 0, "100%", "200%", "300%", "400%")) +
  theme_bw() +
  theme(legend.position = c(0.925, 0.75)) 
dev.off()


## ------------------------------------------------------------
## by J.Kruppa on Tuesday, February 11, 2020 (08:03)
## end
