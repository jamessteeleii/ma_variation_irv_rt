
library(tidyverse)
library(metafor)
library(patchwork)
library(bayestestR)


# variance comparison plots

set.seed(1988)

b_0 <- 0 
b_1 <- 0.5
between <- 0.9
error <- 0.1

sim_data <- tibble(
  id = seq(1:20),
  sigma = rnorm(20, 0, between),
  true_outcome_con = b_0 + sigma,
  true_outcome_int = b_0 + b_1 + sigma,
  obs_outcome_con = true_outcome_con + rnorm(20, 0, error),
  obs_outcome_int = true_outcome_int + rnorm(20, 0, error)
)

var.test(sim_data$true_outcome_con,
         sim_data$true_outcome_int)

var.test(sim_data$obs_outcome_con,
         sim_data$obs_outcome_int)

escalc(measure = "CVR",
       sd1i = sd(sim_data$obs_outcome_con),
       sd2i = sd(sim_data$obs_outcome_int),
       m1i = mean(sim_data$obs_outcome_con),
       m2i = mean(sim_data$obs_outcome_int),
       n1i = 20,
       n2i=20
       )

sd(sim_data$obs_outcome_con)/sd(sim_data$obs_outcome_int)

if((var(sim_data$obs_outcome_con)-var(sim_data$obs_outcome_int))> 0){
  sd_ir = sqrt(var(sim_data$obs_outcome_con)-var(sim_data$obs_outcome_int))
  
} else{
  sd_ir = -1*sqrt(var(sim_data$obs_outcome_int)-var(sim_data$obs_outcome_con))
  
}

(sd(sim_data$obs_outcome_con)^2 + sd(sim_data$obs_outcome_int)^2 - unname(sd_ir)^2) / (2 * sd(sim_data$obs_outcome_int) * sd(sim_data$obs_outcome_con))


sim_data |>
  pivot_longer(3:6, 
               names_to = "outcome",
               values_to = "value") |>
  separate(outcome, into = c("outcome", "x", "group")) |>
  select(-x) |>
  ggplot(aes(x = id, y = value, group = outcome, colour = group)) +
  geom_point(position = position_dodge(width=0.5)) +
  geom_line(aes(linetype = outcome, group=interaction(outcome,id)), color = "black", 
            position = position_dodge(width=0.5)) +
  coord_flip()


dodge_y <- position_dodge(width = 1)

sim_data |>
  pivot_longer(3:6, 
               names_to = "outcome",
               values_to = "value") |>
  separate(outcome, into = c("outcome", "x", "group")) |>
  select(-x) |>
  ggplot(aes(y = id, x = value, colour = group)) +
  geom_point(position = dodge_y) +
  geom_line(
    aes(group = interaction(id, outcome), linetype = outcome),
    position = dodge_y,
    color = "black"
  )


### Meta


data <- read_csv(url("https://github.com/jamessteeleii/Meta-Analysis-of-Variation-in-Resistance-Training/raw/refs/heads/main/data/Polito%20et%20al.%20RT%20Extracted%20Data.csv"))


data_long_m <- data %>%
  select(study, arm, es, outcome, RT_n, CON_n, RT_pre_m, CON_pre_m, RT_post_m, CON_post_m) %>%
  pivot_longer(c(RT_post_m, CON_post_m),
               names_to = "group",
               values_to = "mean")

data_long_m$group <- recode(data_long_m$group, 
                            RT_pre_m = "RT", CON_pre_m = "CON",
                            RT_post_m = "RT", CON_post_m = "CON")

data_long_sd <- data %>%
  select(study, arm, es, outcome, RT_n, CON_n, RT_pre_sd, CON_pre_sd, RT_post_sd, CON_post_sd) %>%
  pivot_longer(c(RT_post_sd, CON_post_sd),
               names_to = "group",
               values_to = "sd")

data_long_sd$group <- recode(data_long_sd$group,  
                                  RT_pre_sd = "RT", CON_pre_sd = "CON",
                                  RT_post_sd = "RT", CON_post_sd = "CON")

data_long <- cbind(data_long_m, sd = data_long_sd$sd) 

data_long_RT <- data_long %>%
  filter(group == "RT") %>%
  mutate(n = RT_n,
         arm = as.factor(unclass(factor(unlist(arm)))),
         es = as.factor(unclass(factor(unlist(es))))) %>%
  select(study, arm, es, outcome, group, group, mean, sd, n)

data_long_CON <- data_long %>%
  filter(group == "CON") %>%
  distinct(study, group, group, mean, sd, .keep_all = TRUE) %>%
  mutate(n = CON_n,
         arm = as.factor(unclass(factor(unlist(arm)))+length(unique(data_long_post_RT$arm))),
         es = as.factor(unclass(factor(unlist(es)))+length(unique(data_long_post_RT$es)))) %>%
  select(study, arm, es, outcome, group, group, mean, sd, n)

data_long <- rbind(data_long_RT, data_long_CON) %>%
  filter(!is.na(mean) |
           !is.na(sd))

# Calculate log SD and variance of log SD
data_long <- escalc(measure = "SDLN",
                         sdi = sd,
                         ni = n,
                         data = data_long)


# Plot raw mean and SD
m_sd_strength <- ggplot(subset(data_long, outcome == "strength"), aes(x=mean, y=sd)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.1) +
  geom_vline(aes(xintercept = 0), alpha = 0.1, lty = "dashed") +
  geom_hline(aes(yintercept = 0), alpha = 0.1, lty = "dashed") +
  geom_point(aes(x=mean, y=sd), alpha = 0.2) +
  labs(x = "Mean", y = "Standard Deviation") +
  guides(fill = "none") +
  theme_classic() +
  ggtitle("Strength outcomes")

m_sd_hypertrophy <- ggplot(subset(data_long, outcome == "hypertrophy"), aes(x=mean, y=sd)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.1) +
  geom_vline(aes(xintercept = 0), alpha = 0.1, lty = "dashed") +
  geom_hline(aes(yintercept = 0), alpha = 0.1, lty = "dashed") +
  geom_point(aes(x=mean, y=sd), alpha = 0.2) +
  labs(x = "Mean", y = "Standard Deviation") +
  guides(fill = "none") +
  theme_classic() +
  ggtitle("Hypertrophy outcomes")

m_sd_strength_log <- ggplot(subset(data_long, outcome == "strength"), aes(x=log(mean), y=yi)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.1) +
  geom_vline(aes(xintercept = 1), alpha = 0.1, lty = "dashed") +
  geom_hline(aes(yintercept = 1), alpha = 0.1, lty = "dashed") +
  geom_point(aes(x=log(mean), y=yi), alpha = 0.2) +
  labs(x = "Log Mean", y = "Log Standard Deviation") +
  guides(fill = "none") +
  theme_classic()  +
  ggtitle("Strength outcomes")

m_sd_hypertrophy_log <- ggplot(subset(data_long, outcome == "hypertrophy"), aes(x=log(mean), y=yi)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.1) +
  geom_vline(aes(xintercept = 1), alpha = 0.1, lty = "dashed") +
  geom_hline(aes(yintercept = 1), alpha = 0.1, lty = "dashed") +
  geom_point(aes(x=log(mean), y=yi), alpha = 0.2) +
  labs(x = "Log Mean", y = "Log Standard Deviation") +
  guides(fill = "none") +
  theme_classic() +
  ggtitle("Hypertrophy outcomes")

# Plot together
mean_variance_raw_plots <- (m_sd_strength | m_sd_hypertrophy) +  plot_layout(axes = "collect")

mean_variance_log_plots <- (m_sd_strength_log | m_sd_hypertrophy_log) +  plot_layout(axes = "collect")


ggsave("mean_variance_raw_plots.tiff", plot = mean_variance_raw_plots, 
       device = "tiff", dpi = 300,
       w = 7.5, h = 3.75)

ggsave("mean_variance_log_plots.tiff", plot = mean_variance_log_plots, 
       device = "tiff", dpi = 300,
       w = 7.5, h = 3.75)



#### Post only modelling ----

# lnVR and lnCVR
data <- escalc(
  measure = "CVR",
  m1i = RT_post_m,
  m2i = CON_post_m,
  sd1i = RT_post_sd,
  sd2i = CON_post_sd,
  n1i = RT_n,
  n2i = CON_n,
  data = data,
  var.names = c("lnCVR_yi", "lnCVR_vi")
)

data <- escalc(
  measure = "VR",
  m1i = RT_post_m,
  m2i = CON_post_m,
  sd1i = RT_post_sd,
  sd2i = CON_post_sd,
  n1i = RT_n,
  n2i = CON_n,
  data = data,
  var.names = c("lnVR_yi", "lnVR_vi")
)

model_strength_VR <- rma.mv(lnVR_yi, V=lnVR_vi, 
                                        data= data %>% 
                                          filter(outcome == "strength"),
                                        slab=paste(label),
                                        random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), method="REML", test="t",
                                        control=list(optimizer="optim", optmethod="Nelder-Mead"))



model_strength_VR <- robust(model_strength_VR, filter(data, outcome == "strength")$study)

model_hypertrophy_VR <- rma.mv(lnVR_yi, V=lnVR_vi, 
                            data= data %>% 
                              filter(outcome == "hypertrophy"),
                            slab=paste(label),
                            random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), method="REML", test="t",
                            control=list(optimizer="optim", optmethod="Nelder-Mead"))



model_hypertrophy_VR <- robust(model_hypertrophy_VR, filter(data, outcome == "hypertrophy")$study)

model_strength_CVR <- rma.mv(lnCVR_yi, V=lnCVR_vi, 
                            data= data %>% 
                              filter(outcome == "strength"),
                            slab=paste(label),
                            random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), method="REML", test="t",
                            control=list(optimizer="optim", optmethod="Nelder-Mead"))



model_strength_CVR <- robust(model_strength_CVR, filter(data, outcome == "strength")$study)

model_hypertrophy_CVR <- rma.mv(lnCVR_yi, V=lnCVR_vi, 
                               data= data %>% 
                                 filter(outcome == "hypertrophy"),
                               slab=paste(label),
                               random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), method="REML", test="t",
                               control=list(optimizer="optim", optmethod="Nelder-Mead"))


model_hypertrophy_CVR <- robust(model_hypertrophy_CVR, filter(data, outcome == "hypertrophy")$study)



model_estimates <- tibble(
  estimate = c(model_strength_VR$b[1], model_hypertrophy_VR$b[1], model_strength_CVR$b[1], model_hypertrophy_CVR$b[1]),
  ci.lb = c(model_strength_VR$ci.lb, model_hypertrophy_VR$ci.lb, model_strength_CVR$ci.lb, model_hypertrophy_CVR$ci.lb),
  ci.ub = c(model_strength_VR$ci.ub, model_hypertrophy_VR$ci.ub, model_strength_CVR$ci.ub, model_hypertrophy_CVR$ci.ub),
  outcome = c("Strength", "Hypertrophy","Strength", "Hypertrophy"),
  effect_size = c("lnVR","lnVR","lnCVR","lnCVR")
) |>
  mutate(effect_size = factor(effect_size, levels = c("lnVR", "lnCVR")))

model_estimates_plot_strength <- ggplot(model_estimates |> filter(outcome=="Strength"), aes(x=effect_size, y=estimate)) +
  geom_hline(aes(yintercept = 0), alpha = 0.1, lty = "dashed") +
  geom_point() +
  geom_linerange(aes(ymin=ci.lb, ymax=ci.ub)) +
  geom_text(
    aes(
      label = glue::glue("{round(estimate,2)}\n[95%CI: {round(ci.lb,2)}, {round(ci.ub,2)}]")
    ),
    position = position_nudge(x=c(-0.3,0.3)),
    size = 2
  ) +
  labs(x = "Effect Size Type", y = "Estimate\n(Positive Values Favour Resistance Training)") +
  theme_classic() +
  ggtitle("Strength Outcomes")

model_estimates_plot_hypertrophy <- ggplot(model_estimates |> filter(outcome=="Hypertrophy"), aes(x=effect_size, y=estimate)) +
  geom_hline(aes(yintercept = 0), alpha = 0.1, lty = "dashed") +
  geom_point() +
  geom_linerange(aes(ymin=ci.lb, ymax=ci.ub)) +
  geom_text(
    aes(
      label = glue::glue("{round(estimate,2)}\n[95%CI: {round(ci.lb,2)}, {round(ci.ub,2)}]")
    ),
    position = position_nudge(x=c(-0.3,0.3)),
    size = 2
  ) +
  labs(x = "Effect Size Type", y = "Estimate\n(Positive Values Favour Resistance Training)") +
  theme_classic() +
  ggtitle("Hypertrophy Outcomes")


model_estimates_plots <- (model_estimates_plot_strength | model_estimates_plot_hypertrophy) +
  plot_layout(guides = 'collect',
              axes = "collect")


model_estimates_plots

ggsave("model_estimates_plots.tiff", plot = model_estimates_plots, 
       device = "tiff", dpi = 300,
       w = 7.5, h = 3.75)


# mean-variance regression


data_long_post_m <- data %>%
  select(study, arm, es, outcome, RT_n, CON_n, RT_post_m, CON_post_m) %>%
  pivot_longer(c(RT_post_m, CON_post_m),
               names_to = "group",
               values_to = "mean")

data_long_post_m$group <- recode(data_long_post_m$group, RT_post_m = "RT", CON_post_m = "CON")

data_long_post_sd <- data %>%
  select(study, arm, es, outcome, RT_n, CON_n, RT_post_sd, CON_post_sd) %>%
  pivot_longer(c(RT_post_sd, CON_post_sd),
               names_to = "group",
               values_to = "sd")

data_long_post_sd$group <- recode(data_long_post_sd$group, RT_post_sd = "RT", CON_post_sd = "CON")

data_long_post <- cbind(data_long_post_m, sd = data_long_post_sd$sd) 

data_long_post_RT <- data_long_post %>%
  filter(group == "RT") %>%
  mutate(n = RT_n,
         arm = as.factor(unclass(factor(unlist(arm)))),
         es = as.factor(unclass(factor(unlist(es))))) %>%
  select(study, arm, es, outcome, group, group, mean, sd, n)

data_long_post_CON <- data_long_post %>%
  filter(group == "CON") %>%
  distinct(study, group, group, mean, sd, .keep_all = TRUE) %>%
  mutate(n = CON_n,
         arm = as.factor(unclass(factor(unlist(arm)))+length(unique(data_long_post_RT$arm))),
         es = as.factor(unclass(factor(unlist(es)))+length(unique(data_long_post_RT$es)))) %>%
  select(study, arm, es, outcome, group, group, mean, sd, n)

data_long_post <- rbind(data_long_post_RT, data_long_post_CON) %>%
  filter(!is.na(mean) |
           !is.na(sd))

# Calculate log SD and variance of log SD
data_long_post <- escalc(measure = "SDLN",
                         sdi = sd,
                         ni = n,
                         data = data_long_post)


# Plot raw mean and SD
m_sd_strength <- ggplot(subset(data_long_post, outcome == "strength"), aes(x=mean, y=sd)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.1) +
  geom_vline(aes(xintercept = 0), alpha = 0.1, lty = "dashed") +
  geom_hline(aes(yintercept = 0), alpha = 0.1, lty = "dashed") +
  geom_point(aes(x=mean, y=sd), alpha = 0.2) +
  labs(x = "Mean", y = "Standard Deviation") +
  guides(fill = "none") +
  theme_classic() +
  ggtitle("Strength groups")

m_sd_hypertrophy <- ggplot(subset(data_long_post, outcome == "hypertrophy"), aes(x=mean, y=sd)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.1) +
  geom_vline(aes(xintercept = 0), alpha = 0.1, lty = "dashed") +
  geom_hline(aes(yintercept = 0), alpha = 0.1, lty = "dashed") +
  geom_point(aes(x=mean, y=sd), alpha = 0.2) +
  labs(x = "Mean", y = "Standard Deviation") +
  guides(fill = "none") +
  theme_classic() +
  ggtitle("Hypertrophy groups")

m_sd_strength_log <- ggplot(subset(data_long_post, outcome == "strength"), aes(x=log(mean), y=yi)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.1) +
  geom_vline(aes(xintercept = 1), alpha = 0.1, lty = "dashed") +
  geom_hline(aes(yintercept = 1), alpha = 0.1, lty = "dashed") +
  geom_point(aes(x=log(mean), y=yi), alpha = 0.2) +
  labs(x = "Log Mean", y = "Log Standard Deviation") +
  guides(fill = "none") +
  theme_classic() 

m_sd_hypertrophy_log <- ggplot(subset(data_long_post, outcome == "hypertrophy"), aes(x=log(mean), y=yi)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.1) +
  geom_vline(aes(xintercept = 1), alpha = 0.1, lty = "dashed") +
  geom_hline(aes(yintercept = 1), alpha = 0.1, lty = "dashed") +
  geom_point(aes(x=log(mean), y=yi), alpha = 0.2) +
  labs(x = "Log Mean", y = "Log Standard Deviation") +
  guides(fill = "none") +
  theme_classic() 

# Plot together
mean_variance_post_plots <- ((m_sd_strength / m_sd_strength_log) | (m_sd_hypertrophy / m_sd_hypertrophy_log)) + 
  plot_layout(guides = 'collect') + plot_annotation(tag_levels = "A") 

mean_variance_post_plots


### Hypertrophy

data_long_post_hyp <- data_long_post |> 
  filter(outcome == "hypertrophy") |>
  filter(!is.na(yi))

MultiLevelModel_ri_only_log_mean_hyp <- rma.mv(yi, V=vi, data=data_long_post_hyp,
                                               random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es),
                                               mods = ~ log(mean) + group,
                                               method="REML", test="t",
                                               # control=list(optimizer="optim", optmethod="Nelder-Mead")
)


### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_ri_only_log_mean_hyp <- robust(MultiLevelModel_ri_only_log_mean_hyp, data_long_post_hyp$study)


# get the predicted log values
data_long_post_hyp <- cbind(data_long_post_hyp, predict(RobuEstMultiLevelModel_ri_only_log_mean_hyp)) %>%
  mutate(wi = 1/sqrt(vi),
         size = 0.5 + 3.0 * (wi - min(wi))/(max(wi) - min(wi)))

model_m_sd_hypertrophy_log <- ggplot(data_long_post_hyp, aes(x=log(mean), y=yi)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.1) +
  geom_vline(aes(xintercept = 0), alpha = 0.1, lty = "dashed") +
  geom_hline(aes(yintercept = 0), alpha = 0.1, lty = "dashed") +
  geom_point(aes(y=yi, color = group, size = size), alpha = 0.2) +
  geom_ribbon(aes(ymax=ci.ub, ymin=ci.lb, fill = group), alpha = 0.2) +
  geom_line(aes(y=pred, color = group)) +
  annotate(x=max(log(data_long_post_hyp$mean))*0.75,y=max(data_long_post_hyp$yi)*0.25,
           geom = "text",
           label = glue::glue("Between Condition Contrast\n{round(RobuEstMultiLevelModel_ri_only_log_mean_hyp$b[3],2)} [95%CI: {round(RobuEstMultiLevelModel_ri_only_log_mean_hyp$ci.lb[3],2)}, {round(RobuEstMultiLevelModel_ri_only_log_mean_hyp$ci.ub[3],2)}]"),
           size = 2) +
  scale_fill_manual("Group", values = alpha(c("Black", "#E69F00"),0.5)) +
  scale_color_manual("Group", values = alpha(c("Black", "#E69F00"),0.5)) +
  labs(x = "Log Mean Post Score", y = "Log Standard Deviation of the Post Score", color = "Group", shape = "", fill = "") +
  theme_classic() +
  ggtitle("Hypertrophy Outcomes") +
  guides(size = "none", fill = "none")



### Strength

data_long_post_str <- data_long_post |> 
  filter(outcome == "strength") |>
  filter(!is.na(yi))

MultiLevelModel_ri_only_log_mean_str <- rma.mv(yi, V=vi, data=data_long_post_str,
                                               random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es),
                                               mods = ~ log(mean) + group,
                                               method="REML", test="t",
                                               # control=list(optimizer="optim", optmethod="Nelder-Mead")
)


### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_ri_only_log_mean_str <- robust(MultiLevelModel_ri_only_log_mean_str, data_long_post_str$study)


# get the predicted log values
data_long_post_str <- cbind(data_long_post_str, predict(RobuEstMultiLevelModel_ri_only_log_mean_str)) %>%
  mutate(wi = 1/sqrt(vi),
         size = 0.5 + 3.0 * (wi - min(wi))/(max(wi) - min(wi)))

model_m_sd_strength_log <- ggplot(data_long_post_str, aes(x=log(mean), y=yi)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.1) +
  geom_vline(aes(xintercept = 0), alpha = 0.1, lty = "dashed") +
  geom_hline(aes(yintercept = 0), alpha = 0.1, lty = "dashed") +
  geom_point(aes(y=yi, color = group, size = size), alpha = 0.2) +
  geom_ribbon(aes(ymax=ci.ub, ymin=ci.lb, fill = group), alpha = 0.2) +
  geom_line(aes(y=pred, color = group)) +
  annotate(x=max(log(data_long_post_str$mean))*0.75,y=max(data_long_post_str$yi)*0.25,
           geom = "text",
           label = glue::glue("Between Condition Contrast\n{round(RobuEstMultiLevelModel_ri_only_log_mean_str$b[3],2)} [95%CI: {round(RobuEstMultiLevelModel_ri_only_log_mean_str$ci.lb[3],2)}, {round(RobuEstMultiLevelModel_ri_only_log_mean_str$ci.ub[3],2)}]"),
           size = 2) +
  scale_fill_manual("Group", values = alpha(c("Black", "#E69F00"),0.5)) +
  scale_color_manual("Group", values = alpha(c("Black", "#E69F00"),0.5)) +
  labs(x = "Log Mean Post Score", y = "Log Standard Deviation of the Post Score", color = "Group", shape = "", fill = "") +
  theme_classic() +
  ggtitle("Strength Outcomes") +
  guides(size = "none", fill = "none")




model_mean_variance_plots <- (model_m_sd_strength_log | model_m_sd_hypertrophy_log) +
  plot_layout(guides = 'collect', 
              axes = "collect")  &
  theme(legend.position = "bottom")


model_mean_variance_plots

ggsave("model_mean_variance_plots.tiff", plot = model_mean_variance_plots, 
       device = "tiff", dpi = 300,
       w = 7.5, h = 3.75)



#### Looking at Gadbury et al. function

p_from_z = function(x,
                    alternative = "two.sided"){
  if(alternative == "two.sided"){
    2*pnorm(-abs(unlist(x)))
  } else  if (alternative  == "greater"){
    pnorm(x, lower.tail = FALSE)
  } else if (alternative  == "less"){
    pnorm(x, lower.tail = TRUE)
  } else{
    stop("alternative must be two.sided, greater, or less")
  }
  
}

# SDir (normal based test from Hopkins variances estimate)
diff_sdir_test = function(sd1, n1,
                          sd2, n2,
                          alternative = c("two.sided",
                                          "greater",
                                          "less")){
  
  alternative = match.arg(alternative)
  var1 = sd1^2
  var1_se <- var1*(sqrt(2/(n1 - 1)))
  var2 = sd2^2
  var2_se <- var2*(sqrt(2/(n2 - 1)))
  df1 = n1 -1
  df2 = n2 - 1
  if((var1-var2)> 0){
    sd_ir = sqrt(var1-var2)
    
  } else{
    sd_ir = -1*sqrt(var2-var1)
    
  }
  est_diff <- var1-var2
  
  est_diff_SE <- sqrt(2 * (sd1^4/df1 + sd2^4/df2)) 
  teststat <- est_diff/est_diff_SE
  pval = p_from_z(teststat, alternative = alternative)
  
  statistic = teststat
  names(statistic) = "z"
  null1 = 0
  names(null1) = "Standard Deviation of the Individual Response"
  estimate = sd_ir
  names(estimate)= "Standard Deviation of the Individual Response"
  std_err = sqrt(est_diff_SE)
  sum_stats = paste0("SD1 = ", sd1, ", SD2 = ", sd2)
  
  rval <- list(statistic = statistic,
               p.value = as.numeric(pval),
               estimate = estimate,
               stderr = std_err,
               null.value = null1,
               alternative = alternative,
               method = "Difference in Variances",
               data.name = sum_stats)
  class(rval) <- "htest"
  return(rval)
}

# Functions for evaluating subject-treatment interaction between two treatments

#' Calculate MLEs and confidence intervals for treatment effect variance
#' @param treatment Numeric vector of observations from treatment group
#' @param control Numeric vector of observations from control group
#' @param rho_xy Correlation between potential outcomes (-1 to 1)
#' @param conf.level Confidence level for intervals (default 0.95)
#' @return List containing estimates and confidence intervals
estimate_sigma_d <- function(treatment, 
                             control, 
                             n1, n2,
                             rho_xy, 
                             conf.level = 0.95) {
  # Sample sizes
  n1 <- n1
  n2 <- n2
  
  # Calculate sample statistics
  sigma_x_hat <- treatment
  sigma_y_hat <- control
  
  # Estimate sigma_d^2
  sigma_d_sq_hat <- sigma_x_hat^2 + sigma_y_hat^2 - 
    2 * sigma_x_hat * sigma_y_hat * rho_xy
  
  # Calculate variance of sigma_d^2 estimate (equation 4 from paper)
  var_sigma_d_sq <- 2 * (
    (sigma_x_hat^2/n1) * (sigma_x_hat - rho_xy * sigma_y_hat)^2 +
      (sigma_y_hat^2/n2) * (sigma_y_hat - rho_xy * sigma_x_hat)^2
  )
  
  # Calculate confidence intervals
  z_crit <- qnorm(1 - (1 - conf.level)/2)
  ci_lower <- sigma_d_sq_hat - z_crit * sqrt(var_sigma_d_sq)
  ci_upper <- sigma_d_sq_hat + z_crit * sqrt(var_sigma_d_sq)
  
  # Take square root for sigma_d
  sigma_d_hat <- sqrt(sigma_d_sq_hat)
  # remove sign errors
  sigma_d_ci <- sqrt(abs(c(ci_lower, ci_upper))) * sign(c(ci_lower,ci_upper))
  
  return(list(
    sigma_d = sigma_d_hat,
    sigma_d_vi = sqrt(var_sigma_d_sq),
    sigma_d_ci = sigma_d_ci,
    sigma_d_sq = sigma_d_sq_hat,
    sigma_d_sq_ci = c(ci_lower, ci_upper)
  ))
}

#' Calculate probability of unfavorable treatment effect
#' @param treatment Numeric vector of observations from treatment group
#' @param control Numeric vector of observations from control group
#' @param ATE Average treatment effect. A mean effect, like those from a ANCOVA, can be supplied in lieu of using the raw data
#' @param ATE_SE Standard error of the ATE.
#' @param rho_xy Correlation between potential outcomes (-1 to 1)
#' @param conf.level Confidence level for intervals (default 0.95)
#' @param lower.tail determines whether the area under the normal distribution curve is returned to the left or right of a specified value.
#'  Default is FALSE which infers higher scores in treatment group are positive (i.e., positivie ATE is good)
#' @return List containing P- estimate and confidence interval
estimate_p_minus <- function(treatment, 
                             control, 
                             n1,
                             n2,
                             ATE = NULL,
                             ATE_SE = NULL,
                             rho_xy, 
                             conf.level = 0.95,
                             lower.tail = FALSE) {
  # Sample sizes
  n1 <- n1
  n2 <- n2
  lt = lower.tail
  # Calculate sample statistics
  sigma_x_hat <- treatment
  sigma_y_hat <- control
  # Set the average effect of treatment
  if(is.null(ATE)){
    mu_d_hat <- mean(treatment) - mean(control)
  } else{
    mu_d_hat = ATE
  }
  
  
  # Get sigma_d estimates
  sigma_est <- estimate_sigma_d(treatment = treatment, 
                                control = control,
                                n1 = n1,
                                n2 = n2,
                                rho_xy = rho_xy, 
                                conf.level = conf.level)
  sigma_d_hat <- sigma_est$sigma_d
  
  # Calculate P- (probability of unfavorable effect)
  p_minus_hat <- pnorm(mu_d_hat/sigma_d_hat, lower.tail = lt)
  
  
  # Calculate variance components (equation 5 from paper)
  if(is.null(ATE_SE)){
    var_mu_d <- sigma_x_hat^2/n1 + sigma_y_hat^2/n2
  } else{
    var_mu_d = ATE_SE^2
  }
  
  phi_term <- dnorm(mu_d_hat/sigma_d_hat)
  
  var_p_minus <- (phi_term^2/sigma_d_hat^2) * (
    var_mu_d + 
      (mu_d_hat^2 * sigma_est$sigma_d_sq_ci[2])/(4 * sigma_d_hat^4)
  )
  
  # Calculate confidence intervals
  z_crit <- qnorm(1 - (1 - conf.level)/2)
  ci_lower <- p_minus_hat - z_crit * sqrt(var_p_minus)
  ci_upper <- p_minus_hat + z_crit * sqrt(var_p_minus)
  
  return(list(
    p_minus = p_minus_hat,
    ci = c(ci_lower, ci_upper)
  ))
}

#' Bootstrap version for small samples
#' @param treatment Numeric vector of observations from treatment group
#' @param control Numeric vector of observations from control group
#' @param rho_xy Correlation between potential outcomes (-1 to 1)
#' @param conf.level Confidence level for intervals (default 0.95)
#' @param B Number of bootstrap samples (default 1999)
#' @return List containing bootstrap estimates and confidence intervals
bootstrap_sigma_d <- function(treatment, 
                              control,
                              rho_xy, 
                              conf.level = 0.95, 
                              B = 1999) {
  n1 <- length(treatment)
  n2 <- length(control)
  
  # Function to compute sigma_d for a bootstrap sample
  boot_sigma_d <- function(treat_sample, ctrl_sample) {
    sigma_x <- sd(treat_sample) * sqrt(n1/(n1-1))
    sigma_y <- sd(ctrl_sample) * sqrt(n2/(n2-1))
    sigma_d_sq <- sigma_x^2 + sigma_y^2 - 2*sigma_x*sigma_y*rho_xy
    return(sqrt(sigma_d_sq))
  }
  
  # Generate bootstrap samples and compute sigma_d
  boot_estimates <- replicate(B, {
    treat_sample <- sample(treatment, size = n1, replace = TRUE)
    ctrl_sample <- sample(control, size = n2, replace = TRUE)
    boot_sigma_d(treat_sample, ctrl_sample)
  })
  
  # Calculate percentile confidence intervals
  ci <- quantile(boot_estimates, probs = c((1-conf.level)/2, 1-(1-conf.level)/2))
  
  return(list(
    sigma_d = mean(boot_estimates),
    ci = ci,
    boot_samples = boot_estimates
  ))
}

#' Sensitivity analysis across range of rho values
#' @param treatment Numeric vector of observations from treatment group
#' @param control Numeric vector of observations from control group
#' @param ATE Average treatment effect. A mean effect, like those from a ANCOVA, can be supplied in lieu of using the raw data
#' @param rho_seq Sequence of rho values to evaluate
#' @param conf.level Confidence level for intervals
#' @param method Either "mle" or "bootstrap"
#' @param B Number of bootstrap samples if method="bootstrap", default is 1999.
#' @param lower.tail determines whether the area under the normal distribution curve is returned to the left or right of a specified value.
#'  Default is FALSE which infers higher scores in treatment group are positive (i.e., positivie ATE is good)
#' @return Data frame of estimates across rho values
sensitivity_analysis <- function(treatment, control, 
                                 n1, n2,
                                 ATE = NULL,
                                 ATE_SE = NULL,
                                 rho_seq = seq(-1, 1, by = 0.1),
                                 conf.level = 0.95,
                                 method = "mle",
                                 B = 1999,
                                 lower.tail = FALSE) {
  
  results <- lapply(rho_seq, function(rho) {
    if(method == "mle") {
      sigma_d <- estimate_sigma_d(treatment = treatment, 
                                  control = control, 
                                  n1 = n1,
                                  n2 = n2,
                                  rho_xy = rho, 
                                  conf.level = conf.level)
      p_minus <- estimate_p_minus(treatment = treatment, 
                                  control = control, 
                                  n1 = n1,
                                  n2 = n2,
                                  ATE = ATE,
                                  ATE_SE = NULL,
                                  rho_xy = rho, 
                                  conf.level = conf.level,
                                  lower.tail = lower.tail)
      
      data.frame(
        rho = rho,
        sigma_d = sigma_d$sigma_d,
        sigma_d_vi = sigma_d$sigma_d_vi,
        sigma_d_lower = sigma_d$sigma_d_ci[1],
        sigma_d_upper = sigma_d$sigma_d_ci[2],
        p_minus = p_minus$p_minus,
        p_minus_lower = p_minus$ci[1],
        p_minus_upper = p_minus$ci[2]
      )
    } else {
      boot_results <- bootstrap_sigma_d(treatment, control, rho, conf.level, B)
      data.frame(
        rho = rho,
        sigma_d = boot_results$sigma_d,
        sigma_d_lower = boot_results$ci[1],
        sigma_d_upper = boot_results$ci[2]
      )
    }
  })
  
  do.call(rbind, results)
}


# How does SDir assume rho
SDs <- tibble(
  ratio = c(1, seq(0.33, 3, length.out=99)),
  sd1i = rep(10, 100),
  sd2i = sd1i * ratio
) |>
  rowwise() |>
  mutate(
    sd_ir_estimate = diff_sdir_test(sd1 = sd1i,
                                    sd2 = sd2i,
                                    n1 = 100,
                                    n2 = 100)$estimate,
    sdir_rho_xy = (sd1i^2 + sd2i^2 - unname(sd_ir_estimate)^2) / (2 * sd1i * sd2i)
  )


rho_assumptions_plot <- SDs |>
  ggplot(aes(x=sd1i/sd2i, y = sdir_rho_xy)) +
  geom_line() +
  labs(x = "Ratio of Standard Deviations Between Int:Con",
       y = expression("Assumed value of " * rho[{"Int:Con"}])) +
  theme_classic()

ggsave("rho_assumptions_plot.tiff", plot = rho_assumptions_plot, 
       device = "tiff", dpi = 300,
       w = 7.5, h = 3.75)


# Calculate pre-post SDs from SEs
data$RT_pre_sd <- ifelse(is.na(data$RT_pre_se), data$RT_pre_sd, data$RT_pre_se * sqrt(data$RT_n))
data$CON_pre_sd <- ifelse(is.na(data$CON_pre_se), data$CON_pre_sd, data$CON_pre_se * sqrt(data$CON_n))
data$RT_post_sd <- ifelse(is.na(data$RT_post_se), data$RT_post_sd, data$RT_post_se * sqrt(data$RT_n))
data$CON_post_sd <- ifelse(is.na(data$CON_post_se), data$CON_post_sd, data$CON_post_se * sqrt(data$CON_n))


data <- data |>
  filter(!is.na(RT_post_sd) & !is.na(CON_post_sd)) |>
  rowwise() |>
  mutate(
    sd_ir_estimate = diff_sdir_test(sd1 = RT_post_sd,
                                    sd2 = CON_post_sd,
                                    n1 = RT_n,
                                    n2 = CON_n)$estimate,
    sdir_rho_xy = (RT_post_sd^2 + CON_post_sd^2 - unname(sd_ir_estimate)^2) / (2 * RT_post_sd * CON_post_sd)
  )

data <- escalc(measure = "MD",
               m1i = RT_post_m,
               m2i = CON_post_m,
               sd1i = RT_post_sd,
               sd2i = CON_post_sd,
               n1i = RT_n,
               n2i = CON_n,
               var.names = c("ATE", "ATE_vi"),
               data = data)  

data <- data |>
  mutate(ATE_SE = sqrt(ATE_vi)) 


data_test <- data |>
  rowwise() |>
  mutate(
    var_d_mle = list(sensitivity_analysis(control = RT_post_sd,
                                          treatment = CON_post_sd,
                                          n1 = RT_n,
                                          n2 = CON_n,
                                          ATE = ATE,
                                          ATE_SE = ATE_SE,
                                          method = "mle",
                                          lower.tail = FALSE))
  )


data_test_unnest <- unnest(data_test) |>
  mutate(sigma_di_se = sqrt(sigma_d_vi)) # NOTE - need to think about what variance for log transformed sigma_d is


data_test_unnest |>
  ggplot(aes(x=log(abs(ATE)), y=log(sigma_d))) +
  geom_hline(yintercept = 0, color="darkred", alpha = .8) +
  geom_vline(xintercept = 0, color="darkred", alpha = .8) +
  geom_abline(intercept = 0, slope = 1, color="darkred", alpha = .8) +
  # geom_ribbon(fill = "grey",
  #             alpha = .2) +
  geom_point(alpha=0.15) +
  facet_wrap("rho", ncol = 7, nrow = 3) +
  labs(x = "log(ATE)",
       y = bquote(log(SD[D])),
       caption = bquote("Facets are different values of " * rho[RT*","*CON])) +
  theme_classic()

hist(data$sdir_rho_xy, breaks = "fd")


rho_hist <- data_test_unnest |> 
  ggplot(aes(x=sdir_rho_xy)) +
  geom_histogram(color="black") +
  labs(x = expression("Assumed value of " * rho[{"Int:Con"}])) +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

ggsave("rho_hist_plot.tiff", plot = rho_hist, 
       device = "tiff", dpi = 300,
       w = 7.5, h = 3.75)








##### pre-post correlations

##### Read csv as data frame into environment - Note: change source address
Data <- data

# Calculate pre-post SDs from SEs
Data$RT_pre_sd <- ifelse(is.na(Data$RT_pre_se), Data$RT_pre_sd, Data$RT_pre_se * sqrt(Data$RT_n))
Data$CON_pre_sd <- ifelse(is.na(Data$CON_pre_se), Data$CON_pre_sd, Data$CON_pre_se * sqrt(Data$CON_n))
Data$RT_post_sd <- ifelse(is.na(Data$RT_post_se), Data$RT_post_sd, Data$RT_post_se * sqrt(Data$RT_n))
Data$CON_post_sd <- ifelse(is.na(Data$CON_post_se), Data$CON_post_sd, Data$CON_post_se * sqrt(Data$CON_n))

# Convert p to t (Change scores)
Data$RT_delta_t_value <- replmiss(Data$RT_delta_t_value, with(Data, qt(RT_delta_p_value/2, df=RT_n-1, lower.tail=FALSE)))
Data$CON_delta_t_value <- replmiss(Data$CON_delta_t_value, with(Data, qt(CON_delta_p_value/2, df=CON_n-1, lower.tail=FALSE)))

# Convert t to SE (Change scores)
Data$RT_delta_se <- replmiss(Data$RT_delta_se, with(Data, ifelse(is.na(RT_delta_m), 
                                                                 (RT_post_m - RT_pre_m)/RT_delta_t_value, RT_delta_m/RT_delta_t_value)))
Data$CON_delta_se <- replmiss(Data$CON_delta_se, with(Data, ifelse(is.na(CON_delta_m), 
                                                                   (CON_post_m - CON_pre_m)/CON_delta_t_value, CON_delta_m/CON_delta_t_value)))
# Make positive
Data$RT_delta_se <- ifelse(Data$RT_delta_se < 0, Data$RT_delta_se * -1, Data$RT_delta_se)
Data$CON_delta_se <- ifelse(Data$CON_delta_se < 0, Data$CON_delta_se * -1, Data$CON_delta_se)

# Convert CI to SE (Change scores)
Data$RT_delta_se <- replmiss(Data$RT_delta_se, with(Data, (RT_delta_CI_upper - RT_delta_CI_lower)/3.92))
Data$CON_delta_se <- replmiss(Data$CON_delta_se, with(Data, (CON_delta_CI_upper - CON_delta_CI_lower)/3.92))

# Convert SE to SD (Change scores)
Data$RT_delta_sd <- replmiss(Data$RT_delta_sd, with(Data, RT_delta_se * sqrt(RT_n)))
Data$CON_delta_sd <- replmiss(Data$CON_delta_sd, with(Data, CON_delta_se * sqrt(CON_n)))

# Calculate pre-post correlation coefficient for those with pre, post, and delta SDs
Data$RT_ri <- (Data$RT_pre_sd^2 + Data$RT_post_sd^2 - Data$RT_delta_sd^2)/(2 * Data$RT_pre_sd * Data$RT_post_sd)
Data$CON_ri <- (Data$CON_pre_sd^2 + Data$CON_post_sd^2 - Data$CON_delta_sd^2)/(2 * Data$CON_pre_sd * Data$CON_post_sd)

# Remove values outside the range of -1 to +1 as they are likely due to misreporting or miscalculations in original studies
Data$RT_ri <- ifelse(between(Data$RT_ri,-1,1) == FALSE, NA, Data$RT_ri)
Data$CON_ri <- ifelse(between(Data$CON_ri,-1,1) == FALSE, NA, Data$CON_ri)

# Then we'll convert using Fishers r to z, calculate a meta-analytic point estimate, and impute that across the studies with missing correlations
Data <- escalc(measure = "ZCOR", ri = RT_ri, ni = RT_n, data = Data)

Meta_RT_ri <- rma.mv(yi, V=vi, data=Data,
                     slab=paste(label),
                     random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), method="REML", test="t",
                     mods = ~ 0 + outcome,
                     control=list(optimizer="optim", optmethod="Nelder-Mead"))

RobuEstMeta_RT_ri <- robust(Meta_RT_ri, Data$study)

z2r_RT <- psych::fisherz2r(RobuEstMeta_RT_ri$b)
z2r_RT_lower <- psych::fisherz2r(RobuEstMeta_RT_ri$ci.lb)
z2r_RT_upper <- psych::fisherz2r(RobuEstMeta_RT_ri$ci.ub)


Data$RT_ri <- ifelse(is.na(Data$RT_ri), z2r_RT, Data$RT_ri)

Data <- escalc(measure = "ZCOR", ri = CON_ri, ni = CON_n, data = Data)

### Note, data is coded with study and arm as having explicit nesting so all random effects are (~ 1 | study, ~ 1 | arm)
Meta_CON_ri <- rma.mv(yi, V=vi, data=Data,
                      slab=paste(label),
                      random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), method="REML", test="t",
                      mods = ~ 0 + outcome,
                      control=list(optimizer="optim", optmethod="Nelder-Mead"))

RobuEstMeta_CON_ri <- robust(Meta_CON_ri, Data$study)

z2r_CON <- psych::fisherz2r(RobuEstMeta_CON_ri$b[1])
z2r_CON_lower <- psych::fisherz2r(RobuEstMeta_CON_ri$ci.lb)
z2r_CON_upper <- psych::fisherz2r(RobuEstMeta_CON_ri$ci.ub)

Data$CON_ri <- ifelse(is.na(Data$CON_ri), z2r_CON, Data$CON_ri)

# Estimate change score difference SD where only pre-post data available
Data$RT_delta_sd <- replmiss(Data$RT_delta_sd, with(Data, sqrt(RT_pre_sd^2 + RT_post_sd^2 - (2*RT_ri*RT_pre_sd*RT_post_sd))))
Data$CON_delta_sd <- replmiss(Data$CON_delta_sd, with(Data, sqrt(CON_pre_sd^2 + CON_post_sd^2 - (2*CON_ri*CON_pre_sd*CON_post_sd))))




pre_post_cors <- tibble(
  group = c("Intervention","Intervention","Control","Control"),
  outcome = c("Hypertrophy", "Strength", "Hypertrophy", "Strength"),
  cor = psych::fisherz2r(c(RobuEstMeta_RT_ri$b, RobuEstMeta_CON_ri$b)),
  ci.lb = psych::fisherz2r(c(RobuEstMeta_RT_ri$ci.lb, RobuEstMeta_CON_ri$ci.lb)),
  ci.ub = psych::fisherz2r(c(RobuEstMeta_RT_ri$ci.ub, RobuEstMeta_CON_ri$ci.ub))
)


pre_post_cors




# PRE-POST EFECT SIZES AS UNBIASED






### Group by design for comparative treatment standardised effect size calculations 

# Between studies 
Data_between_std <- subset(Data, study_design == "between")
Data_between_std$SD_pool <- sqrt(((Data_between_std$RT_n - 1)*Data_between_std$RT_pre_sd^2 + (Data_between_std$CON_n - 1)*Data_between_std$CON_pre_sd^2) / (Data_between_std$RT_n + Data_between_std$CON_n - 2))

Data_between_std_RT <- escalc(measure="SMCR", m1i=RT_post_m, 
                              m2i=RT_pre_m, sd1i=SD_pool, ni=RT_n, ri=RT_ri, data = Data_between_std)
Data_between_std_CON <- escalc(measure="SMCR", m1i=CON_post_m, 
                               m2i=CON_pre_m, sd1i=SD_pool, ni=CON_n, ri=CON_ri, data = Data_between_std)

Data_between_std$yi <- (Data_between_std_RT$yi - Data_between_std_CON$yi)
Data_between_std$vi <- (Data_between_std_RT$vi + Data_between_std_CON$vi)


# Within participant studies
Data_within_std <- subset(Data, study_design == "within")
Data_within_std$SD_pool <- (((Data_within_std$RT_n - 1)*Data_within_std$RT_pre_sd) + ((Data_within_std$CON_n - 1)*Data_within_std$CON_pre_sd)) / (Data_within_std$RT_n + Data_within_std$CON_n - 2)

Data_within_std_RT <- escalc(measure="SMCR", m1i=RT_post_m, 
                             m2i=RT_pre_m, sd1i=SD_pool, ni=RT_n, ri=RT_ri, data = Data_within_std)
Data_within_std_CON <- escalc(measure="SMCR", m1i=CON_post_m, 
                              m2i=CON_pre_m, sd1i=SD_pool, ni=CON_n, ri=CON_ri, data = Data_within_std)

Data_within_std$yi <- (Data_within_std_RT$yi - Data_within_std_CON$yi)
Data_within_std$vi <- (Data_within_std_RT$vi + Data_within_std_CON$vi)

# recombine standardised effect size 
Data <- rbind(Data_between_std,Data_within_std)


Data <- escalc(measure="SMCR", m1i=RT_post_m, 
                         m2i=RT_pre_m, sd1i=RT_pre_sd, ni=RT_n, ri=RT_ri, data = Data,
                         var.names = c("yi_within", "vi_within"))


Data <- Data %>%
  filter(!is.na(yi) & !is.na(yi_within))


### Differences within studies

between_within_effects <- Data |>
  ggplot(aes(x = yi, y = yi_within)) +
  scale_color_viridis_b() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) + 
  geom_point(alpha = 0.5) +
  labs(x = "Standardised Mean Difference\n(Postive values indicate greater treatment response in RT compared to CON)",
       y = "Standardised Mean Change\n(Postive values indicate greater outcome post RT)") +
  theme_classic()


between_within_effects

ggsave("between_within_effects_plot.tiff", plot = between_within_effects, 
       device = "tiff", dpi = 300,
       w = 7.5, h = 3.75)



