prepare_data <- function(data) {
  
  # filter to only untrained participants and between participant RCTs
  data <- data |>
    filter(train_status == "untrained") |>
    filter(study_design == "between")
  
  # get means for post intervention outcomes in long form
  data_long_post_m <- data %>%
    select(study, arm, es, outcome, RT_n, CON_n, RT_post_m, CON_post_m) %>%
    pivot_longer(c(RT_post_m, CON_post_m),
                 names_to = "group",
                 values_to = "mean")
  
  data_long_post_m$group <- recode(data_long_post_m$group, RT_post_m = "RT", CON_post_m = "CON")
  
  # get sds for post intervention outcomes in long form
  data_long_post_sd <- data %>%
    select(study, arm, es, outcome, RT_n, CON_n, RT_post_sd, CON_post_sd) %>%
    pivot_longer(c(RT_post_sd, CON_post_sd),
                 names_to = "group",
                 values_to = "sd")
  
  data_long_post_sd$group <- recode(data_long_post_sd$group, RT_post_sd = "RT", CON_post_sd = "CON")
  
  # recombine all post data in long form
  data_long_post <- cbind(data_long_post_m, sd = data_long_post_sd$sd) 
  
  # recode arms so all are unique for RT groups and CON in each study
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
  
  # recombine and filter any with missing means and sds
  data_long_post <- rbind(data_long_post_RT, data_long_post_CON) %>%
    filter(!is.na(mean) |
             !is.na(sd))
  
  # calculate log SD and variance of log SD
  data_long_post <- escalc(measure = "SDLN",
                           sdi = sd,
                           ni = n,
                           data = data_long_post)
  
}

plot_mean_variance <- function(data_prepared) {
  
  # plot raw means and sds
  m_sd_strength <- ggplot(subset(data_prepared, outcome == "strength"), aes(x=mean, y=sd)) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.1) +
    geom_vline(aes(xintercept = 0), alpha = 0.1, lty = "dashed") +
    geom_hline(aes(yintercept = 0), alpha = 0.1, lty = "dashed") +
    geom_point(aes(x=mean, y=sd), alpha = 0.2) +
    labs(x = "Mean Post Score", y = "Standard Deviation of the Post Score") +
    guides(fill = "none") +
    theme_classic() +
    ggtitle("Strength outcomes")
  
  m_sd_hypertrophy <- ggplot(subset(data_prepared, outcome == "hypertrophy"), aes(x=mean, y=sd)) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.1) +
    geom_vline(aes(xintercept = 0), alpha = 0.1, lty = "dashed") +
    geom_hline(aes(yintercept = 0), alpha = 0.1, lty = "dashed") +
    geom_point(aes(x=mean, y=sd), alpha = 0.2) +
    labs(x = "Mean Post Score", y = "Standard Deviation of the Post Score") +
    guides(fill = "none") +
    theme_classic() +
    ggtitle("Hypertrophy outcomes")
  
  # plot log transformed means and sds
  m_sd_strength_log <- ggplot(subset(data_prepared, outcome == "strength"), aes(x=log(mean), y=yi)) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.1) +
    geom_vline(aes(xintercept = 1), alpha = 0.1, lty = "dashed") +
    geom_hline(aes(yintercept = 1), alpha = 0.1, lty = "dashed") +
    geom_point(aes(x=log(mean), y=yi), alpha = 0.2) +
    labs(x = "Log Mean Post Score", y = "Log Standard Deviation of the Post Score") +
    guides(fill = "none") +
    theme_classic() 
  
  m_sd_hypertrophy_log <- ggplot(subset(data_prepared, outcome == "hypertrophy"), aes(x=log(mean), y=yi)) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.1) +
    geom_vline(aes(xintercept = 1), alpha = 0.1, lty = "dashed") +
    geom_hline(aes(yintercept = 1), alpha = 0.1, lty = "dashed") +
    geom_point(aes(x=log(mean), y=yi), alpha = 0.2) +
    labs(x = "Log Mean Post Score", y = "Log Standard Deviation of the Post Score") +
    guides(fill = "none") +
    theme_classic() 
  
  # Plot together
  mean_variance_post_plots <- ((m_sd_strength / m_sd_strength_log) | (m_sd_hypertrophy / m_sd_hypertrophy_log)) + 
    plot_layout(guides = 'collect') + plot_annotation(tag_levels = "A") 
  
  mean_variance_post_plots
}


fit_mean_var_meta <- function(data_prepared, outcomes) {
  
  # filter to outcome and remove missing rows
  data_prepared <- data_prepared |> 
    filter(outcome == outcomes) |>
    filter(!is.na(yi))
  
  # fit model
  model <- rma.mv(yi, V=vi, data=data_prepared,
                                                 random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es),
                                                 mods = ~ log(mean) + group,
                                                 method="REML", test="t"
  )
  
  # get cluster robust variance estimate
  robust_model <- robust(model, data_prepared$study)
  
  # return robust model
  return(robust_model)
  
}


plot_models <- function(data_prepared, strength_model, hypertrophy_model) {
  
  
  # prepare data getting predicted values from models
  data_prepared_strength <- cbind(data_prepared |> filter(outcome == "strength") |>
                                    filter(!is.na(yi)), predict(strength_model)) %>%
    mutate(wi = 1/sqrt(vi),
           size = 0.5 + 3.0 * (wi - min(wi))/(max(wi) - min(wi)))
  
  data_prepared_hypertrophy <- cbind(data_prepared |> filter(outcome == "hypertrophy") |>
                                       filter(!is.na(yi)), predict(hypertrophy_model)) %>%
    mutate(wi = 1/sqrt(vi),
           size = 0.5 + 3.0 * (wi - min(wi))/(max(wi) - min(wi)))
  
  # make plots
  strength_model_plot <- ggplot(data_prepared_strength, aes(x=log(mean), y=yi)) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.1) +
    geom_vline(aes(xintercept = 0), alpha = 0.1, lty = "dashed") +
    geom_hline(aes(yintercept = 0), alpha = 0.1, lty = "dashed") +
    geom_point(aes(y=yi, color = group, size = size), alpha = 0.2) +
    geom_ribbon(aes(ymax=ci.ub, ymin=ci.lb, fill = group), alpha = 0.2) +
    geom_line(aes(y=pred, color = group)) +
    annotate(x=max(log(data_prepared_strength$mean))*0.75,y=max(data_prepared_strength$yi)*0.25,
             geom = "text",
             label = glue::glue("Between Condition Contrast\n{round(strength_model$b[3],2)} [95%CI: {round(strength_model$ci.lb[3],2)}, {round(strength_model$ci.ub[3],2)}]"),
             size = 2) +
    scale_fill_manual("Group", values = alpha(c("Black", "#E69F00"),0.5)) +
    scale_color_manual("Group", values = alpha(c("Black", "#E69F00"),0.5)) +
    labs(x = "Log Mean Post Score", y = "Log Standard Deviation of the Post Score", color = "Group", shape = "", fill = "") +
    theme_classic() +
    ggtitle("Strength outcomes") +
    guides(size = "none", fill = "none")
  
  hypertrophy_model_plot <- ggplot(data_prepared_hypertrophy, aes(x=log(mean), y=yi)) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.1) +
    geom_vline(aes(xintercept = 0), alpha = 0.1, lty = "dashed") +
    geom_hline(aes(yintercept = 0), alpha = 0.1, lty = "dashed") +
    geom_point(aes(y=yi, color = group, size = size), alpha = 0.2) +
    geom_ribbon(aes(ymax=ci.ub, ymin=ci.lb, fill = group), alpha = 0.2) +
    geom_line(aes(y=pred, color = group)) +
    annotate(x=max(log(data_prepared_hypertrophy$mean))*0.75,y=max(data_prepared_hypertrophy$yi)*0.25,
             geom = "text",
             label = glue::glue("Between Condition Contrast\n{round(hypertrophy_model$b[3],2)} [95%CI: {round(hypertrophy_model$ci.lb[3],2)}, {round(hypertrophy_model$ci.ub[3],2)}]"),
             size = 2) +
    scale_fill_manual("Group", values = alpha(c("Black", "#E69F00"),0.5)) +
    scale_color_manual("Group", values = alpha(c("Black", "#E69F00"),0.5)) +
    labs(x = "Log Mean Post Score", y = "Log Standard Deviation of the Post Score", color = "Group", shape = "", fill = "") +
    theme_classic() +
    ggtitle("Hypertrophy outcomes") +
    guides(size = "none", fill = "none")
  
  model_mean_variance_plots <- (strength_model_plot | hypertrophy_model_plot) +
    plot_layout(guides = 'collect', 
                axes = "collect")  &
    theme(legend.position = "bottom")
  
  
  return(model_mean_variance_plots)
  
}








