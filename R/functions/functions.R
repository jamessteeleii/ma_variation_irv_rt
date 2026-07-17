# Prepare data ----
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

# Plot mean-variance ----
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

# Fit and plot models ----
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
             label = glue::glue("Between Condition Contrast\n{round(strength_model$b[3],2)} [95%CI: {round(strength_model$ci.lb[3],2)}, {round(strength_model$ci.ub[3],3)}]"),
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










# Examine assumptions about rho_Int:Con ----

    # Note, all functions in this indented section are taken from Caldwell et al. code - https://doi.org/10.5281/zenodo.15492100
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
    

plot_rho_assumptions <- function() {
  
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
  
  return(rho_assumptions_plot)
  
}
  
  

## SDir assumptions regarding rho_Int:Con in our dataset ----

check_rho_assumptions_data <- function(data) {
  
  # calculate missing pre and post SDs from SEs
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
}

plot_rho_assumptions_data <- function(rho_assumptions_data) {
  
  rho_hist <- rho_assumptions_data |> 
    ggplot(aes(x=sdir_rho_xy)) +
    geom_histogram(color="black", linewidth =  .5) +
    # geom_vline(aes(xintercept = median(sdir_rho_xy))) +
    # geom_vline(aes(xintercept = quantile(sdir_rho_xy, .25)), linetype = "dashed") +
    # geom_vline(aes(xintercept = quantile(sdir_rho_xy, .75)), linetype = "dashed") +
    
    labs(x = expression("Assumed value of " * rho[{"Int:Con"}])) +
    theme_classic() +
    theme(axis.title.y = element_blank(),
          axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  
  return(rho_hist)
  
}

## Pre-post as unbiased estimator of intervention effect ----
get_between_within_ates <- function(data) {
  
  # filter to only untrained participants and between participant RCTs
  data <- data |>
    filter(train_status == "untrained") |>
    filter(study_design == "between")
  
  # calculate missing pre and post SDs from SEs
  data$RT_pre_sd <- ifelse(is.na(data$RT_pre_se), data$RT_pre_sd, data$RT_pre_se * sqrt(data$RT_n))
  data$CON_pre_sd <- ifelse(is.na(data$CON_pre_se), data$CON_pre_sd, data$CON_pre_se * sqrt(data$CON_n))
  data$RT_post_sd <- ifelse(is.na(data$RT_post_se), data$RT_post_sd, data$RT_post_se * sqrt(data$RT_n))
  data$CON_post_sd <- ifelse(is.na(data$CON_post_se), data$CON_post_sd, data$CON_post_se * sqrt(data$CON_n))
  
  # # convert p to t (Change scores)
  # data$RT_delta_t_value <- replmiss(data$RT_delta_t_value, with(data, qt(RT_delta_p_value/2, df=RT_n-1, lower.tail=FALSE)))
  # data$CON_delta_t_value <- replmiss(data$CON_delta_t_value, with(data, qt(CON_delta_p_value/2, df=CON_n-1, lower.tail=FALSE)))
  # 
  # # convert t to SE (Change scores)
  # data$RT_delta_se <- replmiss(data$RT_delta_se, with(data, ifelse(is.na(RT_delta_m), 
  #                                                                  (RT_post_m - RT_pre_m)/RT_delta_t_value, RT_delta_m/RT_delta_t_value)))
  # data$CON_delta_se <- replmiss(data$CON_delta_se, with(data, ifelse(is.na(CON_delta_m), 
  #                                                                    (CON_post_m - CON_pre_m)/CON_delta_t_value, CON_delta_m/CON_delta_t_value)))
  # make positive
  data$RT_delta_se <- ifelse(data$RT_delta_se < 0, data$RT_delta_se * -1, data$RT_delta_se)
  data$CON_delta_se <- ifelse(data$CON_delta_se < 0, data$CON_delta_se * -1, data$CON_delta_se)
  
  # convert CI to SE (Change scores)
  data$RT_delta_se <- replmiss(data$RT_delta_se, with(data, (RT_delta_CI_upper - RT_delta_CI_lower)/3.92))
  data$CON_delta_se <- replmiss(data$CON_delta_se, with(data, (CON_delta_CI_upper - CON_delta_CI_lower)/3.92))
  
  # convert SE to SD (Change scores)
  data$RT_delta_sd <- replmiss(data$RT_delta_sd, with(data, RT_delta_se * sqrt(RT_n)))
  data$CON_delta_sd <- replmiss(data$CON_delta_sd, with(data, CON_delta_se * sqrt(CON_n)))
  
  # calculate pre-post correlation coefficient for those with pre, post, and delta SDs
  data$RT_ri <- (data$RT_pre_sd^2 + data$RT_post_sd^2 - data$RT_delta_sd^2)/(2 * data$RT_pre_sd * data$RT_post_sd)
  data$CON_ri <- (data$CON_pre_sd^2 + data$CON_post_sd^2 - data$CON_delta_sd^2)/(2 * data$CON_pre_sd * data$CON_post_sd)
  
  # remove values outside the range of -1 to +1 as they are likely due to misreporting or miscalculations in original studies
  data$RT_ri <- ifelse(between(data$RT_ri,-1,1) == FALSE, NA, data$RT_ri)
  data$CON_ri <- ifelse(between(data$CON_ri,-1,1) == FALSE, NA, data$CON_ri)
  
  # convert using Fishers r to z
  data <- escalc(measure = "ZCOR", ri = RT_ri, ni = RT_n, data = data,
                 var.names = c("RT_yi", "RT_vi"))
  
  # convert using Fishers r to z
  data <- escalc(measure = "ZCOR", ri = CON_ri, ni = CON_n, data = data,
                 var.names = c("CON_yi", "CON_vi"))
  
  
  # get meta-analytic estimate for RT arms
  Meta_RT_ri <- rma.mv(RT_yi, V=RT_vi, data=data,
                       slab=paste(label),
                       random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), method="REML", test="t",
                       mods = ~ 0 + outcome,
                       control=list(optimizer="optim", optmethod="Nelder-Mead"))
  
  RobuEstMeta_RT_ri <- robust(Meta_RT_ri, data$study)
  
  z2r_RT <- psych::fisherz2r(RobuEstMeta_RT_ri$b)
  z2r_RT_lower <- psych::fisherz2r(RobuEstMeta_RT_ri$ci.lb)
  z2r_RT_upper <- psych::fisherz2r(RobuEstMeta_RT_ri$ci.ub)
  
  # get meta-analytic estimate for CON arms
  Meta_CON_ri <- rma.mv(CON_yi, V=CON_vi, data=data,
                        slab=paste(label),
                        random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), method="REML", test="t",
                        mods = ~ 0 + outcome,
                        control=list(optimizer="optim", optmethod="Nelder-Mead"))
  
  RobuEstMeta_CON_ri <- robust(Meta_CON_ri, data$study)
  
  z2r_CON <- psych::fisherz2r(RobuEstMeta_CON_ri$b)
  z2r_CON_lower <- psych::fisherz2r(RobuEstMeta_CON_ri$ci.lb)
  z2r_CON_upper <- psych::fisherz2r(RobuEstMeta_CON_ri$ci.ub)
  
  # Impute mising pre-post correlations
  data <- data |>
    mutate(
      RT_ri = case_when(
        is.na(RT_ri) & outcome == "hypertrophy" ~ z2r_RT[1],
        is.na(RT_ri) & outcome == "strength" ~ z2r_RT[2],
        .default = RT_ri
      ),
      CON_ri = case_when(
        is.na(CON_ri) & outcome == "hypertrophy" ~ z2r_CON[1],
        is.na(CON_ri) & outcome == "strength" ~ z2r_CON[2],
        .default = CON_ri
      )
    )
  
  
  
  # Between arms 
  data <- subset(data, study_design == "between")
  data$SD_pool <- sqrt(((data$RT_n - 1)*data$RT_pre_sd^2 + (data$CON_n - 1)*data$CON_pre_sd^2) / (data$RT_n + data$CON_n - 2))
  
  data_RT <- escalc(measure="SMCR", m1i=RT_post_m, 
                                m2i=RT_pre_m, sd1i=SD_pool, ni=RT_n, ri=RT_ri, data = data)
  data_CON <- escalc(measure="SMCR", m1i=CON_post_m, 
                                 m2i=CON_pre_m, sd1i=SD_pool, ni=CON_n, ri=CON_ri, data = data)
  
  data$yi_between <- (data_RT$yi - data_CON$yi)
  data$vi_between <- (data_RT$vi + data_CON$vi)
  
  # Within arms
  data <- escalc(measure="SMCR", m1i=RT_post_m, 
                 m2i=RT_pre_m, sd1i=RT_pre_sd, ni=RT_n, ri=RT_ri, data = data,
                 var.names = c("yi_within", "vi_within"))
  
  
  data <- data %>%
    filter(!is.na(yi_between) & !is.na(yi_within))
  
  return(data)
}

plot_between_within_ates <- function(data) {
  between_within_ates <- data |>
    ggplot(aes(x = yi_between, y = yi_within)) +
    scale_color_viridis_b() +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_abline(intercept = 0, slope = 1, linetype = 2) + 
    geom_point(alpha = 0.5) +
    labs(x = "Standardised Mean Difference\n(Postive values indicate greater treatment response in RT compared to CON)",
         y = "Standardised Mean Change\n(Postive values indicate greater outcome post RT)") +
    theme_classic()
  
  
  return(between_within_ates)
}

## Comparing to rho_Pre:Post ----

get_pre_post_rho <- function(data) {
  
  # calculate missing pre and post SDs from SEs
  data$RT_pre_sd <- ifelse(is.na(data$RT_pre_se), data$RT_pre_sd, data$RT_pre_se * sqrt(data$RT_n))
  data$CON_pre_sd <- ifelse(is.na(data$CON_pre_se), data$CON_pre_sd, data$CON_pre_se * sqrt(data$CON_n))
  data$RT_post_sd <- ifelse(is.na(data$RT_post_se), data$RT_post_sd, data$RT_post_se * sqrt(data$RT_n))
  data$CON_post_sd <- ifelse(is.na(data$CON_post_se), data$CON_post_sd, data$CON_post_se * sqrt(data$CON_n))
  
  # # convert p to t (Change scores)
  # data$RT_delta_t_value <- replmiss(data$RT_delta_t_value, with(data, qt(RT_delta_p_value/2, df=RT_n-1, lower.tail=FALSE)))
  # data$CON_delta_t_value <- replmiss(data$CON_delta_t_value, with(data, qt(CON_delta_p_value/2, df=CON_n-1, lower.tail=FALSE)))
  # 
  # # convert t to SE (Change scores)
  # data$RT_delta_se <- replmiss(data$RT_delta_se, with(data, ifelse(is.na(RT_delta_m), 
  #                                                                  (RT_post_m - RT_pre_m)/RT_delta_t_value, RT_delta_m/RT_delta_t_value)))
  # data$CON_delta_se <- replmiss(data$CON_delta_se, with(data, ifelse(is.na(CON_delta_m), 
  #                                                                    (CON_post_m - CON_pre_m)/CON_delta_t_value, CON_delta_m/CON_delta_t_value)))
  # make positive
  data$RT_delta_se <- ifelse(data$RT_delta_se < 0, data$RT_delta_se * -1, data$RT_delta_se)
  data$CON_delta_se <- ifelse(data$CON_delta_se < 0, data$CON_delta_se * -1, data$CON_delta_se)
  
  # convert CI to SE (Change scores)
  data$RT_delta_se <- replmiss(data$RT_delta_se, with(data, (RT_delta_CI_upper - RT_delta_CI_lower)/3.92))
  data$CON_delta_se <- replmiss(data$CON_delta_se, with(data, (CON_delta_CI_upper - CON_delta_CI_lower)/3.92))
  
  # convert SE to SD (Change scores)
  data$RT_delta_sd <- replmiss(data$RT_delta_sd, with(data, RT_delta_se * sqrt(RT_n)))
  data$CON_delta_sd <- replmiss(data$CON_delta_sd, with(data, CON_delta_se * sqrt(CON_n)))
  
  # calculate pre-post correlation coefficient for those with pre, post, and delta SDs
  data$RT_ri <- (data$RT_pre_sd^2 + data$RT_post_sd^2 - data$RT_delta_sd^2)/(2 * data$RT_pre_sd * data$RT_post_sd)
  data$CON_ri <- (data$CON_pre_sd^2 + data$CON_post_sd^2 - data$CON_delta_sd^2)/(2 * data$CON_pre_sd * data$CON_post_sd)
  
  # remove values outside the range of -1 to +1 as they are likely due to misreporting or miscalculations in original studies
  data$RT_ri <- ifelse(between(data$RT_ri,-1,1) == FALSE, NA, data$RT_ri)
  data$CON_ri <- ifelse(between(data$CON_ri,-1,1) == FALSE, NA, data$CON_ri)
  
  # convert using Fishers r to z
  data <- escalc(measure = "ZCOR", ri = RT_ri, ni = RT_n, data = data,
                 var.names = c("RT_yi", "RT_vi"))
  
  # convert using Fishers r to z
  data <- escalc(measure = "ZCOR", ri = CON_ri, ni = CON_n, data = data,
                 var.names = c("CON_yi", "CON_vi"))
  
  # get meta-analytic estimate for diff in correlations
  data$diff_yi <- data$RT_yi - data$CON_yi
  data$diff_vi <- data$RT_vi + data$CON_vi
  
  # get meta-analytic estimate for diff between arms
  Meta_diff_ri <- rma.mv(diff_yi, V=diff_vi, data=data,
                       slab=paste(label),
                       random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), method="REML", test="t",
                       mods = ~ 0 + outcome,
                       control=list(optimizer="optim", optmethod="Nelder-Mead"))
  
  RobuEstMeta_diff_ri <- robust(Meta_diff_ri, data$study)
  
  RobuEstMeta_diff_ri
  
  
  # get meta-analytic estimate for RT arms
  Meta_RT_ri <- rma.mv(RT_yi, V=RT_vi, data=data,
                       slab=paste(label),
                       random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), method="REML", test="t",
                       mods = ~ 0 + outcome,
                       control=list(optimizer="optim", optmethod="Nelder-Mead"))
  
  RobuEstMeta_RT_ri <- robust(Meta_RT_ri, data$study)
  
  z2r_RT <- psych::fisherz2r(RobuEstMeta_RT_ri$b)
  z2r_RT_lower <- psych::fisherz2r(RobuEstMeta_RT_ri$ci.lb)
  z2r_RT_upper <- psych::fisherz2r(RobuEstMeta_RT_ri$ci.ub)
  
  # get meta-analytic estimate for CON arms
  Meta_CON_ri <- rma.mv(CON_yi, V=CON_vi, data=data,
                        slab=paste(label),
                        random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), method="REML", test="t",
                        mods = ~ 0 + outcome,
                        control=list(optimizer="optim", optmethod="Nelder-Mead"))
  
  RobuEstMeta_CON_ri <- robust(Meta_CON_ri, data$study)
  
  z2r_CON <- psych::fisherz2r(RobuEstMeta_CON_ri$b)
  z2r_CON_lower <- psych::fisherz2r(RobuEstMeta_CON_ri$ci.lb)
  z2r_CON_upper <- psych::fisherz2r(RobuEstMeta_CON_ri$ci.ub)
  
  # combine pre-post cors
  pre_post_cors <- tibble(
    group = c("Intervention","Intervention","Control","Control"),
    outcome = c("Hypertrophy", "Strength", "Hypertrophy", "Strength"),
    cor = psych::fisherz2r(c(RobuEstMeta_RT_ri$b, RobuEstMeta_CON_ri$b)),
    ci.lb = psych::fisherz2r(c(RobuEstMeta_RT_ri$ci.lb, RobuEstMeta_CON_ri$ci.lb)),
    ci.ub = psych::fisherz2r(c(RobuEstMeta_RT_ri$ci.ub, RobuEstMeta_CON_ri$ci.ub))
  )
  
  checks_for_cors <- list(RobuEstMeta_diff_ri,
                          pre_post_cors)
  
  return(checks_for_cors)

}