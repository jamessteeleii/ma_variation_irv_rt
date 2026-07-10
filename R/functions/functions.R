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
    

check_rho_assumptions <- function() {
  
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


