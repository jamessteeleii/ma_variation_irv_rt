# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c(
    "tidyverse",
    "here",
    "metafor",
    "bayestestR",
    "patchwork"
  ), # Packages that your targets need for their tasks.
  memory = "transient",
  format = "qs",  # Optionally set the default storage format. qs is fast.
  garbage_collection = TRUE,
  storage = "worker",
  retrieval = "worker"
  #
  # Pipelines that take a long time to run may benefit from
  # optional distributed computing. To use this capability
  # in tar_make(), supply a {crew} controller
  # as discussed at https://books.ropensci.org/targets/crew.html.
  # Choose a controller that suits your needs. For example, the following
  # sets a controller that scales up to a maximum of two workers
  # which run as local R processes. Each worker launches when there is work
  # to do and exits if 60 seconds pass with no tasks to run.
  #
  #   controller = crew::crew_controller_local(workers = 2, seconds_idle = 60)
  #
  # Alternatively, if you want workers to run on a high-performance computing
  # cluster, select a controller from the {crew.cluster} package.
  # For the cloud, see plugin packages like {crew.aws.batch}.
  # The following example is a controller for Sun Grid Engine (SGE).
  #
  #   controller = crew.cluster::crew_controller_sge(
  #     # Number of workers that the pipeline can scale up to:
  #     workers = 10,
  #     # It is recommended to set an idle time so workers can shut themselves
  #     # down if they are not running tasks.
  #     seconds_idle = 120,
  #     # Many clusters install R as an environment module, and you can load it
  #     # with the script_lines argument. To select a specific verison of R,
  #     # you may need to include a version string, e.g. "module load R/4.3.2".
  #     # Check with your system administrator if you are unsure.
  #     script_lines = "module load R"
  #   )
  #
  # Set other options as needed.
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source("R/functions/.")
# tar_source("other_functions.R") # Source other scripts as needed.

# Replace the target list below with your own:
list(
  
  # Read and prepare data ----
  
  # read data from Steele et al. (2023) repository
  tar_target(
    data,
    readr::read_csv(url("https://github.com/jamessteeleii/Meta-Analysis-of-Variation-in-Resistance-Training/raw/refs/heads/main/data/Polito%20et%20al.%20RT%20Extracted%20Data.csv"))
  ),
  
  # prepare data
  tar_target(
    data_prepared,
    prepare_data(data)
  ),
  
  # Mean variance relationship ----
  tar_target(
    mean_var_plot,
    plot_mean_variance(data_prepared)
  ),
  
  tar_target(
    mean_var_plot_tiff,
    {
      ggsave(
        plot = mean_var_plot,
        filename = "plots/mean_var_plot.tiff",
        device = "tiff",
        dpi = 300,
        width = 10,
        height = 7.5
      )
    }
  ),
  
  # Fit models ----
  tar_target(
    strength_model,
    fit_mean_var_meta(data_prepared, "strength")
  ),
  
  tar_target(
    hypertrophy_model,
    fit_mean_var_meta(data_prepared, "hypertrophy")
  ),
  
  # Plot models ----
  tar_target(
    models_plot,
    plot_models(data_prepared, strength_model, hypertrophy_model)
  ),
  
  tar_target(
    models_plot_tiff,
    {
      ggsave(
        plot = models_plot,
        filename = "plots/models_plot.tiff",
        device = "tiff",
        dpi = 300,
        width = 10,
        height = 5
      )
    }
  ),
  
  # Check rho_Int:Con assumptions ----
  tar_target(
    rho_assumptions_plot,
    plot_rho_assumptions()
  ),
  
  tar_target(
    rho_assumptions_plot_tiff,
    {
      ggsave(
        plot = rho_assumptions_plot,
        filename = "plots/rho_assumptions_plot.tiff",
        device = "tiff",
        dpi = 300,
        width = 8,
        height = 5
      )
    }
  ),
  
  tar_target(
    between_within_ates,
    get_between_within_ates(data)
  ),
  
  tar_target(
    between_within_ates_plot,
    plot_between_within_ates(between_within_ates)
  ),
  
  tar_target(
    between_within_ates_plot_tiff,
    {
      ggsave(
        plot = between_within_ates_plot,
        filename = "plots/between_within_ates_plot.tiff",
        device = "tiff",
        dpi = 300,
        width = 8,
        height = 5
      )
    }
  ),
  
  tar_target(
    rho_assumptions_data,
    check_rho_assumptions_data(data)
  ),
  
  tar_target(
    rho_assumptions_data_plot,
    plot_rho_assumptions_data(rho_assumptions_data)
  ),
  
  tar_target(
    rho_assumptions_data_plot_tiff,
    {
      ggsave(
        plot = rho_assumptions_data_plot,
        filename = "plots/rho_assumptions_data_plot.tiff",
        device = "tiff",
        dpi = 300,
        width = 8,
        height = 5
      )
    }
  ),
  
  tar_target(
    pre_post_rho_checks,
    get_pre_post_rho(data)
  ),
  
  # Add grateful report ----
  
  tar_target(
    grateful_report,
    grateful::cite_packages(out.dir = ".", out.format = "html")
  )
  
)
