#' run_forecast_sims
#' @param pars matrix of parameters
#' @param weeks vector of weeks to simulate for
#' @param configs vector of config files
#' @param return.F boolean to return full output (default F)
#' @param week_offset increments the starting week (default 0)
#' @export
run_forecast_sims <- function(pars,weeks,configs,metadata = NULL,return.full=F, week_offset=0){
  if(is.null(metadata)){
    metadata = data.table::data.table(run_id = 1:nrow(pars),do_chk = F)
  }
  sims <- foreach::foreach(i = 1:nrow(pars),
                           export = c("run_sim","process_vac_data","base_config."),
                           .packages = c("MetaRVM","magrittr"),
                           .combine = rbind,.inorder = T)  %do% {
  withCallingHandlers({
  param_row <- as.numeric(pars[i, 1:(ncol(pars)-1)])
  ts_vals <- param_row[1:6]
  pea_vals <- NULL
  psr_vals <- param_row[7:12]
  config <- configs[as.numeric(pars[i,13])]
  
  # adjust metadata
  metadata_row <- metadata[i,]
  do_chk <- metadata_row[['do_chk']]
  run_id <- metadata_row[['run_id']]
  param_id <- metadata_row[['param_id']]
  if(is.null(param_id)){
    param_id <- i
  }
  sim_res <- run_sim(
    config = config,
    weeks = weeks,
    ts = ts_vals,
    pea = pea_vals,
    psr = psr_vals,
    do_chk = do_chk,
    run_id = run_id,
    week_offset = week_offset,
    return.full = return.full
  )
  
  sim_res[, param_id := i]
  sim_res
  }
  )
 }
  sims
}

make_forecast <- function(
    results_long,
    weeks,
    config_files,
    forecast_horizon = 4,
    nsamples = 100,
    np
){

  # 1. Select the best parameter sets (the posterior)
  posterior <- results_long

  # 2. Weigh them based on error. A lower error gets a higher weight.
  # We use exp(-error) and normalize. Subtracting the max error improves numerical stability.
  weights <- exp(-(posterior$simerr - max(posterior$simerr)))
  normalized_weights <- weights / sum(weights)

  # 3. Sample from the posterior with replacement, using the weights.
  n_trajectories <- nsamples
  sampled_indices <- sample(1:nrow(posterior), size = n_trajectories, replace = TRUE, prob = normalized_weights)

  results_forecast <- posterior[sampled_indices, ]
  pars_rerun <- posterior[sampled_indices, c(3+(1:np), 3), with = F]

  # 4. Rerun 
  forecast_horizon <- 4
  sims_forecast <- run_forecast_sims(
      pars = pars_rerun,
      weeks = c(weeks, max(weeks) + (1:forecast_horizon)),
      configs = config_files,
      metadata = NULL,
      return.full = TRUE,
      week_offset = min(weeks) -1
  )
  return(sims_forecast)
}

