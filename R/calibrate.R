## update: last column of pars should be config_id (1, 2, ..., length(configs)) and configs will be a vector
## of all configs that are being used for the current running window

#' run batch of simulations
#' @param pars matrix of calibration parameter values to simulate
#' @param weeks integer vector denoting weeks to simulate (example: weeks = 1:10)
#' @param configs vector of configuration files to run \code{MetaRBM}
#' @param metadata \code{data.table}
#' @param return.full boolean 
#' @param week_offset integer
#' @param parameter_bounds named list 
#' @param npop number of subpopulations
#' @param parallel boolean for parallel computing
#' @param ncores number of parallel cores to use
#' @param debug optional parameter used for debugging
#' @returns \code{data.table} of \code{MetaRVM} simulation results
#' @export
run_sims <- function(pars,weeks,configs,metadata = NULL,return.full=F,week_offset=0,parameter_bounds,npop,parallel = F,ncores = 2,debug=F){
  
  if(is.null(metadata)){
    metadata = data.table::data.table(run_id = 1:nrow(pars),do_chk = F)
  }
  # initialize do operator
  `%DO%` <- `%do%`
  if(parallel){
    # setup cluster
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    `%DO%` <- `%dopar%`
  }
  
  sims <- foreach::foreach(i = 1:nrow(pars),.export = c("run_sim","process_vac_data")
                  ,.packages = c("MetaRVM","magrittr","data.table"),.combine = rbind,.inorder = T)  %DO% {
  withCallingHandlers({
  .datatable.aware <- TRUE
  param_row <- as.numeric(pars[i, 1:(ncol(pars)-1)])
  np <- as.integer(length(parameter_bounds) * npop)
  # map param_row to parameter idxs
  parameter_idxs = list()
  k <- 1
  for (name in names(parameter_bounds)){
    parameter_idxs[[name]] <- ((k-1)*npop) + (1:npop)
    k <- k + 1
  }
  
  ts_vals <- param_row[parameter_idxs[["ts"]]]
  pea_vals <- param_row[parameter_idxs[["pea"]]]
  psr_vals <- param_row[parameter_idxs[["psr"]]]
  config <- configs[as.numeric(pars[i,np+1])]
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
  sim_res[, param_id := param_id]
  sim_res
  }
  )
 }
  if(parallel){
  parallel::stopCluster(cl)
  }
  sims
}

#' calibrate_metaRVM
#'@description calibrate metaRVM output to data
#'@param ground_truth output data to calibrate to
#'@param base_config configuration file containing the "static" (i.e., non-calibrated) parameters for `\code{MetaRVM}`
#'@param periods named list where each entry is an integer vector of weeks to calibrate
#'@param parameter_bounds names list of parameter:c(minimum val, maximum val), for example list("ts"=c(0.01,0.5))
#'@param npop number of subpopulations (must be same size as mixing matrices)
#'@param n_checkpoints number of checkpoints to save
#'@param n_trajectories number of trajectories to save
#'@param init_nsim size of Latin Hypercube sample for initial simulations
#'@param grid_size size of Latin Hypercube grid during Thompson sampling
#'@param sim_batch_size number of simulations to acquire in each round of Thompson sampling
#'@param tot_nsim total number of simulations to run for calibration
#'@param forecast_horizon number of weeks to forecast for
#'@param n_forecast_samples number of samples for forecast
#'@param checkpoint_dir directory to save checkpoint files
#'@param moving_window_size when \code{calibration_strategy}="moving_window" size of moving window
#'@param calibration_strategy one of "moving_window" (default) or "distinct" 
#'@param ncores number of cores to use for parallel runs (cannot be larger than parallel::detectCores() - 1)
#'@param verbose boolean option to print output during calibration
#'@returns a named list
#'\itemize{
#'  \item \code{posterior}: \code{n_trajectories} of parameter values
#'  \item \code{sims_post}: simulations from \code{posterior}
#'  \item \code{checkpoints}: \code{data.table} with parameter id and stem of checkpoint file name
#'  \item \code{sims_forecast}: simulations sampled from \code{sims_post} for \code{forecast_horizon} weeks
#'  \item \code{config_files}: vector of configuration files with checkpoints for next period of calibration
#'}
#'@export
calibrate_metaRVM <- function(ground_truth,
                                    base_config,
                                    periods,
                                    parameter_bounds,
                                    npop = 6,
                                    param_offset = 0,
                                    n_checkpoints = 10,
                                    n_trajectories = 50, 
                                    init_nsim = 500, 
                                    grid_size = 500,
                                    sim_batch_size = 500, 
                                    tot_nsim = 1500,
                                    checkpoint_dir = './checkpoints',
                                    n_forecast_samples = 100,
                                    forecast_horizon = 4,
                                    moving_window_size = 4,
                                    calibration_strategy = c("moving_window","distinct"),
                                    ncores = 2,
                                    verbose = T){
  
  calibration_strategy <- match.arg(calibration_strategy,c("moving_window","distinct"))
  if (!dir.exists(checkpoint_dir)){
    dir.create(checkpoint_dir)
  }
  checkpoint_dir <- normalizePath(checkpoint_dir)
  # make sure cores aren't overloaded
  max_cores <- parallel::detectCores() - 1
  ncores <- min(ncores,max_cores)
  # make np from parameter_bounds
  np = as.integer(length(parameter_bounds) * npop)
  # make upper and lower bounds
  lower_bound = c()
  upper_bound = c()
  for (name in names(parameter_bounds)){
    l = rep(parameter_bounds[[name]][1],npop)
    u = rep(parameter_bounds[[name]][2],npop)
    lower_bound = c(lower_bound,l)
    upper_bound = c(upper_bound,u)
  }
  if(verbose){
    cat(msg_header_footer(0.75),'\n')
    cat(sprintf('Starting calibration with the following parameters for %d age groups (%d total):',length(parameter_bounds),np),'\n')
    msg_str = c()
    for (name in names(parameter_bounds)){
      msg_str = c(msg_str,sprintf('>> %s: %s to %s',name,parameter_bounds[[name]][1],parameter_bounds[[name]][2]))
    }
    cat(paste(msg_str,collapse = '\n'),'\n')
  }
  analysis_start_time = proc.time()[3]

if(!is.list(periods)){
  periods = list(periods)
}

checkpoint_prefix <- paste0(checkpoint_dir,'/checkpoint_run_id_%08d_week_')

checkpoint_dir <- dirname(checkpoint_prefix)
data_dir <- dirname(base_config[1])

results = list()

config_files <- base_config
week_offset <- min(periods[[1]]) - 1
output_counter <- 1
output = list()




for (weeks in periods){

  if (verbose){
    cat("> Calibrating weeks ", paste0(weeks,collapse = ","),'\n')
  }

  init_nsim_per_config = floor(init_nsim / length(config_files)) # so that total init_nsim stays same as before

  sims_all = data.table::data.table()
  results <- list()
  res <- calibrate_week(np = np,
                        npop = npop,
                        ground_truth = ground_truth,
                        weeks = weeks,
                        parameter_bounds = parameter_bounds,
                        lower_bound = lower_bound, 
                        upper_bound = upper_bound,
                        configs = config_files,
                        init_nsim_per_config = init_nsim_per_config,
                        sim_batch_size = sim_batch_size,
                        grid_size = grid_size,
                        tot_nsim = tot_nsim,
                        verbose=verbose)
  
  if (verbose){
    cat(">> Trajectories obtained\n") 
  }
  
  
  
  sims_long <- res$sims
  results_long <- cbind(res$Y, res$pars[, 1:np])
  results_long[,param_id:= param_id + param_offset]
  sims_long[,param_id:= param_id + param_offset]
  data.table::setkey(results_long, "simerr")

  # Forecasting
  if(verbose){
    cat(">> Forecasting\n")
  }
  
  sims_forecast <- make_forecast(
    results_long = results_long,
    weeks = weeks,
    config_files = config_files,
    forecast_horizon = forecast_horizon,
    nsamples = n_forecast_samples,
    np = np
  )

  # TODO: possibly combine forecast and checkpointing in one go
  # by defining a checkpointing data

  results_to_checkpoint <- head(results_long,n = n_checkpoints)
  
  pars_rerun <- as.matrix(results_to_checkpoint[, c(3+(1:np), 3), with = F])
  metadata = data.table::data.table(param_id=results_to_checkpoint$param_id)
  metadata[,run_id:=sprintf(checkpoint_prefix,param_id)]
  metadata[,do_chk:=T]
  
  sims_rerun <- run_sims(
    pars = pars_rerun,
    weeks = weeks,
    configs = config_files,
    metadata = metadata,
    parameter_bounds = parameter_bounds,
    return.full = T,
    week_offset = week_offset,
    npop = npop,
    parallel = F,
    ncores = ncores
  )
  cat(">> Checkpointed\n")
  
  if(calibration_strategy == 'moving_window'){
    next_week = max(weeks) + 1 - moving_window_size
    week_offset <- next_week - 1
  }
  else{
    next_week = max(weeks)
    week_offset <- next_week - 1
  }
  checkpoint_files <-  sprintf("%s%03d.RDS",metadata[,run_id],next_week)
  new_start_date = sims_rerun[week==(next_week-1),max(date)] + 1
  #cleanup_config_files(base_config,config_files)
  
  config_files <- sapply(checkpoint_files,make_checkpoint_yaml,
                         base_yaml_file=base_config,
                         new_start_date=new_start_date)
  posterior = head(results_long, n = n_trajectories)
  sims_post = data.table::merge.data.table(sims_long,posterior[,c("param_id","simerr")],by=c("param_id"))
  setkey(sims_post,"simerr") # sort
  out = list(posterior = posterior,
             sims_post = sims_post,
             checkpoints = metadata,
             sims_forecast = sims_forecast,
             config_files = config_files
        )
  # update counters
  output[[output_counter]] <- out
  output_counter <- output_counter + 1
  #week_offset <- next_week - 1
  param_offset <- param_offset + tot_nsim
}



analysis_stop_time = proc.time()[3]
analysis_total_time = as.numeric(analysis_stop_time - analysis_start_time)
h = as.integer(analysis_total_time %/% (60*60))
m = as.integer((analysis_total_time %% 3600) %/% 60)
s = as.integer(analysis_total_time %% 60)
analysis_time_hours_minutes_seconds <- sprintf("%02d:%02d:%02d (hh:mm:ss)", h, m, s)

cat('Analysis finished in ',analysis_time_hours_minutes_seconds,'\n')
cat(msg_header_footer(),'\n')
return(output)
}