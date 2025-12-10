#'@title MetaRVM Wrapper
#' @description
#' Wrapper function to run MetaRVM given parameters
#' @param config configuration file path or text
#' @param weeks Vector of weeks to simulate
#' @param ts Vector of transmissibility values (length = number of populations)
#' @param pea Vector of pre-existing immunity values (length = number of populations)
#' @param psr Vector of seasonal reduction values (length = number of populations)
#' @param do_chk checkpoint model
#' @param run_id 
#' @return Data table of weekly hospitalizations by age group
#' @export
run_sim <- function(config, weeks, ts = NULL, pea = NULL, psr = NULL,do_chk = F,run_id = NULL,start_date = NULL,week_offset = 0,return.full=F){
  tmp_file = NULL
  if (is.character(config)){
      yml <- yaml::read_yaml(file=config)
  }
  else{
      yml <- yaml::read_yaml(config$config_data$config_file)
  }
  config_obj <- MetaRVM::MetaRVMConfig$new(config)
  if(!is.null(tmp_file)){
    file.remove(tmp_file)
  }
  nsim <- config_obj$config_data$nsim
  npop <- config_obj$config_data$N_pop
  start_date <- config_obj$config_data$start_date
  # reparse vaccine data
  if(!is.null(weeks)){
    #browser()
    config_obj$config_data$sim_length <- 7*(diff(range(weeks))+1)
    
    vac_file <- file.path(dirname(config_obj$config_file),yml$population_data$vaccination)
    vac_dt <- data.table::fread(vac_file, header = TRUE)
    processed_vac <- process_vac_data(vac_dt,
                                      sim_start_date = start_date,
                                      sim_length = config_obj$config_data$sim_length,
                                      delta_t =  config_obj$config_data$delta_t)
    
    vac_time_id <- processed_vac[, c("t")]
    vac_counts <- as.matrix(processed_vac[, -1])
    config_obj$config_data$vac_counts <- vac_counts
    config_obj$config_data$vac_mat <- cbind(unlist(vac_time_id, use.names = F), vac_counts)
  }
  
  
  if(!is.null(ts) & length(ts)>0){
    config_obj$config_data$ts <- matrix(rep(ts, nsim), nrow = nsim, ncol = npop, byrow = TRUE)
  }
  if(!is.null(pea) & length(pea)>0){
    config_obj$config_data$pea <- matrix(rep(pea, nsim), nrow = nsim, ncol = npop, byrow = TRUE)
  }
  if(!is.null(psr) & length(psr)>0){
    config_obj$config_data$psr <- matrix(rep(psr, nsim), nrow = nsim, ncol = npop, byrow = TRUE)
  }
  if(do_chk){
    config_obj$config_data$do_chk <- T
    
    # create set of timesteps for days
    tsteps <- (1:config_obj$config_data$sim_length) / config_obj$config_data$delta_t
    # subset to weekly output
    tsteps <- tsteps[(tsteps %% (7*config_obj$config_data$delta_t) == 0)]
    config_obj$config_data$chk_time_steps <- tsteps
    # expand run_id to include timesteps
    weeks_as_character <- sprintf('%03d',seq(from=min(weeks),to=max(weeks)))
    chk_file_names <- paste0(run_id,weeks_as_character,'.RDS')
    config_obj$config_data$chk_file_names <- matrix(chk_file_names,nrow = 1)
  }
  results <- MetaRVM::metaRVM(config_obj)
  
  res_H <- results$summarize(group_by = c("age"),
                             disease_states = "n_IsympH",
                             stats = c("mean"))
  if(return.full){
    out <- results$results
    out[, `:=`(week = as.numeric(date - (start_date+1)) %/% 7 + 1)]
    out[,week:= week + week_offset]
  return(out)
  }
  res_H_weekly <- res_H$data[, `:=`(week = as.numeric(date - (start_date+1)) %/% 7 + 1)][, .(weekly_H = sum(mean_value),
                                                                                         date=min(date)), by = .(age, week)]
  res_H_weekly[,week:=week + week_offset]
  return(res_H_weekly)
}

#'@title Process vaccine data
#'@param vac_dt data.table of vaccines over time
#'@param sim_start_date starting date of simulation (such as 2025-10-01 for the first day of flu season)
#'@param sim_length number of timesteps to simulate
#'@param delta_t timestep size, defaults to 0.5 (half day)
#' @export
process_vac_data <- function(vac_dt, sim_start_date, sim_length, delta_t) {
  
  # Ensure the date column is of Date type
  vac_dt[[1]] <- as.Date(vac_dt[[1]], tryFormats = c("%m/%d/%Y"))
  
  # Rename the first column to "date" for easier processing
  data.table::setnames(vac_dt, 1, "date")
  
  date_filtered <- vac_dt %>%
    dplyr::filter(date >= as.Date(sim_start_date)) %>%
    dplyr::mutate(t = (date - as.Date(sim_start_date)) / 0.5) %>%
    dplyr::select(-c(date)) %>%
    dplyr::select(dplyr::last_col(), dplyr::everything())
  
  ## fill in the missing time in vac data
  complete_time <- data.table::data.table(t = seq(0, sim_length / delta_t))
  
  ## merge
  vac_all_dates <- data.table::merge.data.table(complete_time, date_filtered, by = "t", all.x = TRUE)
  vac_all_dates[is.na(vac_all_dates)] <- 0
  
  return(vac_all_dates)
}