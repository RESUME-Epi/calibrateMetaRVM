# acquire new points

## calibrate_week takes configs instead of one config and returns tot_nsim simulations.
## So when we call calibrate_week in a calibration window, we would get tot_nsim trajectories
## out instead of length(config)*tot_nsim trajectories

#' calibrate_week
#' @description
#' Calibrate metaRVM model outputs to data. This function creates a Latin Hypercube sample (LHS) of the calibration parameters
#' and then adaptively explores parameter space using Thompson sampling. The algorith works as follows:
#' \itemize{
#' 1. Create an LHS of size \code{init_nsim_per_config}
#' 2. Run simulations for those parameters
#' 3. Assess the error between model outputs and data
#' 4. Fit a Gaussian process emulator on the (parameters, error) observations
#' 5. Perform Thompson sampling to acquire new parameters and run new simulations, refitting GPs each time
#' 6. Return results of parameters, errors, and simulations
#' }
#' @param np number of dimensions to calibrate
#' @param npop number of populations (same dimension as mixing matrices)
#' @param weeks vector of weeks to calibrate to
#' @param configs vector of yaml file specifying potential configurations
#' @param ground_truth ground truth data to calibrate to
#' @param parameter_bounds named list of parameter -> range (example: ts: (0.01,0.50)) to be passed to run_sims
#' @param lower_bound vector of lower bounds for calibration parameters
#' @param upper_bound vector of upper bounds for calibration parameters
#' @param init_nsim_per_config number of simulations per config file
#' @param sim_batch_size number of sims to run per emulator acquisition
#' @param grid_size size of Latin Hypercube grid to explore per acquisition
#' @param tot_nsim total number of simulations to run
#' @param error_function (data,ground_truth) discrepancy function. Defaults to squared log loss.
#' @param verbose print output during calibration (default F)
#' @returns a named list with 
#' \itemize{
#' \item \code{pars}: parameters acquired during calibration
#' \item \code{Y}: errors between simulations and data
#' \item \code{sims}: simulations corresponding to \code{pars}
#' }
calibrate_week <- function(np,
                           npop,
                           weeks,
                           configs,
                           ground_truth,
                           parameter_bounds,
                           lower_bound,
                           upper_bound,
                           init_nsim_per_config,
                           sim_batch_size,
                           grid_size=100,
                           tot_nsim=100,
                           error_function = NULL,verbose=F){
  if (is.null(error_function)){
    error_function <- squared_log_loss
  }
  nconfig <- length(configs)
  init_nsim <- init_nsim_per_config * nconfig
  # create a design of size (init_nsim x (np+1)), the last one is the config
  x01 <- lhs::randomLHS(init_nsim, np)
  x <- sapply(1:np, function(i) x01[, i] * (upper_bound[i] - lower_bound[i]) + lower_bound[i])

  x01_aug <- cbind(x01[rep(1:init_nsim_per_config, nconfig), ], rep(1:nconfig, each = init_nsim_per_config))
  x_aug <- cbind(x[rep(1:init_nsim_per_config, nconfig), ], rep(1:nconfig, each = init_nsim_per_config))
  
  x01_all <- x01_aug
  x_all <- x_aug
  sims <- run_sims(x_aug,
                    weeks = weeks,
                    configs = configs,
                    metadata = NULL,return.full = F,week_offset = min(weeks) -1,
                    parameter_bounds = parameter_bounds,
                    npop = npop,parallel = F)
  sims_all <- sims
  
  sims_with_ground_truth <- data.table::merge.data.table(sims,ground_truth,by = c("week","age"))
  simerr = sims_with_ground_truth[,.(simerr = error_function(weekly_H,ct_recs_c)),by=param_id]
  simerr[, `:=`(config_idx = x_aug[, ncol(x_aug)])]

  Y_all <- simerr # book-keeping
  ## here we introduce composite GP surrogate 
  
  f <- fit_composite_gp(x01_aug, simerr$simerr, simerr$config_idx, np = np)
  gps <- f$gp_list
  m_sd_mat <- f$m_sd_mat
  nsim <- init_nsim

  
  m <- 2
  ## Thompson sampling
  while (nsim < tot_nsim){
    
    ## keep sampling until the batch size is acheived
    batch_size <- 0
    newx_all <- c()
    while(batch_size < sim_batch_size){
      ## define grid
      xcan <- lhs::randomLHS(grid_size, np)
      xcan_aug <- cbind(xcan[rep(1:grid_size, nconfig), ], rep(1:nconfig, each = grid_size))
      if(verbose){
      cat(">>>> Random grid initialized\n")
      }
      # draw samle from composite GP
      tTS <- draw_sample_composite_gp(gps, xcan_aug, sim_batch_size)     
      best_ids <- unique(apply(tTS, 1, which.min))
      
      newx_aug <- xcan_aug[best_ids, ]
      if(!is.matrix(newx_aug)) newx_aug <- matrix(newx_aug, nrow = 1)
      
      batch_size <- batch_size + length(best_ids)
      newx_all <- rbind(newx_all, newx_aug)
      if(verbose){
      cat(">>>> Batch size = ", batch_size,'\n')
      }
    }
    # trim the newx if the size exceeds sim_batch_size
    newx_all <- newx_all[1:sim_batch_size, ]
        
    ## run new simulations
    newx_native <- cbind(sapply(1:np, function(i) newx_all[, i] * (upper_bound[i] - lower_bound[i]) + lower_bound[i]), newx_all[, np+1])
    
    sims <- run_sims(newx_native,
                      weeks = weeks,
                      configs = configs,
                      metadata = NULL,
                     return.full = F,
                     week_offset = min(weeks) -1,
                     parameter_bounds = parameter_bounds,
                     npop = npop,
                     parallel=F)
    # update param_id
    sims[,param_id:=param_id + max(sims_all$param_id)]
    sims_all <- data.table::rbindlist(list(sims_all,sims))
    
    sims_with_ground_truth <- data.table::merge.data.table(sims,ground_truth,by = c("week","age"))
    simerr = sims_with_ground_truth[,.(simerr = error_function(weekly_H,ct_recs_c)),by=param_id]
    simerr[, `:=`(config_idx = newx_all[, ncol(newx_all)])]

    ## update GPs
    gps <- update_composite_gp(gps, newx_all, simerr$simerr, simerr$config_idx, m_sd_mat)    
    
    nsim <- nsim + nrow(newx_native)
    
    ## store the lasts and hosps
    x01_all <- rbind(x01_all, newx_all)
    x_all <- rbind(x_all, newx_native)
    Y_all <- rbind(Y_all,simerr)
    if(verbose){
      cat(">>>> Total simulations = ", nsim,'\n')
    }
    m <- m + 1
  }
  
  
  return(list(pars = x_all,
              Y = Y_all,
              sims = sims_all))
}