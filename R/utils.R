#'@description
#'Create checkpoint yaml from a checkpoint file and a base yaml file and add a new start date
#'
#'@param checkpoint_file
#'@param base_yaml_file base config file that reads in vaccination, population, mixing matrices
#'@param new_start_date new starting date of simulation (one day after the end of the previous simulation)
#'@returns yaml object to be used as config file for run_sims
#'@note this function returns a yaml object as text so as to not write many, many yaml files to disk (which could)
#'be problematic when running many simulations
#'@export
make_checkpoint_yaml <- function(checkpoint_file,base_yaml_file,new_start_date){
  checkpoint_dir <- dirname(checkpoint_file)
  data_dir <- dirname(base_yaml_file)[1]
  yml <- yaml::read_yaml(base_yaml_file[1])
  yml2 <- copy(yml)
  yml2$simulation_config$start_date <- new_start_date
  yml2$simulation_config$restore_from <- xfun::relative_path(checkpoint_file,dir = data_dir)
  yml2$population_data$initialization <- NULL
  fp_out <- stringr::str_replace(checkpoint_file, 
                                 pattern = checkpoint_dir, 
                                 replacement = data_dir)
  
  fp_out <- stringr::str_replace(fp_out,pattern = '.RDS',replacement = '.yaml')
  if(file.exists(fp_out)){
    warning(sprintf('Overwriting: %s \n Consider removing unused yaml files before calibration',fp_out))
  }
  yaml::write_yaml(yml2,fp_out)
  #fp_out <- normalizePath(fp_out)
  names(fp_out) <- NULL
  return(fp_out)
}

#' @description
#' Removes (temporary) config files that are not the base yaml file
#' @param base_yaml_file file to not remove
#' @param config_files list of files to remove
#' @export
cleanup_config_files <- function(base_yaml_file,config_files){
  files_to_remove <- unlist(config_files)
  files_to_remove <- files_to_remove[files_to_remove!=base_yaml_file]
  if(length(files_to_remove)>0){
    file.remove(files_to_remove)
  }
}


#' Messaging functions
#' 
#' @description set of functions to control display output
#' @param n size of window
#'@export
msg_header_footer <- function(n=getOption("width")){
  WIDTH_MAX = getOption("width")
  n = min(n,WIDTH_MAX)
  if(n<1){
    n = n * WIDTH_MAX
  }
  n = as.integer(n)
  paste0(rep('*',n),collapse = '')
}


