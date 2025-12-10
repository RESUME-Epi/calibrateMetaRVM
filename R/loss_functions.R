#' Loss functions
#' a collection of common loss functions

#' @description sum of squared differences between model output and ground truth
#' @param sim_output simulation output
#' @param ground_truth ground_truth data
#' @export
squared_log_loss <- function(sim_output,ground_truth){
  d = (sim_output-ground_truth)**2 # squared diff
  d = sum(d)
  d = log(d)
  d
}

#' absolute value of difference between model output and ground truth
#' @param sim_output vector of model outputs
#' @param ground_truth vector of data to compare to model output
#' @export
absolute_loss <- function(sim_output,ground_truth){
  d = abs(sim_out - ground_truth)
  d = sum(d)
  d
}

