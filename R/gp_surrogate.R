#' Fit GP surrogate models for multiple configs
#'
#' This function takes a matrix of inputs X, a vector of outputs Y, and a vector of categorical indicators cat_ind.
#' It divides Y and X according to cat_ind, and for each category, it fits a GP surrogate
#' using laGP. It then returns a list of the fitted laGP objects.
#'
#' @param X A matrix of input parameters.
#' @param Y A vector of output values.
#' @param cat_ind A vector of categorical indicators.
#' @return A list of fitted `laGP` objects, with names corresponding to the categories.
fit_composite_gp <- function(X, Y, cat_inds, np) {
  categories <- unique(cat_inds)
  gp_list <- list()
  m_sd_mat <- matrix(NA, nrow = length(categories), ncol = 2)
  for (cat in categories) {
    indices <- which(cat_inds == cat)
    X_cat <- X[indices, -ncol(X), drop = FALSE]
    Y_cat <- Y[indices]

    # Fit a GP model for the current category

    ym <- mean(Y_cat)
    ys <- sd(Y_cat)
    ystd <- (Y_cat - ym) / ys

    m_sd_mat[cat, ] <- c(ym, ys)

    d <- rep(1,np)
    gp <- laGP::newGPsep(X_cat, ystd, d, 1e-6, dK = T)
    laGP::mleGPsep(gp)

    gp_list[[cat]] <- gp
  }
  return(list(gp_list = gp_list, m_sd_mat = m_sd_mat))
}
#' predict_composite_gp
#' @description
#' Make prediction for each GP in \code{gp_list} on \code{newX}
#' @param gp_list list of GPs from \code{fit_composite_gp}
#' @param newX matrix of new locations to predict on
#' @returns list of predictions for each GP
predict_composite_gp <- function(gp_list, newX){

  pred_list <- list()
  cat_inds <- newX[, ncol(newX)]
  categories <- unique(cat_inds)
  for (cat in categories) {
    indices <- which(cat_inds == cat)
    newX_cat <- newX[indices, -ncol(newX), drop = FALSE]
    pred_list[[cat]] <- laGP::predGPsep(gp_list[[cat]], newX_cat)
  }
  return(pred_list)

}

#' update_composite_gp
#' @description Update each GP in \code{gp_list} with \code{newX} and \code{newY}
#' @param gp_list list of GPs
#' @param newX matrix of new locations to fit
#' @param newY vector of new locations to fit
#' @param new_cat_ind vector of integers corresponding to configuration files
#' @param m_sd_mat matrix of mean, standard deviations to normalize newX and newY
#' @param updated \code{gp_list}
update_composite_gp <- function(gp_list, newX, newY, new_cat_ind, m_sd_mat) {
  
  categories <- unique(new_cat_ind)
  for (cat in categories) {
    indices <- which(new_cat_ind == cat)
    newX_cat <- newX[indices, -ncol(newX), drop = FALSE]
    newY_cat <- newY[indices]

    ## normalize Y
    ym <- m_sd_mat[cat, 1]
    ys <- m_sd_mat[cat, 2]
    newY_cat_std <- (newY_cat - ym) / ys

    laGP::updateGPsep(gp_list[[cat]], newX_cat, newY_cat_std)
  }

  return(gp_list)
}

#' draw_sample_composite_gp
#' @description
#' Draw sample from list of GPs
#' 
#' @param gp_list list of GPs
#' @param newX matrix of X locations 
#' @param sim_batch_size size of sample
#' @returns matrix of Thompson sampling results
draw_sample_composite_gp <- function(gp_list, newX, sim_batch_size){

  cat_inds <- newX[, ncol(newX)]
  categories <- unique(cat_inds)
  tTS_all <- matrix(NA, nrow = sim_batch_size, ncol = nrow(newX))
  for (cat in categories){
    indices <- which(cat_inds == cat)
    newX_cat <- newX[indices, -ncol(newX), drop = FALSE]
    pred <- laGP::predGPsep(gp_list[[cat]], newX_cat)

    tTS <- mvtnorm::rmvnorm(sim_batch_size, mean = pred$mean, 
                              sigma = (pred$Sigma + t(pred$Sigma)) / 2)
    
    tTS_all[, indices] <- tTS
  }
  return(tTS_all)
}