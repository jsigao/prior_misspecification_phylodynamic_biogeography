library(stringr)

#' Compute the least common multiple of two integers
#' @param x A integer
#' @param y A integer
#' @return LCM of \code{x} and \code{y}.
lcm <- function(x, y) {
  # choose the greater number
  if (x > y) {
    greater <- x
    increment <- x
  } else {
    greater <- y
    increment <- y
  }
  while (T) {
    if ((greater %% x == 0) && (greater %% y == 0)) {
      lcm <- greater
      break
    }
    greater <- greater + increment
  }
  return(lcm)
}

#' Compute the stationary frequency vector of a given irreducible Q matrix
#' @param Q A irreducible square matrix that characterizes a CTMC
#' @return stationary frequency vector of \code{Q}
get_stationary_freq <- function(Q) {
  
  # first perform spectral decomposition
  Q_eigen <- eigen(Q)
  if (is.complex(Q_eigen$values)) {
    eig_values <- Re(Q_eigen$values)
  } else {
    eig_values <- Q_eigen$values
  }
  
  # find the dominant/largest eigenvalue (which should theoretically be zero)
  eig_value1 <- max(eig_values)
  if (abs(eig_value1) > 1e-6) {
    message("cannot find an eigen value that's (approximately) zero; likely reducible Q matrix")
    return(NULL)
  }
  eig1_idx <- which(eig_values == eig_value1)
  
  # the left eigen vector associated with the zero eigenvalue is a rescale of the stationary frequency vector
  eig_vectors_left <- solve(Q_eigen$vectors)
  if (is.complex(eig_vectors_left)) {
    stat_freq <- Re(eig_vectors_left[eig1_idx, ])
  } else {
    stat_freq <- eig_vectors_left[eig1_idx, ]
  }
  stat_freq <- stat_freq / sum(stat_freq)
  
  if (any(stat_freq < -1e-6)) {
    message("some stationary frequency element(s) is(are) negative")
    return(NULL)
  }
  
  return(stat_freq)
}

#' Extract and store geographic model parameter estimates from a BEAST log-file output
#' @param log_path path to the BEAST log file
#' @param burnin fraction of or absolute number of generations of the log file to be discarded as burnin
#' @return a list (each one correspond to an interval defined by the Q matrix; lacking of this layer if the model is constant)
#' of list of parameter summaries
get_BF_counts <- function(log_path, burnin = 0) {
  
  log_dat <- read.table(log_path, header = T, sep = "\t", check.names = F, stringsAsFactors = F)
  
  # Remove the rows where flatlines are
  if (any(log_dat$posterior == Inf)) {
    log_dat <- log_dat[-(which(log_dat$posterior == Inf)),]
  }
  
  # get the number of states
  # here we assume that the xml file is in the same directory with the log file and have identical name (except the filename extension)
  xml_path <- grep(gsub("_combined\\.log$|\\.log$", "", basename(log_path)), 
                   grep("MLE", list.files(dirname(log_path), recursive = T, full.names = T, pattern = "*.xml$"), value = T, invert = T), value = T)
  if (length(xml_path) == 0) {
    stop("cannot find xml")
  } else if (length(xml_path) > 1 && length(unique(gsub("_run(.*?)\\.xml$", "", basename(xml_path)))) != 1) {
    stop("cannot find xml")
  } else {
    xml_path <- xml_path[1]
  }
  x <- scan(file = xml_path, what = character(), sep = "\n", strip.white = F, blank.lines.skip = F)
  state.names <- gsub("^\\s+|\\s+$", "", gsub("<state code|=|\t|\"|/>", "", x[grep("<state code", x)]))
  state.num <- length(state.names)
  K <- state.num
  
  # compute burnin
  log_nrow <- nrow(log_dat)
  log_freq <- diff(log_dat[1:2, 1])
  burnin_nrow <- 0
  if (burnin > 0 && log_nrow > 1) {
    if (burnin > 1) {
      burnin_nrow <- ceiling(burnin / log_freq)
    } else {
      burnin_nrow <- ceiling(log_nrow * burnin)
    }
    if (burnin_nrow >= log_nrow) burnin_nrow <- log_nrow - 1
    
    log_dat <- log_dat[-c(1:burnin_nrow), ]
    log_nrow <- nrow(log_dat)
    gc()
  }
  burnin_ngens <- ifelse(burnin_nrow == 0, 0, log_dat[1, 1])
  
  # read in the stochastic mapping log file (if exist)
  # this file only won't exist when the "fast stochastic mapping" algorithm is performed,
  # which compute the expected number of events on each branch (instead of simulating the entire history)
  # again, we assume it's named identical to the log file
  txt_nrow <- 0
  txt_path <- gsub("\\.log$", ".txt", log_path)
  if (file.exists(txt_path)) {
    history_log <- read.table(txt_path, header = T, sep = "\t", check.names = F, stringsAsFactors = F)
    
    # assumed the column names obey BEAST convention
    history_colnum <- which(colnames(history_log) == "completeHistory_1")
    if (length(history_colnum) == 0) {
      message("cannot find history column.\n")
      rm(history_log)
      gc()
    } else {
      
      # when there are multiple history columns (therefore multiple sites/characters)
      # we assume the first one is the one we want
      if (length(history_colnum) > 1) { 
        bool_vec <- rep(T, length(history_colnum))
        for (j in 1:length(history_colnum)) {
          if (history_colnum[j] < ncol(history_log) && colnames(history_log)[history_colnum[j] + 1] == "completeHistory_2") {
            bool_vec[j] <- F
          }
        }
        history_colnum <- history_colnum[bool_vec]
      }
      
      if (length(history_colnum) != 1) {
        message("zero or more than one single-site trait columns.\n")
        rm(history_log)
        gc()
      } else {
        
        # find the shared generations between the history log file and the parameter log file
        # and then subsample them to leave the shared generations (if they were sampled at different frequency)
        txt_freq <- diff(history_log[1:2, 1])
        sample_freq <- lcm(txt_freq, log_freq)
        
        if (burnin_ngens > 0) {
          if ((burnin_ngens / txt_freq) %% 1 != 0) {
            burnin_nrow <- ceiling(burnin_ngens / sample_freq)
            burnin_ngens <- burnin_nrow * sample_freq
            
            log_dat <- log_dat[log_dat[, 1] >= burnin_ngens, ]
          }
          history_log <- history_log[history_log[, 1] >= burnin_ngens, ]
        }
        
        log_dat <- log_dat[log_dat[, 1] %% sample_freq == 0, ]
        log_nrow <- nrow(log_dat)
        
        history_log <- history_log[history_log[, 1] %in% log_dat[, 1], ]
        gc()
        
        txt_nrow <- nrow(history_log)
        if (txt_nrow != log_nrow) {
          stop("txt and log rows should match.\n")
        }
      } # end of inner ifelse
    } # end of outer ifelse
  }
  
  # parse the history log as each row contains a collection of events
  # and each event is defined by its age, from and to states
  if (txt_nrow > 0) {
    history_df_all <- vector("list", txt_nrow)
    history_dummy <- data.frame(age = numeric(), from = character(), to = character(), stringsAsFactors = F)
    for (l in 1:txt_nrow) {
      history <- unlist(strsplit(history_log[l, history_colnum], "\\{\\{|\\}\\,\\{|\\}\\}"))
      if (length(history) >= 3) {
        history <- history[-c(1, length(history))]
        history <- data.frame(matrix(unlist(strsplit(history, ",")), ncol = 4, byrow = T), stringsAsFactors = F)
        history <- history[, -1]
        
        colnames(history) <- c("age", "from", "to")
        history$age <- as.numeric(history$age)
        history_df_all[[l]] <- history
      } else {
        history_df_all[[l]] <- history_dummy
      }
    }
  }
  
  # fetch the interval bounds of the average dispersal rate and the Q matrix, respectively (as they can be different)
  mu_bounds_age_linenums <- grep("epoch transitionTime", x)[grep("epoch transitionTime", x) > grep("rateEpochBranchRates id=", x) & 
                                                             grep("epoch transitionTime", x) < grep("</rateEpochBranchRates>", x)]
  Q_bounds_age_linenums <- grep("epoch transitionTime", x)[grep("epoch transitionTime", x) > grep("epochBranchModel id=", x) & 
                                                           grep("epoch transitionTime", x) < grep("</epochBranchModel>", x)]
  mu_epoch_num <- length(mu_bounds_age_linenums) + 1L
  Q_epoch_num <- length(Q_bounds_age_linenums) + 1L

  mu_bounds_age <- NULL
  if (mu_epoch_num > 1) {
    mu_bounds_age <- as.numeric(str_extract(x[mu_bounds_age_linenums], "\\d+"))
  }
  Q_bounds_age <- NULL
  if (Q_epoch_num > 1) {
    Q_bounds_age <- as.numeric(str_extract(x[Q_bounds_age_linenums], "\\d+"))
  }
  
  # loop over the intervals to fetch the geographic parameter estimates for each interval
  BFs_counts_epoch <- vector("list", Q_epoch_num)
  for (i in 1:Q_epoch_num) {
    
    counts_mean_mat <- NULL
    indicator_mat_all <- NULL
    
    # rates
    if (Q_epoch_num == 1) {
      rates_colnums <- grep("*\\.rates\\d*$|*\\.rates\\.*", colnames(log_dat))
    } else {
      rates_colnums <- grep(paste0("*\\.rates.epoch", i, "\\d*$|*\\.rates.epoch", i, "\\.*"), colnames(log_dat))
    }
    
    # set default to symmetric model
    # determine whether it's symmetrical model or asymmetrical model
    rates_sym <- T
    if (choose(K, 2) * 2 == length(rates_colnums)) {
      rates_sym <- F
    } else if (choose(K, 2) != length(rates_colnums)) {
      stop("The number of states extracted from the xml file does not match its counterpart computed from the log file.\n")
    }
    
    # fill out rate matrices
    rates_startcolnum <- min(rates_colnums)
    rate_mat_all <- vector("list", log_nrow)
    for (l in 1:log_nrow) {
      rate_mat <- matrix(0, nrow = K, ncol = K)
      rate_mat[lower.tri(rate_mat)] <- as.numeric(log_dat[l, 1:(choose(K, 2)) + rates_startcolnum - 1])
      rate_mat <- t(rate_mat)
      if (rates_sym) {
        rate_mat[lower.tri(rate_mat)] <- as.numeric(log_dat[l, 1:(choose(K, 2)) + rates_startcolnum - 1])
      } else {
        rate_mat[lower.tri(rate_mat)] <- as.numeric(log_dat[l, (choose(K, 2) + 1):(choose(K, 2) * 2) + rates_startcolnum - 1])
      }
      rate_mat_all[[l]] <- rate_mat
    }
    
    # indicators
    if (Q_epoch_num == 1) {
      indicators_colnums <- grep("*\\.indicators\\d*$|*\\.indicators\\.*", colnames(log_dat))
    } else {
      indicators_colnums <- grep(paste0("*\\.indicators.epoch", i, "\\d*$|*\\.indicators.epoch", i, "\\.*"), colnames(log_dat))
    }
    
    # fill out indicator matrices
    if (length(indicators_colnums) > 0) {
      # set default to symmetric model
      # determine whether it's symmetrical model or asymmetrical model
      indicators_sym <- T
      if (choose(K, 2) * 2 == length(indicators_colnums)) {
        indicators_sym <- F
      } else if (choose(K, 2) != length(indicators_colnums)) {
        stop("The number of states extracted from the xml file does not match its counterpart computed from the log file.\n")
      }
      
      indicators_startcolnum <- min(indicators_colnums)
      indicator_mat_all <- vector("list", log_nrow)
      
      for (l in 1:log_nrow) {
        indicator_mat <- matrix(0, nrow = K, ncol = K)
        indicator_mat[lower.tri(indicator_mat)] <- as.numeric(log_dat[l, 1:(choose(K, 2)) + indicators_startcolnum - 1])
        indicator_mat <- t(indicator_mat)
        if (indicators_sym) {
          indicator_mat[lower.tri(indicator_mat)] <- as.numeric(log_dat[l, 1:(choose(K, 2)) + indicators_startcolnum - 1])
        } else {
          indicator_mat[lower.tri(indicator_mat)] <- as.numeric(log_dat[l, (choose(K, 2) + 1): (choose(K, 2) * 2) + indicators_startcolnum - 1])
        }
        indicator_mat_all[[l]] <- indicator_mat
      }
    }
    
    # calculate the posterior estimates of indicators
    if (indicators_sym) {
      posterior <- numeric(choose(K, 2))
      for (j in 1:choose(K, 2)) {
        posterior[j] <- mean(log_dat[, j + indicators_startcolnum - 1])
      }
    } else {
      posterior <- numeric(choose(K, 2) * 2)
      for (j in 1:(choose(K, 2) * 2)) {
        posterior[j] <- mean(log_dat[, j + indicators_startcolnum - 1])
      }
    }
    
    posterior_mat <- matrix(0, nrow = K, ncol = K)
    if (indicators_sym) {
      posterior_mat[lower.tri(posterior_mat)] <- posterior
      posterior_mat <- t(posterior_mat)
      posterior_mat[lower.tri(posterior_mat)] <- posterior
    } else {
      posterior_mat[lower.tri(posterior_mat)] <- posterior[1:choose(K, 2)]
      posterior_mat <- t(posterior_mat)
      posterior_mat[lower.tri(posterior_mat)] <- posterior[(choose(K, 2) + 1):(choose(K, 2) * 2)]
    }
    
    pois_str <- grep("<poissonPrior", x, value = T)[i] # here we assumes each epoch has its independent prior on mu
    if (length(pois_str) == 0) {
      lambda <- choose(K, 2) * 0.5 * (2 - indicators_sym)
    } else {
      lambda <- sum(as.numeric(gsub("mean", "", unlist(strsplit(gsub(" |<poissonPrior|\t|=|\"|>", "", pois_str), split = "offset")))))
    }
    
    # calculate the so-called analytical prior
    analytical_prior <- lambda / (choose(K, 2) * (2 - indicators_sym))
    
    # calculate the Bayes factors
    BF_original <- (posterior/(1 - posterior))/(analytical_prior/(1 - analytical_prior))
    if (any(BF_original == -Inf)) {
      BF_original[BF_original == -Inf] <- Inf
    }
    
    BF_mat <- matrix(0, nrow = K, ncol = K)
    if (indicators_sym) { # symmetric
      BF_mat[lower.tri(BF_mat)] <- BF_original
      BF_mat <- t(BF_mat)
      BF_mat[lower.tri(BF_mat)] <- BF_original
    } else { # asymmetric
      BF_mat[lower.tri(BF_mat)] <- BF_original[1:choose(K, 2)]
      BF_mat <- t(BF_mat)
      BF_mat[lower.tri(BF_mat)] <- BF_original[(choose(K, 2) + 1):(choose(K, 2) * 2)]
    }
    
    # counts
    if (txt_nrow > 0) { # when we have the simulated complete history
      if (Q_epoch_num > 1) {
        history_df <- vector("list", txt_nrow)
        for (l in 1:txt_nrow) {
          if (i == 1) {
            history_df[[l]] <- history_df_all[[l]][history_df_all[[l]]$age <= Q_bounds_age[i], ]
          } else if (i < Q_epoch_num) {
            history_df[[l]] <- history_df_all[[l]][history_df_all[[l]]$age <= Q_bounds_age[i] & history_df_all[[l]]$age > Q_bounds_age[i - 1], ]
          } else {
            history_df[[l]] <- history_df_all[[l]][history_df_all[[l]]$age > Q_bounds_age[i - 1], ]
          }
        }
      } else {
        history_df <- history_df_all
        rm(history_df_all)
        gc()
      }
      
      counts_all <- matrix(0, nrow = txt_nrow, ncol = K * K)
      for (l in 1:txt_nrow) {
        for (j in 1:K) {
          for (k in 1:K) {
            if (j != k) {
              counts_all[l, (j - 1) * K + k] <- sum((history_df[[l]]$from == state.names[j]) & (history_df[[l]]$to == state.names[k]))
            }
          }
        }
      }
      
      counts_mean_vec <- apply(counts_all, 2, mean)
      counts_HPD_lower_vec <- apply(counts_all, 2, function(x) quantile(x, probs = c(0.025, 0.975))[1])
      counts_HPD_upper_vec <- apply(counts_all, 2, function(x) quantile(x, probs = c(0.025, 0.975))[2])
      
      counts_mean_mat <- matrix(counts_mean_vec, nrow = K, ncol = K, byrow = T)
      counts_HPD_lower_mat <- matrix(counts_HPD_lower_vec, nrow = K, ncol = K, byrow = T)
      counts_HPD_upper_mat <- matrix(counts_HPD_upper_vec, nrow = K, ncol = K, byrow = T)
      
    } else if (length(grep("*.count\\[", colnames(log_dat))) == choose(K, 2) * 2 ||
               length(grep("*.count\\[", colnames(log_dat))) == choose(K, 2) * 2 + 1) { 
      # when we only have the computed expected number of changes
      # todo: probably doesn't work for epochal model
      
      counts_colnums <- grep("c_(.*?)-(.*?)\\.count\\[", colnames(log_dat))
      if (length(counts_colnums) != choose(K, 2) * 2) {
        stop("number of dispersal events columns doesn't match the number of states.\n")
      }
      
      counts_startcolnum <- min(counts_colnums)
      colnum <- counts_startcolnum
      counts_mean_vec <- numeric(K * K)
      counts_HPD_lower_vec <- numeric(K * K)
      counts_HPD_upper_vec <- numeric(K * K)
      
      counts_all <- matrix(0, nrow = log_nrow, ncol = K * K)
      for (counts_num in 1:(K * K)) {
        if ((counts_num - 1) %% (K + 1) != 0) { # exclude diagonal elements
          counts_all[, counts_num] <- log_dat[, colnum]
          counts_mean_vec[counts_num] <- mean(log_dat[, colnum])
          counts_HPD_lower_vec[counts_num] <- quantile(log_dat[, colnum], probs = c(0.025, 0.975))[1]
          counts_HPD_upper_vec[counts_num] <- quantile(log_dat[, colnum], probs = c(0.025, 0.975))[2]
          colnum <- colnum + 1
        }
      }
      
      counts_mean_mat <- matrix(counts_mean_vec, nrow = K, ncol = K, byrow = T)
      counts_HPD_lower_mat <- matrix(counts_HPD_lower_vec, nrow = K, ncol = K, byrow = T)
      counts_HPD_upper_mat <- matrix(counts_HPD_upper_vec, nrow = K, ncol = K, byrow = T)
    }
    
    # mu and qij
    if (Q_epoch_num == 1) {
      mu_colnum <- grep("*\\.clock\\.rate.epoch\\d{1,3}$|*\\.clock\\.rate$", colnames(log_dat))
    } else if (mu_epoch_num == 1) {
      mu_colnum <- grep("*\\.clock\\.rate$", colnames(log_dat))
    } else {
      if (i == 1) {
        mu_bounds_age_this <- mu_bounds_age[mu_bounds_age <= Q_bounds_age[i]]
      } else if (i > 1 && i < Q_epoch_num) {
        mu_bounds_age_this <- mu_bounds_age[mu_bounds_age <= Q_bounds_age[i] & mu_bounds_age > Q_bounds_age[i - 1]]
      } else {
        mu_bounds_age_this <- mu_bounds_age[mu_bounds_age > Q_bounds_age[i - 1]]
      }
      mu_colnum <- grep("*\\.clock\\.rate.epoch\\d{1,3}$", colnames(log_dat))[match(mu_bounds_age_this, mu_bounds_age)]
      if (length(mu_bounds_age_this) == 0 || max(mu_bounds_age_this) == max(mu_bounds_age)) {
        mu_colnum <- c(mu_colnum, max(grep("*\\.clock\\.rate.epoch\\d{1,3}$", colnames(log_dat))))
      }
    }
    
    mu_colnum_all <- grep("*\\.clock\\.rate.epoch\\d{1,3}$|*\\.clock\\.rate$", colnames(log_dat))
    mu_idx <- match(mu_colnum, mu_colnum_all)
    mu_epoch_num_this <- length(mu_idx)
    overall_dispersal_rate_all <- log_dat[, mu_colnum, drop = F]
    
    # calculate qij matrices
    qij_mat_all <- vector("list", log_nrow)
    qijabs_mat_all <- vector("list", mu_epoch_num_this)
    mu_all <- vector("list", mu_epoch_num_this)

    for (j in 1:mu_epoch_num_this) {
      qijabs_mat_all[[j]] <- vector("list", log_nrow)
      mu_all[[j]] <- numeric(log_nrow)
    }
    
    for (l in 1:log_nrow) {
      qij_mat_ori <- rate_mat_all[[l]]
      if (length(indicators_colnums) > 0) {
        qij_mat_ori <- qij_mat_ori * indicator_mat_all[[l]]
      }
      
      # beast way of rescaling
      mu <- sum(apply(qij_mat_ori, 1, sum) / nrow(qij_mat_ori))
      qij_mat_ori <- qij_mat_ori / mu

      for (j in 1:mu_epoch_num_this) {
        qijabs_mat <- qij_mat_ori * overall_dispersal_rate_all[l, j]
        diag(qijabs_mat) <- -apply(qijabs_mat, 1, sum)
        
        # correct way of rescaling
        if (j == 1) {
          stat_freq <- tryCatch(get_stationary_freq(qijabs_mat), error = function(e) {})
        }
        
        if (is.null(stat_freq)) {
          mu <- NaN
        } else {
          mu <- -sum(stat_freq * diag(qijabs_mat))
        }
        
        if (j == 1) {
          if (is.null(stat_freq)) {
            qij_mat_all[[l]] <- NA
          } else {
            qij_mat_all[[l]] <- qijabs_mat / mu
          }
        }
        
        qijabs_mat_all[[j]][[l]] <- qijabs_mat
        mu_all[[j]][l] <- mu
      }
    }
    
    # when there is only one interval for the average rate, then simpify the data structures
    if (mu_epoch_num_this == 1) {
      qijabs_mat_all <- qijabs_mat_all[[1]]
      mu_all <- mu_all[[1]]
      overall_dispersal_rate_all <- overall_dispersal_rate_all[, 1]
    }
    
    # generated the list containing all geographic model parameter summaries
    if (is.null(counts_mean_mat) && is.null(indicator_mat_all)) {
      BFs_counts_epoch[[i]] <- list(state.names = state.names, state.num = state.num, symmetry = as.logical(indicators_sym),
                                    overall_dispersal_rate_all = overall_dispersal_rate_all, rate_mat_all = rate_mat_all,
                                    qijabs_mat_all = qijabs_mat_all, qij_mat_all = qij_mat_all, mu_all = mu_all)
    } else if (is.null(counts_mean_mat)) {
      BFs_counts_epoch[[i]] <- list(state.names = state.names, state.num = state.num, symmetry = as.logical(indicators_sym),
                                    overall_dispersal_rate_all = overall_dispersal_rate_all, rate_mat_all = rate_mat_all,
                                    indicator_mat_all = indicator_mat_all, posterior_mat = posterior_mat, analytical_prior = analytical_prior, BF_mat = BF_mat,
                                    qijabs_mat_all = qijabs_mat_all, qij_mat_all = qij_mat_all, mu_all = mu_all)
    } else if (is.null(indicator_mat_all)) {
      BFs_counts_epoch[[i]] <- list(state.names = state.names, state.num = state.num, symmetry = as.logical(indicators_sym),
                                    overall_dispersal_rate_all = overall_dispersal_rate_all, rate_mat_all = rate_mat_all,
                                    counts_mean_mat = counts_mean_mat, counts_HPD_lower_mat = counts_HPD_lower_mat, 
                                    counts_HPD_upper_mat = counts_HPD_upper_mat, counts_all = counts_all,
                                    qijabs_mat_all = qijabs_mat_all, qij_mat_all = qij_mat_all, mu_all = mu_all)
    } else {
      BFs_counts_epoch[[i]] <- list(state.names = state.names, state.num = state.num, symmetry = as.logical(indicators_sym),
                                    overall_dispersal_rate_all = overall_dispersal_rate_all, rate_mat_all = rate_mat_all,
                                    counts_mean_mat = counts_mean_mat, counts_HPD_lower_mat = counts_HPD_lower_mat, 
                                    counts_HPD_upper_mat = counts_HPD_upper_mat, counts_all = counts_all,
                                    indicator_mat_all = indicator_mat_all, posterior_mat = posterior_mat, analytical_prior = analytical_prior, BF_mat = BF_mat,
                                    qijabs_mat_all = qijabs_mat_all, qij_mat_all = qij_mat_all, mu_all = mu_all)
    }
  }
  
  # when there is only one interval for the Q matrix, then simpify the data structures
  if (Q_epoch_num == 1) {
    BFs_counts_epoch <- BFs_counts_epoch[[1]]
  }

  return(BFs_counts_epoch)
}
