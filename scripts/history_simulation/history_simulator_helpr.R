library(expm)

#' Simulate history of discrete-character changes over a branch (which may comprise one or more pieces delimited by time interval bounds) 
#' conditioning on the start and end states
#' @param node_ages ages of the two nodes defining this branch
#' @param Q a instantaneous-rate matrix (or a list of matrices for a piecewise constant geographic model) characterizes the CTMC
#' @param Q_ages boundaries of time intervals for a piecewise constant geographic model (NULL means a constant model)
#' @param node_state_indices state indices (i.e., the indices to use to fecth the states from the states vector) pf the start and end nodes
#' @param states states of the discrete character
#' @return A map object (as in the simmap data structure) that contains the simulated history (as a series of events) over this branch
sim_history_branch_conditional <- function(node_ages, Q, Q_ages = NULL, node_state_indices, states) {
  
  # we first compute P matrix for each piece (defined by time intervals) of the branch
  # then if there are more than one pieces on this branch, we need to compute the convoluted P matrix vector
  # then if there are more than one pieces on this branch, we draw state at each intermediate node
  # then we sample history conditioning on the end state of each piece using uniformization
  
  # fetch the time intervals this branch overlap with
  indices <- 1L
  if (length(Q) > 1) {
    indices <- sapply(node_ages, function(x) findInterval(-x, -Q_ages, left.open = T)) + 1L
    indices <- indices[1]:indices[2]
  }
  
  nstates <- nrow(Q[[indices[1]]])
  
  npieces <- length(indices)
  if (npieces > 1) {
    node_ages <- append(node_ages, Q_ages[indices[-length(indices)]], after = 1)
  }
  t <- abs(diff(node_ages))
  
  # compute P matrix of each piece using matrix exponentiation
  # todo: we can also use uniformization to compute the P matrices, which should be a faster 
  # as we have already done some part of the computation other where in this function
  P <- lapply(1:npieces, function(x) expm::expm(Q[[indices[x]]] * t[x]))
  
  if (npieces > 1) { # when there are more than one intervals overlapping with this branch
    
    # convolved matrix: first element is P along the entire branch
    # second element is P along the branch without the first piece, so on,
    # last element is P along the last two pieces, so number of elements equals npieces - 1
    # for instance, when there are two pieces, this list should only contain one P 
    # as the P along the entire branch (which should = P1 %*% P2 in this case)
    
    # to save computation, we compute the last elements of this list: P_conv[[npieces - 1]] = P[[npieces - 1]] %*% P[[npieces]]
    # then P_conv[[npieces - 2]] = P[[npieces - 2]] %*% P_conv[[npieces - 1]]
    
    P_conv <- vector("list", npieces - 1)
    P_conv[[npieces - 1]] <- P[[npieces - 1]] %*% P[[npieces]]
    
    if (npieces > 2) {
      for (j in (npieces - 2):1) { # loop over the pieces backwards
        P_conv[[j]] <- P[[j]] %*% P_conv[[npieces - 1]]
      }
    }
    
    # loop over the pieces forwards to draw state of each intermediate node 
    # (defined by the overlap of the boundary of time intervals with this branch) 
    # conditioning on the state of start state of the branch (for the first intermediate node) or state of the previous intermediate node
    # and the end state of this branch
    state_idx_end <- node_state_indices[length(node_state_indices)]
    for (j in 1:(npieces - 1)) { 
      
      if (j == npieces - 1) {
        P_complement <- P[[j + 1]]
      } else {
        P_complement <- P_conv[[j + 1]]
      }
      
      state_idx_start <- node_state_indices[j]
      # compute the probability vector of this intermediate node
      state_probs_inter <- numeric(nstates)
      for (k in 1:nstates) {
        # here prob_i,j = P(this piece)[i, k] * P(complement; i.e., along the remaining pieces of the branch)[k, j] / P(along this + remaining pieces)[i, j]
        # but as the denominator is identical to all ks, so we can omit it
        state_probs_inter[k] <- P[[j]][state_idx_start, k] * P_complement[k, state_idx_end]
      }
      # draw a state of this intermediate node according to the probability vector
      state_idx_inter <- sample.int(nstates, 1, prob = state_probs_inter)
      
      node_state_indices <- append(node_state_indices, state_idx_inter, after = j)
    } # end loop to draw state of each intermediate node (now node_state_indices is a vector of state indices 
      # where the ones between the first and last elements correspond to the intermediate nodes)
  }
  
  # loop over the pieces again to simulate history over each piece conditioning on the associated start and end states 
  to_ages <- lapply(1:npieces, function(i) sim_history_piece_conditional(node_ages[c(i, i + 1)], Q, Q_ages, P = P[[i]], node_state_indices[c(i, i + 1)]))
  to <- unlist(lapply(to_ages, function(x) x$to)) # end state of each event over the branch (we don't need to store the corresponding start state as it must be the previous end state)
  age <- unlist(lapply(to_ages, function(x) x$age)) # age of each event over the branch
  
  # format the return as the map object of the simmap data structure
  ages <- c(node_ages[1], age, node_ages[length(node_ages)])
  map <- abs(diff(ages))
  names(map) <- states[c(node_state_indices[1], to)]
  
  return(map)
}

#' Simulate history of discrete-character changes over a duration of time (i.e., a piece of a branch that the character evolved under a constant model)
#' conditioning on the start and end states (implementing algorithm #5 of Hobolth and Stone, 2009)
#' @param node_ages start and end ages of this duration
#' @param Q a instantaneous-rate matrix (or a list of matrices for a piecewise constant geographic model) characterizes the CTMC
#' @param Q_ages boundaries of time intervals for a piecewise constant geographic model (NULL means a constant model)
#' @param P the transition probability matrix computed from the Q matrix (if not provided then it will be computed in this function)
#' @param node_state_indices state indices (i.e., the indices to use to fecth the states from the states vector) pf the start and end nodes
#' @return a list containing two vectors, one for the end state of each event over this duration and other for the age of each event
sim_history_piece_conditional <- function(node_ages, Q, Q_ages = NULL, P = NULL, node_state_indices) {
  
  # first find interval
  idx <- 1L
  if (length(Q) > 1) {
    idx <- findInterval(-node_ages[2], -Q_ages, left.open = T) + 1L
  }
  
  nstates <- nrow(Q[[idx]])
  
  # compute the P matrix using matrix exponentiation (if it's not provided)
  t <- abs(diff(node_ages))
  if (is.null(P)) {
    P <- expm::expm(Q[[idx]] * t)
  }
  
  # draw number of changes conditioning on the start and end states
  state_idx_start <- node_state_indices[1]
  state_idx_end <- node_state_indices[length(node_state_indices)]
  n <- draw_nevents_conditional(node_ages, Q, Q_ages, node_state_indices, P[state_idx_start, state_idx_end])
  
  to <- integer()
  age <- numeric()
  
  # for each event, draw the time it occurs and the state the character changes to
  if (n > 1 || (n == 1 && node_state_indices[1] != node_state_indices[2])) {
    
    # simulate event times
    ts <- sort(runif(n)) * t
    
    # simulate change per event
    if (n > 1) { # for the first to the second but last event (when there are more than one events)
      for (i in 1:(n - 1)) { # loop over each event (but cannot be parallelized as the next start state is the previous end state)
        
        state_idx_start <- node_state_indices[i]
        state_probs_inter <- numeric(nstates)
        for (j in 1:nstates) {
          state_probs_inter[j] <- R[[idx]][[2]][state_idx_start, j] * R[[idx]][[n - i + 1]][j, state_idx_end] # Remark 7 from Hobolth and Stone, 2009
        }
        state_idx_inter <- sample.int(nstates, 1, prob = state_probs_inter)
        
        node_state_indices <- append(node_state_indices, state_idx_inter, after = i)
        if (state_idx_inter != node_state_indices[i]) {
          to <- c(to, state_idx_inter)
          age <- c(age, node_ages[1] - ts[i])
        }
      }
    }
    
    if (node_state_indices[n] != node_state_indices[n + 1]) { # for the last event (if there is any)
      to <- c(to, node_state_indices[n + 1])
      age <- c(age, node_ages[1] - ts[n])
    }
  }
  
  return(list(to = to, age = age))
}

#' Compute the uniformized matrix (R) of Q and its exponentiations (R^2, R^3, ..., R^n)
#' Here the vector of uniformized matrices in stored globally so we check if the provided n (number of events)
#' exceeds the maximum exponent minus 1 (as the first element is an identity matrix)
#' and only append the additional exponentiation(s) when it exceeds
#' @param Q a instantaneous-rate matrix (or a list of matrices for a piecewise constant geographic model) characterizes the CTMC
#' @param n number of events (which is also the largest exponent to compute in this step)
#' @param indices which interval we need to update (if not provided then update all)
fill_DTMCmats <- function(Q, n, indices = NULL) {
  
  nstates <- nrow(Q[[1]])
  R_this <- R
  
  if (is.null(indices)) {
    indices <- seq_along(R)
  }
  
  for (i in indices) {
    if (length(R_this[[i]]) <= n) {
      
      for (j in length(R_this[[i]]):n) {
        if (j == 0) { # the first element of the vector is an identity matrix as R^0 = I
          R_this[[i]] <- list(diag(nrow = nstates, ncol = nstates))
        } else if (j == 1) {
          
          R_tmp <- Q[[i]] / mu[i] + R_this[[i]][[1]]
          R_this[[i]] <- c(R_this[[i]], list(R_tmp))
          
        } else {
          R_this[[i]] <- c(R_this[[i]], list(R_this[[i]][[j]] %*% R_this[[i]][[2]]))
        }
      }
      
    }
  }
  
  # R is a global variable that is a list (each element corresponds to an interval) of 
  # list (each element corresponds to an exponents) of uniformized matrices
  assign("R", R_this, envir = .GlobalEnv)
}

#' Draw the number of events conditioning on the start and end state given the underlying CTMC 
#' (implementing equation 2.9 from Hobolth and Stone, 2009)
#' @param node_ages start and end ages of this duration
#' @param Q a instantaneous-rate matrix (or a list of matrices for a piecewise constant geographic model) characterizes the CTMC
#' @param Q_ages boundaries of time intervals for a piecewise constant geographic model (NULL means a constant model)
#' @param node_state_indices state indices (i.e., the indices to use to fecth the states from the states vector) pf the start and end nodes
#' @param Pab the transition probability between the start and end states
#' @return an integer that is the sampled number of events
draw_nevents_conditional <- function(node_ages, Q, Q_ages = NULL, node_state_indices, Pab) {
  
  # first find the time interval we are at
  idx <- 1L
  if (length(Q) > 1) {
    idx <- findInterval(-node_ages[2], -Q_ages, left.open = T) + 1L
  }
  
  t <- abs(diff(node_ages))
  
  n <- 0L
  u <- runif(1)
  cum_prob <- 0
  
  repeat {
    fill_DTMCmats(Q, n, indices = idx)
    cum_prob <- cum_prob + dpois(n, mu[idx] * t) * R[[idx]][[n + 1]][node_state_indices[1], node_state_indices[2]] / Pab
    if (cum_prob >= u) {
      break
    }
    n <- n + 1L
  }
  
  return(n)
}
