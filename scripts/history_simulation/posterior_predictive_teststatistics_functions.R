# this script contains functions that compute posterior-predictive summary statistics
# under a constant or piecewise constant geographical model
# specifically, two such statistics are focussed here, including the piecewise constant tipwise multinomial statistic and the piecewise constant parsimony statistic
# it also contains functions that can reconstruct parsimony histories
library(phangorn)
library(stringr)
library(pbapply)

#' Compute the piecewise constant tipwise multinomial statistics for one interval
#' @param nvec a vector of states of all tips
#' @param mvec a vector of states of the tips in a given time interval
#' @return the piecewise constant tipwise multinomial statistic of this interval
multinomiallikelihood_tipwise_calculator <- function(nvec, mvec = NULL) {
  
  nvec <- nvec[nvec != "?"]
  if (is.null(mvec)) { # no interval-specific tips, then all the tips fall in one interval
    return(sum(table(nvec) * log(table(nvec))) - length(nvec) * log(length(nvec)))
  } else {
    if (length(mvec) > 0) {
      mvec <- mvec[mvec != "?"]
      pvec <- table(nvec) / length(nvec)
      nmvec <- table(mvec)
      return(sum(nmvec * log(pvec[names(nmvec)])))
    } else {
      return(0)
    }
  }
}

#' Compute the piecewise constant tipwise multinomial statistics
#' @param tree a bifuracting tree of class phylo
#' @param tipstates the observed tip states
#' @param Q_ages boundaries of time intervals for a piecewise constant geographic model
#' @return a vector (each element correspond to a time interval, from most ancient to the present) of piecewise constant tipwise multinomial statistics
multinomiallikelihood_tipwise_epoch_calculator <- function(tree, tipstates, Q_ages) {
  
  tips <- tree$edge[, 2][!tree$edge[, 2] %in% tree$edge[, 1]]
  node_times <- ape::node.depth.edgelength(tree)
  node_ages <- max(node_times) - node_times
  
  Q_times <- max(node_ages) - sort(Q_ages, decreasing = T)
  epoch_num <- length(Q_times) + 1L
  tipstates <- tipstates[tree$tip.label]
  
  tipepoch_idx <- findInterval(node_times[tips], Q_times) + 1L
  # fetch the tips in each interval and perform the calculation
  return(rev(sapply(1:epoch_num, function(m) multinomiallikelihood_tipwise_calculator(nvec = tipstates, mvec = tipstates[tipepoch_idx == m]))))
}

#' Compute the piecewise constant number of realized states statistics
#' @param tree a bifuracting tree of class phylo
#' @param tipstates the observed tip states
#' @param Q_ages boundaries of time intervals for a piecewise constant geographic model
#' @return a vector (each element correspond to a time interval, from most ancient to the present) of piecewise constant number of realized states statistics 
nstates_epoch_counter <- function(tree, tipstates, Q_ages) {
  
  tips <- tree$edge[, 2][!tree$edge[, 2] %in% tree$edge[, 1]]
  node_times <- ape::node.depth.edgelength(tree)
  node_ages <- max(node_times) - node_times
  
  Q_times <- max(node_ages) - sort(Q_ages, decreasing = T)
  epoch_num <- length(Q_times) + 1L
  tipstates <- tipstates[tree$tip.label]
  
  tipepoch_idx <- findInterval(node_times[tips], Q_times) + 1L
  nstates_epoch <- sapply(1:epoch_num, function(m) {
    tipstates_tmp <- tipstates[tipepoch_idx == m]
    return(length(unique(tipstates_tmp[tipstates_tmp != "?"])))
  })
  nstates_epoch <- rev(nstates_epoch)
  
  return(nstates_epoch)
}

#' Reconstruct a parsimony history of a discrete-character over a bifurcation tree conditioning on the observed tip states and the simulated tip states, respectively
#' @param trees a list of bifuracting trees (each of them is a phylo object with an additional states component for the tip states)
#' @param observed_tipstates a vector of tip states where each element corresponds to a tip; the name attribute contains the tip labels
#' @param states a vector of states of the discrete character
#' @return a tree with two additional components, node.states_parsimony_observed and node.states_parsimony_simulated, attached
sim_par <- function(tree, observed_tipstates, states) { # only work for a single site/character
  
  # parsimony history for the observed tip states
  pr_oberserved <- phangorn::ancestral.pars(tree = tree, 
                                            data = phyDat(t(t(observed_tipstates)), type = "USER", levels = states), 
                                            type = "MPR", return = "phyDat")
  
  # parsimony history for the simulated tip states
  simulated_tipstates <- tree$states
  simulated_tipstates[names(observed_tipstates)[observed_tipstates == "?"]] <- "?"
  pr_simulated <- ancestral.pars(tree = tree, 
                                 data = phyDat(t(t(simulated_tipstates)), type = "USER", levels = states), 
                                 type = "MPR", return = "phyDat")
  
  tree$node.states_parsimony_observed <- matrix(states[as.integer(pr_oberserved[tree$edge])], ncol = 2)
  tree$node.states_parsimony_simulated <- matrix(states[as.integer(pr_simulated[tree$edge])], ncol = 2)
  
  return(tree)
}

#' Reconstruct some number of parsimony histories of a discrete-character over a list of bifurcation trees conditioning on the observed tip states and the simulated tip states, respectively
#' @param trees a list of bifuracting trees (each of them is a phylo object with an additional states component for the tip states)
#' @param observed_tipstates a vector of tip states where each element corresponds to a tip; the name attribute contains the tip labels
#' @param states a vector of states of the discrete character
#' @param indices generation indices that will only used for printing progress to screen
#' @return a list of trees where each has two additional components, node.states_parsimony_observed and node.states_parsimony_simulated, attached
sim_pars <- function(trees, observed_tipstates, states, indices = NULL) {
  
  ntrees <- length(trees)
  verbose <- !is.null(indices)
  for (l in 1:ntrees) {
    if (verbose && l %% 50 == 0) {
      cat(paste0("parsimony reconstruction no. ", indices[l], ".\n"))
    }
    if (!any(is.na(trees[[l]]))) {
      trees[[l]] <- sim_par(tree = trees[[l]], observed_tipstates, states)
    }
  }
  
  return(trees)
}

#' Reconstruct parsimony histories of a discrete-character over a list of bifurcation trees conditioning on the observed tip states and the simulated tip states, respectively
#' @param trees a list of bifuracting trees (each of them is a phylo object with an additional states component for the tip states)
#' @param observed_tipstates a vector of tip states where each element corresponds to a tip; the name attribute contains the tip labels
#' @param states a vector of states of the discrete character
#' @param ncores number of computer cores to use (if more than one then simulations will be parallelized)
#' @return a list of trees where each has two additional components, node.states_parsimony_observed and node.states_parsimony_simulated, attached
simulate_pars <- function(trees, observed_tipstates, states, ncores = 1L) {
  
  ntrees <- length(trees)
  if (ncores > 1L && ntrees > 50L) { # multi-core parallelization
    
    nchunks <- ncores * 4L
    if (ntrees < nchunks) {
      nchunks <- ntrees
    }
    nhis_percore <- floor(ntrees / nchunks)
    his_indices <- vector("list", nchunks)
    k <- 0L
    for (j in 1:nchunks) {
      if (j < nchunks) {
        his_indices[[j]] <- 1L:nhis_percore + k
      } else {
        his_indices[[j]] <- (k + 1L):ntrees
      }
      k <- k + nhis_percore
    }
    
    trees_all <- pblapply(his_indices, function(his_idx) sim_pars(trees = trees[his_idx], observed_tipstates, states, indices = his_idx), cl = ncores)
    trees <- do.call(c, trees_all)
  } else { # single core
    trees <- sim_pars(trees, observed_tipstates, states, indices = seq_along(trees))
  }
  
  return(trees)
}

#' Compute the piecewise constant parsimony statistics
#' @param tree a bifuracting tree of class phylo
#' @param epoch_bound_ages boundaries of time intervals for a piecewise constant geographic model
#' @param node.states identical to the edge matrix of a phylo object except here each element is the node state instead of node index
#' @return a vector (each element correspond to a time interval, from most ancient to the present) of piecewise constant parsimony statistics
parsimonyscore_epoch_calculator <- function(tree, epoch_bound_ages, node.states) {

  epoch_num <- length(epoch_bound_ages) - 1L
  
  nedges <- nrow(tree$edge)
  node_times <- ape::node.depth.edgelength(tree)
  node_ages <- max(node_times) - node_times
  
  pieceduration_mat <- matrix(0, ncol = epoch_num, nrow = nedges, byrow = T)
  for (k in 1:nedges) { # loop over branches to get the interval mapping matrix
    
    node_anc <- tree$edge[k, 1]
    node_dec <- tree$edge[k, 2]
    age_anc <- node_ages[node_anc]
    age_dec <- node_ages[node_dec]
    
    for (m in 1:epoch_num) {
      piece_duration <- min(age_anc, epoch_bound_ages[m]) - max(age_dec, epoch_bound_ages[m + 1])
      if (piece_duration > 0) {
        pieceduration_mat[k, m] <- piece_duration
      }
    }
  }
  
  nchanges_epoch <- integer(epoch_num)
  for (k in 1:nedges) { # loop over branches again to count the number of events
    if (length(unique(node.states[k, ])) > 1 && all(node.states[k, ] != "?")) { # when the start and end states are different (therefore a change must had occurred)
      pieces_duration <- pieceduration_mat[k, ]
      
      if (sum(pieces_duration > 0) > 1) {
        epoch_id <- sample.int(epoch_num, size = 1, prob = pieces_duration)
        nchanges_epoch[epoch_id] <- nchanges_epoch[epoch_id] + 1L
      } else if (sum(pieces_duration > 0) == 1) {
        nchanges_epoch[pieces_duration > 0] <- nchanges_epoch[pieces_duration > 0] + 1L
      } else {
        stop("branch cannot be zero length when there is a change on it.\n")
      }
    }
  }
  
  return(rev(nchanges_epoch))
}


#' Compute the piecewise constant pairwise parsimony statistic for one interval
#' @param tree a bifuracting tree of class phylo
#' @param node.states identical to the edge matrix of a phylo object except here each element is the node state instead of node index
#' @param states a vector of states of the discrete character
#' @return the piecewise constant pairwise parsimony statistic of this interval
parsimonyscore_pairwise_calculator <- function(tree, node.states, states) {
  
  nstates <- length(states)
  nchanges_pairwise <- matrix(0, nrow = nstates, ncol = nstates, byrow = T)
  rownames(nchanges_pairwise) <- states
  colnames(nchanges_pairwise) <- states
  
  node.states <- node.states[node.states[, 1] != node.states[, 2] & node.states[, 1] != "?" & node.states[, 2] != "?", , drop = F]
  
  if (nrow(node.states > 0)) {
    # count the number of events between each pair
    statepairs_tab <- table(apply(node.states, 1, function(x) paste(x, collapse = ", ")))
    statepairs_nchanges <- as.integer(statepairs_tab)
    statepairs_mat <- matrix(unlist(strsplit(names(statepairs_tab), ", ")), ncol = 2, byrow = T)
    
    for (k in 1:nrow(statepairs_mat)) {
      nchanges_pairwise[statepairs_mat[k, 1], statepairs_mat[k, 2]] <- statepairs_nchanges[k]
    }
  }

  return(nchanges_pairwise)
}


#' Compute the piecewise constant pairwise parsimony statistics
#' @param tree a bifuracting tree of class phylo
#' @param epoch_bound_ages boundaries of time intervals for a piecewise constant geographic model
#' @param node.states identical to the edge matrix of a phylo object except here each element is the node state instead of node index
#' @param states a vector of states of the discrete character
#' @return a vector (each element correspond to a time interval, from most ancient to the present) of piecewise constant pairwise parsimony statistics
parsimonyscore_pairwise_epoch_calculator <- function(tree, epoch_bound_ages, node.states, states) {
  
  nstates <- length(states)
  epoch_num <- length(epoch_bound_ages) - 1L
  
  nedges <- nrow(tree$edge)
  node_times <- ape::node.depth.edgelength(tree)
  node_ages <- max(node_times) - node_times
  
  pieceduration_mat <- matrix(0, ncol = epoch_num, nrow = nedges, byrow = T)
  for (k in 1:nedges) { # loop over branches to get the interval mapping matrix
    
    node_anc <- tree$edge[k, 1]
    node_dec <- tree$edge[k, 2]
    age_anc <- node_ages[node_anc]
    age_dec <- node_ages[node_dec]
    
    for (m in 1:epoch_num) {
      piece_duration <- min(age_anc, epoch_bound_ages[m]) - max(age_dec, epoch_bound_ages[m + 1])
      if (piece_duration > 0) {
        pieceduration_mat[k, m] <- piece_duration
      }
    }
  }
  
  nchanges_pairwise_epoch <- vector("list", epoch_num)
  for (m in 1:epoch_num) {
    nchanges_pairwise_epoch[[m]] <- matrix(0, nrow = nstates, ncol = nstates, byrow = T)
    rownames(nchanges_pairwise_epoch[[m]]) <- states
    colnames(nchanges_pairwise_epoch[[m]]) <- states
  }

  for (k in 1:nedges) { # loop over branches again to count the number of events between each pair
    if (length(unique(node.states[k, ])) > 1 && all(node.states[k, ] != "?")) { # when the start and end states are different (therefore a change must had occurred)
      
      pieces_duration <- pieceduration_mat[k, ]
      epoch_id <- which(pieces_duration > 0)
      
      if (length(epoch_id) > 1) {
        epoch_id <- sample.int(epoch_num, size = 1, prob = pieces_duration)
      } else if (length(epoch_id) == 0) {
        stop("branch cannot be zero length when there is a change on it.\n")
      }
      
      nchanges_pairwise_epoch[[epoch_id]][node.states[k, 1], node.states[k, 2]] <- nchanges_pairwise_epoch[[epoch_id]][node.states[k, 1], node.states[k, 2]] + 1L
    }
  }
  
  return(rev(nchanges_pairwise_epoch))
}
