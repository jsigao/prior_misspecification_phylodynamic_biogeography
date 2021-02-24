library(ape)
library(stringr)
source('./history_simulator_helpr.R')

#' Simulate history of discrete-character changes over a bifurcation tree conditional or unconditional on the tip states
#' @param tree a bifurcation tree of class "phylo"
#' @param Q a instantaneous-rate matrix (or a list of matrices for a piecewise constant geographic model) characterizes the CTMC
#' @param Q_ages boundaries of time intervals for a piecewise constant geographic model (NULL means a constant model)
#' @param root_freq the vector state frequencies at the root of the tree (that we will draw the root state from); if NULL then it's uniform
#' @param nsim number of simulations to perform
#' @param conditional whether condition on the observed states at the tip (i.e., stochastic mapping) or not (i.e., forward simulation)
#' @return A phylo and simmap object (or a multiPhylo and multiSimmap object when nsim > 1) that contains the simulated full history
sim_history <- function(tree, Q, Q_ages = NULL, root_freq = NULL, nsim = 1L, conditional = F, trait = NULL) {
  
  # first some sanity checks
  if (!inherits(tree, "phylo")) {
    stop("tree should be an object of class \"phylo\".")
  }
  
  # config Q matrices and the interval time bounds
  if (is.list(Q)) {
    if (is.null(Q_ages) && length(Q) > 1) {
      stop("multiple matrices but no information about how to arrange them chronologically.\n")
    } else if (!is.null(Q_ages)) {
      
      if (all(Q_ages <= 0)) {
        Q_ages <- -Q_ages
      } else if (any(Q_ages < 0)) {
        stop("Q times need to be either all non-negative or non-positive.\n")
      }
      
      if (length(Q) != length(Q_ages[Q_ages != 0]) + 1L) {
        stop("number of matrices does not match the number of epoch time boundaries, which is supposed to be one fewer.\n")
      }
      
      if (all(Q_ages != 0)) {
        if (length(Q_ages) > 1 && identical(Q_ages, sort(Q_ages, decreasing = T))) {
          Q_ages <- c(Q_ages, 0)
        } else {
          Q_ages <- c(0, Q_ages)
        }
      }
      
      Q <- Q[order(Q_ages, decreasing = T)]
      Q_ages <- sort(Q_ages, decreasing = T)
      Q_ages <- Q_ages[-length(Q_ages)]
      
    }
  } else if (!is.matrix(Q)) {
    stop("Q needs to be either a single matrix or a list of matrices.\n")
  } else {
    Q <- list(Q)
  }
  
  nstates <- ncol(Q[[1]])
  if (nstates < 2) {
    stop("there has to be at least two states.\n")
  }
  
  # we assume the state names can be found as either the row or column names of the Q matrix of or the name of the root frequency vector
  states <- rownames(Q[[1]])
  if (is.null(states)) {
    states <- colnames(Q[[1]])
    if (is.null(states)) {
      if (!is.null(root_freq)) {
        states <- names(root_freq)
      }
      if (is.null(states)) {
        stop("state names are not assigned.\n")
      }
    }
  }
  
  # sanity checks for Q matrix
  for (i in 1:length(Q)) {
    if (ncol(Q[[i]]) != nrow(Q[[i]])) {
      stop("all matrices need to be square.\n")
    }
    if (ncol(Q[[i]]) != nstates) {
      stop("all matrices need to have identical dimensions.\n")
    }
    if (!all.equal(as.numeric(apply(Q[[i]], 1, sum)), rep(0, nstates), tolerance = 1e-3)) {
      stop("row of Q needs to sum to 0.\n")
    }
    if (any(Q[[i]][row(Q[[i]]) != col(Q[[i]])] < 0)) {
      stop("off diagonal elements of Q needs to be non-negative.\n")
    }
    if (any(diag(Q[[i]]) > 0)) {
      stop("diagonal elements of Q needs to be non-positive.\n")
    }
    if (all(diag(Q[[i]]) == 0)) {
      stop("some elements in Q need to be positive.\n")
    }
  }
  
  if (conditional) {
    return(sim_history_conditional(tree, Q, Q_ages, states, nsim, trait))
  } else {
    return(sim_history_unconditional(tree, Q, Q_ages, states, root_freq, nsim))
  }
}

#' Simulate history of discrete-character changes over a bifurcation tree without conditioning on the tip states
#' @param tree a bifurcation tree of class "phylo"
#' @param Q a instantaneous-rate matrix (or a list of matrices for a piecewise constant geographic model) characterizes the CTMC
#' @param Q_ages boundaries of time intervals for a piecewise constant geographic model (NULL means a constant model)
#' @param states states of the discrete character
#' @param root_freq the vector state frequencies at the root of the tree (that we will draw the root state from); if NULL then it's uniform
#' @param nsim number of simulations to perform
#' @return A phylo and simmap object (or a multiPhylo and multiSimmap object when nsim > 1) that contains the simulated full history
sim_history_unconditional <- function(tree, Q, Q_ages = NULL, states, root_freq = NULL, nsim = 1L) {
  
  nstates <- length(states)

  # sanity checks for root state freq
  if (is.null(root_freq)) { # if not provided then assume it to be uniform
    root_freq <- rep(1 / nstates, nstates)
    names(root_freq) <- states
  } else {
    
    if (any(root_freq < 0)) {
      stop("root frequencies need to be non-negative.\n")
    } else if (all(root_freq == 0)) {
      stop("at least one root frequency needs to be positive.\n")
    }
    
    if (!is.null(names(root_freq))) {
      if (!identical(states, names(root_freq))) {
        if (!identical(sort(states), sort(names(root_freq)))) {
          stop("names of root frequencies are not identical to Q's.\n")
        } else {
          root_freq <- root_freq[states]
        }
      }
    }
     # normalize the vector
    root_freq <- root_freq / sum(root_freq)
  } # end sanity checks
  
  # tree meta info
  nedges <- nrow(tree$edge)
  nnodes <- nedges + 1L
  ntips <- length(tree$tip.label)
  if (2 * ntips - 1L != nnodes) {
    stop("we assume the tree is fully bifurcating")
  }
   
  root <- unique(tree$edge[, 1][!tree$edge[, 1] %in% tree$edge[, 2]])
  if (length(root) != 1) {
    stop("we should have only one root.\n")
  }
  
  trees <- vector("list", nsim)
  class(trees) <- c("multiSimmap", "multiPhylo")
  
  # start simulation
  for (i in 1:nsim) {
    
    node.states <- matrix(NA, nrow = nedges, ncol = 2, byrow = T)
    maps <- vector("list", nedges)
    
    # helper variables
    node_visited <- rep(F, nnodes)
    node_times <- ape::node.depth.edgelength(tree)
    node_ages <- max(node_times) - node_times
    
    Q_times <- NULL
    if (length(Q) > 1) {
      Q_times <- max(node_ages) - Q_ages
    }

    # sample root state first
    node <- root
    state_current <- states[sample.int(nstates, size = 1, prob = root_freq)]
    node.states[tree$edge[, 1] == node, 1] <- state_current
    
    # now we need to traverse the tree to simulate the history
    while (any(!node_visited)) {
      
      # traversal
      if (node %in% tree$edge[, 1]) { # not tip, so must have two descendants
        descendants <- tree$edge[tree$edge[, 1] == node, 2]
        if (!node_visited[descendants[1]]) { # first visit
          node_next <- descendants[1] # move forward
        } else if (!node_visited[descendants[2]]) { # second visit
          node <- descendants[2] # move forward
          next
        } else { # third visit
          node_visited[node] <- T
          if (node == root) {
            break
          } else {
            ancestor <- tree$edge[tree$edge[, 2] == node, 1]
            node <- ancestor # move backward
            next
          }
        }
      } else { # is a tip
        node_visited[node] <- T
        ancestor <- tree$edge[tree$edge[, 2] == node, 1]
        node_next <- ancestor # move backward
      }
      
      if (node != root) { # do simulation on the subtending branch of this node
        
        map <- numeric()
        edge_idx <- which(tree$edge[, 2] == node)
        if (length(edge_idx) != 1) {
          stop("there should be one and only one subtending branch of each non-root node.\n")
        }
        ancestor <- tree$edge[edge_idx, 1]
        edge_length <- tree$edge.length[edge_idx]
        
        time_last <- node_times[ancestor]
        time_current <- node_times[ancestor]
        time_end <- node_times[node]
        
        state_current <- node.states[edge_idx, 1]
        stateidx_current <- match(state_current, states)
        
        while (time_current < time_end) {
          
          # next matrix change time
          if (length(Q) > 1 && any(Q_times > time_current)) {
            time_wall <- Q_times[Q_times > time_current][1]
          } else {
            time_wall <- time_end
          }
          if (time_wall > time_end) {
            time_wall <- time_end
          }

          # fetch the Q matrix in this interval
          if (length(Q) == 1) {
            Q_current <- Q[[1]]
          } else {
            Q_current <- Q[[findInterval(time_current, Q_times) + 1L]]
          }
          
          # draw the event (if any) in this interval
          rate_dominant <- -Q_current[stateidx_current, stateidx_current]
          if (rate_dominant == 0) { # this state is absorbing
            time_current <- time_wall
          } else {
            piece_duration <- rexp(1, rate_dominant)
            time_current <- time_current + piece_duration
            
            if (time_current < time_wall) {
              map_states <- names(map)
              map <- c(map, time_current - time_last)
              names(map) <- c(map_states, state_current)
              time_last <- time_current
              
              prob_vec <- Q_current[stateidx_current, ]
              prob_vec[stateidx_current] <- 0
              state_current <- states[sample.int(nstates, size = 1, prob = prob_vec)]
              stateidx_current <- match(state_current, states)
              
            } else {
              time_current <- time_wall
            }
          }
          
        } # end while loop
        
        node.states[edge_idx, 2] <- state_current
        if (node %in% tree$edge[, 1]) { # assign start state of the descendant branch (if haven't reached tip)
          node.states[tree$edge[, 1] == node, 1] <- state_current
        }
        
        map_states <- names(map)
        map <- c(map, time_end - time_last)
        names(map) <- c(map_states, state_current)
        maps[[edge_idx]] <- map
      } # end simulation on this branch
      
      node <- node_next
    } # end tree traversal
    
    # format the tree as a simmap object
    # tip states
    tip.states <- node.states[match(1:ntips, tree$edge[, 2]), 2]
    names(tip.states) <- tree$tip.label
    
    # mapped.edge
    mapped.edge <- matrix(data = 0, nrow = nedges, ncol = nstates, 
                          dimnames = list(apply(tree$edge, 1, function(x) paste(x, collapse = ",")), state = states))
    for (j in 1:nedges) {
      for (k in 1:length(maps[[j]])) {
        mapped.edge[j, names(maps[[j]])[k]] <- mapped.edge[j, names(maps[[j]])[k]] + maps[[j]][k]
      }
    }
    
    tree$maps <- maps
    tree$node.states <- node.states
    tree$states <- tip.states
    tree$mapped.edge <- mapped.edge
    class(tree) <- c("simmap", class(tree)[class(tree) != "simmap"])
    
    trees[[i]] <- tree
  } # end simulation
  
  # (as convention established by sim.history of phytools) simplify data struture when there is only one resulted history
  if (nsim == 1) {
    trees <- trees[[1]]
  }
  
  return(trees)
}

#' Simulate history of discrete-character changes over a bifurcation tree conditioning on the tip states (currently on all node states)
#' we can relax this to allow only conditioning on the tip states, which would require an implementation of pruning algorithm (including transition-probablity matrix computation)
#' @param tree a bifurcation tree of class "phylo" 
#' (here we assume it either has a node.states, node.data, or maps component so that we can fetch the state of each node, both internal and tip, in the tree to condition on)
#' @param Q a instantaneous-rate matrix (or a list of matrices for a piecewise constant geographic model) characterizes the CTMC
#' @param Q_ages boundaries of time intervals for a piecewise constant geographic model (NULL means a constant model)
#' @param states states of the discrete character
#' @param nsim number of simulations to perform
#' @return A phylo and simmap object (or a multiPhylo and multiSimmap object when nsim > 1) that contains the simulated full history
sim_history_conditional <- function(tree, Q, Q_ages = NULL, states, nsim = 1L, trait = NULL) {
  
  nstates <- length(states)
  
  # tree meta info
  nedges <- nrow(tree$edge)
  nnodes <- nedges + 1L
  ntips <- length(tree$tip.label)
  if (2 * ntips - 1L != nnodes) {
    stop("we assume the tree is fully bifurcating")
  }

  # we need to fetch the node states when thers is not already node.states
  if (!("node.states" %in% names(tree))) {
    
    if ("maps" %in% names(tree)) {
      
      node.state <- character(nnodes)
      node.state[tree$edge[, 2]] <- sapply(tree$maps, function(x) names(x)[length(x)])
      
      root <- unique(tree$edge[, 1][!tree$edge[, 1] %in% tree$edge[, 2]])
      if (length(root) != 1 || node.state[root] != "") {
        stop ("root ill-formed.\n")
      }
      
      node.state[root] <- names(tree$maps[[which(tree$edge[, 1] == root)[1]]])[1]
      
      if (any(!node.state %in% states)) { # currently we assume internal node states are all resolved
        stop ("node state info ill-formed.\n")
      }
      
    } else if ("node.data" %in% names(tree)) {
      
      newickext_names <- gsub("=\\{\\{(.*?)\\}\\}", "", tree$node.data)
      newickext_names <- gsub("=\\{(.*?)\\}", "", newickext_names)
      newickext_names <- lapply(strsplit(newickext_names, ","), function(x) gsub("=.*$", "", x))
      newickext_name <- unique(unlist(newickext_names))
      
      if (is.null(trait)) {
        trait_nothistory <- grep("history|rate|\\.prob$|\\.set$|length|height|range|posterior|HPD|median", newickext_name, invert = T, value = T)
        if (length(trait_nothistory) == 1 && all(sapply(newickext_names, function(x) trait_nothistory %in% x))) {
          trait <- trait_nothistory
        } else {
          stop("Cannot figure out tag name for the trait (i.e., node state).")
        }
      } else if (any(!sapply(newickext_names, function(x) trait %in% x))) {
        stop("trait not found in (at least) some of the nodes")
      }
      
      node.state <- str_match(tree$node.data, paste0(trait, "=\"(.*?)\""))[, 2]
      if (length(node.state) != nnodes || any(!node.state %in% states)) { # currently we assume internal node states are all resolved
        stop ("node state info ill-formed.\n")
      }

    } else {
      # todo
      stop("currently this function assumes the internal node states are known.\n")
    }
    
    tree$node.states <- matrix(node.state[tree$edge], ncol = 2)
  } # end fetching
  
  nodes_times <- ape::node.depth.edgelength(tree)
  nodes_ages <- max(nodes_times) - nodes_times
  
  # create mu and the uniformized R matrix as global variable as they need to be updated by the fill_DTMCmats function as needed
  assign("mu", sapply(Q, function(x) max(-diag(x))), envir = .GlobalEnv)
  assign("R", vector("list", length(Q)), envir = .GlobalEnv)
  fill_DTMCmats(Q, 1)
  
  trees <- vector("list", nsim)
  # start simulation
  for (i in 1:nsim) {
    maps <- vector("list", nedges)
    for (j in 1:nedges) { # loop over branches
      node_ages <- nodes_ages[tree$edge[j, ]]
      node_state_indices <- match(tree$node.states[j, ], states)
      maps[[j]] <- sim_history_branch_conditional(node_ages, Q, Q_ages, node_state_indices, states)
    }
    
    tree$maps_posthoc <- maps
    trees[[i]] <- tree
  }
  
  if (nsim == 1) {
    trees <- trees[[1]]
  }
  
  return(trees)
  
  # todo: when internal state are not available
  # 1: write pruning algorithm ourselves to obtain the partial likelihoods for each node (based on the precomputed P matrix of each branch)
  # 2: reconstruct the joint ancestral states by doing the up-pass using the computed partial likelihoods, root frequency (if provided), and P matrices
  # 3: stochastic mapping (same as what we do when we have internal state estimates already)
}

