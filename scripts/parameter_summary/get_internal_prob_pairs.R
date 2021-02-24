library(ape)
source('../history_simulation/BeastTree_parser.R')
source("../history_simulation/BeastTree_ExtendedPhylo_Converter.R")

#' Process a pair of MCC trees to summarize the posterior probability of the maximum a posteriori (MAP) state 
#' (determined by the first provided tree) at each of the shared internal nodes between the trees.
#' @param tree_paths A pair of BEAST MCC tree paths
#' @return a matrix (where each row corresponds to a shared internal node) with five columns, 
#' containing the posterior probabilities of the MAP state of the pair of MCC trees (first two columns), 
#' a bool value indicating whether the MAP state are identical between the pair (third column), and
#' the relative time (e.g., root is zero and the most recent tip is 1) of the node in the two trees (last two columns)
getInternalProbPairs <- function(tree_paths) {

  # read in MCC trees
  trees <- lapply(tree_paths, readBeast)
  trees <- lapply(trees, function(x) BeasttreeExtendedPhyloConverter(x, convert = F))
  class(trees) <- "multiPhylo"
  # remove quotation marks wrapping the tip labels
  for (j in 1:length(trees)) {
    trees[[j]]$tip.label <- gsub("\'", "", trees[[j]]$tip.label)
  }
  
  # orient the pair of trees against the first tree
  # comment this out as this may not be necessary as the tip labels are already sorted alphabetically with readBeast
  # todo: make this more robust by performing the reorientation and, importantly, also reorder the node info objects accordingly
  # exemplar <- trees[1]
  # class(exemplar) <- "multiPhylo"
  # these_trees <- c(exemplar, trees)
  # these_trees <- ape::.compressTipLabel(these_trees)
  # trees <- these_trees[-1]
  
  # find the splits that appears in both MCC trees
  # first get the tip descendants of the shared internal nodes/splits
  parts <- lapply(trees, ape::prop.part)
  parts_labels <- lapply(parts, function(part) {
    lapply(part, function(x) sort(attr(part, "labels")[x]))
  })
  parts_labels_intersect <- intersect(parts_labels[[1]], parts_labels[[2]])
  
  # then get the shared internal nodes pairs
  internal_nodes_pairs <- t(sapply(parts_labels_intersect, function(x) {
    sapply(1:2, function(i) ape::getMRCA(trees[[i]], x))
  }))
  
  # traverse the both MCC trees to get the information of each nodes in the extended newick strings
  relative.heights <- lapply(trees, function(x) 1 - node.depth.edgelength(x) / max(node.depth.edgelength(x)))
  node.states <- lapply(trees, function(x) x$node.state.set)
  
  # now go over all the shared internal nodes to find their corresponding ancestral states information
  # here we want the posterior probabilities of the map state of each internal node in the first tree and 
  # the corresponding state (which may not be the map state) of that node in the second tree
  # and the third column of the data frame that would be returned takes that information (whether the map states are shared)
  # the fourth and fifth columns contain the relative heights of both nodes
  internal_nodes_prob_pairs <- matrix(0, nrow = nrow(internal_nodes_pairs), ncol = 5, byrow = T)
  for (j in 1:nrow(internal_nodes_pairs)) {
    
    # obtain the prob and name of each state of that shared node in both trees
    state.probs <- lapply(1:2, function(i) as.numeric(node.states[[i]][[internal_nodes_pairs[j, i]]]))
    state.names <- lapply(1:2, function(i) names(node.states[[i]][[internal_nodes_pairs[j, i]]]))
    
    map_state <- lapply(seq_along(state.probs), function(i) sort(state.names[[i]][state.probs[[i]] == max(state.probs[[i]])]))
    if (length(intersect(map_state[[1]], map_state[[2]])) > 0) {
      internal_nodes_prob_pairs[j, 3] <- 1L
      internal_nodes_prob_pairs[j, 1:2] <- sapply(seq_along(state.probs), function(i) state.probs[[i]][state.names[[i]] == intersect(map_state[[1]], map_state[[2]])[1]])
    } else if (any(map_state[[1]] %in% state.names[[2]])) { # can find map state under default in the alternative state list
      internal_nodes_prob_pairs[j, 1:2] <- sapply(seq_along(state.probs), function(i) state.probs[[i]][state.names[[i]] == intersect(map_state[[1]], state.names[[2]])[1]])
    } else { # cannot find map state under default in the alternative state list
      internal_nodes_prob_pairs[j, 1] <- state.probs[[1]][state.names[[1]] == map_state[[1]][1]]
      internal_nodes_prob_pairs[j, 2] <- 0L
    }
    
    internal_nodes_prob_pairs[j, 4:5] <- sapply(1:2, function(i) relative.heights[[i]][internal_nodes_pairs[j, i]])
  }
  
  return(internal_nodes_prob_pairs)
}
