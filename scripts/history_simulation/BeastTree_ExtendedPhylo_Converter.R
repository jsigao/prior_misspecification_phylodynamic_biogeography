library(ape)
library(stringr)
source("./phylo_ExtendedPhylo_Converter.R")

#' Parse the node.data component of a tree (of class phylo) to get the state of each node
#' @param tree a single tree (of class phylo) with node.data component that contains unparsed node info
#' @param trait tag name of ancestral state reconstruction of the discrete character (here we assume that there is a single character)
#' @param trait_history tag name of simulated history of the discrete character
#' @param convert whether convert the object that will be returned into a customized format (extendedPhylo) or keep it to be phylo and simmap
#' @return a tree with node state information
BeasttreeExtendedPhyloConverter <- function(tree, trait = NULL, trait_history = NULL, convert = T) {
  
  if (!"node.data" %in% names(tree)) {
    stop("cannot find node annotations.\n")
  }
  
  # get all the tag names
  newickext_names <- gsub("=\\{\\{(.*?)\\}\\}", "", tree$node.data)
  newickext_names <- gsub("=\\{(.*?)\\}", "", newickext_names)
  newickext_names <- lapply(strsplit(newickext_names, ","), function(x) gsub("=.*$", "", x))
  newickext_name <- unique(unlist(newickext_names))
  
  # try to figure out trait tag name if it's not specified
  if (is.null(trait)) {
    trait_nothistory <- grep("history|rate|\\.prob$|\\.set$|length|height|range|posterior|HPD|median", newickext_name, invert = T, value = T)
    if (length(trait_nothistory) == 1 && all(sapply(newickext_names, function(x) trait_nothistory %in% x))) {
      trait <- trait_nothistory
    } else {
      stop("cannot figure out tag name for the trait (i.e., node state).")
    }
  } else if (any(!sapply(newickext_names, function(x) trait %in% x))) {
    stop("trait not found in (at least) some of the nodes")
  }
  
  if (is.null(trait_history)) {
    trait_his <- grep("history_all", newickext_name, value = T)
    if (length(trait_his) == 1) {
      trait_history <- trait_his
    } else if (length(trait_his) > 1) {
      stop("cannot figure out tag name for the trait history (i.e., node stochastic maps).")
    }
  } else if (!trait_history %in% newickext_name) {
    stop("trait history not found in the tree")
  }
  
  if (all(paste0(trait, c(".set", ".set.prob")) %in% newickext_name)) { # for the MCC tree
    
    nedges <- nrow(tree$edge)
    nnodes <- nedges + 1L
    
    # get the root
    root <- which(!(1:nnodes %in% tree$edge[, 2]))
    if (length(root) != 1) {
      stop("there should be only one root.\n")
    }
    
    posterior <- NULL
    if ("posterior" %in% newickext_name) {
      posterior <- "posterior"
    } else {
      message("cannot find posterior in this mcc tree.\n")
    }
    posterior.prob <- numeric(length(tree$node.data))
    
    node.state <- character(length(tree$node.data))
    node.state.set <- vector("list", length(tree$node.data))
    
    for (i in 1:length(tree$node.data)) { # loop over all the nodes to parse their associated node info
      
      node.info.split <- unlist(lapply(unlist(strsplit(tree$node.data[i], "\\},|=")), function(x) {
        if (!grepl("\\{", x)) {
          unlist(strsplit(x, ","))
        } else {
          x
        }
      }))
      
      if ((!is.null(posterior)) && i != root && posterior %in% node.info.split) {
        posterior.prob[i] <- as.numeric(node.info.split[which(node.info.split == posterior) + 1])
      }
      
      node.state[i] <- gsub("\"", "", node.info.split[which(node.info.split == trait) + 1])
      node.state.set[[i]] <- as.numeric(gsub("\\{|\"|\\}", "", unlist(strsplit(node.info.split[which(node.info.split == paste0(trait, ".set.prob")) + 1], ","))))
      names(node.state.set[[i]]) <- gsub("\\{|\"|\\}", "", unlist(strsplit(node.info.split[which(node.info.split == paste0(trait, ".set")) + 1], ",")))
      
      if (!(node.state[i] %in% names(node.state.set[[i]]))) {
        state <- unlist(strsplit(node.state[i], "\\+"))
        if (length(state) > 1 && all(state %in% names(node.state.set[[i]]))) {
          node.state[i] <- state[sample.int(length(state), size = 1)]
        } else {
          stop("The state recorded as the map state is inconsistent with the map state in the states set vector")
        }
      }
      
      if (!(node.state[i] %in% names(node.state.set[[i]])[node.state.set[[i]] == max(node.state.set[[i]])])) {
        stop("The state recorded as the map state is inconsistent with the map state in the states set vector")
      }
    }
    
    if (convert) {
      extendedPhylo <- phyloExtendedPhyloConverter(tree) # convert the object from class phylo and simmap into extendedPhylo
      extendedPhylo$state <- node.state[extendedPhylo$node]
      extendedPhylo$state.set <- node.state.set[extendedPhylo$node]
      extendedPhylo$posterior <- posterior.prob[extendedPhylo$node]
    } else {
      tree$node.state <- node.state
      tree$node.state.set <- node.state.set
      tree$node.posterior <- posterior.prob
    }

  } else if ((!is.null(trait_history)) && trait_history %in% newickext_name) { # for a tree (with history) from the posterior distribution
    
    node.heights <- max(node.depth.edgelength(tree)) - node.depth.edgelength(tree)
    edge_maps <- vector("list", nrow(tree$edge)) # prepare simmap format
    
    for (i in 1:nrow(tree$edge)) { # loop over all the branches to parse their associated node info
      
      history_edge <- tree$node.data[tree$edge[i, 2]]
      
      if (!grepl(trait_history, history_edge)) { # no simulated history
        
        state_str <- strsplit(history_edge, "=|,")[[1]]
        state <- gsub("\"", "", state_str[which(state_str == trait) + 1])
        edge_length <- node.heights[tree$edge[i, 1]] - node.heights[tree$edge[i, 2]]
        
        edge_maps[[i]] <- edge_length
        names(edge_maps[[i]]) <- state
        
      } else { # with simulated history

        history_edge_str <- stringr::str_match(history_edge, paste0(trait_history, "=\\{\\{(.*?)\\}\\}"))[1, 2]
        history_edge_df <- data.frame(matrix(unlist(strsplit(history_edge_str, ",")), ncol = 4, byrow = T), stringsAsFactors = F)[, -1]
        
        history_edge_df[, 1] <- as.numeric(history_edge_df[, 1])
        if (any(history_edge_df[, 1] > node.heights[tree$edge[i, 1]] & history_edge_df[, 1] < node.heights[tree$edge[i, 2]])) {
          stop("event not on edge.\n")
        }
        
        edge_length <- -diff(c(node.heights[tree$edge[i, 1]], history_edge_df[, 1], node.heights[tree$edge[i, 2]]))
        state <- c(history_edge_df[, 2], history_edge_df[nrow(history_edge_df), 3])
        
        edge_maps[[i]] <- edge_length
        names(edge_maps[[i]]) <- state
      }
    }
    
    tree$maps <- edge_maps
    
    if (convert) {
      extendedPhylo <- phyloExtendedPhyloConverter(tree)
    }
    
  } else if (!is.null(trait)) { # for a tree (without history) from the posterior distribution
    
    node.state <- character(length(tree$node.data))
    for (i in 1:length(tree$node.data)) {
      state_str <- strsplit(tree$node.data[i], "=|,")[[1]]
      node.state[i] <- gsub("\"", "", state_str[which(state_str == trait) + 1])
    }
    
    if (convert) {
      extendedPhylo <- phyloExtendedPhyloConverter(tree)
      extendedPhylo$state <- node.state[extendedPhylo$node]
    } else {
      tree$node.state <- node.state
    }
  } else {
    stop("either the tree is not formatted as we expect or the provided trait or trait_history argument does not match the tree.\n")
  }
  
  if (convert) {
    return(extendedPhylo)
  } else {
    return(tree)
  }
}

#' Parse the node.data component of a tree (of class phylo) or a list of such trees to get the state of each node
#' @param tree a single tree (of class phylo) with node.data component that contains unparsed node info or a list of such trees
#' @param trait tag name of ancestral state reconstruction of the discrete character (here we assume that there is a single character)
#' @param trait_history tag name of simulated history of the discrete character
#' @param convert whether convert the object that will be returned into a customized format (extendedPhylo) or keep it to be phylo and simmap
#' @return a (list of) tree(s) with node state information
BeastTreeExtendedPhyloConverter <- function(tree, trait = NULL, trait_history = NULL, convert = T) {
  
  if (is.null(tree$edge)) {
    extendedPhylo <- lapply(tree, function(x) BeasttreeExtendedPhyloConverter(x, trait = trait, trait_history = trait_history, convert = convert))
  } else {
    extendedPhylo <- BeasttreeExtendedPhyloConverter(tree, trait = trait, trait_history = trait_history, convert = convert)
  }
  
  return(extendedPhylo)
}
