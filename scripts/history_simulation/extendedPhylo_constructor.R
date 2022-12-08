# this script contains functions that process an extendedPhylo object

#' Remove intermediate nodes (node created by a character-change event) from an extendedPhylo object
#' @param extendedPhylo an extendedPhylo tree object
#' @return the extendedPhylo tree object (without any intermediate nodes)
removeIntermediatenodes <- function(extendedPhylo) {
  
  if (any(is.na(extendedPhylo$descendant.left) + is.na(extendedPhylo$descendant.right) == 1)) { # intermediate nodes exist
    
    for (i in 1:nrow(extendedPhylo)) { # loop over nodes
      
      if ((is.na(extendedPhylo[i, ]$descendant.left) + is.na(extendedPhylo[i, ]$descendant.right) != 1) && 
          (!is.na(extendedPhylo[i, ]$ancestor))) { # if this node is a tree node but not a root node
        
        current <- extendedPhylo[i, ]$node
        subtending.branch.length <- 0
        
        # traverse along that branch until reach another tree node
        repeat{
          subtending.branch.length <- subtending.branch.length + extendedPhylo[extendedPhylo$node == current,]$subtending.edge.length
          
          ancestor <- extendedPhylo[extendedPhylo$node == current, ]$ancestor
          if (is.na(extendedPhylo[extendedPhylo$node == ancestor, ]$descendant.left) + is.na(extendedPhylo[extendedPhylo$node == ancestor, ]$descendant.right) != 1) break
          
          current <- ancestor
        }
        
        # connect the start tree node and the end tree node of each branch
        extendedPhylo[i, ]$ancestor <- ancestor
        if (extendedPhylo[extendedPhylo$node == ancestor, ]$descendant.left == current) {
          extendedPhylo[extendedPhylo$node == ancestor, ]$descendant.left <- extendedPhylo[i, ]$node
        } else if (extendedPhylo[extendedPhylo$node == ancestor, ]$descendant.right == current) {
          extendedPhylo[extendedPhylo$node == ancestor, ]$descendant.right <- extendedPhylo[i, ]$node
        }
        
        extendedPhylo[i, ]$subtending.edge.length <- subtending.branch.length
      }
    }
    
    extendedPhylo <- extendedPhylo[is.na(extendedPhylo$descendant.left) + is.na(extendedPhylo$descendant.right) != 1, ]
  }
  
  return(extendedPhylo)
}

#' Generate the map component like in a simmap object and attach it to the extendedPhylo
#' @param extendedPhylo an extendedPhylo tree object
#' @return the extendedPhylo tree object (with the map component attached)
getBranchmap <- function(extendedPhylo) {
  
  extendedPhylo$branchmap <- NA
  
  if (any(is.na(extendedPhylo$descendant.left) + is.na(extendedPhylo$descendant.right) == 1)) { # intermediate nodes exist
    
    for (i in 1:nrow(extendedPhylo)) { # loop over nodes
      
      if ((is.na(extendedPhylo[i, ]$descendant.left) + is.na(extendedPhylo[i, ]$descendant.right) != 1) &&
          (!is.na(extendedPhylo[i, ]$ancestor))) { # if this node is a tree node but not a root node
        
        current <- extendedPhylo[i, ]$node
        duration_vec <- NULL
        state_vec <- NULL
        
        # traverse along that branch until reach another tree node
        repeat{
          duration_vec <- c(extendedPhylo[extendedPhylo$node == current,]$subtending.edge.length, duration_vec)
          state_vec <- c(extendedPhylo[extendedPhylo$node == current,]$state, state_vec)
          
          ancestor <- extendedPhylo[extendedPhylo$node == current,]$ancestor
          if (is.na(extendedPhylo[extendedPhylo$node == ancestor,]$descendant.left) + is.na(extendedPhylo[extendedPhylo$node == ancestor,]$descendant.right) != 1) break
          
          current <- ancestor
        }
        
        branchmap <- duration_vec
        names(branchmap) <- state_vec
        
        extendedPhylo[i, ]$branchmap <- list(branchmap)
      }
      
    }
  }
  
  return(extendedPhylo)
}

#' Retrieve the set of tip indices 
#' @param node_df an extendedPhylo tree object
#' @param tip tip indices or tip labels
#' @param tip_states the discrete-character state(s) of the desired tips to fetch
#' @param tip_labels the labels of the desired tips to fetch
#' @return an integer vector containing the desired tip indices
get.tip <- function(node_df, tip = NULL, tip_states = NULL, tip_labels = NULL) {
  if (is.null(tip)) {
    if (is.null(tip_states) && is.null(tip_labels)) {
      stop("no feature specified therefore no tip selected")
    } else if (is.null(tip_states)) {
      tip_states <- unique(node_df$state)
    } else if (is.null(tip_labels)) {
      tip_labels <- node_df$label[!is.na(node_df$label)]
    }
    # find the tips according to the provided tip_states and/or tip_labels
    tip <- node_df$node[is.na(node_df$descendant.left) & is.na(node_df$descendant.right) & 
                          node_df$state %in% tip_states & node_df$label %in% tip_labels]
  } else if (is.character(tip)) {
    if (any(!tip %in% node_df$label)) {
      stop("cannot find some tip in this tree.\n")
    }
    tip <- node_df[match(tip, node_df$label), ]$node
  } # could add more sanity checks here
  
  if (length(tip) == 0) stop("cannot find any tips")
  tip <- unique(tip)
  return(tip)
}

#' Retrieve the path from each of the provided tips to the root
#' @param node_df an extendedPhylo tree object
#' @param tip tip labels or tip indices
#' @param include_tip whether the tip itself is included in the path or not
#' @return a list containing path from each of the provide tip to the root
get.path.toroot <- function(node_df, tip = NULL, tip_states = NULL, tip_labels = NULL, include_tip = F) {
  
  tip <- get.tip(node_df, tip, tip_states, tip_labels)
  ntips_clade <- length(tip)
  routes <- vector("list", ntips_clade)
  for (i in 1:ntips_clade) {
    node <- tip[i]
    if (include_tip) routes[[i]] <- tip[i]
    
    repeat { # move from the tip backwards until reaching the root
      if (is.na(node_df[node_df$node == node, ]$ancestor)) break
      routes[[i]] <- c(routes[[i]], node_df[node_df$node == node, ]$ancestor)
      node <- node_df[node_df$node == node, ]$ancestor
    }
  }
  
  if (ntips_clade == 1) {
    routes <- routes[[1]]
  }
  return(routes)
}

#' Retrieve the path from the MRCA of the provided tips to the root
#' @param node_df an extendedPhylo tree object
#' @param tip tip labels or tip indices
#' @return a vector containing the nodes along the way from the MRCA to the root (two ends included)
get.MRCA.path.toroot <- function(node_df, tip) {
  
  if (is.character(tip)) {
    if (any(!tip %in% node_df$label)) {
      stop("cannot find some tip in this tree.\n")
    }
    tip <- node_df[match(tip, node_df$label), ]$node
  } # could add more sanity checks here
  
  tip <- unique(tip)
  ntips_clade <- length(tip)
  routes <- vector("list", ntips_clade)
  for (i in 1:ntips_clade) {
    node <- tip[i]
    
    repeat { # move from the tip backwards until reaching the root
      if (is.na(node_df[node_df$node == node, ]$ancestor)) break
      routes[[i]] <- c(routes[[i]], node_df[node_df$node == node, ]$ancestor)
      node <- node_df[node_df$node == node, ]$ancestor
    }
  }
  
  # find the intersection of the routes
  routes_tab <- table(unlist(routes))
  ancs_inall <- as.integer(names(routes_tab)[routes_tab == ntips_clade])
  
  if (length(ancs_inall) == 0) { # MRCA is the root
    mrca <- node_df[is.na(node_df$ancestor), ]$node
  } else {
    mrca <- ancs_inall[match(ancs_inall, routes[[1]]) == min(match(ancs_inall, routes[[1]]))] # because the first one encounter must be the most tipy one
  }
  
  return(routes[[1]][match(mrca, routes[[1]]):length(routes[[1]])])
}

#' Retrieve the path (as an extendePhylo so that the branch lengths, ages, states associated with the nodes on the path) 
#' from the MRCA of the provided tips to the root
#' @param tree an extendedPhylo tree object
#' @param tip tip indices or tip labels
#' @param tip_states discrete-character state(s) at the tip
#' @param tip_labels tip labels
#' @return a vector containing the nodes along the way from the MRCA to the root (two ends included)
getPathtree <- function(tree, tip = NULL, tip_states = NULL, tip_labels = NULL, include_root = F) {
  
  tip <- get.tip(tree, tip, tip_states, tip_labels)
  ntips <- length(tip)
  # get the path from the MRCA of the provided tips to the root
  if (ntips == 1) {
    mrca_path_toroot <- get.path.toroot(node_df = tree, tip = tip, include_tip = T)
  } else {
    mrca_path_toroot <- get.MRCA.path.toroot(tree, tip)
  }
  if (!include_root) mrca_path_toroot <- mrca_path_toroot[mrca_path_toroot != tree$node[is.na(tree$ancestor)]] # exclude the root from the path
  path_tree <- tree[tree$node %in% mrca_path_toroot, ]
  
  return(path_tree)
}

#' Retrieve the MRCA of the provided tips
#' @param node_df an extendedPhylo tree object
#' @param tip tip labels or tip indices
#' @return the node index (an integer) of the MRCA
get.MRCA <- function(node_df, tip) {
  
  if (is.character(tip)) {
    if (any(!tip %in% node_df$label)) {
      stop("cannot find some tip in this tree.\n")
    }
    tip <- node_df[match(tip, node_df$label), ]$node
  } # could add more sanity checks here
  
  tip <- unique(tip)
  if (length(tip) == 1) { # only a single tip is provided, then return itself
    return(tip[1])
  } else { # more than one tips
    
    ntips_clade <- length(tip)
    routes <- vector("list", ntips_clade)
    for (i in 1:ntips_clade) {
      node <- tip[i]
      
      repeat { # move from the tip backwards until reaching the root
        if (is.na(node_df[node_df$node == node, ]$ancestor)) break
        routes[[i]] <- c(routes[[i]], node_df[node_df$node == node, ]$ancestor)
        node <- node_df[node_df$node == node, ]$ancestor
      }
    }
    
    routes_tab <- table(unlist(routes))
    ancs_inall <- as.integer(names(routes_tab)[routes_tab == ntips_clade])
    
    if (length(ancs_inall) == 0) { # MRCA is the root
      return(node_df[is.na(node_df$ancestor), ]$node)
    } else {
      return(ancs_inall[match(ancs_inall, routes[[1]]) == min(match(ancs_inall, routes[[1]]))]) # because the first one encounter must be the most tipy one
    }
  }
}

#' Retrieve the descendants of a provided node
#' @param node_df an extendedPhylo tree object
#' @param node node index
#' @param full include internal nodes and tips (when true) or only the tip descendants (when false)
#' @return a vector of descendant node indices
#' todo: double check to make sure the recursion work as expected
get.Descendants <- function(node_df, node, full = T) {
  
  descendants_vec <- NULL
  get.Descendants.recursive <- function(node_df, node, descendants_vec, full = T) { # a recursive function to retrieve descendants (only returns when reach a tip)
    
    descendant.left <- node_df[node_df$node == node, ]$descendant.left
    descendant.right <- node_df[node_df$node == node, ]$descendant.right
    
    if (is.na(descendant.left) + is.na(descendant.right) == 2) { # tip node
      descendants_vec <- c(descendants_vec, node)
    } else if (is.na(descendant.left) + is.na(descendant.right) == 0) { # internal node
      
      if (full) {
        descendants_vec <- c(descendants_vec, node)
      }
      descendants_vec <- get.Descendants.recursive(node_df, descendant.left, descendants_vec, full = full)
      descendants_vec <- get.Descendants.recursive(node_df, descendant.right, descendants_vec, full = full)
      
    } else if (is.na(descendant.left) + is.na(descendant.right) == 1) { # intermediate node
      
      if (!is.na(descendant.left)) {
        if (full) {
          descendants_vec <- c(descendants_vec, node)
        }
        descendants_vec <- get.Descendants.recursive(node_df, descendant.left, descendants_vec, full = full)
      } else if (!is.na(descendant.right)) {
        if (full) {
          descendants_vec <- c(descendants_vec, node)
        }
        descendants_vec <- get.Descendants.recursive(node_df, descendant.right, descendants_vec, full = full)
      }
    }
    
    return(descendants_vec)
  }
  
  return(get.Descendants.recursive(node_df, node, descendants_vec, full = full))
}

#' Compute the clade size of each node (the number of its tree nodes descendants, including itself)
#' @param node_df an extendedPhylo tree object
#' @return the extendedPhylo tree object (with the clade.size component attached)
get.clade.size <- function(node_df) {

  # recursively compute the clade size of each node (the number of its tree nodes descendants, including itself) 
  get.clade.size.recursive <- function(descendant.left, descendant.right, node, clade.size) {
    
    descendant.left_this <- descendant.left[node]
    descendant.right_this <- descendant.right[node]
    
    if (is.na(descendant.left_this) + is.na(descendant.right_this) == 2) { # if this node is a tip
      clade.size[node] <- 1L
    } else if (is.na(descendant.left_this) + is.na(descendant.right_this) == 0) { # if this node is an internal node
      
      clade.size <- get.clade.size.recursive(descendant.left, descendant.right, descendant.left_this, clade.size)
      clade.size <- get.clade.size.recursive(descendant.left, descendant.right, descendant.right_this, clade.size)
      clade.size[node] <- sum(clade.size[c(descendant.left_this, descendant.right_this)])
      
    } else if (is.na(descendant.left_this) + is.na(descendant.right_this) == 1) { # if this node is a character-change node
      
      if (!is.na(descendant.left_this)) {
        clade.size <- get.clade.size.recursive(descendant.left, descendant.right, descendant.left_this, clade.size)
        clade.size[node] <- clade.size[descendant.left_this]
      } else if (!is.na(descendant.right_this)) {
        clade.size <- get.clade.size.recursive(descendant.left, descendant.right, descendant.right_this, clade.size)
        clade.size[node] <- clade.size[descendant.right_this]
      }
    }
    
    return(clade.size)
  }

  npieces <- nrow(node_df)
  clade.size <- numeric(npieces)
  # starting the traversal from the root
  root <- which(is.na(node_df$ancestor))
  if (length(root) != 1) {
    stop("there must be one and only one root")
  }
  node_df$clade.size <- get.clade.size.recursive(node_df$descendant.left, node_df$descendant.right, root, clade.size)

  return(node_df)
}

#' Order the tips to ladderize the tree decreasingly (smaller clades is on top) or increasingly
#' @param node_df an extendedPhylo tree object
#' @param order whether ladderize the tree decreasingly (smaller clades is on top) or increasingly
#' @return an ordered extendedPhylo tree object
get.new.tipOrder <- function(node_df, order = "decreasing") {
  
  # sanity check
  if (!order %in% c("decreasing", "increasing")) {
    stop("I don't understand the order; it should be either decreasing or increasing.")
  }
  if (!"clade.size" %in% colnames(node_df)) {
    stop("clade.size has to be included as a column of the provided tree.")
  }
  
  is_orderdec <- order == "decreasing"
  npieces <- nrow(node_df)
  new.tipOrder <- integer(npieces)
  ancestor <- node_df$ancestor
  descendant.left <- node_df$descendant.left
  descendant.right <- node_df$descendant.right
  clade.size <- node_df$clade.size
  descendant.left_exists <- !is.na(descendant.left)
  descendant.right_exists <- !is.na(descendant.right)
  ndescendants_all <- descendant.left_exists + descendant.right_exists
  tipid <- 1L
  
  # now we need to traverse the tree
  node_visited <- rep(F, npieces)
  root <- which(is.na(node_df$ancestor))
  if (length(root) != 1) {
    stop("there must be one and only one root")
  }
  node <- root # starting the traversal from the root
  repeat { # traverse from root to tip then back to root
    
    descendant.left_this <- descendant.left[node]
    descendant.right_this <- descendant.right[node]
    descendants_this <- c(descendant.left_this, descendant.right_this)
    descendants_this <- descendants_this[!is.na(descendants_this)]
    
    # go to its ancestor (or finish if node is root) when node with two descendants both have been visited or with one descendant has been visited
    if ((ndescendants_all[node] == 2 && sum(node_visited[descendants_this]) == 2) || 
        (ndescendants_all[node] == 1 && node_visited[descendants_this])) { 
      node_visited[node] <- T
      if (node == root) {
        break
      } else {
        node <- ancestor[node]
      }
    } else if (ndescendants_all[node] == 1) { # node with one descendant has not been visited, then move to that descendant
      node <- ifelse(descendant.left_exists[node], descendant.left_this, descendant.right_this)
    } else if (ndescendants_all[node] == 2 && sum(node_visited[descendants_this]) == 1) {
      # node with two descendants and only one of them has been visted, then move to the other descendant
      node <- ifelse(node_visited[descendant.left_this], descendant.right_this, descendant.left_this)
    } else if (ndescendants_all[node] == 2) { # node with two descendants that neither have been visited
      # move to the descendant associated with smaller clade while order decreasingly; otherwise to the descendant associated with larger clade
      node <- ifelse(((clade.size[descendant.left_this] <= clade.size[descendant.right_this]) + is_orderdec) != 1, descendant.left_this, descendant.right_this)
    } else { # ndescendants_all[node] must equal 0 now so node must be a tip
      node_visited[node] <- T
      new.tipOrder[node] <- tipid
      tipid <- tipid + 1L
      node <- ancestor[node]
    }
  }
  
  node_df$new.tipOrder <- new.tipOrder
  return(node_df)
}

#' Compute the height (horizontal distance from the root) of each node
#' @param node_df an extendedPhylo tree object
#' @return the extendedPhylo tree object (with the node.heights component attached)
get.node.heights <- function(node_df) {
  
  # recursively traversing the tree to get the height of each node
  get.node.heights.recursive <- function(ancestor, descendant.left, descendant.right, subtending.edge.length, node, node.heights) {
    
    ancestor_this <- ancestor[node]
    descendant.left_this <- descendant.left[node]
    descendant.right_this <- descendant.right[node]
    isna_descendant.left_this <- is.na(descendant.left_this)
    isna_descendant.right_this <- is.na(descendant.right_this)
    
    if (!is.na(ancestor_this)) { # not at the root, compute the height (horizontal distance from the root) of this node
      node.heights[node] <- node.heights[ancestor_this] + subtending.edge.length[node]
    }
    
    # if we are not at a tip, keep traversing the tree until we reach one
    if (isna_descendant.left_this + isna_descendant.right_this != 2) {
      if (!isna_descendant.left_this) {
        node.heights <- get.node.heights.recursive(ancestor, descendant.left, descendant.right, subtending.edge.length, descendant.left_this, node.heights)
      }
      if (!isna_descendant.right_this) {
        node.heights <- get.node.heights.recursive(ancestor, descendant.left, descendant.right, subtending.edge.length, descendant.right_this, node.heights)
      }
    }
    
    return(node.heights)
  }
  
  npieces <- nrow(node_df)
  node.heights <- numeric(npieces)
  # starting the traversal from the root
  root <- which(is.na(node_df$ancestor))
  if (length(root) != 1) {
    stop("there must be one and only one root")
  }
  if (!is.na(node_df$subtending.edge.length[root])) {
    node.heights[root] <- subtending.edge.length[root]
  }
  node_df$node.heights <- get.node.heights.recursive(node_df$ancestor, node_df$descendant.left, node_df$descendant.right, node_df$subtending.edge.length, root, node.heights)
  
  # start from the root
  return(node_df)
}

#' Compute the vertical position of each node
#' @param node_df an extendedPhylo tree object
#' @return the extendedPhylo tree object (with the node.vertical component attached)
get.node.vertical <- function(node_df) {
  
  # recursively compute the vertical position of each node
  get.node.vertical.recursive <- function(descendant.left, descendant.right, new.tipOrder, node, node.vertical) {
    
    descendant.left_this <- descendant.left[node]
    descendant.right_this <- descendant.right[node]
    
    if (is.na(descendant.left_this) + is.na(descendant.right_this) == 2) { # if the node is a tip then its vertical position is the same as its ladderized order index
      node.vertical[node] <- new.tipOrder[node]
      
    } else if (is.na(descendant.left_this) + is.na(descendant.right_this) == 0) { # if the node is an internal tree node its vertical position is the center of its two descendants'
      node.vertical <- get.node.vertical.recursive(descendant.left, descendant.right, new.tipOrder, descendant.left_this, node.vertical)
      node.vertical <- get.node.vertical.recursive(descendant.left, descendant.right, new.tipOrder, descendant.right_this, node.vertical)
      node.vertical[node] <- mean(node.vertical[c(descendant.left_this, descendant.right_this)])

    } else if (is.na(descendant.left_this) + is.na(descendant.right_this) == 1) { # if the node is a character-change node then its vertical position is the same as its descendant's
      
      if (!is.na(descendant.left_this)) {
        node.vertical <- get.node.vertical.recursive(descendant.left, descendant.right, new.tipOrder, descendant.left_this, node.vertical)
        node.vertical[node] <- node.vertical[descendant.left_this]
      } else if (!is.na(descendant.right)) {
        node.vertical <- get.node.vertical.recursive(descendant.left, descendant.right, new.tipOrder, descendant.right_this, node.vertical)
        node.vertical[node] <- node.vertical[descendant.right_this]
      }
    }
    
    return(node.vertical)
  }
  
  npieces <- nrow(node_df)
  node.vertical <- numeric(npieces)
  # starting the traversal from the root
  root <- which(is.na(node_df$ancestor))
  if (length(root) != 1) {
    stop("there must be one and only one root")
  }
  node_df$node.vertical <- get.node.vertical.recursive(node_df$descendant.left, node_df$descendant.right, node_df$new.tipOrder, root, node.vertical)

  return(node_df)
}
