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

#' Retrieve the path from each of the provided tips to the root
#' @param node_df an extendedPhylo tree object
#' @param tip tip labels or tip indices
#' @return a list containing path from each of the provide tip to the root
get.path.toroot <- function(node_df, tip) {
  
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
#' @param node_df an extendedPhylo tree object
#' @param tip_states discrete-character states at the tip
#' @param tip_labels tip labels
#' @return a vector containing the nodes along the way from the MRCA to the root (two ends included)
getPathtree <- function(tree, tip_states = NULL, tip_labels = NULL) {
  
  if (is.null(tip_states) && is.null(tip_labels)) {
    stop("no feature specified therefore no tip selected")
  } else if (is.null(tip_states)) {
    tip_states <- unique(tree$state)
  } else if (is.null(tip_labels)) {
    tip_labels <- tree$label[!is.na(tree$label)]
  }
  
  # find the tips according to the provided tip_states and/or tip_labels
  tip <- tree$node[is.na(tree$descendant.left) & is.na(tree$descendant.right) & 
                       tree$state %in% tip_states & tree$label %in% tip_labels]
  
  mrca_path_toroot <- get.MRCA.path.toroot(tree, tip) # get the path from the MRCA of the provided tips to the root
  mrca_path_toroot <- mrca_path_toroot[mrca_path_toroot != tree$node[is.na(tree$ancestor)]] # exclude the root from the path
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
  
  node_df$clade.size <- 0
  Ntip <- sum(is.na(node_df$descendant.left) + is.na(node_df$descendant.right) == 2)
  
  # starting the traversal from the root
  node <- Ntip + 1L
  
  # recursively compute the clade size of each node (the number of its tree nodes descendants, including itself) 
  get.clade.size.recursive <- function(node_df, node) {
    
    descendant.left <- node_df[node_df$node == node, ]$descendant.left
    descendant.right <- node_df[node_df$node == node, ]$descendant.right
    
    if (is.na(descendant.left) + is.na(descendant.right) == 2) { # if this node is a tip
      node_df[node_df$node == node, ]$clade.size <- 1L
    } else if (is.na(descendant.left) + is.na(descendant.right) == 0) { # if this node is an internal node
      
      node_df <- get.clade.size.recursive(node_df, descendant.left)
      node_df <- get.clade.size.recursive(node_df, descendant.right)
      node_df[node_df$node == node, ]$clade.size <- sum(c(node_df[node_df$node == descendant.left, ]$clade.size, 
                                                          node_df[node_df$node == descendant.right, ]$clade.size))
      
    } else if (is.na(descendant.left) + is.na(descendant.right) == 1) { # if this node is a character-change node
      
      if (!is.na(descendant.left)) {
        node_df <- get.clade.size.recursive(node_df, descendant.left)
        node_df[node_df$node == node, ]$clade.size <- node_df[node_df$node == descendant.left, ]$clade.size
      } else if (!is.na(descendant.right)) {
        node_df <- get.clade.size.recursive(node_df, descendant.right)
        node_df[node_df$node == node, ]$clade.size <- node_df[node_df$node == descendant.right, ]$clade.size
      }
    }
    
    return(node_df)
  }
  
  return(get.clade.size.recursive(node_df, node))
}

#' Order the tips to ladderize the tree decreasingly (smaller clades is on top) or increasingly
#' @param node_df an extendedPhylo tree object
#' @param order whether ladderize the tree decreasingly (smaller clades is on top) or increasingly
#' @return an ordered extendedPhylo tree object
get.new.tipOrder <- function(node_df, order = "decreasing") {
  
  # sanity check
  if (!order %in% c("decreasing", "increasing")) {
    stop("I don't understand the order now; it should be either decreasing or increasing.")
  }
  
  node_df$new.tipOrder <- 0
  
  # register if the left and/or right descendant exist (thus need to be visited)
  node_df$descendant.left.tovisit <- 1 - is.na(node_df$descendant.left)
  node_df$descendant.right.tovisit <- 1 - is.na(node_df$descendant.right)
  
  Ntip <- sum(is.na(node_df$descendant.left) + is.na(node_df$descendant.right) == 2)
  node <- Ntip + 1L
  
  order.number <- 1L
  
  while (order.number <= Ntip) {
    # move back/root wards until reach a node with non-visited descendants
    node <- move.backward(node_df, node)
    
    # move fore/tip wards until reach a tip
    node_df.node <- move.forward(node_df, node, order)
    node <- node_df.node$node
    node_df <- node_df.node$node_df
    
    # (re)numbering tips so that tips in smaller clade get a smaller number
    node_df[node_df$node == node, ]$new.tipOrder <- order.number
    
    # increment order index
    order.number <- order.number + 1L
  }
  
  # remove temporary columns
  node_df$descendant.left.tovisit <- NULL
  node_df$descendant.right.tovisit <- NULL
  
  return(node_df)
}

#' Move back/root wards until reach a node with non-visited descendants
#' @param node_df an extendedPhylo tree object
#' @param node node index
#' @return index of a node with non-visited descendants
move.backward <- function(node_df, node) {
  
  ancestor <- node_df[node_df$node == node, ]$ancestor
  
  if (!is.na(ancestor)) { # not root
    if (node_df[node_df$node == ancestor, ]$descendant.left.tovisit + node_df[node_df$node == ancestor, ]$descendant.right.tovisit == 0) {
      # when no non-visited descendants, keep moving back
      move.backward(node_df, ancestor)
    } else if (node_df[node_df$node == ancestor, ]$descendant.left.tovisit + node_df[node_df$node == ancestor, ]$descendant.right.tovisit > 0) {
      return(ancestor)
    }
  } else { # at the root
    return(node)
  }
}

#' Move fore/tip wards until reach a tip
#' @param node_df an extendedPhylo tree object
#' @param node node index
#' @param order visit the node with smaller descendant clade first (decreasing) or larger first (increasing)
#' @return tip index
move.forward <- function(node_df, node, order) {
  
  descendant.left <- node_df[node_df$node == node, ]$descendant.left
  descendant.right <- node_df[node_df$node == node, ]$descendant.right
  
  if (node_df[node_df$node == node, ]$descendant.left.tovisit + node_df[node_df$node == node, ]$descendant.right.tovisit == 0) { 
    # we reach a tip or a node whose descendants have all been visited (thus effectively a tip as in the pruning)
    return(list(node_df = node_df, node = node))
    
  } else if (node_df[node_df$node == node, ]$descendant.left.tovisit + node_df[node_df$node == node, ]$descendant.right.tovisit == 1) { 
    # we reach a character-change node or an internal tree node one of whose descendants has been visited
    
    if (node_df[node_df$node == node, ]$descendant.left.tovisit == 1) {
      next.node <- descendant.left # move forward
      node_df[node_df$node == node, ]$descendant.left.tovisit <- 1 - node_df[node_df$node == node,]$descendant.left.tovisit # mark this descendant as visited
    } else if (node_df[node_df$node == node, ]$descendant.right.tovisit == 1) {
      next.node <- descendant.right # move forward
      node_df[node_df$node == node, ]$descendant.right.tovisit <- 1 - node_df[node_df$node == node, ]$descendant.right.tovisit # mark this descendant as visited
    }
  } else if (node_df[node_df$node == node, ]$descendant.left.tovisit + node_df[node_df$node == node, ]$descendant.right.tovisit == 2) { 
    # we reach an internal tree node neither of whose descendants have been visited
    
    # visit the node with smaller descendant clade first when order decreasingly, otherwise (i.e., when order increasingly) larger first
    if (((node_df[node_df$node == descendant.left, ]$clade.size <= node_df[node_df$node == descendant.right, ]$clade.size) + (order == "decreasing")) != 1) {
      next.node <- descendant.left # move forward
      node_df[node_df$node == node, ]$descendant.left.tovisit <- 1 - node_df[node_df$node == node, ]$descendant.left.tovisit # mark this descendant as visited
    } else {
      next.node <- descendant.right # move forward
      node_df[node_df$node == node, ]$descendant.right.tovisit <- 1 - node_df[node_df$node == node, ]$descendant.right.tovisit # mark this descendant as visited
    }
  }
  
  move.forward(node_df, next.node, order) # recursively move forward
}

#' Compute the height (horizontal distance from the root) of each node
#' @param node_df an extendedPhylo tree object
#' @return the extendedPhylo tree object (with the node.heights component attached)
get.node.heights <- function(node_df) {
  
  # recursivly traversing the tree to get the height of each node
  get.node.heights.recursive <- function(node_df, node) {
    
    # if we are not at the root, compute the height (horizontal distance from the root) of this node
    if (!is.na(node_df[node_df$node == node, ]$ancestor)) {
      ancestor <- node_df[node_df$node == node, ]$ancestor
      node_df[node_df$node == node, ]$node.heights <- sum(node_df[node_df$node == ancestor, ]$node.heights + node_df[node_df$node == node, ]$subtending.edge.length, na.rm = T)
    }
    
    # if we are not at a tip, keep traversing the tree until we reach one
    if (is.na(node_df[node_df$node == node, ]$descendant.left) + is.na(node_df[node_df$node == node, ]$descendant.right) < 2){
      if (!is.na(node_df[node_df$node == node, ]$descendant.left)) {
        node_df <- get.node.heights.recursive(node_df, node_df[node_df$node == node, ]$descendant.left)
      } 
      if (!is.na(node_df[node_df$node == node, ]$descendant.right)) {
        node_df <- get.node.heights.recursive(node_df, node_df[node_df$node == node, ]$descendant.right)
      }
    }
    
    return(node_df)
  }
  
  Ntip <- sum(is.na(node_df$descendant.left)+is.na(node_df$descendant.right) == 2)
  node <- Ntip + 1L
  
  # in case there is a root edge
  if (!is.na(node_df[node_df$node == node, ]$subtending.edge.length)) {
    node_df$node.heights <- node_df[node_df$node == node, ]$subtending.edge.length
  } else {
    node_df$node.heights <- 0
  }
  
  # start from the root
  return(get.node.heights.recursive(node_df, node))
}

#' Compute the vertical position of each node
#' @param node_df an extendedPhylo tree object
#' @return the extendedPhylo tree object (with the node.vertical component attached)
get.node.vertical <- function(node_df) {
  
  # recusively compute the vertical position of each node
  get.node.vertical.recursive <- function(node_df, node) {
    
    descendant.left <- node_df[node_df$node == node, ]$descendant.left
    descendant.right <- node_df[node_df$node == node, ]$descendant.right
    
    if (is.na(descendant.left) + is.na(descendant.right) == 2) { # if the node is a tip then its vertical position is the same as its ladderized order index
      node_df[node_df$node == node,]$node.vertical <- node_df[node_df$node == node, ]$new.tipOrder
      
    } else if (is.na(descendant.left) + is.na(descendant.right) == 0) { # if the node is an internal tree node its vertical position is the center of its two descendants'
      node_df <- get.node.vertical.recursive(node_df, descendant.left)
      node_df <- get.node.vertical.recursive(node_df, descendant.right)
      node_df[node_df$node == node,]$node.vertical <- mean(c(node_df[node_df$node == descendant.left, ]$node.vertical, 
                                                             node_df[node_df$node == descendant.right, ]$node.vertical))
      
    } else if (is.na(descendant.left) + is.na(descendant.right) == 1) { # if the node is a character-change node then its vertical position is the same as its descendant's
      
      if (!is.na(descendant.left)) {
        node_df <- get.node.vertical.recursive(node_df, descendant.left)
        node_df[node_df$node == node, ]$node.vertical <- node_df[node_df$node == descendant.left, ]$node.vertical
      } else if (!is.na(descendant.right)) {
        node_df <- get.node.vertical.recursive(node_df, descendant.right)
        node_df[node_df$node == node, ]$node.vertical <- node_df[node_df$node == descendant.right, ]$node.vertical
      }
    }
    
    return(node_df)
  }
  
  node_df$node.vertical <- 0
  Ntip <- sum(is.na(node_df$descendant.left) + is.na(node_df$descendant.right) == 2)
  node <- Ntip + 1L
  
  # start from the root
  return(get.node.vertical.recursive(node_df, node))
}
