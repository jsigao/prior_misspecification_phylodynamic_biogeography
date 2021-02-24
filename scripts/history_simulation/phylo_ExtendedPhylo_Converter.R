source("./extendedPhylo_constructor.R") # extendedPhylo helper functions
if (!require(ape)) {
  install.packages("ape", repos = "https://cloud.r-project.org/")
  library(ape)
}
if (!require(phytools)) {
  install.packages("phytools", repos = "https://cloud.r-project.org/")
  library(phytools)
}

#' Take either a phylo or extendedPhylo object and convert it to the other one as the return
#' @param phylo a phylo or extendedPhylo object
#' @return an extendedPhylo or phylo object (reverse of the input)
phyloExtendedPhyloConverter <- function(phylo) {
  
  if (inherits(phylo, "phylo")) { # convert phylo to extendedPhylo
    
    if (!is.null(phylo$maps)) { # with ancestral state reconstruction and/or full history
      node_number_total <- 1L + length(unlist(phylo$maps))
      if (!is.null(phylo$maps.set)) {
        extendedPhylo <- data.frame(node = 1:node_number_total, ancestor = NA, descendant.left = NA, descendant.right = NA, subtending.edge.length = NA, label = NA, state = NA, state.set = NA, stringsAsFactors = F)
      } else {
        extendedPhylo <- data.frame(node = 1:node_number_total, ancestor = NA, descendant.left = NA, descendant.right = NA, subtending.edge.length = NA, label = NA, state = NA, stringsAsFactors = F)
      }
    } else { # without ancestral state reconstruction or full history
      node_number_total <- 1L + nrow(phylo$edge)
      extendedPhylo <- data.frame(node = 1:node_number_total, ancestor = NA, descendant.left = NA, descendant.right = NA, subtending.edge.length = NA, label = NA, stringsAsFactors = F)
    }
    
    if (!is.null(phylo$maps)) { # the phylo object contains stochatic maps as its member
      
      # ancestor
      extendedPhylo[match(phylo$edge[, 2], extendedPhylo$node), ]$ancestor <- phylo$edge[, 1]
      # if the left descendant is empty then fills it by the end node index of that edge, otherwise the other the right descendant
      extendedPhylo[match(phylo$edge[!duplicated(phylo$edge[, 1]), 1], extendedPhylo$node), ]$descendant.left <- phylo$edge[!duplicated(phylo$edge[, 1]), 2]
      extendedPhylo[match(phylo$edge[duplicated(phylo$edge[, 1]), 1], extendedPhylo$node), ]$descendant.right <- phylo$edge[duplicated(phylo$edge[, 1]), 2]
      # edge length
      extendedPhylo[match(phylo$edge[, 2], extendedPhylo$node), ]$subtending.edge.length <- phylo$edge.length
      
      # node state info
      extendedPhylo[match(phylo$edge[, 2], extendedPhylo$node), ]$state <- sapply(phylo$maps, function(x) names(x)[length(x)])
      if (!is.null(phylo$maps.set)) {
        extendedPhylo[match(phylo$edge[, 2], extendedPhylo$node), ]$state.set <- lapply(phylo$maps.set, function(x) x[length(x)])
      }
      
      root <- unique(phylo$edge[, 1][!phylo$edge[, 1] %in% phylo$edge[, 2]])
      if (length(root) != 1) {
        stop("there should be one and only one root.\n")
      }
      root_edge <- which(!phylo$edge[, 1] %in% phylo$edge[, 2])[1]
      
      # root state
      extendedPhylo[extendedPhylo$node == root, ]$state <- names(phylo$maps[[root_edge]])[1]
      if (!is.null(phylo$maps.set)) {
        extendedPhylo[extendedPhylo$node == root, ]$state.set <- phylo$maps.set[[root_edge]][1]
      }
      
      # index of the first character-change node should be 2*n
      intermediate_node <- 2L + nrow(phylo$edge)
      
      # loop over all the branches to get the character-change events
      for (i in 1:nrow(phylo$edge)) {
          
        if (length(phylo$maps[[i]]) >= 2) { # when there is at least one character-change event on this branch
          
          # first segment
          # if the left descendant of this tree node is empty then fills it by the character-change node index, otherwise the other the right descendant
          if (extendedPhylo[extendedPhylo$node == phylo$edge[i, 1], ]$descendant.left == phylo$edge[i, 2]) {
            left <- T
          } else if (extendedPhylo[extendedPhylo$node == phylo$edge[i, 1], ]$descendant.right == phylo$edge[i, 2]) {
            left <- F
          } else {
            stop("neither left nor right descendant.\n")
          }
          if (left) {
            extendedPhylo[extendedPhylo$node == phylo$edge[i, 1], ]$descendant.left <- intermediate_node
          } else {
            extendedPhylo[extendedPhylo$node == phylo$edge[i, 1], ]$descendant.right <- intermediate_node
          }
          
          # ancestor
          extendedPhylo[extendedPhylo$node == intermediate_node, ]$ancestor <- phylo$edge[i, 1]
          
          # state
          # if the node is a root node then there is no subtending branch so start and end state the same
          if (i == root_edge) {
            extendedPhylo[extendedPhylo$node == phylo$edge[i, 1], ]$state <- names(phylo$maps[[i]])[1]
            if (!is.null(phylo$maps.set)) {
              extendedPhylo[extendedPhylo$node == phylo$edge[i, 1], ]$state.set <- phylo$maps.set[[i]][1]
            }
          }
          
          extendedPhylo[extendedPhylo$node == intermediate_node, ]$state <- names(phylo$maps[[i]])[1]
          if (!is.null(phylo$maps.set)) {
            extendedPhylo[extendedPhylo$node == intermediate_node, ]$state.set <- phylo$maps.set[[i]][1]
          }
          
          # edge length
          extendedPhylo[extendedPhylo$node == intermediate_node, ]$subtending.edge.length <- as.numeric(phylo$maps[[i]][1])
          
          if (length(phylo$maps[[i]]) > 2) { # when there is more than one character-change events on this branch
            
            # second to the last but one segment
            for (j in 2:(length(phylo$maps[[i]]) - 1)) {
              
              # there should only be one descendant of this character-change node; let's make it the left one
              if (is.na(extendedPhylo[extendedPhylo$node == intermediate_node, ]$descendant.left)) {
                extendedPhylo[extendedPhylo$node == intermediate_node, ]$descendant.left <- intermediate_node + 1L
              } else {
                stop("there should only be one descendant of this character-change node.\n")
              }
              
              # ancestor
              extendedPhylo[extendedPhylo$node == intermediate_node + 1L, ]$ancestor <- intermediate_node
              
              # state
              extendedPhylo[extendedPhylo$node == intermediate_node + 1L, ]$state <- names(phylo$maps[[i]])[j]
              if (!is.null(phylo$maps.set)) {
                extendedPhylo[extendedPhylo$node == intermediate_node + 1L, ]$state.set <- phylo$maps.set[[i]][j]
              }
              
              # edge length
              extendedPhylo[extendedPhylo$node == intermediate_node + 1L, ]$subtending.edge.length <- as.numeric(phylo$maps[[i]][j])
              
              # increment the character-change node index
              intermediate_node <- intermediate_node + 1L
            }
          }
          
          # the last segment
          # there should only be one descendant of this character-change node, which should be the end node of that branch; let's make it the left one
          if (is.na(extendedPhylo[extendedPhylo$node == intermediate_node, ]$descendant.left)) {
            extendedPhylo[extendedPhylo$node == intermediate_node, ]$descendant.left <- phylo$edge[i, 2]
          } else {
            stop("there should only be one descendant of this character-change node.\n")
          }
          
          # ancestor
          extendedPhylo[extendedPhylo$node == phylo$edge[i, 2], ]$ancestor <- intermediate_node
          
          # state
          extendedPhylo[extendedPhylo$node == phylo$edge[i, 2], ]$state <- names(phylo$maps[[i]])[length(phylo$maps[[i]])]
          if (!is.null(phylo$maps.set)) {
            extendedPhylo[extendedPhylo$node == phylo$edge[i, 2], ]$state.set <- phylo$maps.set[[i]][length(phylo$maps.set[[i]])]
          }
          
          # edge length
          extendedPhylo[extendedPhylo$node == phylo$edge[i, 2], ]$subtending.edge.length <- as.numeric(phylo$maps[[i]][length(phylo$maps[[i]])])
          
          # increment the character-change node index
          intermediate_node <- intermediate_node + 1L
        }
      }
      
      if (!is.null(phylo$maps.set)) {
        extendedPhylo$edge.state.set <- extendedPhylo$state.set
      }
      
    } else { # the phylo object does not contain stochatic maps as its member, so simple phylo
      
      # ancestor
      extendedPhylo[match(phylo$edge[, 2], extendedPhylo$node), ]$ancestor <- phylo$edge[, 1]
      # if the left descendant is empty then fills it by the end node index of that edge, otherwise the other the right descendant
      extendedPhylo[match(phylo$edge[!duplicated(phylo$edge[, 1]), 1], extendedPhylo$node), ]$descendant.left <- phylo$edge[!duplicated(phylo$edge[, 1]), 2]
      extendedPhylo[match(phylo$edge[duplicated(phylo$edge[, 1]), 1], extendedPhylo$node), ]$descendant.right <- phylo$edge[duplicated(phylo$edge[, 1]), 2]
      # edge length
      extendedPhylo[match(phylo$edge[, 2], extendedPhylo$node), ]$subtending.edge.length <- phylo$edge.length
    }
    
    tips <- phylo$edge[, 2][!phylo$edge[, 2] %in% phylo$edge[, 1]]
    extendedPhylo[match(tips, extendedPhylo$node), ]$label <- phylo$tip.label[tips]
    
    # compute the clade size of each node (the number of its tree nodes descendants, including itself)
    # which would be used to order/ladderize the tree in get.new.tipOrder function
    extendedPhylo <- get.clade.size(extendedPhylo)
    
    # ladderize the tree by generating new tip orders (the default is decreasing thus smaller clade on top)
    extendedPhylo <- get.new.tipOrder(extendedPhylo)
    
    # recursivly traversing the tree to compute the height of each node
    extendedPhylo <- get.node.heights(extendedPhylo)
    
    # recursivly traversing the tree to comoute the vertical position of each node
    extendedPhylo <- get.node.vertical(extendedPhylo)
    
    return(extendedPhylo)
    
  } else { # convert extendedPhylo to phylo
    
    # convert to simmap
    # if there is any node with only one descendant (then it's a character-change node) or the state column is not empty
    # then this extendedPhylo object contains stochastic maps
    if (any(is.na(phylo$descendant.left) + is.na(phylo$descendant.right) == 1) || ((!is.null(phylo$state)) && any(!is.na(phylo$state)))) {
      
      num_tips <- sum(is.na(phylo$descendant.left) + is.na(phylo$descendant.right) == 2)
      
      vanillaPhylo_edge <- matrix(0, nrow = 2 * num_tips - 2, ncol = 2)
      vanillaPhylo_edge_length <- numeric(2 * num_tips - 2)
      vanillaPhylo_maps <- vector("list", 2 * num_tips - 2)
      vanillaPhylo_rownum <- 1
      
      for (i in 1:nrow(phylo)) {
        
        # if this node is a tree node
        if (is.na(phylo[i, ]$descendant.left) + is.na(phylo[i,]$descendant.right) != 1) {
          
          current <- phylo[i, ]$node
          
          # if there is an ancestor then it's not a root node
          if (!is.na(phylo[i, ]$ancestor)) {
            
            duration_vec <- NULL
            state_vec <- NULL
            
            # traverse along that branch until reach another tree node
            repeat{
              
              duration_vec <- c(phylo[phylo$node == current,]$subtending.edge.length, duration_vec)
              state_vec <- c(phylo[phylo$node == current,]$state, state_vec)
              
              ancestor <- phylo[phylo$node == current, ]$ancestor
              current <- ancestor
              
              if (is.na(phylo[phylo$node == ancestor, ]$descendant.left) + is.na(phylo[phylo$node == ancestor,]$descendant.right) != 1) break
              
            }
            
            vanillaPhylo_edge[vanillaPhylo_rownum, ] <- c(ancestor, phylo[i, ]$node)
            vanillaPhylo_edge_length[vanillaPhylo_rownum] <- sum(duration_vec)
            vanillaPhylo_maps[[vanillaPhylo_rownum]] <- duration_vec
            names(vanillaPhylo_maps[[vanillaPhylo_rownum]]) <- state_vec
            
            # increment row index
            vanillaPhylo_rownum <- vanillaPhylo_rownum + 1
          }
        }
        
      }
      
      # tip label
      tip_label <- phylo[match((1:num_tips), phylo$node), ]$label
      
      # internal node states
      internal_node_states <- t(sapply(vanillaPhylo_maps, function(x) c(names(x)[1], names(x)[length(x)])))
      
      # tip states
      tip_states <- phylo[match((1:num_tips), phylo$node), ]$state
      
      # mapped.edge
      mapped_edge <- matrix(data = 0, nrow(vanillaPhylo_edge), length(unique(phylo$state)), 
                            dimnames = list(apply(vanillaPhylo_edge, 1, function(x) paste(x, collapse = ",")), state = sort(unique(phylo$state))))
      for (j in 1:length(vanillaPhylo_maps)) {
        for (k in 1:length(vanillaPhylo_maps[[j]])) {
          mapped_edge[j, names(vanillaPhylo_maps[[j]])[k]] <- mapped_edge[j, names(vanillaPhylo_maps[[j]])[k]] + vanillaPhylo_maps[[j]][k]
        }
      }
      
      # make the tree object
      vanillaPhylo <- list(edge        = vanillaPhylo_edge,
                           tip.label   = tip_label,
                           edge.length = vanillaPhylo_edge_length,
                           Nnode       = num_tips - 1,
                           maps        = vanillaPhylo_maps,
                           node.states = internal_node_states,
                           states      = tip_states,
                           mapped.edge = mapped_edge)
      
      class(vanillaPhylo) <- "phylo"
      vanillaPhylo <- ladderize.simmap(vanillaPhylo, right = FALSE)
      class(vanillaPhylo) <- c("simmap", class(vanillaPhylo))
      
      return(vanillaPhylo)
      
    } else { # convert to plain phylo
      
      num_tips <- sum(is.na(phylo$descendant.left) + is.na(phylo$descendant.right) == 2)
      
      vanillaPhylo_edge <- matrix(0, nrow = 2 * num_tips - 2, ncol = 2)
      vanillaPhylo_edge_length <- numeric(2 * num_tips - 2)
      vanillaPhylo_rownum <- 1
      
      for (i in 1:nrow(phylo)) {
        if (!is.na(phylo[i,]$ancestor)) {
          vanillaPhylo_edge[vanillaPhylo_rownum,] <- c(phylo[i,]$ancestor, phylo[i,]$node)
          vanillaPhylo_edge_length[vanillaPhylo_rownum] <- phylo[i,]$subtending.edge.length
          vanillaPhylo_rownum <- vanillaPhylo_rownum + 1
        }
      }
      
      # tip label
      tip_label <- phylo[match((1:num_tips), phylo$node), ]$label
      
      # make the tree object
      vanillaPhylo <- list(edge        = vanillaPhylo_edge,
                           tip.label   = tip_label,
                           edge.length = vanillaPhylo_edge_length,
                           Nnode       = num_tips - 1)
      class(vanillaPhylo) <- c("phylo")
      
      vanillaPhylo <- ladderize(vanillaPhylo, right = FALSE)
      
      # todo: add root edge if that's not NA in the extendedPhylo data frame;
      # take care of the case when the root node has its associated property (root edge length) -> phylo$root.edge != 0
      return(vanillaPhylo)
    }
  }
}
