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
    
    nedges <- nrow(phylo$edge)
    nnodes <- nedges + 1L
    ancestors <- phylo$edge[, 1]
    descendants <- phylo$edge[, 2]
    has_maps <- "maps" %in% names(phylo)
    has_maps.set <- "maps.set" %in% names(phylo)
    
    if (has_maps) { # with ancestral state reconstruction and/or full history
      npieces_all <- lengths(phylo$maps)
      node_number_total <- 1L + sum(npieces_all)
      state <- rep(NA_character_, node_number_total)
      if (has_maps.set) {
        state.set <- rep(NA_character_, node_number_total)
      }
    } else {
      node_number_total <- 1L + nedges
    }
    
    node <- 1:node_number_total
    ancestor <- descendant.left <- descendant.right <- rep(NA_integer_, node_number_total)
    subtending.edge.length <- rep(NA_real_, node_number_total)
    label <- rep(NA_character_, node_number_total)
    
    ancestor[descendants] <- ancestors
    subtending.edge.length[descendants] <- phylo$edge.length
    
    second_appearance <- duplicated(ancestors)
    first_appearance <- !second_appearance
    descendant.left[ancestors[first_appearance]] <- descendants[first_appearance]
    descendant.right[ancestors[second_appearance]] <- descendants[second_appearance]
    
    tips <- descendants[!descendants %in% ancestors]
    label[tips] <- phylo$tip.label[tips]
    
    if (has_maps) { # the phylo object contains stochatic maps as its member
      # node state info
      state[descendants] <- vapply(phylo$maps, function(x) names(x)[length(x)], FUN.VALUE = character(1L))
      if (has_maps.set) {
        state.set[descendants] <- lapply(phylo$maps.set, function(x) x[length(x)])
      }
      
      # root state
      # if the node is a root node then there is no subtending branch so start and end state the same
      root_edges <- which(!ancestors %in% descendants)
      root <- unique(phylo$edge[root_edges, 1])
      if (length(root) != 1) {
        stop("there should be one and only one root.\n")
      }
      root_edge <- root_edges[1]
      
      state[root] <- names(phylo$maps[[root_edge]])[1]
      if (!is.null(phylo$maps.set)) {
        state.set[root] <- phylo$maps.set[[root_edge]][1]
      }
      
      # loop over all the branches that have at least one character-change event to get the events
      event_edge_idx <- which(npieces_all > 1)
      events_edge_bools <- npieces_all > 2
      if (length(event_edge_idx) > 0) {
        pieces_all <- unlist(phylo$maps[event_edge_idx])
        pieces_all_state <- names(pieces_all)
        pieces_all_length <- as.numeric(pieces_all)
        piece_id <- 0L
        
        # index of the first character-change node should be 2*n
        intermediate_node <- nnodes + 1L
        
        for (i in event_edge_idx) {
          
          piece_id <- piece_id + 1L
          
          # first segment
          # if the left descendant of this tree node is empty then fills it by the character-change node index, otherwise the right descendant
          if (descendant.left[ancestors[i]] == descendants[i]) {
            descendant.left[ancestors[i]] <- intermediate_node
          } else {
            descendant.right[ancestors[i]] <- intermediate_node
          }
          
          # ancestor
          ancestor[intermediate_node] <- ancestors[i]
          # edge length
          subtending.edge.length[intermediate_node] <- pieces_all_length[piece_id]
          # state
          state[intermediate_node] <- pieces_all_state[piece_id]
          if (has_maps.set) {
            state.set[intermediate_node] <- phylo$maps.set[[i]][1]
          }
          
          if (events_edge_bools[i]) { # when there is more than one character-change events on this branch
            
            # second to the last but one segment
            for (j in 2:(npieces_all[i] - 1L)) {
              
              intermediate_node_next <- intermediate_node + 1L
              piece_id <- piece_id + 1L
              
              # there should only be one descendant of this character-change node; let's make it the left one
              descendant.left[intermediate_node] <- intermediate_node_next
              
              # ancestor
              ancestor[intermediate_node_next] <- intermediate_node
              # edge length
              subtending.edge.length[intermediate_node_next] <- pieces_all_length[piece_id]
              # state
              state[intermediate_node_next] <- pieces_all_state[piece_id]
              if (has_maps.set) {
                state.set[intermediate_node_next] <- phylo$maps.set[[i]][j]
              }
              
              # increment the character-change node index
              intermediate_node <- intermediate_node_next
            }
          }
          
          piece_id <- piece_id + 1L
          
          # the last segment
          # there should only be one descendant of this character-change node, which should be the end node of that branch; let's make it the left one
          descendant.left[intermediate_node] <- descendants[i]
          
          # ancestor
          ancestor[descendants[i]] <- intermediate_node
          # edge length
          subtending.edge.length[descendants[i]] <- pieces_all_length[piece_id]
          # state
          state[descendants[i]] <- pieces_all_state[piece_id]
          if (has_maps.set) {
            state.set[descendants[i]] <- phylo$maps.set[[i]][length(phylo$maps.set[[i]])]
          }
          
          # increment the character-change node index
          intermediate_node <- intermediate_node + 1L
        }
      }
    }
    
    extendedPhylo <- data.frame(node = node, ancestor = ancestor, descendant.left = descendant.left, descendant.right = descendant.right, 
                                subtending.edge.length = subtending.edge.length, label = label, stringsAsFactors = F)
    if (has_maps) { # with ancestral state reconstruction and/or full history
      extendedPhylo$state <- state
      if (has_maps.set) {
        extendedPhylo$edge.state.set <- extendedPhylo$state.set <- state.set
      }
    }
    
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
        if (is.na(phylo[i, ]$descendant.left) + is.na(phylo[i, ]$descendant.right) != 1) {
          
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
