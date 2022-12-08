# this script contains functions that parse BEAST tree(s) log files into R (multi)phylo object

#' Read in BEAST tree(s) log files and then parse each extended newick string into a phylo object
#' @param file_path path to the BEAST tree(s) log file
#' @param tree_text text content of the BEAST tree(s) log file
#' @return a (list of) phylo object(s)
readBeast <- function(file_path = NULL, tree_text = NULL) {
  
  # read in the tree file content as vector of strings
  if (!is.null(file_path) && !is.null(tree_text)) {
    message("either file_path or tree_text needs to be provided, but not both; we will only use file_path then.\n")
    tree_text <- NULL
  }
  
  if (!is.null(file_path) && length(file_path) == 1 && file.exists(file_path)) {
    x_raw <- scan(file_path, quiet = T, sep = "\n", what = character())
  } else if (!is.null(tree_text)) {
    x_raw <- tree_text
  } else {
    stop("don't know what to do.\n")
  }
  
  # only retain the relevant tree text
  start_linenum <- grep("^begin trees;$", tolower(x_raw))
  if (length(start_linenum) != 1) {
    stop("cannot find the tree starting line")
  }
  end_linenum <- grep("^end;$", tolower(x_raw))
  end_linenum <- min(end_linenum[end_linenum > start_linenum])
  if (length(end_linenum) != 1) {
    stop("cannot find the tree ending line")
  }
  x_raw <- x_raw[start_linenum:end_linenum]
  
  # new just take the tree lines and parse them
  treeline_num <- grep('\\[&R\\].*?;', x_raw)
  x <- x_raw[treeline_num]
  x <- gsub("^.*\\[&R\\] |^.*\\[&R\\]", "", x)
  my.data <- lapply(x, parseBeastNewick) # parsing
  
  # decide whether tip label translation is needed
  tl_int <- suppressWarnings(as.integer(my.data[[1]]$tip.label))
  if (all(is.na(tl_int))) {
    translate <- F
  } else if (all(!is.na(tl_int))) {
    translate <- T
  } else {
    stop("the tip labels are mixed of pure numbers and characters")
  }
  
  # order all the tips (and of course node.data) the same way 1:tip.num
  if (translate) {
    ntips <- length(my.data[[1]]$tip.label)
    for (i in 1:length(my.data)) {
      tip.num <- as.integer(my.data[[i]]$tip.label)
      my.data[[i]]$edge[, 2][my.data[[i]]$edge[, 2] <= ntips] <- tip.num[my.data[[i]]$edge[, 2][my.data[[i]]$edge[, 2] <= ntips]]
      
      my.data[[i]]$node.data[tip.num] <- my.data[[i]]$node.data[1:ntips]
      my.data[[i]]$tip.label <- 1:ntips
    }
  }
  
  if (!translate) { # Trees have names in them
    for (i in 1:length(my.data)) {
      my.data[[i]]$tip.label <- gsub("\'", "", my.data[[i]]$tip.label)
    }
    
    # Trees use numbers not names, we will use taxon name block of format "number name" to 
    # link the numbers used in the tree file to the names
  } else {
    
    x_raw <- x_raw[1:(treeline_num[1] - 1)]
    first <- grep('translate', tolower(x_raw)) + 1
    last <- grep(';', x_raw)
    last <- min(last[last > first]) - 1
    x <- x_raw[first:last]
    x <- gsub("\t|,|\'", "", x)
    x <- strsplit(x, " ")
    
    tip_ids <- as.integer(sapply(x, "[[", 1))
    tip_names <- sapply(x, "[[", 2)
    
    for (i in 1:length(my.data)) {
      my.data[[i]]$tip.label <- tip_names[match(my.data[[i]]$tip.label, tip_ids)]
    }
  }
  
  for (i in 1:length(my.data)) {
    class(my.data[[i]]) <- "phylo"
  }
  
  if (length(my.data) > 1) {
    return(my.data)
  } else {
    return(my.data[[1]])
  }
}


#' Read in an extended newick string, parse it to convert it into a phylo object
#' @param string the extended newick string
#' @return a phylo object
parseBeastNewick <- function(string) {
  
  # parse around front brackets and before back brackets while preserve them
  string <- gsub("\\(", "++(++", string)
  string <- gsub("\\)", "++)", string)
  string <- gsub("++++", "++", string, fixed = T)
  tree.vec <- strsplit(string, "++", fixed = T)[[1]][-1]
  
  n.tips <- length(grep("\\(", tree.vec)) + 1L # assuming at least two tips
  n.nodes <- as.integer(2 * n.tips - 1) # assuming fully bifurcating rooted tree
  n.edges <- as.integer(2 * n.tips - 2)
  edges <- matrix(0L, ncol = 2, nrow = n.edges, byrow = T)

  max.node <- n.tips
  at.node <- n.tips
  tip.node <- 0L # assuming fully bifurcating
  edges_rownum <- 1L
  
  tip.names <- character(n.tips)
  node_bls <- numeric(n.nodes)
  node_annotations <- character(n.nodes)
  
  # strip out extended newick part
  newicks <- gsub("\\[(.*?)\\]", "", tree.vec)
  newickexts <- gsub("]:[&", ",", tree.vec, fixed = T)
  extnewick_exists <- grepl("\\[&", tree.vec)
  
  is_firstvisits <- newicks == "("
  is_finalvisits <- grepl("\\)", newicks)
  has_subtendingbranchs <- grepl(":", newicks)
  
  for (i in 1:length(tree.vec)) {
    
    newick <- newicks[i]
    newickext <- newickexts[i]
    
    if (is_firstvisits[i]) { # adding a new node we've never seen before, guaranteed to be internal
      
      if (at.node != n.tips) {
        edges[edges_rownum, ] <- c(at.node, max.node + 1L)
        edges_rownum <- edges_rownum + 1L
      }
      max.node <- max.node + 1L
      at.node <- max.node

    } else if (is_finalvisits[i]) { # we're going back through a previously visited internal node
      
      old.node <- at.node
      
      if (has_subtendingbranchs[i]) { # with a subtending branch (i.e., non-root node)
        
        at.node <- edges[edges[, 2] == old.node, 1] # traverse back one node to their parent

        coloncommasplit <- strsplit(newick, ":|,|;")[[1]][-1]
        node_bls[old.node] <- as.numeric(coloncommasplit[1])

        if (length(coloncommasplit) == 3) { # there is a tip that's sister to this internal node, so this internal node is not root
          
          tip.node <- tip.node + 1L
          tip.names[tip.node] <- coloncommasplit[2]
          node_bls[tip.node] <- as.numeric(coloncommasplit[3])
          
          edges[edges_rownum, ] <- c(at.node, tip.node)
          edges_rownum <- edges_rownum + 1L
          
          if (extnewick_exists[i]) { # parse node annotation, the extended newick part, if there exists
            # when there are multiple extended newick parts
            newickext <- strsplit(newickext, ",(?![^[]*])", perl = T)[[1]]

            # first node in the string is the node we just passed through
            # second node in the string is the newly found tip
            for (j in 1:2) {
              if (grepl("\\[&", newickext[j])) {
                node_annotations[c(old.node, tip.node)[j]] <- gsub(".*\\[&|\\].*", "", newickext[j])
              }
            }
          }
          
        } else if (length(coloncommasplit) == 1) { # sister to another internal node or is the root
          
          if (extnewick_exists[i]) {
            node_annotations[old.node] <- gsub(".*\\[&|\\].*", "", newickext)
          }
        } else {
          stop("the format of this internal node string is not what we expect")
        }
      } else { # root
        node_bls[old.node] <- 0
        if (extnewick_exists[i]) {
          node_annotations[old.node] <- gsub(".*\\[&|\\].*", "", newickext)
        }
      }

    } else if (has_subtendingbranchs[i]) { # there must be one or two tips descending from the current node
      
      coloncommasplit <- strsplit(newick, ":|,")[[1]]
      if (!length(coloncommasplit) %in% c(2, 4)) {
        stop("the format of this tip node string is not what we expect")
      }
      
      tip.node_old <- tip.node
      for (j in 1:(length(coloncommasplit) / 2)) {
        tip.node <- tip.node + 1L
        tip.names[tip.node] <- coloncommasplit[j * 2 - 1]
        node_bls[tip.node] <- as.numeric(coloncommasplit[j * 2])
        
        edges[edges_rownum, ] <- c(at.node, tip.node)
        edges_rownum <- edges_rownum + 1L
      }
      
      if (extnewick_exists[i]) { # associated node info exists
        newickext <- strsplit(newickext, ",(?![^[]*])", perl = T)[[1]]
        
        for (j in 1:length(newickext)) {
          if (grepl("\\[&", newickext[j])) {
            node_annotations[tip.node_old + j] <- gsub(".*\\[&|\\].*", "", newickext[j])
          }
        }
      }
      
    } else {
      stop("cannot understand some components of this newick string")
    }
  }
  
  if (any(is.na(node_bls))) {
    stop("some branch lengths cannot be found")
  }
  
  # format the phylo object
  edge.lengths <- numeric(n.edges)
  node_indices <- (1:n.nodes)[-(n.tips + 1)]
  edge.lengths[match(node_indices, edges[, 2])] <- node_bls[node_indices]

  my.tree <- list(edge = edges, tip.label = tip.names, edge.length = edge.lengths, Nnode = n.tips - 1, node.data = node_annotations)
  if (node_bls[n.tips + 1] > 0) {
    my.tree$root.edge <- node_bls[n.tips + 1]
  }
  
  return(my.tree)
}
