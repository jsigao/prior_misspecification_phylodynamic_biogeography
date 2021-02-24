# this script contains functions that simulate full history over a bifurcating phylogeny
# under a constant or piecewise constant geographical model
# specifically, it performs stochastic mapping (conditioning on the tree, node states, and geographic model parameter estimates)
# and posterior-predictive simulation (conditioning on the tree and geographic model parameter estimates)
library(stringr)
library(pbapply)
source("./history_simulator.R") # history simulation
source("./BeastTree_parser.R") # tree parsing
source("./BeastTree_ExtendedPhylo_Converter.R") # node info parsing

#' Compute the least common multiple of two integers
#' @param x A integer
#' @param y A integer
#' @return LCM of \code{x} and \code{y}.
lcm <- function(x, y) {
  # choose the greater number
  if (x > y) {
    greater <- x
    increment <- x
  } else {
    greater <- y
    increment <- y
  }
  repeat {
    if ((greater %% x == 0) && (greater %% y == 0)) {
      lcm <- greater
      break
    }
    greater <- greater + increment
  }
  return(lcm)
}

#' Simulate a single history of discrete-character changes over a bifurcation tree conditional or unconditional on the tip states
#' @param tree a bifurcation tree of class "phylo"
#' @param Q a instantaneous-rate matrix (or a list of matrices for a piecewise constant geographic model) characterizes the CTMC
#' @param Q_ages boundaries of time intervals for a piecewise constant geographic model (NULL means a constant model)
#' @param root_freq the vector state frequencies at the root of the tree (that we will draw the root state from); if NULL then it's uniform
#' @param numstates_conditioned whether to require the number of tip states to reach the number of states in the geographic model when perform forward simulation (i.e., not stochastic mapping)
#' @param nrejections_max number of rejections (due to not realizing all the tip states in the simulated history, so only will be used when numstates_conditioned is true) to reach before skip that round of simulation
#' @param conditional whether condition on the observed states at the tip (i.e., stochastic mapping) or not (i.e., forward simulation)
#' @return A phylo and simmap object (or a multiPhylo and multiSimmap object when nsim > 1) that contains the simulated full history
history_sim <- function(tree, Q, Q_ages = NULL, root_freq = NULL, states, numstates_conditioned = T, nrejections_max = 100L, conditional = F) {
  
  history <- sim_history(tree = tree, Q = Q, Q_ages = Q_ages, root_freq = root_freq, nsim = 1L, conditional = conditional)

  # require the simulated history to have the same number of states as in the observation
  # (only applicable to posterior-predictive simulation but not stochastic mapping)
  if (numstates_conditioned && !conditional) {
    num_rejections <- 0
    while (length(unique(history$states)) != length(states) && num_rejections <= nrejections_max) {
      num_rejections <- num_rejections + 1
      history <- sim_history(tree = tree, Q = Q, Q_ages = Q_ages, root_freq = root_freq, nsim = 1L)
    }
    if (num_rejections > nrejections_max) {
      history <- NA
      cat("too many rejections, skip this sample.\n")
    }
  }
  
  return(history)
}

#' Simulate a number of histories of discrete-character changes over a list of bifurcation trees conditioning on the tip states or not
#' @param trees a list of bifurcation tree of class "phylo"
#' @param Qs a list of instantaneous-rate matrices (or a list of lists of matrices for a piecewise constant geographic model) characterizes the CTMC
#' @param Q_ages boundaries of time intervals for a piecewise constant geographic model (NULL means a constant model)
#' @param root_freqs a matrix of state frequencies (each row corresponds to a generation) at the root of the tree (that we will draw the root state from); if NULL then it's uniform
#' @param numstates_conditioned whether to require the number of tip states to reach the number of states in the geographic model when perform forward simulation (i.e., not stochastic mapping)
#' @param nrejections_max number of rejections (due to not realizing all the tip states in the simulated history, 
#' so only will be used when numstates_conditioned is true) to reach before skip that round of simulation
#' @param indices generation indices that will only used for printing progress to screen
#' @param conditional whether condition on the observed states at the tip (i.e., stochastic mapping) or not (i.e., forward simulation)
#' @return a list of phylo and simmap objects (or a multiPhylo and multiSimmap object when nsim > 1) that each contains the simulated full history
histories_sim <- function(trees, Qs, Q_ages = NULL, root_freqs = NULL, states, numstates_conditioned = T, nrejections_max = 100L, indices = NULL, conditional = F) {
  
  nhis <- length(trees)
  if (nhis != length(Qs)) {
    stop("number of iterations don't match.\n")
  }
  if ((!is.null(root_freqs)) && nhis != nrow(root_freqs)) {
    stop("number of iterations don't match.\n")
  }
  
  his_all <- vector("list", nhis)
  verbose <- !is.null(indices)
  for (l in 1:nhis) {
    if (verbose) {
      if (l %% 10 == 0) {
        cat(paste0("simulated history no.", indices[l], ".\n"))
      }
    }
    his_all[[l]] <- history_sim(tree = trees[[l]], Q = Qs[[l]], Q_ages = Q_ages, root_freq = unlist(root_freqs[l, ]), states = states,
                                numstates_conditioned = numstates_conditioned, nrejections_max = nrejections_max, conditional = conditional)
  }
  
  return(his_all)
}

#' Simulate history of discrete-character changes over a bifurcation tree conditioning on the node states or not
#' @param file_path path to the BEAST log file that contains all the parameter estimates
#' @param numstates_conditioned whether to require the number of tip states to reach the number of states in the geographic model when perform forward simulation (i.e., not stochastic mapping) 
#' @param ncores number of computer cores to use (if more than one then simulations will be parallelized)
#' @param nrejections_max number of rejections (due to not realizing all the tip states in the simulated history, so only will be used when numstates_conditioned is true) to reach before skip that round of simulation
#' @param conditional whether condition on the observed states at the tip (i.e., stochastic mapping) or not (i.e., forward simulation)
#' @param folder_path path to the directory that contains the BEAST trees file output (associated with the provided log file)
#' @param file_paths paths to all the log files in the \code{folder_path} to fetch the tree file associated with the provided \code{file_path}
#' @param sample_size number of full histories to simulate
#' @return A multiPhylo and multiSimmap object that contains the simulated full histories
history_simulator <- function(file_path, numstates_conditioned = T, ncores = 1L, nrejections_max = 100L, conditional = F, 
                              folder_path = NULL, file_paths = NULL, sample_size = NULL) {
  
  cat(file_path, sep = "\n")
  
  #####################################
  # read in parameter values and data #
  #####################################
  # read in discrete trait information
  xml_path <- grep(gsub("_combined\\.log$|\\.log$", "", basename(file_path)), 
                   grep("MLE", list.files(dirname(file_path), recursive = T, full.names = T, pattern = "*.xml$"), value = T, invert = T), value = T)
  if (length(xml_path) == 0) {
    stop("cannot find xml")
  } else if (length(xml_path) > 1 && length(unique(gsub("_run(.*?)\\.xml$", "", basename(xml_path)))) != 1) {
    stop("cannot find xml")
  } else {
    xml_path <- xml_path[1]
  }
  x <- scan(file = xml_path, what = character(), sep = "\n", strip.white = F, blank.lines.skip = F)
  states <- gsub("^\\s+|\\s+$", "", gsub("<state code|=|\t|\"|/>", "", x[grep("<state code", x)]))
  K <- length(states)
  
  tree_fixed <- F
  if (any(grepl("<newick id=", x))) {
    tree_fixed <- T
  }
  
  # read in log file to get the Q matrix parameter estimates
  log_dat <- read.table(file_path, header = T, sep = "\t", check.names = F, stringsAsFactors = F)
  
  ######################################
  # determine number of time intervals #
  ######################################
  # here we are assume the number of epochal rates is the same as the number of epocal matrices
  epoch_num <- length(grep("*\\.clock\\.rate.epoch\\d*$|*\\.clock\\.rate$", colnames(log_dat)))
  Q_ages <- NULL
  if (epoch_num > 1) {
    bounds_age_linenums <- grep("epoch transitionTime", x)[grep("epoch transitionTime", x) > grep("epochBranchModel id=", x) & 
                                                             grep("epoch transitionTime", x) < grep("</epochBranchModel>", x)]
    if (length(bounds_age_linenums) != epoch_num - 1L) {
      stop("incorrect number of bounds")
    }
    Q_ages <- as.numeric(str_extract(x[bounds_age_linenums], "\\d+"))
  }
  
  # deciding how many samples to simulate
  # here considering the running time, we sample a given number of samples from the approximated joint posterior distribution
  if (is.null(sample_size)) {
    if ((length(states) > 10 && numstates_conditioned) || conditional) {
      sample_size <- 1000L
    } else {
      sample_size <- 2500L
    }
  }
  
  ##############################################
  # read in the trees to simulate histories on #
  ##############################################
  tree_path <- gsub("\\.log$", ".trees", file_path)
  currenttree_colnum <- grep("Current Tree", colnames(log_dat))
  
  # depending on the analyses we did, there could be three scenarios on the tree input
  # 1. if a single tree was fixed in the geographic analysis, we need to read that tree in
  # 2. if a distribution of trees was sampled over in the geographic analysis, we need to read in that distribution 
  # as well as the log file to fetch the specific tree that was sampled at each generation
  # 3. if the tree was inferred jointly in the geographic analysis, we need to read in the inferred distribution of trees
  # for the last two scenarios, an extra misc. complication is that the log file and tree file might have been sampled at different frequency
  # so we would need to match the sampled generations between them
  # one other variation to this combinatorics is that if we are to perform stochastic mapping (instead of posterior-prective simulation), 
  # no matter which of the above three scenarios was applicable, we will treat it as scenario 3 as we need to know the sampled node states 
  # (which is likely to only have been logged in the tree file)
  if (conditional || (!tree_fixed && length(currenttree_colnum) == 0 && file.exists(tree_path))) { # scenario 3
    
    tree_gens_combined <- as.integer(gsub("STATE_", "", system(paste("grep 'tree STATE_'", tree_path, "| cut -f2 -d' '"), intern = T)))
    if (identical(tree_gens_combined, log_dat[, 1])) { # tree file and log file were sampled at the same frequency
      
      # subsample log file according to the specified number of simulations
      if (sample_size > nrow(log_dat)) sample_size <- nrow(log_dat)
      sample_indices <- sample(1:nrow(log_dat), size = sample_size)
      log_dat <- log_dat[sample_indices, ]
      
      if (conditional) { # conditioning on the node states
        
        trees <- NULL
        failed <- F
        if (!is.null(folder_path) && dir.exists(folder_path)) { # if folder_path is provided then we assumes that an rds file summarizing the tree already exists
          trees_augmented_rdspath <- list.files(folder_path, full.names = T, recursive = T, pattern = "trees_augmented_simmap.*\\.rds$")

          if (length(trees_augmented_rdspath) == 1) { # one augmented (i.e., with node states fetched) tree rds file found
            
            trees <- readRDS(trees_augmented_rdspath)
            if (!is.null(file_paths)) {
              if (length(trees) != length(file_paths)) {
                message("number of elements in trees_augmented list doesn't match number of tree files.\n")
                failed <- T
              } else {
                trees <- trees[[match(file_path, file_paths)]] # there could be multiple tree rds files in this folder, we use the provided file_paths to fetch the matched one
              }
            } else if (length(trees) == 1) {
              trees <- trees[[1]]
            } else {
              message("multiple elements in trees_augmented list so don't know which one to take.\n")
              failed <- T
            }
            
            trees <- trees[sample_indices]
            
          } else { # no augmented (i.e., with node states fetched) tree rds file found
            
            trees_rdspath <- grep("augmented|simmap|dataframe", list.files(folder_path, full.names = T, recursive = T, pattern = "trees_.*\\.rds$"), value = T, invert = T)
            if (length(trees_rdspath) == 1) { # one tree (i.e., no known node states) rds file found
              trees <- readRDS(trees_rdspath)
              
              if (!is.null(file_paths)) {
                if (length(trees) != length(file_paths)) {
                  message("number of elements in trees list doesn't match number of tree files.\n")
                  failed <- T
                } else {
                  trees <- trees[[match(file_path, file_paths)]]
                }
              } else if (length(trees) == 1) {
                trees <- trees[[1]]
              } else {
                message("multiple elements in trees_augmented list so don't know which one to take.\n")
                failed <- T
              }
              
              trees <- trees[sample_indices]
              trees <- lapply(trees, function(x) BeastTreeExtendedPhyloConverter(x, convert = F)) # parse the node.data component of each tree to get the node states
            }
          }
        }
        
        if (failed || is.null(trees)) { # when no tree rds can be found, we need to parse the raw BEAST tree output file
          
          # subsample trees
          # first get the meta (including the tip label) info above the extended newick strings
          tree_other_txt <- scan(tree_path, quiet = T, sep = "\n", what = character(), blank.lines.skip = F)
          tree_linenum <- grep("tree STATE_", tree_other_txt)
          tree_head_txt <- ""
          if (tree_linenum[1] >= 2) {
            tree_head_txt <- tree_other_txt[1:(tree_linenum[1] - 1)]
          }
          tree_tail_txt <- ""
          if (length(tree_other_txt) - tree_linenum[length(tree_linenum)] >= 1) {
            tree_tail_txt <- tree_other_txt[(tree_linenum[length(tree_linenum)] + 1):length(tree_other_txt)]
          }
          
          tree_txt <- system(paste("grep 'tree STATE_'", tree_path), intern = T)
          tree_txt <- c(tree_head_txt, tree_txt[sample_indices], tree_tail_txt) # subsample here to reduce the amount of parsing needed
          trees <- readBeast(tree_text = tree_txt) # parse the extended newick strings to get trees (of class phylo) and their associated node.data
          rm(tree_other_txt, tree_txt)
          trees <- lapply(trees, function(x) BeastTreeExtendedPhyloConverter(x, convert = F)) # parse the node.data component of each tree to get the node states
        }
        
      } else { # not conditioning on the node states, then we simply read the trees in as newick strings (i.e., node info can be discarded)
        trees <- read.nexus(tree_path)
        trees <- trees[sample_indices]
      }
      
    } else { # tree file and log file were sampled at different frequencies, then we need to subsample them to only keep their shared generations
      
      # find the raw BEAST tree files
      tree_paths <- grep(gsub("_combined.*", "", basename(tree_path)), 
                         grep("powerposterior|prior|combined", list.files(dirname(tree_path), full.names = T, recursive = T, pattern = ".trees$"), 
                              value = T, invert = T), value = T)
      tree_paths <- tree_paths[sapply(tree_paths, function(x) grepl("^End;$", system(paste("tail -n 1", x), intern = T)))]
      nreps <- length(tree_paths) # number of replicated runs
      
      tree_gens <- as.integer(gsub("STATE_", "", system(paste("grep 'tree STATE_'", tree_paths[1], "| cut -f2 -d' '"), intern = T)))
      tree_logby <- diff(tree_gens[1:2])
      
      # sampling frequency in the combined log file
      log_logby_combined <- diff(log_dat[1:2, 1])
      log_nrow_combined <- nrow(log_dat)
      log_firstgen_combined <- log_dat[1, 1]
      
      # sampling frequency in each replicated log file
      log_paths <- gsub("\\.trees$", ".log", tree_paths)
      log_dats <- lapply(log_paths, function(x) read.table(x, header = T, sep = "\t", check.names = F, stringsAsFactors = F))
      log_logby <- diff(log_dats[[1]][1:2, 1])
      log_nrow <- nrow(log_dats[[1]])
      log_firstgen <- log_dats[[1]][1, 1]
      
      # computing the burnin generations (total gens minus the gens remained) for the log files
      log_burnin_nrow <- log_nrow - ((log_nrow_combined/nreps - 1) * log_logby_combined / log_logby + 1)
      gens_remain <- log_dats[[1]][(log_burnin_nrow + 1):log_nrow, 1]
      
      # get the LCM sampling frequency between the log file and the tree file
      sample_freq <- lcm(log_logby, tree_logby)
      gens_remain <- gens_remain[gens_remain %% sample_freq == 0]
      for (j in 1:nreps) {
        log_dats[[j]] <- log_dats[[j]][log_dats[[j]][, 1] %in% gens_remain, ]
      }
      
      # combining the subsampled log files and the tree files
      log_dat <- do.call(rbind, log_dats)
      rm(log_dats)
      
      tree_remain <- tree_gens %in% gens_remain
      tree_txt <- lapply(tree_paths, function(x) system(paste("grep 'tree STATE_'", x), intern = T)[tree_remain])
      tree_txt <- do.call(c, tree_txt)
      
      if (nrow(log_dat) != length(tree_txt)) {
        stop("log and tree generations don't match.\n")
      }
      
      # subsample the combined log again according to the specified number of simulations
      if (sample_size > nrow(log_dat)) sample_size <- nrow(log_dat)
      if (sample_size == nrow(log_dat)) {
        sample_indices <- 1:sample_size
      } else {
        sample_indices <- sample(1:nrow(log_dat), size = sample_size)
      }
      log_dat <- log_dat[sample_indices, ]
      
      # subsample the combined trees again according to the specified number of simulations
      tree_other_txt <- scan(tree_paths[1], quiet = T, sep = "\n", what = character(), blank.lines.skip = F)
      tree_linenum <- grep("tree STATE_", tree_other_txt)
      tree_head_txt <- ""
      if (tree_linenum[1] >= 2) {
        tree_head_txt <- tree_other_txt[1:(tree_linenum[1] - 1)]
      }
      tree_tail_txt <- ""
      if (length(tree_other_txt) - tree_linenum[length(tree_linenum)] >= 1) {
        tree_tail_txt <- tree_other_txt[(tree_linenum[length(tree_linenum)] + 1):length(tree_other_txt)]
      }
      
      tree_txt <- c(tree_head_txt, tree_txt[sample_indices], tree_tail_txt)
      
      if (conditional) { # fetch the node states when conditioning on them
        trees <- readBeast(tree_text = tree_txt)
        trees <- lapply(trees, function(x) BeastTreeExtendedPhyloConverter(x, convert = F))
      } else { # simply read the tree topology and branch lengths in when not conditioning on the node states
        tree_path_tmp <- tempfile(pattern = "combined", fileext = ".trees")
        cat(tree_txt, file = tree_path_tmp, sep = "\n")
        trees <- read.nexus(tree_path_tmp)
        file.remove(tree_path_tmp)
      }
      
      rm(tree_other_txt, tree_txt)
    }

    tree_indices <- 1:sample_size # as we already did the subsampling
    
  } else { # scenarios 1 or 2
    
    if (sample_size > nrow(log_dat)) sample_size <- nrow(log_dat)
    if (sample_size == nrow(log_dat)) {
      sample_indices <- 1:sample_size
    } else {
      sample_indices <- sample(1:nrow(log_dat), size = sample_size)
    }
    log_dat <- log_dat[sample_indices, ] # subsample the log file according to the specified number of simulations
    
    if (length(currenttree_colnum) == 1) { # scenario 2
      # read in the tree file that contains the distribution of trees we sampled over in the geographic analysis
      # todo: file finding needs to be more robust
      tree_path <- list.files(paste0(unlist(strsplit(file_path, "/phylogeography"))[1], "/phylogeny"), recursive = T, full.names = T, pattern = "*sample.trees")
      if (length(tree_path) != 1) {
        stop("cannot find tree")
      }
      trees <- read.nexus(tree_path)
      tree_indices <- as.integer(log_dat[, currenttree_colnum] + 1L)
      
      trees <- trees[unique(tree_indices)] # only store the unique trees to save space
      tree_indices <- match(tree_indices, unique(tree_indices))
    } else if (tree_fixed) { # scenario 1
      
      if (file.exists(tree_path)) {
        trees <- read.nexus(tree_path)
        trees <- trees[1]
        tree_indices <- rep(1L, nrow(log_dat))
      } else {
        # todo: get the mcc tree path
        stop("cannot find tree")
      }
      
    } else {
      stop("don't know what to do")
    }
  }
  
  gc()
  
  ########################
  # prepare the Q matrix #
  ########################
  log_nrow <- nrow(log_dat)
  Q_all <- vector("list", log_nrow)
  for (l in 1:log_nrow) {
    Q_all[[l]] <- vector("list", epoch_num)
  }
  
  # using the estimated root freq vector as the prior probability at the root
  # following the default behavior of beast
  root_freqs_mat <- NULL
  rootfreq_colnum <- grep("\\.root\\.frequencies", colnames(log_dat))
  if (length(rootfreq_colnum) == K) {
    root_freqs_mat <- log_dat[, rootfreq_colnum]
    colnames(root_freqs_mat) <- states
  } else if (length(rootfreq_colnum) > 0 && length(rootfreq_colnum) != K) {
    stop ("cannot correctly identify root freqs.\n")
  }
  
  for (i in 1:epoch_num) {
    
    # overall rate
    mu_colnum <- grep(paste0("*\\.clock\\.rate.epoch", i, "|*\\.clock\\.rate$"), colnames(log_dat))
    averagerate_all <- log_dat[, mu_colnum]
    
    # rates
    if (epoch_num == 1) {
      rates_colnums <- grep("*\\.rates\\d*$|*\\.rates\\.*", colnames(log_dat))
    } else {
      rates_colnums <- grep(paste0("*\\.rates.epoch", i, "\\d*$|*\\.rates.epoch", i, "\\.*"), colnames(log_dat))
    }
    
    # set default to symmetric model
    # determine whether it's symmetrical model or asymmetrical model
    rates_sym <- T
    if (choose(K, 2) * 2 == length(rates_colnums)) {
      rates_sym <- F
    } else if (choose(K, 2) != length(rates_colnums)) {
      stop("The number of states extracted from the xml file does not match its counterpart computed from the log file.\n")
    }
    
    rates_startcolnum <- min(rates_colnums)
    rate_mat_all <- vector("list", log_nrow)
    for (l in 1:log_nrow) {
      rate_mat <- matrix(0, nrow = K, ncol = K)
      rate_mat[lower.tri(rate_mat)] <- as.numeric(log_dat[l, 1:(choose(K, 2)) + rates_startcolnum - 1])
      rate_mat <- t(rate_mat)
      if (rates_sym) {
        rate_mat[lower.tri(rate_mat)] <- as.numeric(log_dat[l, 1:(choose(K, 2)) + rates_startcolnum - 1])
      } else {
        rate_mat[lower.tri(rate_mat)] <- as.numeric(log_dat[l, (choose(K, 2) + 1):(choose(K, 2) * 2) + rates_startcolnum - 1])
      }
      rate_mat_all[[l]] <- rate_mat
    }
    
    # indicators
    if (epoch_num == 1) {
      indicators_colnums <- grep("*\\.indicators\\d*$|*\\.indicators\\.*", colnames(log_dat))
    } else {
      indicators_colnums <- grep(paste0("*\\.indicators.epoch", i, "\\d*$|*\\.indicators.epoch", i, "\\.*"), colnames(log_dat))
    }
    
    if (length(indicators_colnums) > 0) {
      # set default to symmetric model
      # determine whether it's symmetrical model or asymmetrical model
      indicators_sym <- T
      if (choose(K, 2) * 2 == length(indicators_colnums)) {
        indicators_sym <- F
      } else if (choose(K, 2) != length(indicators_colnums)) {
        stop("The number of states extracted from the xml file does not match its counterpart computed from the log file.\n")
      }
      
      indicators_startcolnum <- min(indicators_colnums)
      indicator_mat_all <- vector("list", log_nrow)
      
      for (l in 1:log_nrow) {
        indicator_mat <- matrix(0, nrow = K, ncol = K)
        indicator_mat[lower.tri(indicator_mat)] <- as.numeric(log_dat[l, 1:(choose(K, 2)) + indicators_startcolnum - 1])
        indicator_mat <- t(indicator_mat)
        if (indicators_sym) {
          indicator_mat[lower.tri(indicator_mat)] <- as.numeric(log_dat[l, 1:(choose(K, 2)) + indicators_startcolnum - 1])
        } else {
          indicator_mat[lower.tri(indicator_mat)] <- as.numeric(log_dat[l, (choose(K, 2) + 1): (choose(K, 2) * 2) + indicators_startcolnum - 1])
        }
        indicator_mat_all[[l]] <- indicator_mat
      }
    }
    
    # construct Q matrix
    for (l in 1:log_nrow) {
      
      Q_all[[l]][[i]] <- rate_mat_all[[l]]
      if (length(indicators_colnums) > 0) {
        Q_all[[l]][[i]] <- Q_all[[l]][[i]] * indicator_mat_all[[l]]
      }
      diag(Q_all[[l]][[i]]) <- -apply(Q_all[[l]][[i]], 1, sum)
      
      # note here the Q matrix would be correctly rescaled if it's symmetric
      # but would be incorret if it's asymmetric (as pi won't be flat)
      # however, this is how beast has been rescaling internally so far, so we follow it for now
      # todo: providing a toggle to choose to compute stationary freq and then using it to rescale correctly
      # this shouldn't matter though as the sim.history function takes the unrescaled Q matrix as its input
      pi <- 1 / K
      
      mu <- -sum(diag(Q_all[[l]][[i]]) * pi)
      Q_all[[l]][[i]] <- (Q_all[[l]][[i]] / mu) * averagerate_all[l]
      
      dimnames(Q_all[[l]][[i]]) <- list(states, states)
      
    } # end Q matrix construction
  } # end epoch loop
  
  # take care of memory issue
  rm(rate_mat_all)
  if (length(indicators_colnums) > 0) {
    rm(indicator_mat_all)
  }
  rm(averagerate_all)
  rm(log_dat)
  gc()
  
  ######################
  # perform simulation #
  ######################
  cat("start simulation.\n")
  
  if (ncores == 1 || sample_size < 100L) { # single core
    histories <- vector("list", sample_size)
    for (l in 1:log_nrow) {
      histories[[l]] <- history_sim(tree = trees[[tree_indices[l]]], Q = Q_all[[l]], Q_ages = Q_ages, root_freq = unlist(root_freqs_mat[l, ]), states = states,
                                    numstates_conditioned = numstates_conditioned, nrejections_max = nrejections_max, conditional = conditional)
      if (l %% 10 == 0) {
        cat(paste0(l, " histories has been simulated.\n"))
      }
    }
  } else { # multi-core parallelization
    
    # partition the simulations for each core
    nchunks <- ncores * 5L
    if (sample_size < nchunks) {
      nchunks <- sample_size
    }
    nhis_percore <- floor(sample_size / nchunks)
    his_indices <- vector("list", nchunks)
    k <- 0L
    for (j in 1:nchunks) {
      if (j < nchunks) {
        his_indices[[j]] <- 1L:nhis_percore + k
      } else {
        his_indices[[j]] <- (k + 1L):sample_size
      }
      k <- k + nhis_percore
    }
    
    # parse the parameters to simulate a number of histories
    histories_all <- pblapply(his_indices, function(his_idx) {
      histories_sim(trees = trees[tree_indices[his_idx]], Qs = Q_all[his_idx], Q_ages = Q_ages, root_freqs = root_freqs_mat[his_idx, ], states = states,
                    numstates_conditioned = numstates_conditioned, nrejections_max = nrejections_max, indices = his_idx, conditional = conditional)
    }, cl = ncores)
    
    histories <- do.call(c, histories_all)
  }
  
  rm(Q_all)
  gc()
  
  names(histories) <- sample_indices
  return(histories)
}
