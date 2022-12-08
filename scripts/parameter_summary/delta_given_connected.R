# this script contains functions that compute the probability of number of dispersal routes conditioning on connectivity
library(gmp)
library(PolynomF)
library(partitions)

#' Compute the number of graphs (or digraphs when sym is true) given n labeled vertices
#' @param n number of vertices
#' @param sym whether count the graph is undirected (i.e., symmetric) or directed (i.e., asymmetric)
#' @return a vector where the ith element is the number of connected (di)graphs given n vertices and i - 1 edges
nGraphsPerDelta <- function(n, sym = TRUE) {
  k <- choose(n, 2) * (2 - sym)
  total_num_graphs <- chooseZ(k, 0:k)
  return(total_num_graphs)
}

#' Compute the number of connected graphs given n labeled vertices
#' @param n number of vertices
#' @param sym whether count the graph is undirected (i.e., symmetric) or directed (i.e., asymmetric)
#' @param precomp whether obtain the value by fetching the precomputed table (much faster for large n)
#' @param precomp_filepath path of the file containing the precomputed table for the number of connected graphs
#' or computing on the fly (support up to a certain number of vertices)
#' @return a numeric vector where the ith element is the number of connected graphs given n vertices and i - 1 edges
nConnectedGraphsPerDelta <- function(n, sym = TRUE, precomp = TRUE, precomp_filepath = NULL) {
  if (precomp) {
    if (is.null(precomp_filepath)) {
      precomp_filepath_basename <- "num_connected_graphs.txt"
      if (!sym) {
        precomp_filepath_basename <- gsub("_graphs\\.txt$", "_digraphs.txt", precomp_filepath_basename)
      }
      precomp_filepath <- list.files(recursive = T, all.files = T, full.names = T, pattern = precomp_filepath_basename)
      if (length(precomp_filepath) != 1) {
        stop("The precomputed table for the number of connected graphs is not found.")
      }
    }

    graphs_text <- system(paste0("head -n", n, " ", precomp_filepath, " | tail -n1"), intern = T)
    graphs_text <- unlist(strsplit(gsub(" |\\[|\\]", "", graphs_text), ","))
    return(as.bigz(graphs_text))
  } else {
    if (sym) {
      return(nConnectedGraphsPerDeltaSym(n))
    } else {
      return(nConnectedGraphsPerDeltaAsym(n))
    }
  }
}

#' Compute the number of connected graphs given n labeled vertices
#' implementation of P29, Exercise 1.5(a), Harary and Palmer 2014
#' @param n number of vertices
#' @return a vector where the ith element is the number of connected graphs given n vertices and i - 1 edges
nConnectedGraphsPerDeltaSym <- function(n) {
  
  # obtain all partitions of n
  if (n <= 83) {
    parts <- as.matrix(parts(n))
    nparts <- ncol(parts)
    parts_nozero <- vector("list", nparts)
    for (i in 1:nparts) {
      t <- parts[, i]
      t <- t[t > 0]
      parts_nozero[[i]] <- t
    }
    rm(parts)
  } else {
    nparts <- P(n)
    parts_nozero <- vector("list", nparts)
    part_this <- firstpart(n)
    i <- 1
    while(i <= nparts) {
      parts_nozero[[i]] <- part_this[part_this > 0]
      part_this <- nextpart(part_this)
      i <- i + 1
    }
  }

  # gc(verbose = F)
  nops <- lengths(parts_nozero)
  
  first_term   <- (-1)^(nops + 1) / nops
  
  second_term <- as.bigq(numeric(nparts))
  nfact <- factorialZ(n)
  for (i in 1:nparts) {
    second_term[i] <- nfact / prod(factorialZ(parts_nozero[[i]]))
  }
  
  third_term <- as.bigq(numeric(nparts))
  nopsfact <- factorialZ(nops)
  for (i in 1:nparts) {
    third_term[i] <- nopsfact[i] / prod(factorialZ(table(parts_nozero[[i]])))
  }
  
  firstthree_term <- first_term * second_term * third_term
  rm(first_term, second_term, third_term)
  # gc(verbose = F)
  
  num_connected_graphs <- as.bigq(numeric(choose(n, 2) + 1))
  pts <- vapply(parts_nozero, function(x) sum(choose(x, 2)), FUN.VALUE = numeric(1), USE.NAMES = F)
  for (j in (n - 1):choose(n, 2)) {
    num_connected_graphs[j + 1] <- sum(chooseZ(pts, j) * firstthree_term)
  }
  num_connected_graphs <- as.bigz.bigq(num_connected_graphs)
  
  return(num_connected_graphs)
}

#' Compute the number of connected digraphs given n labeled vertices
#' implementation of the recurrence algorithm provided in Corollary 7 of Archer et al. 2020
#' @param n number of vertices
#' @return a numeric vector where the ith element is the number of connected graphs given n vertices and i - 1 edges
nConnectedGraphsPerDeltaAsym <- function(n) {
  return(coef(esPoly(n))[[n]])
}

#' the Ita recurrence in Corollary 7 of Archer et al. 2020 
#' (helper function for computing the number of connected digraphs given n labeled vertices)
#' @param n number of vertices
#' @return a list of polynomials
itaPoly <- function(n) {
  v <- as_polylist(lapply(1:n, function(i) return(1)))
  if (n <= 1) return(v)
  
  y <- polynomial()
  for (i in 2:n) {
    v[[i]] <- (1 + y)^(i * (i - 1))
    for (k in 1:(i - 1)) {
      p2 <- (1 + y)^((i - 1) * (i - k))
      v[[i]] <- v[[i]] - choose(i, k) * p2 * v[[k]]
    }
  }
  
  return(v)
}

#' the S recurrence in Corollary 7 of Archer et al. 2020 
#' (helper function for computing the number of connected digraphs given n labeled vertices)
#' @param n number of vertices
#' @return a list of polynomials
esPoly <- function(n) {
  u <- itaPoly(n)
  v <- as_polylist(lapply(1:n, function(i) return(1)))
  
  for (i in 2:n) {
    v[[i]] <- u[[i]]
    for (k in 1:(i - 1)) {
      v[[i]] <- v[[i]] + choose(i - 1, k - 1) * u[[i - k]] * v[[k]]
    }
  }
  
  return(v)
}

#' Compute the fraction of strongly connected graphs given n labeled vertices
#' @param n number of vertices
#' @param sym whether count the graph is undirected (i.e., symmetric) or directed (i.e., asymmetric)
#' @param precomp whether obtain the value by fetching the precomputed table (much faster for large n)
#' @param precomp_filepath path of the file containing the precomputed table for the number of connected graphs
#' or computing on the fly (support up to a certain number of vertices)
#' @return a numeric vector where the ith element is the fraction of strongly connected graphs given Delta = i - 1 and n vertices
pConnectedPerDelta <- function(n, sym = TRUE, precomp = TRUE, precomp_filepath = NULL) {
  if (n > 100) {
    stop("Currently the largest supported number of vertices is 100.")
  } else if ((!precomp) && ((sym == TRUE && n > 45) || (sym == FALSE && n > 32))) {
    message("Currently the largest supported number of vertices for computing on the fly is ", ifelse(sym, 45, 32), 
            " for the", ifelse(sym, "", "a"), "symmetric model; switching to use the precomputed value.")
    precomp <- TRUE
  }
  n_connected_graphs <- nConnectedGraphsPerDelta(n, sym = sym, precomp = precomp, precomp_filepath = precomp_filepath)
  n_total_graphs <- nGraphsPerDelta(n, sym = sym)
  p_connected_givenDelta <- as.numeric(n_connected_graphs / n_total_graphs)
  return(p_connected_givenDelta)
}

#' Compute the fraction of strongly connected graphs given n labeled vertices and the prior on the number of edges, Delta
#' @param n number of vertices
#' @param sym whether count the graph is undirected (i.e., symmetric) or directed (i.e., asymmetric)
#' @return the fraction of connected graphs given n labeled vertices (a numeric value between 0 and 1)
pConnectedAnalytical <- function(n, f = NULL, sym = TRUE) {
  ms <- 0:(choose(n, 2) * (2 - sym))
  if (is.function(f)) {
    p_Delta <- sapply(ms, function(x) f(x = x, n = n, sym = sym))
  } else if (is.null(f)) {
    p_Delta <- sapply(ms, function(x) pDelta(x = x, n = n, sym = sym))
  } else {
    p_Delta <- f
  }
  
  p_connected_givenDelta <- pConnectedPerDelta(n, sym = sym)
  p_connected <- sum(p_connected_givenDelta * p_Delta)
  
  return(p_connected)
}

#' Compute the probability of Delta (number of edges) given all the corresponding graphs are strongly connected
#' @param n number of vertices
#' @param sym whether count the graph is undirected (i.e., symmetric) or directed (i.e., asymmetric)
#' @param pDelta the prior probability density function on Delta without conditioning on connectivity
#' @return a numeric vector where the ith element is the probability of Delta = i - 1 given n vertices
pDeltaGivenConnected <- function(n, f = NULL, sym = TRUE) {
  ms <- 0:(choose(n, 2) * (2 - sym))
  if (is.function(f)) {
    p_Delta <- sapply(ms, function(x) f(x = x, n = n, sym = sym))
  } else if (is.null(f)) {
    p_Delta <- sapply(ms, function(x) pDelta(x = x, n = n, sym = sym))
  } else {
    p_Delta <- f
  }
  
  p_connected_givenDelta <- pConnectedPerDelta(n, sym = sym)
  p_connected <- sum(p_connected_givenDelta * p_Delta)
  p_Delta_given_connected <- p_connected_givenDelta * p_Delta / p_connected
  
  return(p_Delta_given_connected)
}

#' Compute the expected Delta (number of edges) given all the corresponding graphs are strongly connected
#' @param n number of vertices
#' @param sym whether count the graph is undirected (i.e., symmetric) or directed (i.e., asymmetric)
#' @param pDelta the prior probability density function on Delta without conditioning on connectivity
#' @return expected number of edges given the prior on Delta conditioning on strong connectivity (a positive numeric value)
expectedDeltaGivenConnected <- function(n, pDelta, sym = T) {
  ms <- 0:(choose(n, 2) * (2 - sym))
  p_Delta_given_connected <- pDeltaGivenConnected(n = n, pDelta = pDelta, sym = sym)
  expected_Delta_given_connected <- sum(p_Delta_given_connected * ms)
  return(expected_Delta_given_connected)
}

#' Compute prior probability that each edge exists given all the corresponding graphs are strongly connected
#' @param n number of vertices
#' @param sym whether count the graph is undirected (i.e., symmetric) or directed (i.e., asymmetric)
#' @param pDelta the prior probability density function on Delta without conditioning on connectivity
#' @return prior probability that each edge exists given the prior on Delta conditioning on strong connectivity (a numeric value between 0 and 1)
pEdgeGivenConnected <- function(n, pDelta, sym = T) {
  expected_Delta_given_connected <- expectedDeltaGivenConnected(n = n, pDelta = pDelta, sym = sym)
  k <- choose(n, 2) * (2 - sym)
  p_edge_given_connected <- expected_Delta_given_connected / k
  return(p_edge_given_connected)
}

#' Compute the probability of a truncated Poisson random variable.
#' @param x vector of quantiles.
#' @param lambda vector of (non-negative) means.
#' @param x_max maximum value that the random variable can take.
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @return The (log) density of the random variable(s).
dtpois <- function(x, lambda, x_max = NULL, log = FALSE, upper_trunc = TRUE) {
  if (upper_trunc == FALSE) {
    return(dpois(x = x, lambda = lambda, log = log))
  }
  
  if (x > x_max) {
    p <- 0
  } else {
    normalizing_constant <- ppois(x_max, lambda)
    numerator <- dpois(x, lambda)
    p <- numerator / normalizing_constant
  }
  if (log) p <- log(p)
  return(p)
}

#' Prior probability density function on Delta (number of edges) following the BEAST way
#' @param x number of edges
#' @param n number of vertices
#' @param sym whether count the graph is undirected (i.e., symmetric) or directed (i.e., asymmetric)
#' @param default using the default or alternative prior on Delta
#' @param upper_trunc whether truncate at the upper limit of the number of edges (n(n - 1)/2 for the symmetric graph and n(n-1) for the asymmetric graph)
#' @param lambda rate of the Poisson distribution
#' @param offset offset of the Poisson distribution
#' @return prior probability that each edge exists (a numeric value between 0 and 1)
pDelta <- function(x, n, sym = TRUE, default = TRUE, upper_trunc = TRUE, lambda = NULL, offset = NULL) {
  if (is.null(lambda)) {
    if (default) {
      lambda <- ifelse(sym, log(2), n - 1)
    } else {
      lambda <- ifelse(sym, ceiling(choose(n, 2) / 2) - n + 1, choose(n, 2))
      if (lambda < 1) {
        lambda <- log(2)
      }
    }
  }
  if (is.null(offset)) {
    offset <- ifelse(sym, n - 1, 0)
  }
  
  k <- choose(n, 2) * (2 - sym)
  return(dtpois(x - offset, lambda = lambda, x_max = k - offset, upper_trunc = upper_trunc))
}

#' Compute prior probability that each edge exists following the BEAST+Spread way
#' @param n number of vertices
#' @param sym whether count the graph is undirected (i.e., symmetric) or directed (i.e., asymmetric)
#' @param default using the default or alternative prior on Delta
#' @return prior probability that each edge exists (a numeric value between 0 and 1)
pEdgeBeast <- function(n, sym = T, default = T) {
  if (default) {
    lambda <- ifelse(sym, log(2), n - 1)
  } else {
    lambda <- ifelse(sym, ceiling(choose(n, 2) / 2) - n + 1, choose(n, 2))
    if (lambda < 1) {
      lambda <- log(2)
    }
  }
  
  expected_Delta <- lambda
  if (sym) expected_Delta <- expected_Delta + n - 1 
  
  k <- choose(n, 2) * (2 - sym)
  p_edge <- expected_Delta / k
  return(p_edge)
}
