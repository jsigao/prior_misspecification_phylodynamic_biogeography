library(ape)

#' Generate BEAST XML scripts for performing phylodynamic analysis to infer biogeographic history 
#' (should be compatible with beast 1.8.2-1.10.5; no gurantee for other versions)
#' @param folder_path path to read in relevant supporting file (assuming a certain folder structure)
#' @param writeto_folder_path path to write out the generated XML script (if not provided then \code{folder_path} will be used)
#' @param fixed_tree_path path to the summary tree that will be conditioned on (if not provided then not condition on a single tree)
#' @param discrete_trait_path path to the spreadsheet that contains the sampling time and geographic location (or other discrete trait) information
#' @param date_direction whether the sampling time provided in \code{discrete_trait_path} is time ("forwards") or age ("backwards")
#' @param date_units units (either days or years) of the sampling time provided in \code{discrete_trait_path} 
#' @param lheat number of copies to duplicate the data for data-cloning analyses (default 1 means no data cloning)
#' @param mcmc_chainlength number of generations the MCMC simulation to run
#' @param mcmc_samplingfreq frequency (number of generations) of the MCMC simulation to write parameter estimates to the log file
#' @param mcmc_screen_samplingfreq frequency (number of generations) of the MCMC simulation to write to stdout
#' @param mcmc_tree_samplingfreq frequency (number of generations) of the MCMC simulation to write the tree file
#' @param mcmc_adaptationDelay number of generations at the beginning of the MCMC before the proposals are tunned continuously
#' @param ml_chainlength number of generations of the MCMC simulation at each power of the power-posterior analysis
#' @param ml_burnin number of generations that no sample would be taken (at the beginning of the MCMC simulation) at each power of the power-posterior analysis
#' @param stones_num number of powers to use for the power-posterior analysis
#' @param ml_samplingfreq frequency (number of generations) of the power-posterior MCMC simulation to write the log-likelihood log file
#' @param ml_param_samplingfreq frequency (number of generations) of the power-posterior MCMC simulation to write parameter estimates to the parameter log file
#' @param ml_tree_samplingfreq frequency (number of generations) of the power-posterior MCMC simulation to write the tree file
#' @param posterior whether to infer the joint posterior distribution (true) or the joint prior distribution
#' @param complete_history whether to perform the simulation-free fast stochastic mapping (false) or simulate full history (true)
#' @param symmetry instantaneous-rate matrix to be symmetric (true) or asymmetric (false); in asymmetric matrix the forward and backward rates can be different
#' @param ctmc using the BEAST default ctmc rate-reference prior (true) or an alternative diffuse prior (false) on the average dispersal rate
#' @param poisson_default using the BEAST default poisson prior (true) or an alternative diffuse poisson prior (false) on the number of dispersal routes
#' @param bssvs using BSSVS (Bayesian Stochastic Search Variable Selection) to average over rate matrices (true) or not (false)
#' @param poisson specifiying a poisson prior on the number of dispersal routes (true) or not (false)
#' @param empiricaltree_mh treating the tree sampled from a provided distribution of trees as a state of the Markov chain (true) or simply averging over the distribution (false)
#' @param powerposterior_vanilla whether doing stochastic mapping for the power-posterior analysis (false) or not (true)
#' @param tree_operator_weight weight of the proposal that samples tree from the provided distribution of trees
#' @param rootfreq_proposal_weightb weight of the proposal on the root frequency vector
#' @param clockrate_proposal_weight weight of the proposal on the average dispersal rate
#' @param clockratemean_proposal_weight weight of the proposal on hyper-parameter of the average dispersal rate (when applicable when the alternative diffuse prior is used)
#' @param rates_proposal_weight weight of the proposal on the relative rates in the instantaneous-rate matrix
#' @param indicators_proposal_weight weight of the proposal on the indicator variables that enable BSSVS

generate_phylogeography_xml <- function(folder_path, writeto_folder_path = NULL, fixed_tree_path = NULL, discrete_trait_path = NULL, 
                                        date_direction = "forwards", date_units = "years", 
                                        lheat = 1, mcmc_chainlength, mcmc_samplingfreq, mcmc_screen_samplingfreq = 0, 
                                        mcmc_tree_samplingfreq = 0, mcmc_adaptationDelay = 0,
                                        ml_chainlength = 0, ml_burnin = 0, stones_num = 0, 
                                        ml_samplingfreq = 0, ml_param_samplingfreq = 0, ml_tree_samplingfreq = 0,
                                        posterior = T, complete_history = F, symmetry = T, ctmc = T, poisson_default = T, 
                                        bssvs = T, poisson = T, empiricaltree_mh = T, powerposterior_vanilla = F,
                                        tree_operator_weight = 25, rootfreq_proposal_weight = 0.5,
                                        clockrate_proposal_weight = 10, clockratemean_proposal_weight = 5,
                                        rates_proposal_weight = 15, indicators_proposal_weight = 30)
{
  
  if (is.null(writeto_folder_path)) {
    writeto_folder_path <- folder_path
  }
  
  if (!is.null(discrete_trait_path)) {
    states_dat <- read.table(discrete_trait_path, header = T, sep = ",", stringsAsFactors = F)
  } else { # try to find the discrete trait file
    # stop("cannot find discrete trait information")
    discrete_trait_path <- list.files(unlist(strsplit(folder_path, "/ctmc|/exp_hyper"))[1], full.names = T, pattern = "*discrete_trait.txt$")
    if (length(discrete_trait_path) != 1) {
      stop("cannot find discrete trait information")
    } else {
      states_dat <- read.table(discrete_trait_path, header = T, sep = ",", stringsAsFactors = F)
    }
  }
  
  # assume there is only one discrete trait for now, could be easily relaxed
  # also assumes the name of that given discrete trait is either geography or host; could be easily relaxed as well
  discrete_trait_name <- colnames(states_dat)[colnames(states_dat) == "geography" | colnames(states_dat) == "host"]
  colnames(states_dat)[colnames(states_dat) == "geography" | colnames(states_dat) == "host"] <- "trait"
  colnames(states_dat) <- gsub("taxon_id", "taxon", colnames(states_dat))
  
  states <- sort(as.vector(unique(states_dat$trait)))
  if (any(states == "?")) {
  	states <- states[-which(states == "?")]
  }
  states_num <- length(states)
  
  file_name <- unlist(strsplit(unlist(strsplit(folder_path, "/phylogeography/"))[1], "runs/current/|analyses/"))
  file_name <- gsub("/", "_", file_name[length(file_name)])
  if (complete_history) {
    file_name <- paste0(file_name, "_completeHistory")
  } else if (posterior && powerposterior_vanilla) {
    file_name <- paste0(file_name, "_powerposterior")
  }
  
  # retrieve the name of the marginal distribution of trees
  treesfile_name <- basename(list.files(paste0(unlist(strsplit(folder_path, "/phylogeography")), "/phylogeny"), recursive = T, full.names = T, pattern = "*sample.trees"))
  relativedir_layernum <- length(unlist(strsplit(unlist(strsplit(folder_path, "/phylogeography"))[2], "/")))
  treesfile_name <- paste0(paste(rep("../", relativedir_layernum), collapse = ""), "phylogeny/", treesfile_name)
  
  x <- c("<?xml version=\"1.0\" standalone=\"yes\"?>", "<beast>", "</beast>")
  x <- append(x, c("\t<taxa id=\"taxa\">", "\t</taxa>\n"), after = 2)
  
  # insert discrete data
  for (i in 1:nrow(states_dat)) {
    
    taxon <- paste0("\t\t<taxon id=\"", states_dat$taxon[i], "\">\n")
    if ("date" %in% colnames(states_dat)) {
      taxon <- paste0(taxon, "\t\t\t<date value=\"", states_dat$date[i], "\" direction=\"", date_direction, "\" units=\"", date_units, "\"/>\n")
    }
    if (lheat == 1) {
      taxon <- paste0(taxon, "\t\t\t<attr name=\"", discrete_trait_name, "\">", ifelse(posterior, states_dat$trait[i], "?"), "</attr>\n")
    }
    x <- append(x, paste0(taxon, "\t\t</taxon>"), after = length(x) - 2)
  }
  
  # instert list of discrete trait
  x <- append(x, c(paste0("\t<generalDataType id=\"", discrete_trait_name, ".dataType\">"), "\t</generalDataType>"), after = length(x) - 1)
  for (i in 1:states_num) {
    if (lheat == 1) {
      state_line <- paste0("\t\t<state code=\"", states[i], "\"/>")
    } else if (lheat > 1) {
      state_line <- paste0("\t\t<state code=\"", c(LETTERS, letters, 0:9)[i], "\"/>")
    }
    
    x <- append(x, state_line, after = length(x) - 2)
  }
  
  # insert data pattern chunk for discrete trait
  if (lheat == 1) {
    x <- append(x, paste0("\t<attributePatterns id=\"", discrete_trait_name, ".pattern\" attribute=\"", discrete_trait_name, 
                          "\">\n\t\t<taxa idref=\"taxa\"/>\n",
                          "\t\t<generalDataType idref=\"", discrete_trait_name, ".dataType\"/>\n",
                          "\t</attributePatterns>\n"), after = length(x) - 1)
  } else if (lheat > 1) {
    x <- append(x, c(paste0("\n\t<alignment id=\"", discrete_trait_name, "\">"), "\t</alignment>"), after = length(x) - 1)
    x <- append(x, paste0("\t<generalDataType idref=\"", discrete_trait_name, ".dataType\"/>"), after = length(x) - 2)
    
    for (i in 1:nrow(states_dat)) {
      
      if (!(states_dat$trait[i] %in% states)) {
        state_dc <- states_dat$trait[i]
      } else {
        state_dc <- c(LETTERS, letters, 0:9)[match(states_dat$trait[i], states)]
      }
      
      taxon_stateseq <- paste0("\t\t<sequence>\t<taxon idref=\"", states_dat$taxon[i], "\"/>\t", 
                               paste(rep(state_dc, lheat), collapse = ""), "\t</sequence>")
      x <- append(x, taxon_stateseq, after = length(x) - 2)
    }
    
    # add site compression chunk
    x <- append(x, paste0("\t<patterns id=\"patterns\" from=\"1\" strip=\"false\">\n",
                          "\t\t<alignment idref=\"", discrete_trait_name, "\"/>\n", "\t</patterns>\n"), after = length(x) - 1)
  }
  
  # insert marginal distribution of trees chunk if not fixed tree
  if (!is.null(fixed_tree_path)) {
    
    tree <- ape::read.nexus(fixed_tree_path)
    tree$tip.label <- gsub("\'", "", tree$tip.label)
    tree_newick <- write.tree(tree)
    
    x <- append(x, paste0("\t<newick id=\"startingTree\" units=\"", date_units, "\" usingDates=\"true\">\n",
                          "\t\t", tree_newick, "\n",
                          "\t</newick>"), after = length(x) - 1)
    x <- append(x, paste0("\t<treeModel id=\"treeModel\">\n",
                       "\t\t<tree idref=\"startingTree\"/>\n",
                       "\t\t<rootHeight>\n",
                       "\t\t\t<parameter id=\"treeModel.rootHeight\"/>\n",
                       "\t\t</rootHeight>\n",
                       "\t\t<nodeHeights internalNodes=\"true\">\n",
                       "\t\t\t<parameter id=\"treeModel.internalNodeHeights\"/>\n",
                       "\t\t</nodeHeights>\n",
                       "\t\t<nodeHeights internalNodes=\"true\" rootNode=\"true\">\n",
                       "\t\t\t<parameter id=\"treeModel.allInternalNodeHeights\"/>\n",
                       "\t\t</nodeHeights>\n",
                       "\t</treeModel>\n"), after = length(x) - 1)
  } else {
    x <- append(x, paste0("\t<empiricalTreeDistributionModel id=\"treeModel\" fileName=\"", treesfile_name, "\">\n",
                          "\t\t<taxa idref=\"taxa\"/>\n",
                          "\t</empiricalTreeDistributionModel>"), after = length(x) - 1)
    x <- append(x, paste0("\t<statistic id=\"treeModel.currentTree\" name=\"Current Tree\">\n", 
                          "\t\t<empiricalTreeDistributionModel idref=\"treeModel\"/>\n", 
                          "\t</statistic>\n"), after = length(x) - 1)
  }
  
  # insert the chunk for clock rate of discrete trait
  discrete_trait_clock_rate <- paste0("\t<strictClockBranchRates id=\"", discrete_trait_name, ".branchRates\">\n",
                                    "\t\t<rate>\n\t\t\t<parameter id=\"", discrete_trait_name, ".clock.rate\" value=\"0.01\" lower=\"0.0\"/>\n",
                                    "\t\t</rate>\n",
                                    "\t</strictClockBranchRates>\n")
  if (!ctmc) { # not the default prior
    discrete_trait_clock_rate <- paste0(discrete_trait_clock_rate, "\t<distributionLikelihood id=\"", discrete_trait_name, ".clock.rate.exp\">\n",
                                      "\t\t<distribution>\n\t\t\t<exponentialDistributionModel>\n",
                                      "\t\t\t\t<mean>\n",
                                      "\t\t\t\t\t<parameter id=\"", discrete_trait_name, ".clock.rate.exp.mean\" value=\"1\"/>\n",
                                      "\t\t\t\t</mean>\n",
                                      "\t\t\t</exponentialDistributionModel>\n",
                                      "\t\t</distribution>\n",
                                      "\t\t<data>\n",
                                      "\t\t\t<parameter idref=\"", discrete_trait_name, ".clock.rate\"/>\n",
                                      "\t\t</data>\n",
                                      "\t</distributionLikelihood>\n")
  }
  x <- append(x, discrete_trait_clock_rate, after = length(x) - 1)
  
  # insert the chunk for the substitution model of discrete trait
  substitution_model <- paste0("\t<generalSubstitutionModel id=\"", discrete_trait_name, ".model\">\n",
                               "\t\t<generalDataType idref=\"", discrete_trait_name, ".dataType\"/>\n")
  frequencies <- paste0("\t\t<frequencies>\n",
                        "\t\t\t<frequencyModel id=\"", discrete_trait_name, ".frequencyModel\" normalize=\"true\">\n",
                        "\t\t\t\t<generalDataType idref=\"", discrete_trait_name, ".dataType\"/>\n",
                        "\t\t\t\t<frequencies>\n",
                        "\t\t\t\t\t<parameter id=\"", discrete_trait_name, ".frequencies\" dimension=\"", states_num, "\"/>\n",
                        "\t\t\t\t</frequencies>\n",
                        "\t\t\t</frequencyModel>\n",
                        "\t\t</frequencies>\n")
  
  if (symmetry) {
    edges_max <- choose(states_num, 2)
  } else {
    edges_max <- choose(states_num, 2) * 2
  }
  rates <- paste0("\t\t<rates>\n",
                  "\t\t\t<parameter id=\"", discrete_trait_name, ".rates\" dimension=\"", edges_max, "\" value=\"1.0\"/>\n",
                  "\t\t</rates>\n")
  indicators <- paste0("\t\t<rateIndicator>\n",
                       "\t\t\t<parameter id=\"", discrete_trait_name, ".indicators\" dimension=\"", edges_max, "\" value=\"1.0\"/>\n",
                       "\t\t</rateIndicator>\n")
  
  if (!bssvs) {
    indicators <- ""
  }
  substitution_model <- paste0(substitution_model, frequencies, rates, indicators, "\t</generalSubstitutionModel>")
  x <- append(x, substitution_model, after = length(x) - 1)
  
  # insert the non-zero rates chunk
  if (bssvs) {
    nonZeroRates <- paste0("\t<sumStatistic id=\"", discrete_trait_name, ".nonZeroRates\" elementwise=\"true\">\n",
                           "\t\t<parameter idref=\"",
                           discrete_trait_name, ".indicators\"/>\n",
                           "\t</sumStatistic>")
    x <- append(x, nonZeroRates, after = length(x) - 1)
  }
  
  # insert the site model chunk
  site_model <- paste0("\t<siteModel id=\"", discrete_trait_name, ".siteModel\">\n",
                       "\t\t<substitutionModel>\n",
                       "\t\t\t<generalSubstitutionModel idref=\"", discrete_trait_name, ".model\"/>\n",
                       "\t\t</substitutionModel>\n",
                       "\t</siteModel>\n")
  x <- append(x, site_model, after = length(x) - 1)
  
  # insert the markov jumps chunk
  if (lheat == 1) {
    markov_jumps <- paste0("\t\t<attributePatterns idref=\"", discrete_trait_name, ".pattern\"/>\n")
  } else if (lheat > 1) {
    markov_jumps <- paste0("\t\t<patterns idref=\"patterns\"/>\n")
  }
  markov_jumps <- paste0(markov_jumps, "\t\t<treeModel idref=\"treeModel\"/>\n",
                         "\t\t<siteModel idref=\"", discrete_trait_name,".siteModel\"/>\n",
                         "\t\t<generalSubstitutionModel idref=\"", discrete_trait_name, ".model\"/>\n",
                         "\t\t<strictClockBranchRates idref=\"", discrete_trait_name, ".branchRates\"/>\n")
  
  if (complete_history) {
    markov_jumps <- paste0("\t<markovJumpsTreeLikelihood id=\"", discrete_trait_name,
                         ".treeLikelihood\" useUniformization=\"true\" saveCompleteHistory=\"true\" logCompleteHistory=\"true\" compactHistory=\"true\">\n", markov_jumps)
  } else if (powerposterior_vanilla) {
    markov_jumps <- paste0("\t<markovJumpsTreeLikelihood id=\"", discrete_trait_name, ".treeLikelihood\">\n", markov_jumps)
  } else {
    markov_jumps <- paste0("\t<markovJumpsTreeLikelihood id=\"", discrete_trait_name, ".treeLikelihood\" stateTagName=\"", 
                           discrete_trait_name, ".states\">\n", markov_jumps)
  }
  
  # add root frequencies for asymmetric model
  if (!symmetry) {
    root_frequncies <- paste0("\t\t<frequencyModel id=\"root.frequencyModel\" normalize=\"true\">\n",
                              "\t\t\t<generalDataType idref=\"", discrete_trait_name, ".dataType\"/>\n",
                              "\t\t\t<frequencies>\n",
                              "\t\t\t\t<parameter id=\"", discrete_trait_name, ".root.frequencies\" dimension=\"", states_num,"\"/>\n",
                              "\t\t\t</frequencies>\n",
                              "\t\t</frequencyModel>\n")
    markov_jumps <- paste0(markov_jumps, root_frequncies)
  }
  x <- append(x, markov_jumps, after = length(x) - 1)
  
  if ((!complete_history) && lheat == 1 && (!powerposterior_vanilla)) {
    # total counts
    total_count_name <- paste0("\t\t<parameter id=\"", discrete_trait_name, ".count\" value=\" ")
    total_count <- numeric(length = states_num^2)
    for(i in 1:states_num) {
      for(j in 1:states_num) {
        if(j != i) {
          total_count[(i - 1) * states_num + j] <- 1
        }
      }
    }
    
    total_count_line <- paste0(total_count_name, paste0(total_count, collapse = " "), "\"/>")
    x <- append(x, total_count_line, after = length(x) - 1)
    
    # pairwise counts
    for (i in 1:states_num) {
      for (j in 1:states_num) {
        if(j != i) {
          count_name <- paste0("\t\t<parameter id=\"", i, "-", j, ".count", "\" value=\" ")
          counts <- numeric(length = states_num^2)
          counts[(i - 1) * states_num + j] <- 1
          count_line <- paste0(count_name, paste0(counts, collapse = " "), "\"/>")
          x <- append(x, count_line, after = length(x) - 1)
        }
      }
    }
  }
  
  x <- append(x, "\t</markovJumpsTreeLikelihood>\n", after = length(x) - 1)
  
  # add the operators chunk
  if (!is.null(fixed_tree_path)) {
    trees_operator <- ""
  } else {
    trees_operator <- paste0("\t\t<empiricalTreeDistributionOperator weight=\"", tree_operator_weight, 
                             "\" metropolisHastings=\"", ifelse(empiricaltree_mh, "true", "false"), "\">\n",
                             "\t\t\t<empiricalTreeDistributionModel idref=\"treeModel\"/>\n",
                             "\t\t</empiricalTreeDistributionOperator>\n")
  }
  
  clock_rate_operator <- paste0("\t\t<scaleOperator scaleFactor=\"0.75\" weight=\"", clockrate_proposal_weight, "\">\n",
                                "\t\t\t<parameter idref=\"", discrete_trait_name, ".clock.rate\"/>\n",
                                "\t\t</scaleOperator>\n")
  
  if (!ctmc) {
    clock_rate_operator <- paste0(clock_rate_operator, "\t\t<scaleOperator scaleFactor=\"0.75\" weight=\"", clockratemean_proposal_weight, "\">\n",
                                  "\t\t\t<parameter idref=\"", discrete_trait_name, ".clock.rate.exp.mean\"/>\n",
                                  "\t\t</scaleOperator>\n")
  }
  
  rates_operator <- paste0("\t\t<scaleOperator scaleFactor=\"0.75\" weight=\"", rates_proposal_weight, "\" scaleAllIndependently=\"true\">\n",
                           "\t\t\t<parameter idref=\"", discrete_trait_name, ".rates\"/>\n",
                           "\t\t</scaleOperator>\n")
  bitflip_operator <- ""
  if (bssvs) {
    bitflip_operator <- paste0("\t\t<bitFlipOperator weight=\"", indicators_proposal_weight, "\">\n",
                               "\t\t\t<parameter idref=\"", discrete_trait_name, ".indicators\"/>\n",
                               "\t\t</bitFlipOperator>\n")
  }
  
  root_frequencies_operator <- ""
  if (!symmetry) {
    root_frequencies_operator <- paste0("\t\t<deltaExchange delta=\"0.01\" weight=\"", rootfreq_proposal_weight, "\">\n",
                                        "\t\t\t<parameter idref=\"", discrete_trait_name, ".root.frequencies\"/>\n",
                                        "\t\t</deltaExchange>\n")
  }
  
  operators <- paste0("\t<operators id=\"operators\" optimizationSchedule=\"log\">\n", trees_operator, clock_rate_operator, 
                      rates_operator, bitflip_operator, root_frequencies_operator, "\t</operators>\n")
  
  x <- append(x, operators, after = length(x) - 1)
  
  # add the mcmc chunk
  # first get the prior chunk
  if (ctmc) {
    clock_rate_prior <- paste0("\t\t\t\t<ctmcScalePrior>\n",
                               "\t\t\t\t\t<ctmcScale>\n",
                               "\t\t\t\t\t\t<parameter idref=\"", discrete_trait_name, ".clock.rate\"/>\n",
                               "\t\t\t\t\t</ctmcScale>\n",
                               "\t\t\t\t\t<treeModel idref=\"treeModel\"/>\n",
                               "\t\t\t\t</ctmcScalePrior>\n")
  } else {
    clock_rate_prior <- paste0("\t\t\t\t<distributionLikelihood idref=\"", discrete_trait_name, ".clock.rate.exp\"/>\n",
                               "\t\t\t\t<gammaPrior shape=\"0.5\" scale=\"2.0\" offset=\"0.0\">\n",
                               "\t\t\t\t\t<parameter idref=\"", discrete_trait_name, ".clock.rate.exp.mean\"/>\n",
                               "\t\t\t\t</gammaPrior>\n")
  }
  
  nonZeroRates_prior <- ""
  if (bssvs && poisson) {
    
    if (poisson_default) {
      if (symmetry) {
        poisson_mean <- 0.6931471805599453
      } else {
        poisson_mean <- states_num - 1
      }
    } else {
      if (symmetry) {
        poisson_mean <- ceiling(choose(states_num, 2)/2) - states_num + 1
        if (poisson_mean < 1) {
          poisson_mean <- 0.6931471805599453
        }
      } else {
        poisson_mean <- choose(states_num, 2)
      }
    }
    
    poisson_offset <- states_num - 1
    if (!symmetry) {
      poisson_offset <- 0
    }
    nonZeroRates_prior <- paste0("\t\t\t\t<poissonPrior mean=\"", poisson_mean, "\" offset=\"", poisson_offset, "\">\n",
                                 "\t\t\t\t\t<statistic idref=\"", discrete_trait_name, ".nonZeroRates\"/>\n",
                                 "\t\t\t\t</poissonPrior>\n")
  }
  
  rates_prior <- paste0("\t\t\t\t<cachedPrior>\n",
                        "\t\t\t\t\t<gammaPrior shape=\"1.0\" scale=\"1.0\" offset=\"0.0\">\n",
                        "\t\t\t\t\t\t<parameter idref=\"", discrete_trait_name,".rates\"/>\n",
                        "\t\t\t\t\t</gammaPrior>\n",
                        "\t\t\t\t\t<parameter idref=\"", discrete_trait_name, ".rates\"/>\n",
                        "\t\t\t\t</cachedPrior>\n")
  
  root_frequencies_prior <- ""
  if (!symmetry) {
    root_frequencies_prior <- paste0("\t\t\t\t<uniformPrior lower=\"0.0\" upper=\"1.0\">\n",
                                     "\t\t\t\t\t<parameter idref=\"", discrete_trait_name, ".root.frequencies\"/>\n",
                                     "\t\t\t\t</uniformPrior>\n")
  }
  prior <- paste0("\t\t\t<prior id=\"prior\">\n",
                  clock_rate_prior, nonZeroRates_prior, root_frequencies_prior, rates_prior, 
                  paste0("\t\t\t\t<generalSubstitutionModel idref=\"", discrete_trait_name, ".model\"/>\n"), 
                  "\t\t\t</prior>\n")
  
  likelihood <- paste0("\t\t\t<likelihood id=\"likelihood\">\n",
                       "\t\t\t\t<markovJumpsTreeLikelihood idref=\"", discrete_trait_name, ".treeLikelihood\"/>\n",
                       "\t\t\t</likelihood>\n")
  
  posterior_log <- paste0("\t\t<posterior id=\"posterior\">\n", prior, likelihood, "\t\t</posterior>\n")
  
  if (mcmc_screen_samplingfreq == 0) {
    mcmc_screen_samplingfreq <- mcmc_samplingfreq
  }
  if (mcmc_tree_samplingfreq == 0) {
    mcmc_tree_samplingfreq <- mcmc_samplingfreq
  }
  
  log_screen <- paste0("\t\t<log id=\"screenLog\" logEvery=\"", as.integer(mcmc_screen_samplingfreq), "\">\n",
                       "\t\t\t<posterior idref=\"posterior\"/>\n",
                       "\t\t\t<prior idref=\"prior\"/>\n",
                       "\t\t\t<likelihood idref=\"likelihood\"/>\n",
                       "\t\t\t<parameter idref=\"", discrete_trait_name, ".clock.rate\"/>\n",
                       ifelse(bssvs, paste0("\t\t\t<sumStatistic idref=\"", discrete_trait_name, ".nonZeroRates\"/>\n"), ""),
                       "\t\t</log>\n\n")
  
  log_file <- paste0("\t\t<log id=\"fileLog\" logEvery=\"", as.integer(mcmc_samplingfreq),"\" fileName=\"", file_name, ".log\" overwrite=\"false\">\n",
                   "\t\t\t<posterior idref=\"posterior\"/>\n",
                   "\t\t\t<prior idref=\"prior\"/>\n",
                   "\t\t\t<likelihood idref=\"likelihood\"/>\n")
  
  parameters_log <- paste0("\t\t\t<parameter idref=\"",discrete_trait_name,".clock.rate\"/>\n",
                           "\t\t\t<parameter idref=\"", discrete_trait_name, ".rates\"/>\n")
  if (!ctmc) {
    parameters_log <- paste0("\t\t\t<parameter idref=\"", discrete_trait_name, ".clock.rate.exp.mean\"/>\n", parameters_log)
  }
  if (!symmetry) {
    parameters_log <- paste0("\t\t\t<parameter idref=\"", discrete_trait_name, ".root.frequencies\"/>\n", parameters_log)
  }
  if (bssvs) {
    parameters_log <- paste0(parameters_log, "\t\t\t<parameter idref=\"", discrete_trait_name, ".indicators\"/>\n",
                             "\t\t\t<sumStatistic idref=\"", discrete_trait_name, ".nonZeroRates\"/>\n")
  }
  
  if (is.null(fixed_tree_path)) {
    parameters_log <- paste0(parameters_log, "\t\t\t<statistic idref=\"treeModel.currentTree\"/>\n")
  }
  
  parameters_log <- paste0(parameters_log, "\t\t\t<markovJumpsTreeLikelihood idref=\"", discrete_trait_name, ".treeLikelihood\"/>\n")
  
  log_file <- paste0(log_file, parameters_log, "\t\t</log>\n\n")
  
  log_history <- ""
  if (complete_history) {
    tree_log <- paste0("\t\t\t<treeModel idref=\"treeModel\"/>\n",
                       "\t\t\t<markovJumpstreelikelihood idref=\"", discrete_trait_name, ".treeLikelihood\"/>\n")
    
    log_history <- paste0("\t\t<log id=\"historyLogger\" logEvery=\"", as.integer(mcmc_samplingfreq), "\" fileName=\"", file_name, ".txt\">\n",
                        "\t\t\t<completeHistoryLogger>\n",
                        "\t\t\t\t<markovJumpstreelikelihood idref=\"", discrete_trait_name, ".treeLikelihood\"/>\n",
                        "\t\t\t</completeHistoryLogger>\n",
                        "\t\t</log>\n")
  } else {
    tree_log <- paste0("\t\t\t<treeModel idref=\"treeModel\"/>\n",
                       "\t\t\t<trait name=\"", discrete_trait_name, ".states\" tag=\"", discrete_trait_name, "\">\n",
                       "\t\t\t\t<ancestralTreeLikelihood idref=\"", discrete_trait_name, ".treeLikelihood\"/>\n",
                       "\t\t\t</trait>\n")
  }
  log_tree <- paste0("\t\t<logTree id=\"treeFileLog\" logEvery=\"", as.integer(mcmc_tree_samplingfreq), "\" nexusFormat=\"true\" fileName=\"",
                     file_name, ".trees\" sortTranslationTable=\"true\">\n", tree_log, "\t\t</logTree>\n")
  
  if (powerposterior_vanilla) {
    log_file <- ""
    log_tree <- ""
  }
  
  if (mcmc_adaptationDelay == 0) {
    mcmc_adaptationDelay <- as.integer(mcmc_chainlength) * 0.01
  }
  mcmc <- paste0("\t<mcmc id=\"mcmc\" chainLength=\"", as.integer(mcmc_chainlength), 
                 "\" autoOptimize=\"true\" adaptationDelay=\"", as.integer(mcmc_adaptationDelay), "\">\n",
                 posterior_log,"\t\t<operators idref=\"operators\"/>\n\n",
                 log_screen, log_file, log_tree, log_history, "\t</mcmc>\n")
  x <- append(x, mcmc, after = length(x) - 1)
  
  if (ml_burnin == 0) {
    ml_burnin <- as.integer(ml_chainlength) * 0.1
  }
  
  if (posterior && (!complete_history) && (lheat == 1) && ml_chainlength > 0) {
    mle <- paste0("\t<marginalLikelihoodEstimator chainLength=\"", as.integer(ml_chainlength), "\" burnin=\"", 
                  as.integer(ml_burnin), "\" pathSteps=\"", stones_num,
                "\" pathScheme=\"betaquantile\" alpha=\"0.30\" printOperatorAnalysis=\"true\">\n",
                "\t\t<samplers>\n",
                "\t\t\t<mcmc idref=\"mcmc\"/>\n",
                "\t\t</samplers>\n",
                "\t\t<pathLikelihood id=\"pathLikelihood\">\n",
                "\t\t\t<source>\n",
                "\t\t\t\t<posterior idref=\"posterior\"/>\n",
                "\t\t\t</source>\n\t\t\t<destination>\n",
                "\t\t\t\t<prior idref=\"prior\"/>\n",
                "\t\t\t</destination>\n",
                "\t\t</pathLikelihood>\n")
    mle_log <- paste0("\t\t<log id=\"MLELog\" logEvery=\"", as.integer(ml_samplingfreq), "\" fileName=\"MLE.log\">\n",
                      "\t\t\t<pathLikelihood idref=\"pathLikelihood\"/>\n",
                      "\t\t</log>\n")
    mle_filelog <- paste0("\t\t<log id=\"MLEFileLog\" logEvery=\"", as.integer(ml_param_samplingfreq), "\" fileName=\"powerposterior.log\">\n",
                          parameters_log, "\t\t</log>\n")
    
    if (powerposterior_vanilla) {
      tree_log <- paste0("\t\t\t<treeModel idref=\"treeModel\"/>\n")
    }
    mle_treefilelog <- paste0("\t\t<logTree id=\"MLETreeFileLog\" logEvery=\"", as.integer(ml_tree_samplingfreq), 
                              "\" nexusFormat=\"true\" fileName=\"powerposterior.trees\" sortTranslationTable=\"true\">\n", 
                              tree_log, "\t\t</logTree>\n")
    
    mle <- paste0(mle, mle_log, mle_filelog, mle_treefilelog, "\t</marginalLikelihoodEstimator>\n")
    x <- append(x, mle, after = length(x) - 1)
  }
  
  report <- paste0("\t<report>\n\t\t<property name=\"timer\">\n",
                   "\t\t\t<mcmc idref=\"mcmc\"/>\n",
                   "\t\t</property>\n",
                   "\t</report>")
  x <- append(x, report, after = length(x) - 1)
  
  cat(x, file = paste0(writeto_folder_path, "/", file_name, ".xml"), sep = "\n")
}
