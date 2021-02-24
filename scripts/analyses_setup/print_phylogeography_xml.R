# this script serves as the interface to automatically generate the BEAST XML scripts to set up all the phylodynamic analyses
source('./BEAUtj.R')

folder_path <- "../../data"
discrete_trait_paths <- list.files(folder_path, recursive = T, full.names = T, pattern = "*discrete_trait.txt")
all_folder <- grep("/phylogeography$", list.dirs("../../analyses", recursive = T, full.names = T), value = T)

for (i in 1:length(all_folder)) {
  
  scalar_prior <- c("ctmc", "exp_hyper")
  mat_sym <- c("asymmetric", "symmetric")
  mat_aver <- c("poisson_default", "poisson_intermediate")
  analysis_type <- c("posterior", "prior")
  
  terminal_dir <- paste0("/", apply(expand.grid(scalar_prior, mat_sym, mat_aver, analysis_type, KEEP.OUT.ATTRS = F, 
                                                stringsAsFactors = F), 1, function(x) paste(x, collapse = "/")))
  
  for (j in seq_along(terminal_dir)) {
    
    folder_path <- paste0(all_folder[i], terminal_dir[j])
    if (!dir.exists(folder_path)) { # don't overwrite
      dir.create(folder_path, recursive = T) # create dir first
    }
    
    if (length(list.files(folder_path, pattern = ".xml$")) == 0) { # don't overwrite
      
      fixed_tree_path <- list.files(paste0(unlist(strsplit(output_dir, "/phylogeography")), "/phylogeny"), recursive = T, full.names = T, pattern = "*MCC.tree")
      if (length(fixed_tree_path) != 1) {
        stop ("can not find the fixed tree path")
      }
      
      generate.phylogeography.xml(folder_path = folder_path, writeto_folder_path = folder_path, 
                                  discrete_trait_path = discrete_trait_paths[i],
                                  fixed_tree_path = ifelse(length(fixed_tree_path) > 0 && grepl("fixedMCC", output_dir), fixed_tree_path, NULL),
                                  date_direction = ifelse(grepl("sarscov2", output_dir), "backwards", "forwards"), 
                                  date_units = ifelse(grepl("sarscov2", output_dir), "days", "years"), 
                                  posterior = !grepl("/prior", folder_path), 
                                  complete_history = grepl("/complete_history|/completeHistory_fixedMCC", folder_path),
                                  symmetry = grepl("/symmetric", folder_path), ctmc = grepl("/ctmc", folder_path), 
                                  poisson_default = grepl("/poisson_default", folder_path), 
                                  bssvs = !grepl("/nobssvs", folder_path), poisson = !grepl("/nopoisson", folder_path), lheat = 1,
                                  mcmc_chainlength = 10000000 + 10000000 * grepl("/posterior", folder_path), mcmc_samplingfreq = 2000, 
                                  mcmc_screen_samplingfreq = 1000, mcmc_tree_samplingfreq = 2000, mcmc_adaptationDelay = 200000, 
                                  powerposterior_vanilla = grepl("/posterior_fixedMCC", output_dir),
                                  ml_chainlength = 250000, ml_burnin = 100000, ml_samplingfreq = 100, ml_param_samplingfreq = 1000, 
                                  ml_tree_samplingfreq = 10000, stones_num = 63, 
                                  tree_operator_weight = 35, rootfreq_proposal_weight = 1.25, clockrate_proposal_weight = 16, 
                                  clockratemean_proposal_weight = 0.9, rates_proposal_weight = 30, indicators_proposal_weight = 24)
      
    }
  }
}
