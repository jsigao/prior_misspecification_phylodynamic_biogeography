# Supplemental Archive for: The Impact of Prior Misspecification on Bayesian Phylodynamic Inference of Biogeographic History
This archive contains the materials that are necessary and sufficient to replicate this study; it is divided into three subdirectories: [`data`](#data) (containing the sequence alignment and the associated sampling time and location data), [`analyses`](#analyses) (containing the `BEAST` phylodynamic analyses XML scripts), and [`scripts`](#scripts) (containing the analyses-setup and analyses-postprocessing `R` scripts), each described in detail below.
This archive is also available as an [Dryad repository](link here).

## <a name="data"></a>Data
The `data` subdirectory contains the sequence alignment and the associated sampling time and location information.
We organized the 11 datasets by virus first (corresponding to six immediate subdirectories), and then by subdatasets when there is more than one such for a given virus.
Each (sub)dataset subdirectory contains a `discrete_trait.txt` spreadsheet, which stores the sampling time and location information.
We also include the sequence alignment (as well as the GenBank accession numbers) for the datasets that we don't have access to the marginal distribution of phylogenies or the sequence alignment from the original study.

## <a name="analyses"></a> Analyses
The `analyses` subdirectory contains the `BEAST` XML scripts we used in the phylodynamic analyses.
This subdirectory is structured similarly with the `data` subdirectory in terms of the organization of the 11 datasets.
Each (sub)dataset directory contains two immediate subdirectories: `phylogeny` and `phylogeography`.
Each `phylogeny` subdirectory contains a `dataset_sample.trees` file, <!-- (not available in the GitHub repository due to size limit, but can be found in the Dryad repository) --> which is the corresponding marginal posterior distribution of phylogenies inferred using the sequence data and their associated sampling time (*i.e.*, result of the first step of the sequential phylodynamic inference), and a `dataset_MCC.tree` file, a summary phylogeny computed from the posterior distribution.
`dataset_sample.trees` is used in the second step of the `BEAST` sequential phylodynamic inference where we marginalized over the distribution of phylogenies to estimate the geographic model parameters, to infer the biogeographic history, and to assess geographic (prior)model fit.
`dataset_MCC.tree` is used in the data cloning analyses where we conditioned on the summary tree to evaluate informativeness of the prior relative to the data.
We include the `BEAST` XML scripts for the datasets that we don't have access to the marginal distribution of phylogenies from the original study so that we inferred it in this study.

The `phylogeography` subdirectory contains the `BEAST` XML scripts for the core model-evaluation phylodynamic analyses of this study.
This subdirectory is further structured into three nested layers: (1) `ctmc` and `exp_hyper`, (2) `asymmetric` and `symmetric`, and (3) `poisson_default` and `poisson_intermediate`, corresponding to the eight prior model combinations of: (1) default and alternative priors on the overall dispersal rate, (2) a symmetric and asymmetric rate matrix, and (3) default and alternative priors on the number of dispersal routes.

Then each of the eight prior model combination directories contains two subdirectories: `posterior` and `prior`.
The former contains `BEAST` XML scripts for the analyses estimating the geographic model parameters and biogeographic history first and then running power-posterior inferences to compute the marginal likelihood.
The latter contains `BEAST` XML scripts that set up analyses to infer the joint prior distribution.
In addition, `ctmc/symmetric/poisson_default` and `exp_hyper/symmetric/poisson_intermediate` directories have four other subdirectories&mdash;`completeHistory_fixedMCC`, `datacloning_fixedMCC_lheat5`, `datacloning_fixedMCC_lheat10`, and `datacloning_fixedMCC_lheat20`&mdash;containing `BEAST` XML scripts for performing the data-cloning analyses under 1, 5, 10, and 20 copies of geographic data, respectively.

## <a name="scripts"></a>Scripts
The `scripts` subdirectory contains the `R` scripts we used in this study; it is further divided into three subdirectories, corresponding to the three major types of processing we did (either before or after running the `BEAST` analyses), including:
* the scripts that automatically generate the `BEAST` XML scripts to set up all the analyses (`scripts/analyses_setup`),
* the scripts that perform posterior-predictive simulations using the inferred geographic model parameters and dated phylogenies (`scripts/history_simulation`), and
* the scripts that process the `BEAST` output files to produce summaries of geographic model parameters (`scripts/parameter_summary`).

### <a name="analyses_setup_scripts"></a>`BEAST` phylodynamic analyses setup
To thoroughly explore the statistical behavior of various prior model combinations for each of the 11 datasets, we had to set up hundreds of `BEAST` phylodynamic analyses.
We automated this process by developing a function in R that generates `BEAST` XML scripts under various model and inference configurations.
`scripts/history_simulation/BEAUtj.R` contains this function and `scripts/history_simulation/print_phylogeography_xml.R` serves as an example of an interface (that evokes the function) that allows user to configure the `BEAST` phylodynamic analyses one would like to set up (and as an example it can be executed immediately to recreate the [`analyses`](#analyses) subdirectory of this repository).

### <a name="history_simulation_scripts"></a>History Simulations
We simulated full dispersal history over the sampled dated phylogeny to assess the adequacy of geographic models.
`scripts/history_simulation/history_simulator_functions.R` is the interface script that a user can execute the `history_simulator` routine within it to perform posterior-predictive simulation (by setting the argument `conditional` to false).
Other scripts in this subdirectory contains subroutines that will be evoked by the interface routine.
One exception is `scripts/history_simulation/posterior_predictive_teststatistics_functions`, which provides functions to calculate posterior-predictive summary statistics (including the tipwise multinomial statistic and the parsimony statistic).

### <a name="parameter_summary_scripts"></a>Geographic Model Parameter Summary
`scripts/parameter_summary/get_BFs_counts.R` processes the `BEAST` parameter log file (and simulated history log file if available) to summarize the geographic model parameter estimates for each (prior)model, including pairwise dispersal routes, average dispersal rate, and number of pairwise dispersal events.
`scripts/parameter_summary/get_internal_prob_pairs.R` processes a pair of MCC trees to summarize the posterior probability of the maximum a posteriori state (determined by the first provided tree) at each of the shared internal nodes between the trees.
