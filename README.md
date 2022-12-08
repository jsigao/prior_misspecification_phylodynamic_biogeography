# Supplementary Archive for: Model Misspecification Misleads Inference of the Spatial Dynamics of Disease Outbreaks
This archive contains the materials that are necessary and sufficient to reproduce this study; it is divided into four directories: 
(1) the [`data`](#data) directory contains the sequence alignment and the associated sampling time and location data;
(2) the [`analyses`](#analyses) directory contains the `BEAST` phylodynamic analyses XML scripts;
(3) the [`scripts`](#scripts) directory contains the `R` scripts for setting up and summarizing analyses, and;
(4) the [`program`](#program) directory contains a modified version of the `BEAST` software package that implements a robust function for checking the irreducibility of a geographic model, required to perform the analyses under the asymmetric geographic model in our study.
This archive is also available as a [Dryad repository](https://datadryad.org/stash/share/vTbeDwLq2uSL9rL4NCe_Cocp2bY7BgWTI2tUgoNrLDA).

## <a name="data"></a>Data
The `data` subdirectory contains the sequence alignment and the associated sampling time and location information.
We organized the 14 datasets by the source study first (corresponding to eight immediate subdirectories), and then by data(sub)sets when there is more than one such for a given study.
Each (sub)dataset subdirectory contains a `discrete_trait.txt` spreadsheet, which stores the sampling time and location information.
We also include the sequence alignment (as well as the GenBank accession numbers) for the datasets that we don't have access to the marginal distribution of phylogenies or the sequence alignment from the original study.

## <a name="analyses"></a> Analyses
The `analyses` subdirectory contains the `BEAST` XML scripts we used in the phylodynamic analyses.
This subdirectory is structured similarly with the `data` subdirectory in terms of the organization of the 14 datasets.
Each (sub)dataset directory contains two immediate subdirectories: `phylogeny` and `phylogeography`.
Each `phylogeny` subdirectory contains a `dataset_sample.trees` file, <!-- (not available in the GitHub repository due to size limit, but can be found in the Dryad repository) --> which includes the corresponding marginal posterior distribution of phylogenies inferred using the sequence data and their associated sampling time (*i.e.*, result of the first step of the sequential phylodynamic inference), and a `dataset_MCC.tree` file, a summary phylogeny computed from the posterior distribution.
For all datasets (except the SARS-CoV-2 Global and B.1.1.7 US datasets), `dataset_sample.trees` is used in the second step of the `BEAST` sequential phylodynamic inference where we marginalized over the distribution of phylogenies to estimate the geographic model parameters, to infer the biogeographic history, and to assess geographic (prior) model fit.
`dataset_MCC.tree` is used in the data cloning analyses where we conditioned on the summary tree to evaluate informativeness of the prior relative to the data, and in all the analyses for the SARS-CoV-2 Global and B.1.1.7 US datasets where we conditioned on the MCC summary tree to ensure numerical stability.
We also include the `BEAST` XML scripts for the phylogenetic analyses using the sequence data and sampling time to infer the marginal distribution of phylogenies in the `phylogeny` subdirectory of a study for which we don't have access to that tree distribution.

The `phylogeography` subdirectory contains the `BEAST` XML scripts for the core model-evaluation phylodynamic analyses of this study.
This subdirectory is further structured into three nested layers: (1) `ctmc` and `exp_hyper`, (2) `asymmetric` and `symmetric`, and (3) `poisson_default` and `poisson_intermediate`, corresponding to the eight prior model combinations of: (1) default and alternative priors on the overall dispersal rate, (2) a symmetric and asymmetric rate matrix, and (3) default and alternative priors on the number of dispersal routes.

Then each of the eight prior model combination directories contains two subdirectories: `posterior` and `prior`.
The former contains `BEAST` XML scripts for the analyses estimating the geographic model parameters and biogeographic history first and then running power-posterior inferences to compute the marginal likelihood.
The latter contains `BEAST` XML scripts that set up analyses to infer the joint prior distribution.
In addition, `ctmc/symmetric/poisson_default` and `exp_hyper/symmetric/poisson_intermediate` directories have four other subdirectories&mdash;`completeHistory_fixedMCC`, `datacloning_fixedMCC_lheat5`, `datacloning_fixedMCC_lheat10`, and `datacloning_fixedMCC_lheat20`&mdash;containing `BEAST` XML scripts for performing the data-cloning analyses under 1, 5, 10, and 20 copies of geographic data, respectively.

## <a name="scripts"></a>Scripts
The `scripts` subdirectory contains the `R` scripts we used in this study; it is further divided into three subdirectories, corresponding to the three major types of processing we did (either before or after running the `BEAST` analyses), including:
* the scripts that can be used to generate the `BEAST` XML scripts to set up all the analyses (`scripts/analyses_setup`),
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
`scripts/parameter_summary/get_BFs_counts.R` processes the `BEAST` parameter log file (and simulated history log file if available) to summarize the geographic model parameter estimates for each (prior) model, including the Bayes factor support for each pairwise dispersal route, average dispersal rate, and number of pairwise dispersal events.
The functions for computing the prior probability that each dispersal route exist conditioning on the geographic model being irreducible (which is in turn used in computing the Bayes factor support) are available in `scripts/parameter_summary/delta_given_connected.R`.
By default, these functions fetch the number of strongly connected graphs from the pre-computed tables (`scripts/parameter_summary/num_connected_graphs.txt` and `scripts/parameter_summary/num_connected_digraphs.txt` for the symmetric and asymmetric models, respectively); these pre-computed tables were generated by `scripts/parameter_summary/ngraphs_connected_sym.gp` (symmetric model) and `scripts/parameter_summary/ngraphs_connected_asym.gp` (asymmetric model), which are [`PARI/GP`](https://pari.math.u-bordeaux.fr/) scripts that can be used to compute the number of strongly connected graphs with arbitrary number of areas (if the desired number of areas is greater than 100, the maximum number the pre-computed table currently provides).
`scripts/parameter_summary/get_internal_prob_pairs.R` processes a pair of MCC trees to summarize the posterior probability of the maximum a posteriori state (determined by the first provided tree) at each of the shared internal nodes between the trees.

## <a name="program"></a>Modified BEAST Program
We provide an executable for a modified version of the `BEAST` program, which implements a robust function for checking the irreducibility of a geographic model, required to perform the analyses under the asymmetric geographic model in our study.
The modified source code (from which we compiled the executable) is located in [the `master` branch, commit `1561762`](https://github.com/jsigao/beast-mcmc/commit/1561762d8c14d17ec2fdc4b1547bc562d6af2658) of a forked `BEAST` source-code repository.
Briefly, this executable can be run by invoking:
```
java -jar beast.jar -working ./analyses/sarscov2/phylogeography/ctmc/asymmetric/poisson_default/completeHistory_fixedMCC/sarscov2_041920_030820_coalExp_ucln_fine_completeHistory.xml
```
See [`BEAST` official tutorial website](http://beast.community/index.html) for further details about running `BEAST` analyses using a `Java` executable.