This folder contains all files and code used to generate estimates for thermodynamic constants (free energies of reactions and formation) mapped to reaction and metabolites in the model. Our procedure builds on the component-contribution method from [1] and the data from the EQuilibrator database [2]. However, rather than using the automated Equilibrator API to produce our estimates, we followed a manual procedure enabling us to generate estimates for reactions involving conserved protein side-groups (which are not convered in the EQuilibrator database).
#### `./free_energy_estimation`
contains code and files to reproduce our custom thermodynamic constants estimation pipeline
#### `free_energy_estimates`
Contains all the produced estimates of the thermodynamic constants (see readmes inside the subfolder for more information on the different quantities)
## Note
In order to run the pipeline, we produced a custom compounds cache to use with equilibrator, where protein side-groups were replaced with chemical groups approximating the chemical environment for the compound (see our paper for more details). To run the estimation pipeline, a cust9om-made EQuilibrator local cache file [custom_compounds_cache.sqlite](`custom_compounds_cache.sqlite`) must first be found by the script in the `./free_energy_estimation/data` folder. Due to the large size of the script, this file is only available in the ZENODO record of this repository. Further, running this script requires the `equilibrator-assets` package to be installed (see https://gitlab.com/equilibrator/equilibrator-assets)
#### References
1. Noor, E., Haraldsdóttir, H. S., Milo, R. & Fleming, R. M. T. Consistent Estimation of Gibbs Energy Using Component Contributions. PLOS Computational Biology 9, e1003098 (2013).

2. Beber, M. E. et al. eQuilibrator 3.0: a database solution for thermodynamic constant estimation. Nucleic Acids Research 50, D603–D609 (2022).
