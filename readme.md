This repository contains code and data required to reproduce all results in: 

_A compact model of *Escherichia Coli* core and biosynthetic metabolism_

available at

https://arxiv.org/abs/2406.16596

## Using *i*CH360

This repository is only intended to provide the files and tools to reproduce all results in the paper. If you wish to use *i*CH360 for your own work, please go to the model repo

https://github.com/marco-corrao/iCH360

where you'll find the most up-to-date version of the model and its variants.

## Navigating the repository
#### [./Model](./Model)
Contains the metabolic models (in `JSON` and `SBML` formats) mentioned in the paper, namely:
- The main stoichiometric model, *i*CH360
- The enzyme-constrained model variant, EC-*i*CH360
- The reduced model variant *i*CH360red

#### [./Visualisation](./Visualisation)
Contains all the relevant metabolic maps, ready to be loaded in Escher for visualisation [1]:
- The full model map
- The compressed model map
- The maps for each metabolic subsystem
- The maps for the pathways not included in the model, but used to compute the equivalent biomass reaction used in the model.
#### [./Annotation](./Visualisation)
Contains annotation maps to the EcoCyc database [2].
#### [./Knowledge_Graph](./Knowledge_Graph)
Contains the computational pipeline used to build the knowledge graph complementing the stoichiometric model, as well as the final graph structure in GML (.gml) and cytoscape (.cyjs) formats
#### [./Analysis](./Analysis)
Contains the python scripts required to reproduce all analysis mentioned in the paper. More specific details are provided in each subfolder
#### [./EFM](./EFM)
Contains the pipeline for creating the reduced model variant *i*CH360red, as well as counting and enumerating its elementary flux modes (EFMs).
#### [./Enzyme_Constraints](./Enzyme_Constraints)
Contains the data and scripts used to construct the enzyme constrained model EC-*i*CH360 and fit its turnover numbers to measured enzyme abundances.
#### [./Thermodynamics](./Thermodynamics)
Contains the file and script required to compute the estimates of thermodynamic constants for the reactions and metaboites in the model
#### [./Experimental_data](./Experimental_data)
Contains experimental data (proteomics, metabolomics, and fluxomics) from other works, mapped to the model.
#### [./External_database_data](./External_database_data)
Contains mappings between genes and polypeptides retrieved from the EcoCyc database [2].
#### [./Manuscript_Figures](./Manuscript_Figures)
Contains all the notebooks (in python and R) required to reproduce the figures in the paper.


## Dependencies
The following packages are used throughout the repo:
```
# General dependencies (used throughout)
cobra==0.29.0
numpy==1.24.0
scipy==1.10.1
pandas==1.5.3
matplotlib==3.7.1
seaborn==0.12.2
networkx==3.0
tqdm==4.65.0
requests==2.28.2

Additional dependencies are required to reproduce some analyses:

# EFM enumeration
efmtool==0.2.1

# turnover number fitting procedure:
gurobipy==11.0.1 #requires a valid GUROBI licence
casadi==3.6.3

# MDF analysis
gurobipy==11.0.1 #requires a valid GUROBI licence

# Thermodynamic constant estimation
equilibrator-api==0.4.7
equilibrator-assets==0.4.1
cvxpy==1.5.2
```
**Note**


The following steps may be needed to correctly run enkie and eQuilibrator to reproduce the thermodynamic analysis performed on iCH360.

1. For first time use of enkie, used in  [pta.ipynb](./Analysis/PTA/pta.ipynb), it may be necessary to create the folder ```~/.cache/enkie``` in your home directory
(see https://gitlab.com/csb.ethz/enkie/-/issues/1)

2. If issues  are encountered running eQuilibrator in [drg0_estimation.ipynb](./Thermodynamics/free_energy_estimation/drg0_estimation.ipynb), it may be necessary to manually save the files from the following repos: 
    - https://zenodo.org/records/4128543 
    - https://zenodo.org/records/4013789
    -  https://zenodo.org/records/4010930

    to ```~/.cache/equilibrator```

We kindly thank Benjamin Luke Coltman for suggesting these fixes.


## References
 1. King, Z. A. et al. Escher: A Web Application for Building, Sharing, and Embedding Data-Rich Visualizations of Biological Pathways. PLOS Computational Biology 11, e1004321 (2015).
 2. Keseler, I. M. et al. The EcoCyc database: reflecting new knowledge about _Escherichia coli_ K-12. Nucleic Acids Res 45, D543–D550 (2017).