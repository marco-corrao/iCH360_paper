This repository contains code and data required to reproduce all results in: 

_*A compact model of *Escherichia Coli* core and biosynthetic metabolism_

availble at

https://arxiv.org/abs/2406.16596


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
#### [./Knowledge_graph](./Knowledge_graph)
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



## References
 1. King, Z. A. et al. Escher: A Web Application for Building, Sharing, and Embedding Data-Rich Visualizations of Biological Pathways. PLOS Computational Biology 11, e1004321 (2015).
 2. Keseler, I. M. et al. The EcoCyc database: reflecting new knowledge about _Escherichia coli_ K-12. Nucleic Acids Res 45, D543–D550 (2017).