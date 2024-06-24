This folder contains a collection of custom code used throughout the repo.
#### `./model_assembly_utils.py`
contains useful functions for model building
#### `./visualisation_utils.py`
contains useful functions for visualsiation of the model in Escher [1]
#### `./fba_utils.py`
useful functions for stoichiometric modelling, including the generation of production envelopes.
#### `./biocyc_query_utils.py`
contains a number of functions automating HTTP querying of the EcoCyc database [2]
#### `./graph_utils.py`
contains a number of functions that were used to assemble and annotate the knowledge graph
#### `./EC_utils.py`
contains a number of functions useful to construct, handle, and optimize the enzyme-constrained model, including code for saturation FBA (satFBA) simulations [3].
#### `./EQuilibrator_utils.py`
contains code required to implement manually the thermodynamic constant estimaton pimplemented by EQuilibrator (see Supplementary Information in [4])

### References
1. King, Z. A. et al. Escher: A Web Application for Building, Sharing, and Embedding Data-Rich Visualizations of Biological Pathways. PLOS Computational Biology 11, e1004321 (2015).
2. Keseler, I. M. et al. The EcoCyc database: reflecting new knowledge about Escherichia coli K-12. Nucleic Acids Res 45, D543–D550 (2017).
3. Müller, S., Regensburger, G. & Steuer, R. Resource allocation in metabolic networks: kinetic optimization and approximations by FBA. Biochemical Society Transactions 43, 1195–1200 (2015).
4. Beber, M. E. et al. eQuilibrator 3.0: a database solution for thermodynamic constant estimation. Nucleic Acids Research 50, D603–D609 (2022).
