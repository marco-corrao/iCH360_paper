This folder contains all the computational pipeline used to build the knowledge graph complementing the stoichiometric model, as well as the final graph structure in GML (.gml) and cytoscape (.cyjs) formats.

The graph was assembled using information retrieved from the EcoCyc [1] database and by manual curation. In the graph, nodes represent biological entities (reactions, proteins, genes or compounds) and edges represent (potentially quantitative) functional relationships between them, including catalysis, protein subunit composition, protein modification, and regulatory interactions. See our paper for a complete description of nodes and edges type in the graph.

The graph is provided in GML (.gml) and Cytoscape (.cyjs) formats, both of which can be read by various graph-handling packages. A popular example is NetworkX [2], which is the one we use to load, manipulate, and investigate the graph everywhere in this repository.

#### `./ich360_graph.gml`
The knowledge graph in GML format
#### `./ich360_graph.cyjs`
The knowledge graph in Cytoscape format
#### `./graph_assembly`
Contains all relevant files and code used to parse and assemble the graph

# Note
The pipeline used to construct the graph (`./graph_assembly/build_graph.ipynb`) uses extensive HTTP querying to the EcoCyc Database to retrieve information. To guarantee reproducibility, the output of all biocyc queries was serialised and cached, enabling the script to be run without querying the databse, and guaranteing reproducibility of results in this paper.

However, it is possible to run the pipeline from scratch, retrieving required data from the online database. To this end, one must:
1. create a file  with valid credentials (username and password) of an EcoCyc account. To this end, simply use `./graph_assembly/biocyc_username_and_password_template.csv` as a template, replace the placeholder username and password with yours, and save the modified file in the same directory as `./graph_assembly/biocyc_username_and_password.csv` (this file, not the template, will be called by the graph assembly script.)
2. Run the graph-assembly pipeline (`./graph_assembly/build_graph.ipynb`),changing `USE_CACHE=True` into `USE_CACHE=False` at the start of the script. More details on this can be found on the script itself.

### References
1. Keseler, I. M. et al. The EcoCyc database: reflecting new knowledge about Escherichia coli K-12. Nucleic Acids Res 45, D543–D550 (2017).
2. Hagberg, A. A., Schult, D. A. & Swart, P. J. Exploring Network Structure, Dynamics, and Function using NetworkX. (2008).
