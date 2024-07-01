This folder contains the code required to reproduce all analysis presented in the paper, namely:
- statistics on the available annotation in *i*CH360 with that of its parent *i*ML1515 [1] (`./compare_annotation_w_iml1515`)
- Catalytic disruption analysis and validation with mutant fitness data from [2] (`./catalytic_disruption_analysis`)
- Screening of the EFMs of *i*CH360red for growth and yield (`./EFM_growth_yield_screening`). Note that the output files for this analysis are quite large and are not included in this repository, but can be generated running `/EFM_growth_yield_screening/compute_efm_cost_yield.ipynb`. This required the EFM enumeration files to have also been generated in advance (see `../EFM/efm_enumeration`)
- saturation FBA analysis [3] (`./satFBA`)
- Probabilistic Max-min driving force (MDF) analysis (adapted from [4]) (`./MDF`)
- Probabilistic thermodynamic analysis [5] (`./PTA`).

**Note**

To run the Probabilistic thermodynamic analysis script, you will need to install the `pta` package. To reproduce results in this repositry, it suffices to install the package without the sampler dependencies (whose einstallation can be more involved) by running
```
PTA_NO_SAMPLERS=1 pip install pta             # Linux, OSX
set "PTA_NO_SAMPLERS=1" && pip install pta    # Windows
```


#### References
1. Monk, J. M. et al. iML1515, a knowledgebase that computes Escherichia coli traits. Nat Biotechnol 35, 904–908 (2017).
2. Price, M. N. et al. Mutant phenotypes for thousands of bacterial genes of unknown function. Nature 557, 503–509 (2018).
3.Noor, E. et al. Pathway Thermodynamics Highlights Kinetic Obstacles in Central Metabolism. PLOS Computational Biology 10, e1003483 (2014).
4. Müller, S., Regensburger, G. & Steuer, R. Resource allocation in metabolic networks: kinetic optimization and approximations by FBA. Biochemical Society Transactions 43, 1195–1200 (2015).
5. Gollub, M. G., Kaltenbach, H.-M. & Stelling, J. Probabilistic thermodynamic analysis of metabolic networks. Bioinformatics 37, 2938–2945 (2021).


