This folder contains the various models used in the paper

Subfolders and content:
### `./iCH360`
contains the main model (`Escherichia_coli_iCH360`). In addition, the folder contains a variant of the  model where only catalytic relationships annotated as primary were used to construct gene-protein-reaction (GPR) rules (`Escherichia_coli_iCH360_primary_catalysis_only`)

### `./iCH360red`
contains the iCH360red submodel. iCH360red is a curated reduction of iCH360 containing 18 fewer metabolic reactions (most of which are related to fatty acid biosynthesis) and is designed to be amenable for EFM enumeration and analysis. see `../EFM/ich360red_construction` for the code required to build this model.

### `./EC_iCH360`
Contain the Enzyme-Constrained (EC) variant of the model, constructed in sMOMENT fomat [1]. The folder contains two versions of the enzyme contrained model. The first (`EC_iCH360_unadjusted_kcats.xml`) is parametrised with the original set of turnover numbers from [2] (see `../Enzyme_Constraints/EC_model_building_pipeline` for the code used to build this model ). The second (`EC_iCH360_fitted_kapps.xml`)is parameterised with the adjusted apparent turnover numbers, which were fitted to experimental measurements of enzyme abundance (see paper for more information and `../Enzyme_Constraints/kcat_adjustment_and_enzyme_allocation_predictions` for the code used to produce this model.)

### model_tables
A collection of miscellaneous model data (e.g. reaction, genes, reaction-enzyme-subunit mappings, etc.), in tabular form, useful for data parsing/mapping pipelines e.g. in R

#### References
1. Bekiaris, P. S. & Klamt, S. Automatic construction of metabolic models with enzyme constraints. BMC Bioinformatics 21, 19 (2020).
2. Heckmann, D. et al. Kinetic profiling of metabolic specialists demonstrates stability and consistency of in vivo enzyme turnover numbers. Proceedings of the National Academy of Sciences 117, 23182–23190 (2020).