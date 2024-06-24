In this folder, you can find the data and pipeline used to parse kinetic and proteomic data from other works, construct an enzyme constrained model in sMOMENT format [1], and fit its turnover number to experimental measurements of enzyme abundance.
#### `./data`
contains raw and parsed data from other works, including:
1. *In vivo* turnover number estimates from [1]
2. Measurements of polypeptide abundances across conditions from [2]
3. Polypeptide abundances from the integrated PAX database [3]
#### `./EC_model_building_pipeline`
contains the code to: parse kinetic data from [1]; create unique mappings between reactions and enzymes; compute enzyme-costs for each reaction/direction; use such table to construct the EC-model in sMOMENT [4] format.
#### `./enzyme_abundance_estimation`
Contains the code used to estimate enzyme counts (in copies/cell) and mass abundances (in g/gDW) from polypeptides counts in [2], using non-negative least-square (NNLS) estimation.
#### `./kcat_adjustment_and_enzyme_allocation_predictions`
contains the code required to implement our custom turnover number adjustment procedure (see paper for more information). The code also produces enzyme-allocation predictions pre- and post- adjustment, which are used to generate the relevant figures in the paper, and export a version of EC-*i*CH360 with the fitted turnover numbers.
#### References
1. Heckmann, D. et al. Kinetic profiling of metabolic specialists demonstrates stability and consistency of in vivo enzyme turnover numbers. Proceedings of the National Academy of Sciences 117, 23182–23190 (2020).
2. Schmidt, A. et al. The quantitative and condition-dependent Escherichia coli proteome. Nat Biotechnol 34, 104–110 (2016).
3. Huang, Q., Szklarczyk, D., Wang, M., Simonovic, M. & Mering, C. von. PaxDb 5.0: Curated Protein Quantification Data Suggests Adaptive Proteome Changes in Yeasts. Molecular & Cellular Proteomics 22, (2023).
4. Bekiaris, P. S. & Klamt, S. Automatic construction of metabolic models with enzyme constraints. BMC Bioinformatics 21, 19 (2020).

