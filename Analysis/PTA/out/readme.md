# Output of Probabilistic Thermodynamic Analysis (PTA)
## [pta_fluxes.csv](./pta_fluxes.csv)
Fluxes and thermodynamic state computed by PTA:
- **v**: flux in mmol/gDW/h
- **drg0**: standard Gibbs free energy of reaction, in kJ/mol
- **drg**: Gibbs free energy of reaction, in kJ/mol
- **z_drg0**: z-score for the deviation of the standard Gibbs free energy of reaction from its prior value
- **z_drg0**: z-score for the deviation of the  Gibbs free energy of reaction from its prior value
- **sp_drg**: shadow price of the thermodynamic constraint in the PTA solutiojn

## [pta_metabolite_concentrations](./pta_metabolites_concentrations.csv)
Metabolite concentration computed by PTA
- **conc**: predicted concentration, in mM
- **z_log_c**: z-score for the deviation of the predicted log concentration from the prior

