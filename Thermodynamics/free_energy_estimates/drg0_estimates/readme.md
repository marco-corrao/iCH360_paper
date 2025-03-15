### Standard Gibbs free energies of reactions ($\Delta_rG^{\circ}$)
In the component contribution framework, the standard Gibbs free energies of reactions ($\Delta_rG^{\circ}$) are estimated probabilistically as a multivariate gaussian random vector:
$$\Delta_rG^{\circ} \sim N(\bar{\Delta_rG}^{\circ}, \Sigma)$$
 where $\bar{\Delta_rG}^{\circ}$ is a the vector of mean estimates and  $\Sigma=QQ^T$, is the covariance matrix of the estimates. Here, $Q$ is a full-column-rank square root of $\Sigma$ (i.e. the number of columns of $Q$ corresponds the the rank of $\Sigma$)

### Transformed Standard Gibbs free energies of reactions ($\Delta_rG'^{\circ}$)
The *transformed* standard Gibbs free energies of reactions ($\Delta_rG'^{\circ}$) account for compartment-specific pH, pMG, ionic strength, potential and temperature. These are computed as a linear transformation untransformed estimates (\Delta_rG^{\circ}):
$$\Delta_rG'^{\circ}=\Delta_rG^{\circ}+(S^\top \epsilon_l)$$
Here, $S$ is the stoichiometric matrix of the network, $\epsilon_l$ is a vector of legendre transform (one per metabolite) that corrects for compartment-specific parameters. See [1] for details on the computation of these this correction vector. Note that this linear transform only affects the mean (but not the covariance) of the estimates, giving:
$$\Delta_rG'^{\circ}\sim N(\bar{\Delta_rG'^{\circ}}, \Sigma)$$
where $\bar{\Delta_rG'^{\circ}} \equiv \bar{\Delta_rG^{\circ}}+(S^\top \epsilon_l)$
### Multi-Compartment corrections
To account for reactions occurring across compartments, we further add a constant correction the the random vector of estimates. Hence, we compute the multi-compartment-corrected transformed Gibbs free energies of reactions ($\Delta_rG'^{\circ}_{MC}$) as: $$\Delta_rG'^{\circ}_{MC}=\Delta_rG^{\circ}+(S^\top \epsilon_{MC})$$
where $\epsilon_{MC}$ is another vector of corrections (again, one per metabolite), used to correct for reactions occurring across compartments. Again, see [1] for details on the computation of these this correction vector.
Similarly to before, the multi-compartment correction only affects the mean estimate, yielding:$$\Delta_rG'^{\circ}_{MC} \sim N(\bar{\Delta_rG'^{\circ}}_{MC},~\Sigma)$$
where $\bar{\Delta_rG'^{\circ}}_{MC} \equiv \bar{\Delta_rG'^{\circ}}+(S^\top \epsilon_{MC})$

The compartment parameters used for Legendre transform and multi-compartment corrections are

|Compartment|pH  |pMg|I  |phi   |T     |
|-----------|----|---|---|------|------|
|e          |6.90|3.0|0.1|0.000 |310.15|
|p          |6.91|3.0|0.1|-0.001|310.15|
|c          |7.80|3.0|0.1|-0.086|310.15|

### File description
- `drg0_mean_df.csv` : The mean of the standard Gibbs free energies of reaction ($\bar{\Delta_rG^{\circ}}$)
- `drg0_cov_df.csv` : the covariance matrix $\Sigma$
- `drg0_cov_sqrt.csv` : the square root of the covariance matrix, $Q$
- `eps_l.csv`: the vector of legendre corrections ($$\epsilon_l)
- `legendre_corrections.csv` : the correction term $(S^\top \epsilon_l)$
- `drg0_prime_mean.csv` : the mean of the transformed standard free energies of reactions ($\bar{\Delta_rG'^{\circ}}$)
- `eps_mc.csv` : the vector of multi-compartment corrections ($\epsilon_{MC}$)
- `multicompartment_corrections.csv` the correction term $(S^\top \epsilon_{MC})$
- `drg0_prime_mean_mcc.csv` the mean of multi-compartment-corrected (MCC) transformed standard Gibbs free energies of reaction ($\bar{\Delta_rG'^{\circ}}_{MC}$)
