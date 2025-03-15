import cobra
from cobra.flux_analysis import phenotype_phase_plane, pfba
from cobra.core import model
import pandas as pd
from typing import List, Tuple, Iterable, Any
import pandas as pd
import matplotlib.pyplot as plt

# modified version of plotting.plot_production_envelope
# let's us put multiple planes on the same plot. 
def production_envelope(
        wt_model: model,
        knockouts: Iterable[str],
        carbon_sources: Iterable[str],
        target_reaction: str,
        points: int = 20
) -> pd.DataFrame:
    """Calculates the data of the production envelope for a KO.

    :param wt_model: the Wild-Type model without any knockouts
    :param knockouts: the list of KOs
    :param carbon_sources: the list of carbon sources
    :param target_reaction: the reaction to place on the x-axis of the PPP
    :return: A DataFrame with the biomass yield ranges for each x-value
    """
    ko_model = wt_model.copy()
    for ko in knockouts:
        for k in ko.split('|'):
            ko_model.reactions.get_by_id(k).knock_out()
    for cs in carbon_sources:
        ko_model.reactions.get_by_id("EX_" + cs + "_e").lower_bound = -10.0 / len(carbon_sources)
    return phenotype_phase_plane.production_envelope(
        ko_model, target_reaction, points=points)


def plot_envelope(
        wt_model: model,
        knockouts: Iterable[str],
        carbon_sources: Iterable[str],
        target_reaction: str,
        ax: plt.Axes,
        label: str = None,
        color: Any = 'b',
        fill: bool = True,
        ls: Any = '-'):
    prod_env_df = production_envelope(
        wt_model=wt_model,
        knockouts=knockouts,
        carbon_sources=carbon_sources,
        target_reaction=target_reaction,
        points=200)

    label = label #or ','.join(carbon_sources)
    ax.plot(prod_env_df.flux_minimum,
            prod_env_df[target_reaction],
            c=color, ls=ls, lw=2, label=label)
    ax.plot(prod_env_df.flux_maximum,
            prod_env_df[target_reaction],
            c=color, ls=ls, lw=2, label='')

    if fill:
        ax.fill_betweenx(
            prod_env_df[target_reaction].values,
            prod_env_df.flux_minimum.values,
            prod_env_df.flux_maximum.values,
            linewidth=0,
            alpha=0.3,
            facecolor=color,
            label='')
    ax.grid()