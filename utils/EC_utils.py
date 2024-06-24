'''
A collection of utilities to add/harness enzyme constraints in the metabolic models
'''

import cobra
import numpy as np
import pandas as pd
import re
import tqdm
def split_reversible_reactions(model):
    '''
    Splits reversible reactions into two irreversible reactions
    '''
    model_splitted=model.copy()
    r_ids=[r.id for r in model_splitted.reactions]
    for reaction_id in r_ids:
        reaction=model_splitted.reactions.get_by_id(reaction_id)
        if reaction.reversibility or (reaction in model_splitted.boundary):
            new_reaction = reaction.copy()
            new_reaction.id = reaction.id + '_bw'

            

            new_reaction.upper_bound = -reaction.lower_bound
            new_reaction.lower_bound = 0
            new_reaction.add_metabolites({metabolite: -stoichiometry for metabolite, stoichiometry in reaction.metabolites.items()}, combine=False)
            model_splitted.add_reactions([new_reaction])

            #Change the original
            reaction.lower_bound = 0
            reaction.id = reaction.id + '_fw'
        else:
            reaction.id = reaction.id + '_fw'
    return model_splitted

def backward_rxn(rxn_id):
    return rxn_id + '_bw'
def forward_rxn(rxn_id):
    return rxn_id+'_fw'
def is_reversible(model,rxn_id):
    all_rxns=[rxn.id for rxn in model.reactions]
    if forward_rxn(rxn_id) in all_rxns and backward_rxn(rxn_id) in all_rxns:
        return True
    else:
        return False
def strip_direction(rxn_id):
    return re.sub('_fw$|_bw$','',rxn_id)
def create_ec_model_from_table(model,ec_table,method='sMOMENT',verbose=False):
    '''
    Creates an enzyme-constrained model from a metabolic model and an enzyme table
    '''

    model_ec=model.copy()
    if method=='sMOMENT':
        # Split reversible reactions
        model_ec=split_reversible_reactions(model_ec)
        #create enzyme pool pseudometabolite
        enzyme_pool_metabolite=cobra.Metabolite(id='enzyme_pool',name='enzyme pool pseudometabolite',compartment='e')
        model_ec.add_metabolites([enzyme_pool_metabolite])
        #Loop through the EC table, and add the cost to each reaction
        for i,row in ec_table.iterrows():
            rxn_id=row['reaction_id']
            direction=row['direction']
            cost=row['ec_cost']
            enzyme=row['enzyme']

            if direction in ['fwd','fw','forward']:
                try:
                    r=model_ec.reactions.get_by_id(forward_rxn(rxn_id))
                except:
                    if verbose:
                        print(f'Warning: EC table quotes backward value for {rxn_id}, but no backward reaction found in model. Skipping.')
                    r=None
            elif direction in ['bwd','bw','backward']:
                try:
                    r=model_ec.reactions.get_by_id(backward_rxn(rxn_id))
                except:
                    if verbose:
                        print(f'Warning: EC table quotes backward value for {rxn_id}, but no backward reaction found in model. Skipping.')
                    r=None
            else:
                raise ValueError('direction must be forward or backward')
            if r is not None:
                r.add_metabolites({enzyme_pool_metabolite:-cost},combine=True)
                #Note down which enzyme the cost is based on
                r.annotation['smoment_enzyme']=enzyme
                r.annotation['smoment_mw']=row['MW']
                r.annotation['smoment_kcat_per_s']=row['EC_kcat_value']
        #If any reaction is left without cost, set thgat to 0
        for r in model_ec.reactions:
            if enzyme_pool_metabolite not in r.metabolites.keys():
                model_ec.reactions.get_by_id(r.id).add_metabolites({enzyme_pool_metabolite:0},combine=True)


        # Create enzyme pool pseudoreaction
        v_pool_ub=np.max(ec_table['enzyme_pool_bound_g_gDw']) #note max is redundant, all rows should be the same value
        R_pool=cobra.Reaction(id='enzyme_pool_supply',name='enzyme pool pseudoreaction',lower_bound=0,upper_bound=v_pool_ub)
        R_pool.add_metabolites({enzyme_pool_metabolite:1})
        model_ec.add_reactions([R_pool])

        return model_ec

def get_reaction_costs(ec_model,sigma=1):
    '''
    Returns a dataframe with the reaction costs
    '''
    out={}
    for r in ec_model.reactions:
        if ec_model.metabolites.enzyme_pool in r.metabolites:
            out[r.id]=-r.metabolites[ec_model.metabolites.enzyme_pool]/sigma
        else:
            out[r.id]=0
    return out

def optimize_ec_model(ec_model,sigma=None,proteome_bound=None):
    '''
    Optimize ec_model and Post-process the solution
    '''

    with ec_model as model:
        enzyme_pool=model.metabolites.get_by_id('enzyme_pool')
        enzyme_pool_supply=model.reactions.get_by_id('enzyme_pool_supply')
        if sigma is not None:
            if isinstance(sigma,(int,float)):
                # In this case, we assume that the provided sigma is an average saturation to apply to all reactions
                for r in model.reactions:
                    if r.id != enzyme_pool_supply.id and enzyme_pool in r.metabolites.keys():
                        r.add_metabolites({enzyme_pool:(r.metabolites[enzyme_pool])/sigma},combine=False)
            elif isinstance(sigma,dict):
                # In this case, we assume that the provided sigma is a dictionary of saturation values for a subset of specified reactions
                for r_id in sigma.keys():
                    try:
                        r=model.reactions.get_by_id(r_id)
                    except:
                        print(f'Warning: Reaction {r_id} not found in model. Skipping.')
                        continue
                    if r.id != enzyme_pool_supply.id and enzyme_pool in r.metabolites.keys():
                        r.add_metabolites({enzyme_pool:(r.metabolites[enzyme_pool])/sigma[r.id]},combine=False)
                
        
        if proteome_bound is not None:
            
            enzyme_pool_supply.upper_bound=proteome_bound

        #SOLVE
        sol=model.optimize()
        if sol.status!='optimal':
            print(f'Warning: Solution status is {sol.status}.')
            return None
        #Post-process the solution
        net_rxn_ids=list(set([strip_direction(r.id) for r in ec_model.reactions if r.id != 'enzyme_pool_supply']))
        out=[]

        enzyme_costs=get_reaction_costs(ec_model)
        for rxn_id in net_rxn_ids:
            rxn_fw=forward_rxn(rxn_id)
            rxn_bw=backward_rxn(rxn_id)
        
            if 'smoment_enzyme' in ec_model.reactions.get_by_id(rxn_fw).annotation: 
                enzyme=ec_model.reactions.get_by_id(rxn_fw).annotation['smoment_enzyme'] 
            elif rxn_bw in ec_model.reactions and 'smoment_enzyme' in ec_model.reactions.get_by_id(rxn_bw).annotation:
                enzyme=ec_model.reactions.get_by_id(rxn_bw).annotation['smoment_enzyme']
            else:
                enzyme='NA'
            if is_reversible(ec_model,rxn_id):
                flux_fw=sol.fluxes[rxn_fw]
                flux_bw=sol.fluxes[rxn_bw]
                net_flux=flux_fw-flux_bw

    
                enzyme_cost_fw=enzyme_costs[rxn_fw]*flux_fw
                enzyme_cost_bw=enzyme_costs[rxn_bw]*flux_bw
                enzyme_cost=enzyme_cost_fw+enzyme_cost_bw

                
            else:
                rxn_fw=forward_rxn(rxn_id)
                flux_fw=sol.fluxes[rxn_fw]
                net_flux=flux_fw

                
                enzyme_cost=enzyme_costs[rxn_fw]*flux_fw

        
            out.append({'reaction_id':rxn_id,'flux':net_flux,'enzyme':enzyme,'enzyme_abundance_g_gDW':enzyme_cost})
        out.append({'reaction_id':'enzyme_pool_supply','flux':sol.fluxes['enzyme_pool_supply'],'enzyme':'NA','enzyme_abundance_g_gDW':0})
    out_df=pd.DataFrame([pd.Series(row) for row in out])
    return out_df.set_index('reaction_id')
def compute_net_fluxes(flux_series):
    directional_reaction_ids=(flux_series.index.tolist())
    undirectional_reaction_ids=list(set([strip_direction(x) for x in directional_reaction_ids]))
    out=pd.Series(index=undirectional_reaction_ids,data=0.)
    for rxn_id in out.index:
        if forward_rxn(rxn_id) in flux_series.index:
            out[rxn_id]+=flux_series[forward_rxn(rxn_id)]
        if backward_rxn(rxn_id) in flux_series.index:
            out[rxn_id]-=flux_series[backward_rxn(rxn_id)]
    return out
                              


def compute_efm_yield_cost(efm:pd.Series,
                           ec_cost:dict=None,
                           biomass_id='Biomass',
                           carbon_uptake_id='EX_glc__D_e',
                           flux_supporting=[],
                           verbose=False):
    #Normalise efm by biomass
    efm_norm=efm/efm[biomass_id]

    #compute yield
    efm_yield=np.abs(1/efm_norm[carbon_uptake_id])

    #compute cost
    efm_cost=0

    enzyme_costs=ec_cost

    for r_id in efm.index:
        if efm[r_id]==0 or r_id==biomass_id:
            pass
        elif 'EX_' in r_id or r_id=='ATPM':
            pass
        elif efm[r_id]>0:
            cur_r_cost=enzyme_costs[forward_rxn(r_id)]
            efm_cost+=cur_r_cost*efm_norm[r_id]
        elif efm[r_id]<0:
            cur_r_cost=enzyme_costs[backward_rxn(r_id)]
            efm_cost+=cur_r_cost*(-efm_norm[r_id])

    out={'yield':efm_yield,'cost':efm_cost}
    for flux in flux_supporting:
        out[flux]=efm_norm[flux]
    
    return out

def compute_efm_df_yield_cost(efms_df:pd.DataFrame,
                           ec_cost:dict=None,
                           biomass_id='Biomass',
                           carbon_uptake_id='EX_glc__D_e',
                           flux_supporting=[],
                           verbose=False):
    out_df=pd.DataFrame(index=efms_df.index,
                        columns=['yield','cost']+flux_supporting)
    for i,row in efms_df.iterrows():
        out_df.loc[i]=compute_efm_yield_cost(row,
                                             ec_cost=ec_cost,
                                             biomass_id=biomass_id,
                                             carbon_uptake_id=carbon_uptake_id,
                                             flux_supporting=flux_supporting,
                                             verbose=verbose)
    

    return out_df

def run_satFBA(ec_model, 
               saturable_reaction_id,
               n=100,
               log_spaced=False,
               saturation_min=1e-3,
               saturation_max=0.99,
               sigma=1,
               proteome_bound=None,
               carbon_uptake_reaction_id='EX_glc__D_e',
               biomass_reaction_id='Biomass',
               ):
    '''
    Run saturation FBA by simulating the an ensamble of EC models under varying saturation levels
    '''
    
    if log_spaced:
        saturation_values=np.logspace(start=np.log10(saturation_min),
                                      stop=np.log10(saturation_max),
                                      num=n)
    else:
        saturation_values=np.linspace(start=saturation_min,
                                      stop=saturation_max,
                                      num=n)
    reactions_saturation={r.id:sigma for r in ec_model.reactions if r.id !=saturable_reaction_id}
    reactions_saturation[saturable_reaction_id]=1 #Placeholder, will be changed in the loop

    out={s:None for s in saturation_values}
    for s in tqdm.tqdm(saturation_values):
        reactions_saturation[saturable_reaction_id]=s
        with ec_model as model:
            sol=optimize_ec_model(model,
                                  sigma=reactions_saturation,
                                  proteome_bound=proteome_bound)
            out[s]=sol
    
    #Parse into a tall dataframe
    merged_df=pd.DataFrame()
    for s in out.keys():
        cur_s_results=out[s]
        cur_s_results['saturation']=s
        cur_s_results['yield']=-cur_s_results.loc[biomass_reaction_id,'flux']/cur_s_results.loc[carbon_uptake_reaction_id,'flux'] #gDW per mmol glucose
        cur_s_results['saturable_reaction_proteome_fraction']=cur_s_results.loc[strip_direction(saturable_reaction_id),'enzyme_abundance_g_gDW']/cur_s_results.loc['enzyme_pool_supply','flux']
        merged_df=pd.concat([merged_df,cur_s_results],axis=0)

    return merged_df

