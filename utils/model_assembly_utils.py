import json
import cobra

def get_reactions_from_excher_map(map_path):
    with open(map_path,'r') as file:
        map=json.load(file)
    reactions_dict=map[1]['reactions']
    out=[]
    for k,v in reactions_dict.items():
        out.append(v['bigg_id'])
    return out

def import_reactions_from_other_model(reaction_ids:list,model_origin:cobra.Model,model_destination:cobra.Model):
    updated_model_destination=model_destination.copy()
    for r in reaction_ids:
        cur_r=model_origin.reactions.get_by_id(r)
        #Check if a reaction with the same id already exists in model_destination
        if cur_r.id in model_destination.reactions:
            print(f'Unable to import {cur_r.id}. A reaction with the same id already exist in the destination model')
        rxn_to_import=cur_r.copy()
        updated_model_destination.add_reactions([rxn_to_import])
    return updated_model_destination