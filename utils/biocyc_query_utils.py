'''
A collection of functions to query and retrieve information from biocyc.
'''

import requests
import xml
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np
import pickle
import os
import json
from tqdm import tqdm




def establish_biocyc_session(email, password):
    """
    Establishes a session with BioCyc and returns the session id.
    """
    s= requests.Session()
    url = 'https://websvc.biocyc.org/credentials/login/'
    data = {'email': email, 'password': password}
    s.post(url, data=data)
    return s

def split_db_id(input):
    '''
    split database identifier and object id (e.g. ECOLI:TRANS-RXN0-603 will return ECOLI, TRAN`S-RXN0-603)
    '''
    
    db_id=input.split(':')
    if len(db_id)==2:
        return db_id[0],db_id[1]
    else:
        return None,None
def get_biocyc_object(session, object_id,db='ECOLI'):
    """
    Returns the XML object from BioCyc.
    """
    url = f'https://websvc.biocyc.org/getxml?id={db}:{object_id}'
    response = session.get(url)
    #build xml element tree
    et=ET.fromstring(response.content)
    return et

def get_biocyc_id_from_foreign(session,foreign_id,foreing_db='Kegg',db='ECOLI'):

    """
    Returns the XML object from BioCyc.
    """
    url = f'https://websvc.biocyc.org/{db}/foreignid?ids={foreing_db}:{foreign_id}&fmt=json'
    response = session.get(url)
    data=json.loads(response.content)[0]
    if data['STATUS']==1:
        return [data['RESULTS'][i]['ID'] for i in range(len(data['RESULTS']))]
    else:
        return []
    
def rxn2enzyme(session,rxn_id,db='ECOLI'):
    '''
    Find all enzymes that catalyze a given reaction. Returns protein b IDs
    '''
    #get the reaction object
    rxn=get_biocyc_object(session,rxn_id,db=db)
    #get the enzymes
    path='.//Reaction/enzymatic-reaction/Enzymatic-Reaction/enzyme/Protein'
    enzymes=rxn.findall(path)
    enzymes_ids=list(map(lambda x: x.attrib['frameid'],enzymes))

    return enzymes_ids
def protein_components(session,protein_id,db='ECOLI'):
    '''
    Find all protein components of a given protein. Returns proteins as xml nodes
    '''
    prt=get_biocyc_object(session,protein_id,db=db)
    #get the protein components
    path='.//Protein/component'
    components=prt.findall(path)
    out={}
    for c in components:
        try:
            c_id=c.find('.//Protein').attrib['frameid']
            try: c_stoich=int(c.find('.//coefficient').text) 
            except:c_stoich=1
            out[c_id]=c_stoich
        except:
            print(f'Warning: failed to parse component {c.text}')
    return out
def pp_mw(session,pp_id,db='ECOLI'):
    '''
    Find the molecular weight (in kDa of a polypeptide
    '''
    pp=get_biocyc_object(session,pp_id,db=db)
    try:
        mw=float(pp.find('.//Protein/molecular-weight-seq').text)
    except:
        print(f'Warning: failed to parse molecular weight for {pp_id}')
        mw=np.nan
    return mw
def find_components_recursively(session,protein_id,objects,nodes,edges,db='ECOLI',cache=None,ignore=[]):
    '''
    Recursively find all protein components of a given protein. 
    '''
    protein_biocyc_id=f'{db}:{protein_id}'
    if cache is not None and protein_biocyc_id in cache.keys():
        obj=cache[protein_biocyc_id]
    else:
        obj=get_biocyc_object(session,protein_id,db=db)

    
    objects[protein_biocyc_id]=obj

    
    if len(obj.findall('.//Protein/component'))>0:
        #find the components
        components=protein_components(session,protein_id,db=db)
        for c_id,c_stoich in components.items():
            edges.append({'source':protein_id,
                         'target':c_id,
                         'weight':c_stoich,
                         'type':'subunit_composition',
                         'subtype':"requirement",
                         'notes':'',
                         'references':''})
            objects,edges,nodes,child_type=find_components_recursively(session,c_id,objects,nodes,edges,db=db,cache=cache)
        protein_subtype='multimeric_protein'
    # Some proteins have modified/unmodified forms. In this case often only a single form of the two has the components. 
    # In this case, we simply add an edge with weight one to the other modified form
    elif len(obj.findall('.//Protein/unmodified-form/Protein'))==1:
        unmodified_protein_id=obj.findall('.//Protein/unmodified-form/Protein')[0].attrib['frameid']
        edges.append({'source':protein_id,
                    'target':unmodified_protein_id,
                    'weight':'NA',
                    'type':'protein_modification',
                    'subtype':"NA",
                    'notes':'',
                    'references':''})
        objects,edges,nodes,child_type=find_components_recursively(session,unmodified_protein_id,objects,nodes,edges,db=db,cache=cache)
        protein_subtype='modified_protein'
    else:
        protein_subtype='polypeptide'
        
    nodes.append({'id':protein_id,
                  'type':"protein",
                  'subtype':protein_subtype,
                  'biocyc_id':protein_biocyc_id}
                  )
    return objects,edges,nodes,protein_subtype

def find_all_rxns(session,protein_id,db='ECOLI'):
    '''
    Find all reactions catalyzed by a given protein
    '''
    #get the protein object
    prt=get_biocyc_object(session,protein_id,db=db)
    #navigate up to the uppermost complex
    while len(prt.findall('.//Protein/component-of'))>0:
        prt=get_biocyc_object(session,prt.find('.//Protein/component-of/Protein').attrib['frameid'],db=db)
    #get the reactions
    path='.//Protein/catalyzes/Enzymatic-Reaction/reaction/Reaction'
    rxns=prt.findall(path)
    rxns_ids=list(map(lambda x: x.attrib['frameid'],rxns))
    return rxns_ids



def pp2gene(session,pp_id) :
    try:
        obj=get_biocyc_object(session,pp_id)
    except:
        print(f'Unable to retrieve object for polypeptide {pp_id} from biocyc')
    genes=obj.findall('.//')
def get_dg0(session,cpd_id,type):
    try:
        obj=get_biocyc_object(session,cpd_id)
    except:
        print(f'Unable to retrieve object {cpd_id} from biocyc')
        return []
    conversion={'kcal/mol':4.184,'kj/mol':1}
    if type=='formation':
        out=[float(x.text) for  x in obj.findall('.//Compound/gibbs-0')]
    elif type=='reaction':
        out=[float(x.text)*conversion[x.attrib['units']] for  x in obj.findall('.//Reaction/gibbs-0')]
    return out

def flatten_list(x):
    """
    flatten a nested list of arbitrary depth
    """
    if isinstance(x,list):
        return [a for i in x for a in flatten_list(i)]
    else:
        return [x]
def find_all_parent_proteins(session,protein_id,db='ECOLI'):
    """
    Given a protein, find proteins in EcoCyc of which this protein is a component of.
    """
    biocyc_object=get_biocyc_object(session=session,object_id=protein_id,db=db)
    parent_proteins=[x.attrib['frameid'] for x in biocyc_object.findall('.//Protein/component-of/Protein')]
    if len(parent_proteins)==0:
        return protein_id
    else:
        return flatten_list([find_all_parent_proteins(session=session,protein_id=parent_protein,db=db) for parent_protein in parent_proteins])
    
def parse_regulation_info(biocyc_rxn_ids,biocyc_session,cache=None):
    regulation_entry_cache=dict()
    #Loop across provided biocyc reaction IDs
    for biocyc_rxn_id in tqdm(biocyc_rxn_ids):
        try:
            if cache is not None and biocyc_rxn_id in cache.keys():
                biocyc_rxn_object=cache[biocyc_rxn_id]
            else:
                object_db,object_id=split_db_id(biocyc_rxn_id)
                biocyc_rxn_object=get_biocyc_object(session=biocyc_session,
                                                                object_id=object_id,
                                                                db=object_db)
            #Now identify all enzymatic reactions associated with this object
            enzymatic_rxns=biocyc_rxn_object.findall('./Reaction/enzymatic-reaction/Enzymatic-Reaction')
            enzymatic_rxns_ids=[x.attrib['ID'] for x in enzymatic_rxns]
            #For each enzymatic reaction, identify the relevant regulatory entries
            for enzymatic_rxn_id in enzymatic_rxns_ids:
                enzymatic_rxn_object_db,enzymatic_rxn_object_id=split_db_id(enzymatic_rxn_id)
                enzymatic_rxn_object=get_biocyc_object(session=biocyc_session,
                                                                object_id=enzymatic_rxn_object_id,
                                                                db=enzymatic_rxn_object_db)
                #Now, find all regulatory entries associated with this
                regulation_entries=enzymatic_rxn_object.findall('./Enzymatic-Reaction/regulated-by/Regulation')
                regulation_entries_ids=[x.attrib['ID'] for x in regulation_entries]
                for regulation_id in regulation_entries_ids:
                    regulation_object_db,regulation_object_id=split_db_id(regulation_id)
                    regulation_object=get_biocyc_object(session=biocyc_session,
                                                        object_id=regulation_object_id,
                                                        db=regulation_object_db)
                    regulation_entry_cache[regulation_id]=regulation_object
        except:
            print(f'Unable to parse regulatory information from {biocyc_rxn_id}')
    return regulation_entry_cache