'''
A collection of function used to build the underlying graph structure of ich360 reaction-enzyme-gene mapping
'''
import requests
import xml
import biocyc_query_utils
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np
import pickle
import os
from tqdm import tqdm
import networkx as nx
from pyvis.network import Network
import re



def add_rxns_to_graph(session,bigg_rxns,bigg2biocyc_map,graph,db='ECOLI',cache=None):
    objects=graph['objects']
    nodes=graph['nodes']
    edges=graph['edges']
    for rxn_id in tqdm(bigg_rxns):
        if rxn_id in bigg2biocyc_map.keys():
            biocyc_rxn_id=bigg2biocyc_map[rxn_id]
        else:
            biocyc_rxn_id=None
        node_id=f'bigg:{rxn_id}'
        nodes.append({'id':node_id,
                      'type':'reaction',
                      'subtype':"NA",
                      'biocyc_id':biocyc_rxn_id}
                      )
        #If we do have a biocyc ID, we proceed with the search for components (enzymes, proteins etc)
        if biocyc_rxn_id is not None:
            cur_db,object_id=biocyc_query_utils.split_db_id(biocyc_rxn_id)
            if cur_db!=db:
                print(f'Skipping {rxn_id} as the provided Biocyc ID is not in the specified database ({db})')
            else:
                try:
                    if cache is not None and biocyc_rxn_id in cache.keys():
                        rxn_obj=cache[biocyc_rxn_id]
                    else:
                        rxn_obj=biocyc_query_utils.get_biocyc_object(session,object_id,db=db)
                    objects[biocyc_rxn_id]=rxn_obj

                    #Recursive components search
                    try:
                        enzymes=biocyc_query_utils.rxn2enzyme(session,object_id,db=db)
                        for enzyme_id in enzymes:
                            objects,edges,nodes,type=biocyc_query_utils.find_components_recursively(session,enzyme_id,objects,nodes,edges,db=db,cache=cache)
                            edges.append({'source':node_id,
                                          'target':enzyme_id,
                                          'weight':'NA',
                                          'type':'catalysis',
                                          'subtype':"NA",
                                          'notes':'',
                                          'references':''
                                          })
                    except:
                        print(f'Recursive components search failed for {biocyc_rxn_id}')

                except:
                    print(f'Failed to get {rxn_id} from biocyc or cache')
        
    return {'objects':objects,'nodes':nodes,'edges':edges}

def add_edge(graph_dict,
             parent_node,
             child_node,
             weight=1,
             type=None,
             metadata={},
             parse_from_biocyc=False,
             session=None,
             overwrite=False):
    '''
    Add an edge to the graph, potentially parsing recursive components from biocyc if a node doesn't exist
    '''
    if "subtype" not in metadata.keys():
        metadata['subtype']='NA'
    objects=graph_dict['objects']
    nodes=graph_dict['nodes']
    nodes_ids=[node['id'] for node in nodes]
    edges=graph_dict['edges']

    if parent_node['id'] not in nodes_ids:
        new_parent_node=True
    else:
        new_parent_node=False
    if child_node['id'] not in nodes_ids:
        new_child_node=True
    else:
        new_child_node=False
    # Does an edge between these two nodes already exist?
    if (parent_node['id'],child_node['id']) in [(edge['source'],edge['target']) for edge in edges]:
        edge_already_exists=True
        edge_ix=[(edge['source'],edge['target']) for edge in edges].index((parent_node['id'],child_node['id']))
    else:   
        edge_already_exists=False
    if edge_already_exists and not overwrite:
        print(f'Edge between {parent_node["id"]} and {child_node["id"]} already exists. Skipping')
    elif edge_already_exists and overwrite:
        print(f'Overwriting edge between {parent_node["id"]} and {child_node["id"]}')
        edges[edge_ix]=\
                {'source':parent_node['id'],
                  'target':child_node['id'],
                  'weight':weight,
                  'type':type,
                  }
        for key,value in metadata.items():
            edges[edge_ix][key]=value
    else:
        print(f'Adding new edge between {parent_node["id"]} and {child_node["id"]}')
        new_edge=({'source':parent_node['id'],
                    'target':child_node['id'],
                    'weight':weight,
                    'type':type,
                    })
        for key,value in metadata.items():
            new_edge[key]=value
        edges.append(new_edge)
    
    if new_child_node:
        if parse_from_biocyc:
            if session is None:
                print(f'No session provided. Unable to parse components of {child_node["id"]} from biocyc')
            else:
                print(f'Parsing components of {child_node["id"]} from biocyc as this node was not found in the graph')
                objects,edges,nodes,type=biocyc_query_utils.find_components_recursively(session,child_node['id'],objects,nodes,edges)
        else:
            print(f"Adding child node {child_node} to the graph dict")
            nodes.append(child_node)

    if new_parent_node:
        if parse_from_biocyc:
            if session is None:
                print(f'No session provided. Unable to parse components of {child_node["id"]} from biocyc')
            else:
                print(f'Parsing components of {child_node["id"]} from biocyc as this node was not found in the graph')
                objects,edges,nodes,type=biocyc_query_utils.find_components_recursively(session,child_node['id'],objects,nodes,edges)
        else:
            print(f"Adding parent node {parent_node} to the graph dict")
            nodes.append(parent_node)
    return {'objects':objects,'nodes':nodes,'edges':edges}
    

def build_graph_tables(session,bigg_rxns,bigg2biocyc_map,db='ECOLI',out_path='./',cache=None,export_nodes_and_edges=False,export_biocyc_objects=True):
    graph={'objects':{},'nodes':[],'edges':[]}
    graph=add_rxns_to_graph(session,bigg_rxns,bigg2biocyc_map,graph,db=db,cache=cache)
    if export_nodes_and_edges:
        nodes_df=pd.DataFrame([pd.Series(node) for node in graph['nodes']])
        edges_df=pd.DataFrame([pd.Series(edge) for edge in graph['edges']])
        
        if not os.path.isdir(out_path):
            os.mkdir(out_path)

        nodes_df.to_csv(out_path+'nodes.csv')
        edges_df.to_csv(out_path+'edges.csv')
    if export_biocyc_objects:
        if not os.path.isdir(out_path):
            os.mkdir(out_path)
        with open(out_path+'all_objects.pkl', 'wb') as handle:
            pickle.dump(graph['objects'], handle)
    return graph


def build_graph(nodes_df,edges_df,objects=None):
    '''
    Makes a directed, weighted graph with NetworkX from a pandas dataframe of edges. The DF should have three columns, namely:
    source: source node ID
    target: target node ID
    weight: edge weight, i.e. stoichiometry of the connection (for non stoichiometric connections, e.g. reaction to enzyme, 1 is used)
    '''
    G=nx.DiGraph()
    #Now populate nodes with attributes
    nodes_attributes=list(nodes_df.columns)
    nodes_list=[(row['id'],
                 {attribute:row[attribute] for attribute in nodes_attributes}
                 ) for i,row in nodes_df.iterrows()]
    G.add_nodes_from(nodes_list)
    if objects is not None:
        for node in G.nodes():
            try:
                G.nodes[node]['object']=objects[node]['biocyc_id']
            except:
                pass
    #Now populate edges with attributes
    edge_list=[(row['source'],
                row['target'],
                {'weight':row['weight'],'type':row['type'],'notes':row['notes'],'references':row['references']}
                ) 
                for i,row in edges_df.iterrows()
                ]
    G.add_edges_from(edge_list)
    return G

def add_pp_annotation(session,graph,biocyc_objects=None,db='ECOLI'):
    '''
    Extends the graph with gene nodes and edges to the polypeptides according to a specified map. Each polypeptide should be uniquely mapped to a gene.
    '''
    # Find all polypeptides nodes
    pp_nodes=[node for node in graph.nodes() if graph.nodes[node]['subtype']=='polypeptide']
    for pp in tqdm(pp_nodes):
        if biocyc_objects is None:
            _,object_id=biocyc_query_utils.split_db_id(graph.nodes[pp]['biocyc_id'])
            try:
                pp_object=biocyc_query_utils.get_biocyc_object(session,object_id,db=db)
            except:
                pp_object=None
                print(f'Unable to get {pp} from biocyc')
        else:
            _,object_id=biocyc_query_utils.split_db_id(graph.nodes[pp]['biocyc_id'])
            if graph.nodes[pp]["biocyc_id"] in biocyc_objects.keys():
                pp_object=biocyc_objects[graph.nodes[pp]['biocyc_id']]
            else:
                try:
                    pp_object=biocyc_query_utils.get_biocyc_object(session,object_id,db=db)
                except:
                    pp_object=None
                    print(f'Unable to get {pp} from biocyc')
        if pp_object is not None:
            pp_data=pp_gene_and_annotation(session,pp_object)
            graph.nodes[pp]['annotation']=pp_data['polypeptide_annotation']
            graph.nodes[pp]['gene']=pp_data['gene']
        # graph.add_node(gene,type='gene')
        # graph.add_edge(pp,gene,{'weight':1})
    return graph



def pp_gene_and_annotation(session,biocyc_pp):
    '''
    Given a biocyc polypeptide object, returns associated genes and other relevant annotation 
    '''
    pp_id=biocyc_pp.findall('.//Protein')[0].attrib['frameid']
    out={'gene':{},'polypeptide_annotation':{}}
    genes=biocyc_pp.findall('Protein/gene/Gene')
    genes_ids=[g.attrib['frameid'] for g in genes]

    if len(genes_ids)>1:
        print(f'Warning: Multiple genes found for polypeptide {pp_id}. Selection first one ({genes_ids[0]})')
    #Retrieve gene object
    if len(genes_ids)>0:
        biocyc_gene=biocyc_query_utils.get_biocyc_object(session=session,object_id=genes_ids[0])
        out['gene']=gene_annotation(biocyc_gene=biocyc_gene)
    else:
        out['gene']=None

    annotation=biocyc_pp.findall('Protein/dblink')
    for ann in annotation:
        cur_db=ann.find('.//dblink-db').text
        cur_entry=ann.find('.//dblink-oid').text
        out['polypeptide_annotation'][cur_db]=cur_entry
    return out


def gene_annotation(biocyc_gene):
    out={}
    Gene=biocyc_gene.find('Gene')
    gene_id=Gene.attrib['frameid']
    out['id']=gene_id
    #bnum
    try:
        out['bnum']=Gene.find('.//accession-1').text
    except:
        print(f'No bnum found for gene {gene_id}')
    #common name
    try:
        out['name']=Gene.find('.//common-name').text
    except:
        print(f'No common name found for gene {gene_id}')#
    return out
    

def add_node_attributes_from_df(graph,df):
    '''
    Given a df indexed by node IDs, add the available attributes in the table to the corresponding nodes.
    The attributes will take the name of the df columns
    '''
    for id,row in df.iterrows():
        if id in graph.nodes():
            for col in df.columns:
                graph.nodes[id][col]=row[col]
    return graph







def graph_visualisation(graph,rxn=None):
    '''
    Visualise a graph using PyVis
    '''
    if rxn is not None:
        g=nx.subgraph(graph,
                   nx.node_connected_component(nx.Graph(graph),rxn)
                   )
    else:
        g=graph
    plot=Network(notebook=True,cdn_resources='remote',directed=True)
    attributes=g.nodes.data()
    for node in g.nodes():
        plot.add_node(node,label=node,group=attributes[node]['type'])
    for edge in g.edges.data():
        plot.add_edge(edge[0],edge[1],label=str(edge[2]['weight']))
    return plot

def create_subgraph_for_visualisation(graph,rxn=None,omit_putative=True):
    '''
    Visualise a the subgraph for rxn graph using Gravis
    '''
    if rxn is not None:
        g=graph.subgraph(
                   nx.node_connected_component(nx.Graph(graph),rxn)
                   ).copy()
        rxn_enzymes=list(g.successors(rxn))
    else:
        g=graph.copy()
    
    #Prune all other reaction nodes
   
    for node in list(g.nodes):
        for edge in list(g.out_edges(node)):
            if g.edges[edge]['type']=='regulation' and edge[1] not in rxn_enzymes:
                g.remove_edge(edge[0],edge[1])
            elif g.edges[edge]['type'] in ['putative_gpr','putative gpr'] and omit_putative:
                g.remove_edge(edge[0],edge[1])
    for edge in g.edges:
        if g.edges[edge]['weight']=='NA':
            g.edges[edge]['weight']=''

    g=g.subgraph(
            nx.node_connected_component(nx.Graph(g),rxn)
            ).copy()       
    #g.remove_nodes(to_remove)

    cmap={'reaction':'red',
          'protein':{'multimeric_protein':'blue','modified_protein':'yellow','polypeptide':'lightblue'},
          'gene':'green',
          'compound':'gray',
          'logical_OR':'pink'}
    
    cmap_edges={'catalysis':'red',
                'secondary_catalysis':'pink',
                'regulation':'#0BBCFF',
                
                }

    
    plot=Network(notebook=True,cdn_resources='remote',directed=True)
    attributes=g.nodes.data()
    for node in g.nodes():
        if attributes[node]['type'] in cmap.keys():
            if attributes[node]['type']=='protein':
                protein_subtype=attributes[node]['subtype']
                g.nodes[node]['color']=cmap['protein'][protein_subtype]
            else:
                g.nodes[node]['color']=cmap[attributes[node]['type']]
        else:
            g.nodes[node]['color']='black'
    for edge in g.edges:
        if g.edges[edge]['type'] in cmap_edges.keys():
            g.edges[edge]['color']=cmap_edges[g.edges[edge]['type']]
        else:
            g.edges[edge]['color']='black'
    return g

def genes_in_gpr(gpr_rule):
    s=gpr_rule.replace('(','')
    s=s.replace(')','')
    s=s.replace(' and ',',')
    s=s.replace(' or ',',')
    return s.split(',')

    
def gpr_join(gpr_list,operator,brackets=True):
    if '' in gpr_list:
        gpr_list.remove('')
    elif '()' in gpr_list:
        gpr_list.remove('()')
    
    if len(gpr_list)==1:
        return gpr_list[0]
    elif len(gpr_list)==0:
        return ''
    else:
        joint= f' {operator} '.join(gpr_list)
        if brackets:
            return f'({joint})'
        else:   
            return joint
def compute_node_gpr(graph,node,spontaneous_gpr='s0001',
                     catalysis_types=['catalysis','secondary_catalysis','spontaneous_reaction'],
                     protein_requirement_types=['non_catalytic_requirement'],
                     ):
    node_type=graph.nodes[node]['type']
    node_subtype=graph.nodes[node]['subtype']
    if node_type=='gene':
        node_gpr= graph.nodes[node]['annotation']['bnum'] if 'bnum' in graph.nodes[node]['annotation']['bnum'] else ''
    elif node_type=='reaction':
        #The GPR of a reaction is the OR of the GPRs of its children connected via a catalysis node.
        # In addition each catalytic child is ANDed with the gpr of any protein_metabolite
        catalytic_edge_types=catalysis_types
        catalytic_children=[child for child in graph.successors(node) if graph.edges[node,child]['type'] in catalytic_edge_types]
        catalytic_children_gprs=[compute_node_gpr(graph,child) for child in catalytic_children]

        protein_requirements_edge_types=protein_requirement_types
        protein_requirement_children=[child for child in graph.successors(node) if graph.edges[node,child]['type']in protein_requirements_edge_types]
        protein_requirement_children_gprs=[compute_node_gpr(graph,child) for child in protein_requirement_children]
        protein_requirements_gpr=gpr_join(protein_requirement_children_gprs,'and')
        
        reaction_gpr=gpr_join([
                                gpr_join([gpr,protein_requirements_gpr],'and') 
                               for gpr in catalytic_children_gprs
                               ],
                              'or',
                              brackets=False)
        
        #Rem
        node_gpr= reaction_gpr
    elif node_type=='logical_OR':
        children_gprs=[compute_node_gpr(graph,node) for node in graph.successors(node)]
        node_gpr= gpr_join(children_gprs,'or')
    elif node_type=='logical_AND':
        children_gprs=[compute_node_gpr(graph,node) for node in graph.successors(node)]
        node_gpr= gpr_join(children_gprs,'and')
    elif node_type=='protein' :
        if node_subtype=='polypeptide':
            children_genes=[child for child in graph.successors(node) if graph.nodes[child]['type']=='gene']
            node_gpr= gpr_join(children_genes,'or') #in practice, each polypeptide should map to a unique gene
        elif node_subtype=='modified_protein':
            edges_to_include=['protein_modification','protein_modification_requirement']
            children_gprs=[compute_node_gpr(graph,child) for child in  graph.successors(node) if graph.edges[node,child]['type'] in edges_to_include]
            node_gpr= gpr_join(children_gprs,'and')
        else:
            edges_to_include=['subunit_composition']
            children_gprs=[compute_node_gpr(graph,child) for child in  graph.successors(node) if graph.edges[node,child]['type']  in edges_to_include]
            node_gpr= gpr_join(children_gprs,'and')
    elif node_type=='spontaneous':
        node_gpr= spontaneous_gpr
    else:
        raise ValueError(f'undefined GPR parsing for node of type {node_type}')

    
    return node_gpr
    
def parse_gpr_in_graph(graph):
    '''
    Compute a GPR from the graph
    '''

    rxns= [ n for n in graph.nodes if graph.nodes[n]['type']=='reaction']

    for r in rxns:
        if r not in graph.nodes:
            print(f'Warning: {r} not found in the graph. Skipping')
            continue
        
def bnum2pp(graph,gene):
    '''
    Find the polypeptide associated with a gene
    '''
    for node in graph.nodes:
        if graph.nodes[node]['subtype']=='polypeptide':
            if graph.nodes[node]['gene']['bnum']==gene:
                return graph.nodes[node]['annotation']['UNIPROT'] if 'UNIPROT' in graph.nodes[node]['annotation'].keys() else None
    return None

def compute_node_mw(graph,node,pp_mw_map):
    '''
    Compute the MW of each protein or complex in the graph recursively, based on known MW of the polypeptides in the graph.
    (All other proteins or complexes can be reduced to stoichiometric compositions of polypeptides)
    '''
    assert node['type']!='reaction'
    #Now recursively compute the MW of the node
    if node['subtype']=='polypeptide':
        mw=pp_mw_map[node['id']]
    else: 
        mw=np.sum([compute_node_mw(graph,graph.nodes[child],pp_mw_map)*graph.edges[node['id'],child]['weight'] for child in graph.successors(node['id'])if graph.edges[node['id'],child]['type']=='subunit_composition'])

    return mw
def compute_graph_mw(graph,pp_mw_map):
    '''
    Compute the MW of each protein or complex in the graph recursively, based on known MW of the polypeptides in the graph.
    (All other proteins or complexes can be reduced to stoichiometric compositions of polypeptides)
    '''
    #Find all nodes that are not reactions
    nodes=[node for node in graph.nodes() if graph.nodes[node]['type']!='reaction']
    for i,node in enumerate(nodes):
        #print(node)
        graph.nodes[node]['mw']=compute_node_mw(graph,graph.nodes[node],pp_mw_map)
    return graph

def bnum2uniprot(graph,gene):
    '''
    Find the polypeptide associated with a gene
    '''
    matched_nodes=[node for node in graph.nodes if graph.nodes[node]['type']=='gene' and graph.nodes[node]['annotation']['bnum']==gene]
    if len(matched_nodes)==0:
        print("No node in the graoh could be matched to the provided b-number")
        gene_node=None
    elif len(matched_nodes)>1:
        print(f"Warning: Multiple gene nodes mapped to this bnum {matched_nodes}. Using the first one ")
        gene_node=matched_nodes[0]
    else:
        gene_node=matched_nodes[0]
    if gene_node is None:
        return None
    else:
        associated_pps=[node for node in graph.predecessors(gene_node) if graph.edges[node,gene_node]['type']=='coding_relation']
        if len(associated_pps)==0:
            print("Warning. Unable to find any polypeptide in the graph associated with provided gene")
            return None
        else:
            try:
                return graph.nodes[associated_pps[0]]['annotation']['UNIPROT']
            except:
                print(f"Unable to retrive polypeptie=de UNIPROT annotation for {associated_pps[0]}")
                return None

    

def uniprot_to_pp(graph,uniprot_id):
    '''
    Find the polypeptide associated with a UNIUPROT ID, if exits
    '''
    for node in graph.nodes:
        if graph.nodes[node]['subtype']=='polypeptide':
            if 'UNIPROT' in graph.nodes[node]['annotation'].keys():
                if graph.nodes[node]['annotation']['UNIPROT']==uniprot_id:
                    return node
    return None
def weighted_number_of_children(graph,node):
    '''
    Compute the number of children of a node, weighted by the stoichiometry of the connection
    '''
    return np.sum([graph.edges[node,child]['weight'] for child in graph.successors(node)])
def number_of_pp_subunits(graph,node):
    '''
    recursively Compute the number of polypeptide subunits of a node, weighted by the stoichiometry of the connection
    '''
    if graph.nodes[node]['subtype']=='polypeptide':
        return 1
    else:
        return np.sum([number_of_pp_subunits(graph,child)*graph.edges[node,child]['weight'] for child in graph.successors(node)])

def compute_node_pp_composition(graph,node):
    """
    Given a protein node, compute its stoichiometric composition in terms of polypeptides.
    Returns a dictionary with the polypeptide IDs as keys and the stoichiometry of each polypeptide as values
    """
    if graph.nodes[node]['subtype']=='polypeptide':
        return {node:1}
    else:
        #Recursively search all path to polypeptides
        pp_composition={}
        for child in graph.successors(node):
            child_composition=compute_node_pp_composition(graph,child)
            for pp in child_composition.keys():
                if pp in pp_composition.keys():
                    pp_composition[pp]+=child_composition[pp]*graph.edges[node,child]['weight']
                else:
                    pp_composition[pp]=child_composition[pp]*graph.edges[node,child]['weight']
    return pp_composition
def compute_catalytic_nodes(graph,catalysis_types=['catalysis','secondary_catalysis']):
    '''
    Compute all enzyme nodes in the graph, i.e. all nodes directly connected to a reaction node via a catalysis or secondary catalysis role
    '''
    enzyme_nodes=[]
    for node in graph.nodes:
        if graph.nodes[node]['type']=='reaction':
            for child in graph.successors(node):
                if graph.edges[node,child]['type'] in catalysis_types:
                    enzyme_nodes.append(child)
    return list(set(enzyme_nodes))