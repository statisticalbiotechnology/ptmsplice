'''
Giving a percolator .tab file, and a qvalue threshold, get the min list of 
proteins covering all the peptides identified at that threshold. A random 
protein is chosen whenever ties.  
'''

import data_handling.percolator as percolator 
import sys 
import random
from pygraph.classes.graph import graph
from pygraph.algorithms.accessibility import connected_components
import itertools

def get_dict(peptides, qval_threshold):
    peptide_dict = {}
    protein_dict = {}
    for p in peptides:
        if p.is_decoy==False and p.qvalue<qval_threshold:
            target_proteins = [prot for prot in p.proteins if prot.find("random") == -1]        
            peptide_dict[p.sequence] = target_proteins
            for prot in target_proteins:
                if prot in protein_dict:
                    protein_dict[prot] = protein_dict[prot].union(set([p.sequence]))
                else:
                    protein_dict[prot] = set([p.sequence])
            
    
    return peptide_dict, protein_dict
            

def print_message(msg, verbosity, verbosity_threshold):
    if verbosity > verbosity_threshold:
        print msg
    
    
def get_parsimony_list(peptides, qval_threshold, verbosity=2):
    '''
    The input is a list of PercolatorPeptide objects
    '''
    # create dictionaries
    peptide_dict, protein_dict = get_dict(peptides, qval_threshold)        
    print_message("\n{} peptides, mapping to {} proteins are used".format(\
            len(peptide_dict), len(protein_dict)), verbosity, 2)
    
    parsimony_list = []      
    explained = set([])    
    
    # include the proteins that have at least one unique peptide
    print_message("\nFind the proteins with unique peptides...", verbosity, 2)
    for pr, pe in protein_dict.iteritems():
        has_unique_pe = any(len(peptide_dict[p])==1 for p in pe) 
        if has_unique_pe:
            parsimony_list.append(pr)
            explained = explained.union(pe)
    print_message("{} such proteins found, explaining {} peptides".format(\
        len(parsimony_list), len(explained)), verbosity, 2)
    
    # update the peptide and protein dictionary
    tmp = set(peptide_dict.keys())
    unexp_peptides = tmp.difference(explained) 
    for pr in protein_dict.keys():
        protein_dict[pr] = protein_dict[pr].intersection(unexp_peptides)
        
    if len(unexp_peptides) == 0:
        print_message("\n------\n{} proteins in the final list\n".format(\
                len(parsimony_list)), verbosity, 2)
        return parsimony_list
        
    # build a graph
    print_message("\nBuild a graph...", verbosity, 2)
    gr = graph()
    for pe in unexp_peptides: 
        proteins = peptide_dict[pe]
        for my_protein in proteins:
            if not gr.has_node(my_protein):
                gr.add_node(my_protein)
        # add edges betwen all these proteins 
        for i in range(len(proteins) - 1):
            for j in range(i + 1, len(proteins)):
                if not gr.has_edge((proteins[i], proteins[j])):
                    gr.add_edge((proteins[i], proteins[j]))
    print_message("{} unexplained peptides, {} nodes, with {} edges ".format(\
            len(unexp_peptides), len(gr.nodes()), len(gr.edges())), \
            verbosity, 2)
    # get the connected components
    print_message("\nGet the connected components ...", verbosity, 2)
    con_components = connected_components(gr)
    subgraphs_dict = {}
    for pr, comp in con_components.iteritems():
        if comp in subgraphs_dict:
            subgraphs_dict[comp] = subgraphs_dict[comp].union(set([pr]))
        else:
            subgraphs_dict[comp] = set([pr]) 
    subgraphs = subgraphs_dict.values()
    tmp = [len(s) for s in subgraphs]
    print_message("{} subgraphs, with {} to {} proteins".format(\
            len(subgraphs_dict), min(tmp), max(tmp)), verbosity, 2)
    # get the exact solution in each subgraphs
    print_message("\nPerform exhaustive search in subgraphs", verbosity, 2)
    for sub in subgraphs:
        min_subset = get_exact_solution(sub, protein_dict)
        parsimony_list += min_subset
    print_message("\n------\n{} proteins in the final list\n".format(\
        len(parsimony_list)), verbosity, 2)
    
    return parsimony_list
  

def get_exact_solution(subgraph, protein_dict):
    '''
    given a subgraph of connected proteins, get the minimum subset that 
    explains all peptides; the first one is chosen when more than one  
    '''
    to_cover = get_peptides(subgraph, protein_dict)
    min_subsets = []
    found = False
    for k in range(1, len(subgraph) + 1):
        for sub in itertools.combinations(subgraph, k):
            covered = get_peptides(sub, protein_dict)
            if len(covered) == len(to_cover):                
                min_subsets.append(sub)
                found = True
                break
        if found:
            break
    return min_subsets[0]
    
"""    
def get_exact_solution(subgraph, protein_dict):
    '''
    given a subgraph of connected proteins, get the minimum subset that 
    explains all peptides; 
    '''
    to_cover = get_peptides(subgraph, protein_dict)
    min_subsets = []
    found = False
    for k in range(1, len(subgraph) + 1):
        # get all the protein subsets of length k 
        for sub in itertools.combinations(subgraph, k):
            covered = get_peptides(sub, protein_dict)
            if len(covered) == len(to_cover):
                found = True
                min_subsets.append(sub)
        # we stop here since all the other sets are larger than the current 
        if found:
            break 
    # if more than one minimal subset, choose a random one 
    idx = 0 
    if len(min_subsets) > 1:
        idx = random.randint(0,  len(min_subsets) - 1)
    return min_subsets[idx]
"""

def get_peptides(protein_set, protein_dict):
  peptides = itertools.chain(*[protein_dict[p] for p in protein_set])
  return set(peptides)


def main():
    perc_tab_file = "/scratch/lumi_work/projects/various/oliver/data/yeast/103111-Yeast-2hr-01.percolator.tab"
    qval_threshold = 0.01    
    get_parsimony_list(perc_tab_file, qval_threshold, verbosity=3)
    
    
if __name__ == '__main__':
    main()
    
