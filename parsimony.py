#!/usr/bin/env python

# parsimony.py - takes percolator xml output file and outputs parsimonous
#                protein identifications with peptide matches
#
# Author: Yrin Eldfjell
# Version: 20140318.1

from lxml import etree
from pygraph.classes.graph import graph
from pygraph.algorithms.accessibility import connected_components
from itertools import combinations, chain
import sys

# Outputs the identified proteins identified using parsimony. If a pair of 
# equally sized sets of proteins both cover all peptides in their group, the 
# set that matched first is returned.
#
# This program assumes that the input percolator results XML validates.
#
# Tab-delimited output format: 
# <protein name>\t<peptide 1>\t...\t<peptide N>

def parse_percolator_xml(filename, q_value_cutoff=0.01):
    """
    parse_percolator_xml - takes a file name and q value cutoff and returns 
        a dictionary keyed on identified peptides (give the specified cutoff) 
        with a list of potential protein matches as values.
    """
    NS = '{http://per-colator.com/percolator_out/14}'
    with open(filename, 'r') as f:
        perc_output = etree.parse(f).getroot()
    peptides = {}
    for pep in perc_output.iter('{}peptide'.format(NS)):
        pep_seq = pep.attrib.get('{}peptide_id'.format(NS))
        q_value_str = pep.find('{}q_value'.format(NS)).text
        try:
            q_value = float(q_value_str)
        except ValueError:
            raise ValueError("Invalid q value {}".format(q_value_str))
        if pep_seq in peptides:
            raise RuntimeError("Duplicate peptide in percolator xml output")
        if q_value > q_value_cutoff:
            continue
        protein_ids = [p.text for p in pep.iter('{}protein_id'.format(NS))]
        peptides[pep_seq] = protein_ids
    return peptides

def parsimonous_protein_identification(peptides):
    """
    parsimonous_protein_identification - takes a dict of the form
        {<peptide_seq>: <protein_name>, [<protein_na,e> ...] } and returns the 
        proteins identified using parsimony.
    """
    detected_proteins = {}
    protein2peptides = {}

    # start with the uniquely determined proteins
    for peptide, proteins in peptides.items():
        if len(proteins) == 1:
            detected_proteins[proteins[0]] = [peptide]
            peptides.pop(peptide)
        else:
            for p in proteins:
                if not p in protein2peptides:
                    protein2peptides[p] = [peptide]
                else:
                    protein2peptides[p].append(peptide)

    # remaining peptides have multiple potential proteins, use parsimony
    g = graph()

    # identify protein clusters
    for peptide, proteins in peptides.items():
        for protein in proteins:
            if not g.has_node(protein):
                g.add_node(protein)
        for p1, p2 in combinations(proteins, 2):
            if not g.has_edge((p1, p2)):
                g.add_edge((p1, p2))
    connected = connected_components(g).items()
    clusters = {subgraph: set() for protein, subgraph in connected}
    for protein, subgraph in connected:
        clusters[subgraph] = clusters[subgraph].union(set((protein,)))

    def find_covering(proteins):
        peptides = set(chain(*(tuple(protein2peptides[p]) for p in proteins)))
        for k in range(1, len(proteins) + 1):
            for covering in combinations(proteins, k):
                covered = set(chain(*(tuple(protein2peptides[p]) for p in 
                    covering)))
                if len(covered) == len(peptides):
                    return covering
        return None

    # find the minimal protein covering of each cluster
    for cluster in clusters.values():
        covering = find_covering(cluster)
        if covering is None:
            print "Error, failed to cover " + str(subgraph)
            sys.exit(1)
        else:
            for protein in covering:
                detected_proteins[protein] = protein2peptides[protein]

    return detected_proteins

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Usage: parsimony.py <name of percolator xml output file>"
        sys.exit(1)
    per_xml_filename = sys.argv[1]
    peptides = parse_percolator_xml(per_xml_filename)
    proteins = parsimonous_protein_identification(peptides)
    for protein, peptides in proteins.items():
        print "{}\t{}".format(protein, "\t".join(peptides))
