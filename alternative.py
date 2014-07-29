"""
Parse an Ensembl proteome FASTA file and return which genes have 
multiple splice variants. Count the number of splice variants. Map 
phosphosites from tissue data to a parsed proteome.

per.warholm@scilifelab.se
2014-07-29
"""

from Bio import SeqIO
import numpy as np
import re, glob
import sys, os

class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value 

def parse_proteome(fn_fasta):
    """Parse an Ensembl proteome FASTA file, return nested dictionary."""
    
    ## create class instance for storing parsed proteome
    proteome = AutoVivification()
    translation_table = {}
    
    with open(fn_fasta, 'r') as fh_fasta:
        for record in SeqIO.parse(fh_fasta, 'fasta'):
            
            ## parse header
            header = record.description.split()
            protein_name = header[0]
            gene_name = header[3].split(':')[1]
            transcript_name = header[4].split(':')[1]
            
            ## store sequence
            proteome[gene_name][protein_name] = record.seq
            translation_table[protein_name] = gene_name

    return proteome, translation_table

def count_splice_variants(proteome):
    """Count the number of splice variants for each gene."""

    ## print gene name, # splice variants, protein names to stdout
    for k, v in proteome.iteritems():
        average.append(len(v))
        if len(v) > 1:
            print k, len(v), [l for l in v.iterkeys()]

def count_splice_variants_average(proteome):
    """Count the average number of splice variants in a proteome."""
    
    ## count the number of splice variants for each gene
    splice_variants = []
    for k, v in proteome.iteritems():
        splice_variants.append(len(v))
    
    ## return average number of splice variants
    average = np.mean(splice_variants)
    return average

def map_phosposites(tissues, proteome, translation):
    """Map detected phosphosites from tissue data to a reference
    proteome."""

    phosphosites = AutoVivification()
    ## parse tissue files
    for tissue in tissues:
        tissue_name = os.path.basename(tissue).split('-')[0]
        with open(tissue, 'r') as fh_tissue:
            for line in fh_tissue:
                if not line.strip():
                    continue
                line_parts = line.split()
                protein_name = line_parts[0]
                peptides = line_parts[1:]
                
                ## map peptides to proteome
                splice_variants = proteome[translation[protein_name]]
                n = len(splice_variants)
                hits = [0]*len(peptides)
                for idx, peptide in enumerate(peptides):
                    for splice_name, sequence in splice_variants.iteritems():                        

                        ## remove e.g. [UNIMOD:21] from peptide sequence
                        if re.sub(r'\[.*?\]', '', peptide) in sequence:
                            hits[idx]+=1
                
                ## output to stdout
                # print protein_name, hits, n
                
                ## save result in dictionary
                phosphosites[tissue_name][protein_name] = [float(i)/n for i in hits]

    return phosphosites

def determine_tissue_specificity(phosphosites):
    """Determine which phosphorylation sites are tissue specific."""
    
    tissue_specificity = {}
    for tissue, proteins in phosphosites.iteritems():
        for protein, expressed_phosphosites in proteins.iteritems():
            ## check if expressed protein is available in other tissues
            tissue_specific = True
            for query in phosphosites.iterkeys():
                if query == tissue:
                    continue
                elif protein in phosphosites[query]:
                    tissue_specific = False
            tissue_specificity[protein] = tissue_specific
    return tissue_specificity


## main pipeline
fn_fasta = '../data/mouse/Mus_musculus.GRCm38.75.pep.all.fa'
# fn_tissues = ['../data/mouse/tissues/brain-identified-proteins.txt']
fn_tissues = sorted(glob.glob('../data/mouse/tissues/*-identified-proteins.txt'))

mouse_proteome, translation_table = parse_proteome(fn_fasta)
# count_splice_variants(mouse_proteome)
# count_splice_variants_average(mouse_proteome)
phosphosites = map_phosposites(fn_tissues, mouse_proteome, translation_table)
tissue_specificity = determine_tissue_specificity(phosphosites)
