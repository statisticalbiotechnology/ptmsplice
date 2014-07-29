"""
Parse a proteome FASTA file from the Ensembl database and return which
genes have multiple splice variants. Count number of splice variants.
Map phosphosites to parsed proteome.

per.warholm@scilifelab.se
2014-07-29
"""

from Bio import SeqIO
import numpy as np
import re, sys

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

    ## parse tissue file
    for tissue in tissues:
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
                
                ## output fraction of peptides that are found in each splice variant
                # print protein_name, [float(i)/n for i in hits]
                print protein_name, hits, n

## main pipeline
fn_fasta = '../data/mouse/Mus_musculus.GRCm38.75.pep.all.fa'
fn_tissues = ['../data/mouse/tissues/brain-identified-proteins.txt']

mouse_proteome, translation_table = parse_proteome(fn_fasta)
# count_splice_variants(mouse_proteome)
# count_splice_variants_average(mouse_proteome)
map_phosposites(fn_tissues, mouse_proteome, translation_table)
