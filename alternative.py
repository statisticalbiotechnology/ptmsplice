"""
Parse an Ensembl proteome FASTA file and return which genes have 
multiple splice variants. Count the number of splice variants. Map 
phosphosites from tissue data to a parsed proteome.

per.warholm@scilifelab.se
2014-07-30
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
                # phosphosites[tissue_name][protein_name] = [float(i)/n for i in hits]
                phosphosites[tissue_name][protein_name] = hits, n

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

def probabilities(proteome, translation_table, phosphosites, tissue_specificity):
    """Calculate probabilities according to project notes."""
    
    ## First probability:
    ## Pr(A | "Gene is alt.spliced") > Pr(A | "Gene is not alt.spliced")
    ## A = "Gene has a phosphosite"
    
    ## count alternative slicing
    n_genes_AS = 0
    n_genes_nonAS = 0
    for k, v in proteome.iteritems():
        if len(v) > 1:
            n_genes_AS+=1
        elif len(v) == 1:
            n_genes_nonAS+=1
    
    ## conditional count of phosphosites
    count_phosphogenes = []
    A_given_AS = []
    A_given_nonAS = []
    for tissue, data in phosphosites.iteritems():
        for protein, phosphosite_expression in data.iteritems():
            count_phosphogenes.append(protein)
            if len(proteome[translation_table[protein]]) > 1:
                A_given_AS.append(protein)
            elif len(proteome[translation_table[protein]]) == 1:
                A_given_nonAS.append(protein)
    
    ## counts
    n_genes = len(proteome)
    n_phosphogenes = len(set(count_phosphogenes))
    n_A_given_AS = len(set(A_given_AS))
    n_A_given_nonAS = len(set(A_given_nonAS))
    
    ## probabilities
    P_phosposite = float(n_phosphogenes)/n_genes
    P_A_given_AS = float(n_A_given_AS)/n_genes_AS
    P_A_given_nonAS = float(n_A_given_nonAS)/n_genes_nonAS
    
    ## lambda pseudo function for rounding off results
    rnd = lambda x: '{0:.3f}'.format(x)
    
    print 'First probability:'
    print 'Pr(A | "Gene is alt.spliced")', rnd(P_A_given_AS)
    print 'Pr(A | "Gene is not alt.spliced")', rnd(P_A_given_nonAS)
    print 'Ratio:', rnd(P_A_given_AS/P_A_given_nonAS)
    
    
    ## Second Probability:
    ## Pr(B | "Gene is alt.spliced") > Pr(B | "Gene is not alt.spliced")
    ## B: The genes have phosphosites that are unique for one of the tissues
    
    ## filter away proteins that aren't tissue specific
    n_B_given_AS = 0
    for protein in set(A_given_AS):
        if tissue_specificity[protein]:
            n_B_given_AS+=1
    n_B_given_nonAS = 0
    for protein in set(A_given_nonAS):
        if tissue_specificity[protein]:
            n_B_given_nonAS+=1
    
    ## probabilities
    P_B_given_AS = float(n_B_given_AS)/n_genes_AS
    P_B_given_nonAS = float(n_B_given_nonAS)/n_genes_nonAS
    
    print '\nSecond probability:'
    print 'Pr(B | "Gene is alt.spliced")', rnd(P_B_given_AS)
    print 'Pr(B | "Gene is not alt.spliced")', rnd(P_B_given_nonAS)
    print 'Ratio:', rnd(P_B_given_AS/P_B_given_nonAS)

## locations of input data
fn_fasta = '../../data/mouse/Mus_musculus.GRCm38.75.pep.all.fa'
# fn_tissues = ['../data/mouse/tissues/brain-identified-proteins.txt']
fn_tissues = sorted(glob.glob('../../data/mouse/tissues/*-identified-proteins.txt'))

## main pipeline
mouse_proteome, translation_table = parse_proteome(fn_fasta)
# count_splice_variants(mouse_proteome)
n = count_splice_variants_average(mouse_proteome)
phosphosites = map_phosposites(fn_tissues, mouse_proteome, translation_table)
tissue_specificity = determine_tissue_specificity(phosphosites)
probabilities(mouse_proteome, translation_table, phosphosites, tissue_specificity)
