from Bio import SeqIO
import sys
import itertools
from math import sqrt
import random

TISSUES = ('brain', 'lung', 'pancreas')

protein2gene = {}
gene2protein = {}
gene_expr = {}
gene_expr_global = {}

#
# parse genome and protein/peptide expression data
#

fasta_filename = '../msgf/Mus_musculus.GRCm38.74.pep.all.fa'
with open(fasta_filename, 'r') as f_fa:
    for record in SeqIO.parse(f_fa, 'fasta'):
        attribs = dict(
            a.split(':', 1) for a in record.description.split() if ':' in a)
        if not 'gene' in attribs:
            print "ERROR: No gene found for protein ", record.id
            sys.exit(1)
        gene = attribs['gene']
        protein2gene[record.id] = gene
        if not gene in gene2protein:
            gene2protein[gene] = []
        gene2protein[gene].append(record.id)

for tissue in TISSUES:
    gene_expr[tissue] = {}
    protein_filename = '{}-identified-proteins'.format(tissue)
    with open(protein_filename, 'r') as f_prot:
        for line in f_prot:
            if not line.strip():
                continue
            line_parts = line.split()
            protein_name, peptides = line_parts[0], line_parts[1:]
            if not protein_name in protein2gene:
                print "ERROR: Protein", protein_name, "not found in fasta file!"
                sys.exit(1)
            gene = protein2gene[protein_name]
            if gene not in gene_expr[tissue]:
                gene_expr[tissue][gene] = [(protein_name, peptides)]
            else:
                gene_expr[tissue][gene].append((protein_name, peptides))
            if not gene in gene_expr_global:
                gene_expr_global[gene] = [] 
            if not protein_name in gene_expr_global[gene]:
                gene_expr_global[gene].append(protein_name)

#
# parse panther gene onthology data
#

gene_function_db = {}
with open('panther-filtered', 'r') as f:
    for line in f:
        gene, functions_cs = line.split('\t')
        if functions_cs.strip():
            functions = [fn.strip() for fn in functions_cs.split(';')]
            for function in functions:
                if not function in gene_function_db:
                    gene_function_db[function] = []
                gene_function_db[function].append(gene)

#
# perform analysis
#

def calc_avg_num_splice_forms(genes):
    iso_form_count = 0
    for gene in genes:
        iso_form_count += len(gene2protein[gene])
    return iso_form_count / float(len(genes))

def calc_variance(genes):
    avg = calc_avg_num_splice_forms(genes)
    var_sum = 0
    for gene in genes:
        var_sum += (len(gene2protein[gene]) - avg)**2
    return var_sum / (len(genes) - 1)

def intersection(set1, set2):
    return set(set1).intersection(set(set2))

def same(sets):
    for i in range(0, len(sets) - 1):
        if sets[i] != sets[i+1]:
            return False
    return True

def random_subset_averages(genes):
    genes = list(genes)
    num_boxes = 20
    boxes = {i:[] for i in range(0, num_boxes)}
    random.shuffle(genes)
    i = 0
    for gene in genes:
        boxes[i].append(gene)
        i+=1
        if i >= num_boxes:
            i = 0
    return [calc_avg_num_splice_forms(boxes[i]) for i in range(0, num_boxes)]

# identify tissue specific genes (genes expressed in exactly 1 or 2 tissues)
tissue_specific_genes = []
tissue_non_specific_genes = []
for gene in gene_expr_global.keys():
    num_found = 0
    for tissue in TISSUES:
        if gene in gene_expr[tissue]:
            num_found += 1
    if num_found < len(TISSUES):
        tissue_specific_genes.append(gene)
    else:
        tissue_non_specific_genes.append(gene)

# identify genes with tissue specific splicing 
# (different proteins in different tissues)
tissue_specific_spliced_genes = []
tissue_non_specific_spliced_genes = []
for gene in gene_expr_global.keys():
    expressed_proteins = [{ p[0] for p in gene_expr[t][gene] }
            for t in TISSUES if gene in gene_expr[t]]
    if same(expressed_proteins):
        tissue_non_specific_spliced_genes.append(gene)
    else:
        tissue_specific_spliced_genes.append(gene)

# identify genes with tissue specific phospho-sites
tissue_specific_phospho_site_genes = []
for gene in gene_expr_global:
    # algo: for each protein expressed in more than one tissue - compare peptides
    peptide_sets = {}
    for tissue in TISSUES:
        if not gene in gene_expr[tissue]:
            continue
        for prot, peps in gene_expr[tissue][gene]:
            if not prot in peptide_sets:
                peptide_sets[prot] = []
            peptide_sets[prot].append(set(peps))
    for prot, pep_sets in peptide_sets.iteritems():
        if len(pep_sets) <= 1:
            continue
        # see blog post Thu 20 Mar 2014:
        if not same(pep_sets) and all(len(s) > 2 for s in pep_sets):
            tissue_specific_phospho_site_genes.append(gene)
tissue_specific_phospho_site_genes = list(set(tissue_specific_phospho_site_genes))

# average number of iso-forms for all genes
average_num_iso_forms_all_genes = calc_avg_num_splice_forms(gene2protein.keys())
print "standard deviation for all genes", sqrt(calc_variance(gene2protein.keys()))

# average number of iso-forms for phospho-protein coding genes
average_num_iso_forms_all_phospho_genes = calc_avg_num_splice_forms(gene_expr_global.keys())
print "standard deviation for phospho genes", sqrt(calc_variance(gene_expr_global.keys()))
print "random subset averages for phospho genes\n",\
    '\n'.join(str(f) for f in random_subset_averages(gene_expr_global.keys()))

# average number of iso-forms for genes NOT coding for phospho genes
average_num_iso_forms_non_phospho_genes = calc_avg_num_splice_forms(
    set(gene2protein.keys()) - set(gene_expr_global.keys()))
print "standard deviation for non-phospho genes", sqrt(calc_variance(
    set(gene2protein.keys()) - set(gene_expr_global.keys())))

# calculate averages for randomized subsets of the non-phospho genes
print "random subset averages for non phospho genes\n", \
    '\n'.join(str(f) for f in random_subset_averages(
    set(gene2protein.keys()) - set(gene_expr_global.keys())))

# average number of iso-forms for phospho-protein genes (per tissue)
average_num_iso_forms_all_phospho_genes_per_tissue = {}
for tissue in TISSUES:
    average_num_iso_forms_all_phospho_genes_per_tissue[tissue] = calc_avg_num_splice_forms(
        gene_expr[tissue].keys())

# average number of iso-forms for tissue-specific genes (phospho-)
average_num_iso_forms_tissue_specific_phospho_genes = calc_avg_num_splice_forms(
        tissue_specific_genes)

# average number of iso-forms for non tissue-specific genes (phospho-)
average_num_iso_forms_tissue_non_specific_phospho_genes = calc_avg_num_splice_forms(
        tissue_non_specific_genes)

# avgerage number of iso-forms for genes with tissue specific phospho sites
average_num_iso_forms_tissue_specific_phospho_site_genes = calc_avg_num_splice_forms(
        tissue_specific_phospho_site_genes)

# avgerage number of iso-forms for phospho genes without tissue specific phospho sites
average_num_iso_forms_tissue_non_specific_phospho_site_genes = calc_avg_num_splice_forms(
    set(gene_expr_global.keys()) - set(tissue_specific_phospho_site_genes))

# average number of iso-forms for genes with tissue specific splicing
# (same gene, different proteins across tissues)
# (including genes expressed in 2-3 tissues)
average_num_iso_forms_tissue_specific_spliced_phospho_genes = calc_avg_num_splice_forms(
    tissue_specific_spliced_genes)

# average number of iso-forms for genes without tissue specific splicing
# (same gene, same proteins)
# (including genes expressed in 1-3 tissues)
average_num_iso_forms_tissue_non_specific_spliced_phospho_genes = calc_avg_num_splice_forms(
    tissue_non_specific_spliced_genes)

# average number of iso-forms for genes with tissue specific splicing (per tissue)
average_num_iso_forms_tissue_specific_spliced_phospho_genes_per_tissue = {}
for tissue in TISSUES:
    average_num_iso_forms_tissue_specific_spliced_phospho_genes_per_tissue[tissue] = \
            calc_avg_num_splice_forms(intersection(gene_expr[tissue].keys(), 
                tissue_specific_spliced_genes))

# average number of iso-forms for genes with non tissue specific splicing (per tissue)
average_num_iso_forms_tissue_non_specific_spliced_phospho_genes_per_tissue = {}
for tissue in TISSUES:
    average_num_iso_forms_tissue_non_specific_spliced_phospho_genes_per_tissue[tissue] = \
            calc_avg_num_splice_forms(intersection(gene_expr[tissue].keys(), 
                tissue_non_specific_spliced_genes))

# exact number of iso-forms for each of our identified genes
# => This can easily be calculated from gene_expr_global

# exact number of iso-forms for each of our identified genes (per tissue)
# => This can easily be calculated from gene_expr

def printvar(var, val):
    if type(val) == dict:
        for key in val:
            printvar("{} ({})".format(var, key), val[key])
    else:
        if type(val) == float:
            val = '{:.5f}'.format(val)
        print("{:85} {}".format(var, val))

for var in sorted(globals().keys()):
    if 'average' in var:
        printvar(var, globals()[var])
    if 'genes' in var and 'average' not in var:
        printvar(var, len(globals()[var]))
for tissue in TISSUES:
    printvar("genes expressed in {}".format(tissue), len(gene_expr[tissue]))
printvar("genes expressed in any tissue", len(gene_expr_global))
printvar("genes expressed in all tissues", len(set(gene_expr['brain']).intersection(
    set(gene_expr['lung']).intersection(set(gene_expr['pancreas'])))))

with open('splice_form_counts.m', 'w') as fh:
    for tissue in TISSUES:
        fh.write('splice_count_%s = [ ' % tissue)
        for gene in gene_expr[tissue]:
            fh.write('{}\n'.format(str(len(gene2protein[gene]))))
        fh.write(']\n ')
    fh.write('splice_count_all = [ ')
    for gene in (set(gene2protein.keys()) - set(gene_expr_global.keys())):
        fh.write('{}\n'.format(str(len(gene2protein[gene]))))
    fh.write(']\n ')
    

# TODO: generate graphs
