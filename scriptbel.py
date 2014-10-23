"""Python implementation of sciptbel.R."""

import matplotlib.pyplot as plt
from Bio import SeqIO

## Location of data files 
fn_ubi = '../data/ACETYLATION data/Ubiquitylation/proteins.txt'
fn_genome = '../data/mouse/Mus_musculus.GRCm38.74.pep.all.fa'

## Read ubiquitylation data
with open(fn_ubi, 'r') as fh_ubi:
    lines = fh_ubi.readlines()
    
## Remove newlines, the first line (header) and the last line (empty)
lines = [line.strip() for line in lines[1:-1]]

## Count number of proteins per line in ubiquitylation dataset
count_ubi = [line.count(';') for line in lines]

## Set variables needed for plot
n_max = max(count_ubi)
bs = [float(i)+0.5 for i in range(n_max)]

## Plot histogram
f = plt.figure(figsize=(12,12))
ax = plt.subplot(1,1,1)
histogram = ax.hist(count_ubi, bins=bs,color='#4D4D4D', log=False)
ax.set_xlim(0.5, n_max+0.5)
ax.set_xlabel('Number of isoforms')
ax.set_ylabel('Frequency')
plt.suptitle('Number of protein isoforms per ubiquitylated peptide')
plt.savefig('../plots/histogram_isoforms.png')

## Count the number of proteins in mouse genome
with open(fn_genome, 'r') as fh_genome:
    records = SeqIO.parse(fh_genome, 'fasta')
    n_sequences = len(list(records))

## Count number of alternatively and non-alternatively spliced proteins
nonalternative = histogram[0][0]
alternative = sum(histogram[0][1:])

## Calculating probabilities
## Protein is ubiquitylated and is alternatively spliced/number of proteins
probability1 = float(alternative)/n_sequences 
## Protein is ubiquitylated but is not alternatively spliced/number of proteins
probability2 = float(nonalternative)/n_sequences 

## Output results
print 'Highest number of isoforms is {0}'.format(n_max)
print 'Probability 1:', probability1
print 'Probability 2:', probability2
