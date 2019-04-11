# modified on 2015-11-17 count polyHis presence

# parse FASTA data set into a Pandas DataFrame

import pandas as pd; from pandas import Series, DataFrame

# parse FASTA data set into a Pandas DataFrame
def FASTA_DataFrame(file_name):
    proteins_parsed, title = {}, ''
    for line in open(file_name, 'rt'):
        line = line.rstrip()
        if line[0] == '>':
            title = line
            proteins_parsed[title] = ''
        else: 
            proteins_parsed[title] += line                                            
    return DataFrame(proteins_parsed.items(), columns=['ref','sequence'])

polyhis_proteins = dict()
# count k-mers
def polyHis_count(x):
    k = 0
    polyhis_proteins[x] = set()
    for i in proteins.index:
        seq = proteins.sequence[i]
        if (seq.count('H' * x) >= 1) & (seq.count('H' * (x + 1)) == 0 ): 
            name = proteins.ix[i, 0].split('| ')[-1].split(' [')[0]
            print name          
            polyhis_proteins[x].add(name.split('isoform')[0])
    return len(polyhis_proteins[x])


proteins = FASTA_DataFrame('human.protein.faa.new')
n = arange(5, 18)
p = array([polyHis_count(i) for i in n])
#  array([13, 16, 20, 10,  8,  7,  5,  4,  2,  1,  1,  0,  0])
bar(n, p, color='k')

# polyHis = array([sum(proteins.sequence.apply(lambda seq: seq.count('H' * i))) for i in n])
# len_polyHis = polyHis[:-1] - polyHis[1:]

