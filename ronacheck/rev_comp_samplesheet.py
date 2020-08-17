
alt_map = {'ins': '0'}
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

def reverse_complement(seq):
    for k, v in alt_map.items():
        seq = seq.replace(k, v)
    bases = list(seq)
    bases = reversed([complement.get(base, base) for base in bases])
    bases = ''.join(bases)
    for k, v in alt_map.items():
        bases = bases.replace(v, k)
    return bases

dir = '/home/ubuntu/transfer/incoming/QIB_Sequencing/Nextseq_1_runs/200612_NB501819_0141_AH5WJ2AFX2'

import os 
from shutil import copy
samplesheet = os.path.join(dir, 'SampleSheet.csv')
old_samplesheet = os.path.join(dir, 'old.SampleSheet.csv')
new_samplesheet = os.path.join(dir, 'new.SampleSheet.csv')



copy(samplesheet, old_samplesheet)

out = open(new_samplesheet, 'w')

with open(old_samplesheet) as f:
    for x in f.readlines():
        if x.startswith('PID'):
            lineArray = x.split(',')
            lineArray[5] = reverse_complement(lineArray[5])
            print(','.join(lineArray).strip())
            out.write(','.join(lineArray).strip() + '\n')
        else:
            out.write(x.strip()+ '\n')
            print(x.strip())
