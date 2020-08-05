rerun = '/home/ubuntu//transfer/Incoming/QIB_Sequencing/Nextseq_1_runs/200507_NB501819_0134_AH5773AFX2/deplexed-zero-mismatch/result.ivar.pipeline/ncovIllumina_sequenceAnalysis_callVariants'
oldrun = '/home/ubuntu//transfer/Incoming/QIB_Sequencing/Covid-19_Seq/result.illumina.20200429/ncovIllumina_sequenceAnalysis_callVariants'

oldqc ='/home/ubuntu//transfer/Incoming/QIB_Sequencing/Covid-19_Seq/result.illumina.20200429/NORW-20200429.qc.csv'
newqc = '/home/ubuntu//transfer/Incoming/QIB_Sequencing/Nextseq_1_runs/200507_NB501819_0134_AH5773AFX2/deplexed-zero-mismatch/result.ivar.pipeline/test-blank.qc.csv'

import os 
import subprocess
import csv

def format_csv(path_file):
    new_dict = {}
    for x in csv.DictReader(open(path_file), dialect=csv.excel_tab):
        new_dict[x['POS']] = x
    return new_dict

snps = {} 

qcdict = csv.DictReader(open(newqc))
oldqcdict = [x for x in csv.DictReader(open(oldqc))]
print('Sample_name\tSIGMA_pct_covered_bases\tIDT_pct_covered_bases\tdiff')
for x in qcdict:
    old_value = [y for y in oldqcdict if y['sample_name'].split('_')[0] == x['sample_name'].split('_')[0]]
    if old_value:
        real_old_value = old_value[0]['pct_covered_bases']
        print(f"{x['sample_name']}\t{x['pct_covered_bases']}\t{real_old_value}\t{round(float(x['pct_covered_bases']) -  float(real_old_value)) }")

for x in os.listdir(rerun):
    if x.endswith('tsv'):
        newrunpath = os.path.join(rerun , x)
        old_file_name = [y for y in os.listdir(oldrun) if y.startswith(os.path.basename(x).split('_')[0])]
        if old_file_name:
            print('\n' + old_file_name[0])
            if not snps.get(old_file_name[0]):
                snps[old_file_name[0]] = 0
            oldrunpath = os.path.join(oldrun, old_file_name[0])
            new_run = format_csv(newrunpath)
            old_run = format_csv(oldrunpath)
            for it, val in old_run.items():
                if not new_run.get(it):
                    snps[old_file_name[0]] +=1 

                    print(val)
for x, y  in snps.items():
    print(f"NORW-{x.split('_')[0]}\t{y}")