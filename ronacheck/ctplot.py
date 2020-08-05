#!/usr/bin/env python3
from oauth2client.service_account import ServiceAccountCredentials
import gspread
import collections
import logging
from marshmallow import Schema, fields, EXCLUDE, pre_load, validate
import csv 
import datetime
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO
import csv
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import pysam
from collections import Counter
import os 
import numpy as np

"""
This script can incorporate as many QC checks as required
as long as it outputs a csv file containing a final column
headed with 'qc_pass' and rows for each sample indcating
'TRUE' if the overall QC check has passed or 'FALSE' if not.
"""

def make_qc_plot(depth_pos, n_density, samplename, window=200):
    depth_df = pd.DataFrame( { 'position' : [pos[1] for pos in depth_pos], 'depth' : [dep[2] for dep in depth_pos] } )
    depth_df['depth_moving_average'] = depth_df.iloc[:,1].rolling(window=window).mean()

    n_df = pd.DataFrame( { 'position' : [pos[0] for pos in n_density], 'n_density' : [dens[1] for dens in n_density] } )

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    ax1.set_xlabel('Position')

    ax1.set_ylabel('Depth', color = 'g')
    ax1.set_ylim(top=10**5, bottom=1)
    ax1.set_yscale('log')
    ax1.plot(depth_df['depth_moving_average'], color = 'g')

    ax2.set_ylabel('N density', color = 'r')  
    ax2.plot(n_df['n_density'], color = 'r')
    ax2.set_ylim(top=1)

    plt.title(samplename)
    plt.savefig(samplename + '.depth.png')

def read_depth_file(bamfile):
    p = subprocess.Popen(['samtools', 'depth', '-a', '-d', '0', bamfile],
                       stdout=subprocess.PIPE)
    out, err = p.communicate()
    counter = 0

    pos_depth = []
    for ln in out.decode('utf-8').split("\n"):
       if ln:
          pos_depth.append(ln.split("\t"))
    
    return pos_depth


def count_full_mapped_reads(bam_file, read_length = 148):
    bam = pysam.AlignmentFile(bam_file, "rb")
    full_reads = []
    total_reads = 0 
    any_mapped_reads  = 0 
    total_bases = 0 
    mapped_read_50 = 0 
    for read in bam:
        total_reads += 1 
        if len(read.cigar) == 1:
            if (read.cigar[0][0] == 0 and read.cigar[0][1] >= read_length):
                full_reads.append(read.cigar[0][1])
        if not read.is_unmapped:
            any_mapped_reads += 1
            total_bases += len(read.seq)
            if read.alen > 50:
                mapped_read_50 += 1 
    count_ = dict(Counter(full_reads))
    sum_mapped_reads = sum(count_.values())
    
    return sum_mapped_reads, total_reads, any_mapped_reads, total_bases, mapped_read_50
    

def get_covered_pos(pos_depth, min_depth):
    counter = 0
    for contig, pos, depth in pos_depth:
        if int(depth) >= min_depth:
            counter = counter + 1
    
    return counter

def get_N_positions(fasta):
    n_pos =  [i for i, letter in enumerate(fasta.seq.lower()) if letter == 'n']

    return n_pos

def get_pct_N_bases(fasta):
    
    count_N = len(get_N_positions(fasta))

    pct_N_bases = count_N / len(fasta.seq) * 100

    return pct_N_bases

def get_largest_N_gap(fasta):
    n_pos = get_N_positions(fasta)

    n_pos = [0] + n_pos + [len(fasta.seq)]

    n_gaps = [j-i for i, j in zip(n_pos[:-1], n_pos[1:])]

    return sorted(n_gaps)[-1]

def get_ref_length(ref):
    record = SeqIO.read(ref, "fasta")
    return len(record.seq)

def get_coverage(depth_pos):
    return str(round(sum([int(x[2]) for x in depth_pos]) / len(depth_pos), 4)) + 'X'


def sliding_window_N_density(sequence, window=10):

    sliding_window_n_density = []
    for i in range(0, len(sequence.seq), 1):
        window_mid = i + ( window / 2)
        window_seq = sequence.seq[i:i+window]
        n_count = window_seq.lower().count('n')
        n_density = n_count / window

        sliding_window_n_density.append( [ window_mid, n_density ] )

    return sliding_window_n_density


import random

if os.path.exists('qc_database'):
    all_qc_info = {x['sample_name']: x for x in csv.DictReader(open('qc_database'))}
    all_dates_sequenced = list(set([x['date_sequenced']for x in all_qc_info.values()]))
    for date_sequenced in sorted(all_dates_sequenced):
        high_ct  = {x['sample_name'] : x['mapped_read_50'] for x in all_qc_info.values() if x['date_sequenced'] == date_sequenced and float(x['max_ct']) >= 35 and float(x['max_ct']) < 40}
        blanks = {x['sample_name'] : int(x['mapped_read_50']) for x in all_qc_info.values() if x['date_sequenced'] == date_sequenced and float(x['max_ct']) == 40}
        top_blank_name = [k for k,v in sorted(blanks.items(),key=lambda item: item[1],reverse=True)][0]
        if len(high_ct) > 0 and len(blanks) > 0 :
            random_sample = random.choice(list(high_ct.keys()))
            blank_depth = read_depth_file(all_qc_info[top_blank_name]['bam'])
            sample_depth = read_depth_file(all_qc_info[random_sample]['bam'])

            blank_depth_x = np.array([ int(x[2]) for x in blank_depth  ])
            sample_depth_x = np.array([int(x[2]) for x in sample_depth ])
            plt.plot(sample_depth_x, label='Sample')
            plt.plot(blank_depth_x, label='Blank')
            plt.legend()
            plt.savefig("blankvssample." + str(date_sequenced) + '.png')
        # fig, ax = plt.subplots()
        # ax.bar(blank_depth)
        # ax.bar(sample_depth)
        # ax.set_xlabel('Run')
        # ax.set_ylabel('No. of mapped reads (>50bp)')
        # plt.savefig(all_dates_sequenced + '.blankvssample.png')



scope = ['https://spreadsheets.google.com/feeds','https://www.googleapis.com/auth/drive']
creds = ServiceAccountCredentials.from_json_keyfile_name('cogsub/credentials.json', scope)
client = gspread.authorize(creds)
logging.basicConfig(level=logging.DEBUG)

sheet_name='SARCOV2-Metadata'
sheet = client.open(sheet_name).sheet1

ct = {}
all_values = sheet.get_all_records()
for x in all_values:
    ct_1 = 0 
    if len(str(x['ct_1_ct_value'])) > 1: 
        ct_1 = float(x['ct_1_ct_value'])
    ct_2 = 0 
    if len(str(x['ct_2_ct_value'])) > 1: 
        ct_2 = float(x['ct_2_ct_value'])        
    ct[x["central_sample_id"]] = max([ct_1, ct_2])


import shutil
results_dir = '/home/ubuntu/transfer/incoming/QIB_Sequencing/Covid-19_Seq'
depth = 10
ref_length = 29903
# Update all QC values
all_qc_info = {}
all_dates_sequenced = []
if os.path.exists('qc_database'):
    all_qc_info = {x['sample_name']: x for x in csv.DictReader(open('qc_database'))}
for result in sorted(os.listdir(results_dir)):
    if result.startswith('result.illumina'):
        date_sequenced = result.replace('result.illumina.', "")
        all_dates_sequenced.append(date_sequenced)

for result in sorted(os.listdir(results_dir)):
    if result.startswith('result.illumina') and False:
        read_mapping_dir = os.path.join(results_dir, result, "ncovIllumina_sequenceAnalysis_readMapping")
        read_consensus_dir = os.path.join(results_dir, result, "ncovIllumina_sequenceAnalysis_makeConsensus")
        date_sequenced = result.replace('result.illumina.', "")
        for bam_file in os.listdir(read_mapping_dir):
            if bam_file.endswith('.bam'):
                if bam_file.split('_')[0].startswith('E'): 
                    sample_name = 'NORW-' + bam_file.split('_')[0]
                else:
                    sample_name = date_sequenced + '.' + bam_file.split('.')[0]
                # Depth calcs
                if sample_name not in all_qc_info:
                    bam_file_path = os.path.join(read_mapping_dir, bam_file)
                    depth_pos = read_depth_file(bam_file_path)                
                    full_mapped_reads, total_reads, any_mapped_reads, total_bases, mapped_read_50 = count_full_mapped_reads(bam_file_path)
                    coverage = get_coverage(depth_pos)
                    depth_covered_bases = get_covered_pos(depth_pos, depth)
                    pct_covered_bases = depth_covered_bases / ref_length * 100
                    
                    all_qc_info[sample_name] = {
                            'sample_name': sample_name,
                            'date_sequenced': date_sequenced,
                            'full_mapped_reads' : full_mapped_reads,
                            'total_reads' : total_reads,
                            'coverage' : coverage,
                            'max_ct': ct.get(sample_name, 40),
                            'total_bases' : total_bases,
                            'mapped_reads' : any_mapped_reads,
                            'mapped_read_50': mapped_read_50,
                            'pct_covered_bases' : "{:.2f}".format(pct_covered_bases),
                            'bam' : bam_file_path,
                            'pct_N_bases': "0",
                            "longest_no_N_run": "0",
                            'fasta' : "",
                            'probable_false_positive': "",
                            'ambiguous_bases': "0",
                            'qc_pass': ""                            
                            }                    
                    # Unknown base calcs
                    fasta_added = False
                    for fasta_file in os.listdir(read_consensus_dir):
                        if fasta_file.endswith('fa') and  fasta_file.startswith(bam_file.split('_')[0]):
                            fasta_file_path = os.path.join(read_consensus_dir, fasta_file)
                            fasta = SeqIO.read(fasta_file_path, "fasta")
                            pct_N_bases   = 0
                            largest_N_gap = 0
                            qc_pass       = "FALSE"
                            probable_false_positive = "FALSE"
                            ambig_bases = 0
                            for base in fasta.seq:
                                if base.upper() not in ['A', 'C', 'T', 'G', 'N']:
                                    ambig_bases += 1
                            if len(fasta.seq) != 0:
                                pct_N_bases = get_pct_N_bases(fasta)
                                largest_N_gap = get_largest_N_gap(fasta)
                                # QC PASS / FAIL
                                pct_covered_bases = float(all_qc_info[sample_name]["pct_covered_bases"])
                                if largest_N_gap >= 10000 or pct_covered_bases >= 50.0:
                                    if pct_N_bases < 50.0:
                                        qc_pass = "TRUE"
                                # Flag FP
                                full_mapped_reads = float(all_qc_info[sample_name]["full_mapped_reads"])
                                if pct_N_bases >= 98.0 and full_mapped_reads <= 100:
                                    probable_false_positive = "TRUE"
                            all_qc_info[sample_name]['pct_N_bases'] = "{:.2f}".format(pct_N_bases)
                            all_qc_info[sample_name]['longest_no_N_run'] = largest_N_gap
                            all_qc_info[sample_name]['fasta'] = fasta_file_path
                            all_qc_info[sample_name]['probable_false_positive'] = probable_false_positive
                            all_qc_info[sample_name]['ambiguous_bases'] = ambig_bases
                            all_qc_info[sample_name]['qc_pass'] = qc_pass
                            fasta_added = True
                    if not fasta_added: 
                        all_qc_info.pop(sample_name)
                        print('removed ' + sample_name)
                if all_qc_info.get(sample_name):                        
                    shutil.copy('qc_database', 'old_qc_database')                            
                    with open('qc_database', 'w') as db: 
                        backup = csv.DictWriter(db, fieldnames=list(all_qc_info[sample_name].keys()) + ["max_ct"]) 
                        backup.writeheader()
                        backup.writerows(all_qc_info.values())
# a) CT 35/36/37 samples have more reads than their respective blanks
# print(f'sample_name\trun_name\tct')
# for date_sequenced in all_dates_sequenced:
#     mapped_blanks = [int(x["mapped_read_50"]) for x in all_qc_info.values() if x['date_sequenced'] == date_sequenced and str(x["max_ct"]) == "40"]
#     ave_mapped_blanks = 0
#     if len(mapped_blanks) > 0:
#         ave_mapped_blanks = round(sum(mapped_blanks) / len(mapped_blanks), 2)
#     less_than_blank = {x['sample_name'] : x['max_ct'] for x in all_qc_info.values() if x['date_sequenced'] == date_sequenced and float(x['mapped_read_50']) >= ave_mapped_blanks and float(x['max_ct']) != 40}
#     #print(f'Run: {date_sequenced}')
#     for x,y  in less_than_blank.items():
#         if float(y) > 0:
#             print(f'{x}\t{date_sequenced}\t{y}')

data_plot = []

for date_sequenced in all_dates_sequenced:
    total_cov = [0] * (ref_length +1)
    mapped_blanks = [int(x["mapped_read_50"]) for x in all_qc_info.values() if x['date_sequenced'] == date_sequenced and str(x["max_ct"]) == "40"]
    ave_mapped_blanks = 0
    if len(mapped_blanks) > 0:
        ave_mapped_blanks = round(sum(mapped_blanks) / len(mapped_blanks), 2)
    high_ct  = {x['sample_name'] : x['mapped_read_50'] for x in all_qc_info.values() if x['date_sequenced'] == date_sequenced and float(x['max_ct']) >= 35 and float(x['max_ct']) < 40}
    data_this = [int(x['mapped_read_50']) - int(ave_mapped_blanks) for x in all_qc_info.values() if x['date_sequenced'] == date_sequenced and float(x['max_ct']) >= 35 and float(x['max_ct']) < 40]
    data_plot.append(data_this)
    #print(f'Run: {date_sequenced}')
    for x,y  in high_ct.items():
        if float(y) > 0:
            print(f'{x}\t{date_sequenced}\t{y}\t{ave_mapped_blanks}')
    # b) Are there any particular regions that we always detect in samples that are later than CT 35, that give us good coverage in these late CT samples
    if not  os.path.exists(date_sequenced + '.png'):
        for values in [x for x in all_qc_info.values() if x['date_sequenced'] == date_sequenced]:
            if float(values['max_ct']) > 35 :  
                depth_pos = read_depth_file(values['bam'])
                for pos in depth_pos:
                    if int(pos[2]) > 10:
                        total_cov[int(pos[1])] += 1
        with open(date_sequenced + '.csv', 'w') as zzz :
            zzz.write('pos\tcount\n')
            count = 0 
            for jj in total_cov:
                zzz.write(f'{count}\t{jj}\n')
                count += 1
        fig = plt.figure()
        plt.bar(range(len(total_cov)), total_cov)
        #plt.hist(total_cov)
        plt.savefig( date_sequenced + '.png')

fig, ax = plt.subplots()
bp = ax.violinplot(data_plot)
ax.set_xlabel('Run')
ax.set_ylabel('No. of mapped reads (>50bp)')
plt.savefig('this.png')

scope = ['https://spreadsheets.google.com/feeds','https://www.googleapis.com/auth/drive']
creds = ServiceAccountCredentials.from_json_keyfile_name('cogsub/credentials.json', scope)
client = gspread.authorize(creds)
logging.basicConfig(level=logging.DEBUG)

sheet_name='SARCOV2-Metadata'
sheet = client.open(sheet_name).sheet1

