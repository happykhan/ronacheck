import csv
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import os 
import numpy as np
import random

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


if os.path.exists('qc_database'):
    all_qc_info = {x['sample_name']: x for x in csv.DictReader(open('qc_database'))}
    all_dates_sequenced = list(set([x['date_sequenced']for x in all_qc_info.values()]))
    for date_sequenced in sorted(all_dates_sequenced):
        high_ct  = {x['sample_name'] : int(x['mapped_read_50']) for x in all_qc_info.values() if x['date_sequenced'] == date_sequenced and float(x['max_ct']) >= 35 and float(x['max_ct']) < 40}
        blanks = {x['sample_name'] : int(x['mapped_read_50']) for x in all_qc_info.values() if x['date_sequenced'] == date_sequenced and float(x['max_ct']) == 40}
        if len(high_ct) > 0 and len(blanks) > 0 :
            top_blank_name = [k for k,v in sorted(blanks.items(),key=lambda item: item[1],reverse=True)][0]
            # random_sample = random.choice(list(high_ct.keys()))
            top_sample_name = [k for k,v in sorted(high_ct.items(),key=lambda item: item[1],reverse=True)][0]
            blank_depth = read_depth_file(all_qc_info[top_blank_name]['bam'])
            sample_depth = read_depth_file(all_qc_info[top_sample_name]['bam'])

            blank_depth_x = np.array([ int(x[2]) for x in blank_depth  ])
            sample_depth_x = np.array([int(x[2]) for x in sample_depth ])
            plt.plot(sample_depth_x, label='Sample')
            plt.plot(blank_depth_x, label='Blank')
            plt.xlabel("Position")
            plt.ylabel("No. of mapped reads (>50bp)")
            plt.title(f"Run {date_sequenced}: {top_blank_name} vs {top_sample_name}")
            plt.legend()
            plt.savefig("blankvssample." + str(date_sequenced) + '.png')
            plt.close()
        # fig, ax = plt.subplots()
        # ax.bar(blank_depth)
        # ax.bar(sample_depth)
        # ax.set_xlabel('Run')
        # ax.set_ylabel('No. of mapped reads (>50bp)')
        # plt.savefig(all_dates_sequenced + '.blankvssample.png')
