#!/usr/bin/env python3
"""
truecrown is a tool for testing low yield samples for real signs of SARScov2

Script requires BAM file mapped to SARScov2 reference, and the blank from the same run. 
Script will try to detect if this sample is a true presence of SARSCOV2. It will try to detect errors 
from a contaminated blank, and primer dimer. 

### CHANGE LOG ### 
2020-08-05 Nabil-Fareed Alikhan <nabil@happykhan.com>
    * Initial build
"""
import ronalib
import meta
import sys
import time
import logging
import argparse
import pysam 
import os
import subprocess
import numpy as np 
import matplotlib.pyplot as plt


epi = "Licence: "+meta.__licence__ +  " by " +meta.__author__ + " <" +meta.__author_email__ + ">"
logging.basicConfig()
log = logging.getLogger()

# Fetch reads that full map
def get_full_mapped_reads_depth(bam_file, read_len=150):
    # Produces coverage pileup from fully mapped reads (only)
    p = subprocess.Popen("samtools view -F 4 -h %s | awk '{if($0 ~ /^@/ || $6 ~ %dM) {print $0}}'  | samtools view -Sb   | samtools depth -a -d 0  -" %(bam_file, read_len),
                       stdout=subprocess.PIPE, shell=True)
    out, err = p.communicate()
    if err:
        log.error(err)
    pos_depth = []
    for ln in out.decode('utf-8').split("\n"):
       if ln:
          pos_depth.append(ln.split("\t"))
        # Produces coverage pileup from fully mapped reads (only)
    return pos_depth
    
# Fetch reads that partially match, fetch clip region and re-align. 
def get_part_match():
    pass

# Draw some plots 
def plot_full_mapped_depth(sample_depth, blank_depth, prefix, sample_name):
    if blank_depth:
        blank_depth_x = np.array([ int(x[2]) for x in blank_depth  ])
    else:
        blank_depth_x = np.zeros(len(sample_depth))
    sample_depth_x = np.array([int(x[2]) for x in sample_depth ])
    plt.plot(sample_depth_x, label='Sample')
    plt.plot(blank_depth_x, label='Blank')
    plt.xlabel("Position")
    plt.ylabel("No. of mapped reads (>50bp)")
    plt.title(f"Fully mapped read coverage: {sample_name} vs Blank")
    plt.legend()
    plt.savefig(prefix + ".png")
    plt.close()




# CLI stuff to launch other tools 


def main(args):
    if os.path.exists(args.samplebam) and os.path.exists(args.blankbam):
        full_depth = get_full_mapped_reads_depth(args.samplebam)
        sample_name = os.path.basename(args.samplebam).replace('.sorted.bam', '')
        full_depth_count = 0 
        for base in full_depth:
            if int(base[2]) > args.coverage:
                full_depth_count +=1
        if full_depth_count < args.readcount:
            log.info('FALSE POSITIVE: Not enough fully mapped reads')
        log.info(f'Genome coverage: {round(full_depth_count / len(full_depth) * 100 , 2)}%')
        blank_depth = get_full_mapped_reads_depth(args.blankbam)
        plot_full_mapped_depth(full_depth, blank_depth, args.output, sample_name)
    else: 
        log.error(f'BAM file does not exist: {args.samplebam}')


if __name__ == '__main__':
    start_time = time.time()
    log.setLevel(logging.INFO)
    desc = __doc__.split('\n\n')[1].strip()
    parser = argparse.ArgumentParser(description=desc,epilog=epi)
    parser.add_argument ('-v', '--verbose', action='store_true', default=False, help='verbose output')
    parser.add_argument('--version', action='version', version='%(prog)s ' + meta.__version__)
    parser.add_argument('-o','--output',action='store',help='output prefix', default='truecrown')
    parser.add_argument('-c','--coverage',action='store',help='Minimum coverage for fully mapped reads', default=1)
    parser.add_argument('-n','--readcount',action='store',help='Minimum number of fully mapped reads', default=10)

    parser.add_argument('samplebam', action='store', help='SARSCOV2 BAM file')
    parser.add_argument('blankbam', action='store', help='Negative control BAM file')
    args = parser.parse_args()
    if args.verbose: 
        log.setLevel(logging.DEBUG)
        log.debug( "Executing @ %s\n"  %time.asctime())
    main(args)
    if args.verbose: 
        log.debug("Ended @ %s\n"  %time.asctime())
        log.debug('total time in minutes: %d\n' %((time.time() - start_time) / 60.0))
    sys.exit(0)