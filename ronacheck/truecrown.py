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
import sys
import meta
import time
import logging
import argparse
import pysam 
import os
import subprocess
import numpy as np 
import matplotlib.pyplot as plt
import re 
import csv

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
def get_read_stats(bam_file, read_length=150, min_read_length=30, ref_len=29903):
    bam = pysam.AlignmentFile(bam_file, "rb")
    total_reads = 0 
    any_mapped_reads  = 0 
    total_bases = 0 
    match_len = []
    full_mapped_depth = np.zeros(ref_len) 
    full_match_reads = 0
    part_mapped_depth = np.zeros(ref_len) 
    part_match_reads = 0
    primer_mapped_depth = np.zeros(ref_len)
    primer_reads = 0
    primer_cigar = []
    part_seq_path = 'temp.part_clipped.fasta'
    primer_seq_out = open('temp.primer_clipped.fasta', 'w')
    part_seq_out = open(part_seq_path, 'w')
    part_match_out = open('temp.part_match.txt', 'w')
    part_cigar = [] 
    for read in bam:
        total_reads += 1 
        # Ignore unmapped reads
        if not read.is_unmapped:
            any_mapped_reads += 1
            total_bases += len(read.seq)
            read_len = read.inferred_length
            # Fully mapped reads
            if read_len > min_read_length and read.cigarstring == f'{read_len}M':
                full_match_reads += 1
                for pos in read.positions:
                    full_mapped_depth[pos] += 1 
            elif read_len > min_read_length:
            # Partial mapped reads 
                order = re.findall("([A-Z])", read.cigarstring)
                values = re.findall("(\d+)", read.cigarstring)
                read_position = 0 
                matched_bases = 0 
                for match_type, match_len in zip(order, [int(x) for x in values]):
                    if match_type == 'M':
                        matched_bases += match_len
                    read_position  += match_len
                if matched_bases > min_read_length: 
                    # Detect short matches, indels, deletions
                    # log.debug(f'Partial match detected: {read.qname} ({read.cigarstring}) at pos {read.pos}' )
                    part_cigar.append(read.cigarstring)
                    part_match_reads += 1
                    part_match_out.write(f'{matched_bases}\n')
                    for pos in read.positions:
                        part_mapped_depth[pos] += 1                    
                    # Fetch out all the soft clipped sequences 
                    read_position = 0
                    inner_count = 1
                    for match_type, match_len in zip(order, [int(x) for x in values]):
                        if match_type == 'S':
                            soft_clip_seq = read.seq[read_position:read_position + match_len]
                            part_seq_out.write(f'>{total_reads}-{inner_count}-{read.pos}\n{soft_clip_seq}\n')    
                            inner_count += 1                           
                        read_position += match_len                
                else:
                    # 1 - Detect primer dimer. 
                    # log.debug(f'Primer dimer detected: {read.qname} ({read.cigarstring}) at pos {read.pos}' )
                    primer_reads += 1
                    primer_cigar.append(read.cigarstring)
                    for pos in read.positions:
                        primer_mapped_depth[pos] += 1
                    read_position = 0
                    inner_count = 1                        
                    for match_type, match_len in zip(order, [int(x) for x in values]):
                        if match_type == 'S':
                            soft_clip_seq = read.seq[read_position:read_position + match_len]
                            primer_seq_out.write(f'>{total_reads}-{inner_count}-{read.pos}\n{soft_clip_seq}\n')
                            inner_count += 1                              
                        read_position += match_len
            # 2 - Detect short matches, indels, deletions
                    #match_len.append(match_group[0])

    part_out = open('temp.part_cigar.txt', 'w')
    part_out.write('\n'.join(part_cigar))
    part_out.close()

    primer_out = open('temp.primer_cigar.txt', 'w')
    primer_out.write('\n'.join(primer_cigar))
    primer_out.close()

    primer_seq_out.close()   
    part_seq_out.close()
    part_match_out.close()

    sample_depths = dict(full_match_reads = full_mapped_depth, part_match_reads = part_mapped_depth,  primer_match_reads = primer_mapped_depth)
    
    # Check content of clipped reads 

    check_clipped(part_seq_path)


    return full_match_reads, part_match_reads, primer_reads,  total_reads, any_mapped_reads, sample_depths
    
    

# Draw some plots 
def plot_mapped_series(sample_depths, blank_depth, prefix, sample_name, ref_len = 29903):
    plt.figure(figsize=(8, 6))
    if blank_depth:
        blank_depth_x = np.array([ int(x[2]) for x in blank_depth  ])
    else:
        blank_depth_x = np.zeros(ref_len)
    for series_name, series in sample_depths.items():
        plt.plot(series, label=series_name)
    plt.plot(blank_depth_x, label='Blank')
    plt.xlabel("Position")
    plt.yscale('log')
    plt.ylabel("No. of mapped reads (log)")
    plt.title(f"Mapped read pileup: {sample_name} vs Blank")
    plt.legend()
    plt.savefig(prefix + '.' + sample_name +  ".png")
    plt.close()

def check_clipped(seq_file, primer_file='dat/ncov_primer.tsv'):
    primer_db_file = 'dat/ncov_primer.fasta'
    genome_db_file = 'dat/ncov_ref.fasta'
    blast_headers = ['qseqid', "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

    with open(primer_db_file, 'w') as f:
        for primer_seq in csv.DictReader(open(primer_file), dialect=csv.excel_tab):
            f.write(f">{primer_seq.get('name')}\n{primer_seq.get('seq')}\n")
    subprocess.Popen(f'makeblastdb -in {primer_db_file} -dbtype nucl', shell=True)
    primer_output = subprocess.check_output(f'blastn -db {primer_db_file} -query {seq_file} -outfmt 6', shell=True)
    primer_output_dict = [dict(zip(blast_headers, x.split('\t'))) for x in primer_output.decode("utf-8").split('\n')]
    if primer_output_dict:
        log.info('Primer detected')
    subprocess.Popen(f'makeblastdb -in {genome_db_file} -dbtype nucl', shell=True)
    genome_output = subprocess.check_output(f'blastn -db {genome_db_file} -query {seq_file} -outfmt 6', shell=True)
    blast_headers = ['qseqid', "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    genome_output_dict = [dict(zip(blast_headers, x.split('\t'))) for x in genome_output.decode("utf-8").split('\n')]
    deltas = [] 
    for res in genome_output_dict:
        if len(res['qseqid'].split('-')) > 2: 
            qseq_loc = int(res['qseqid'].split('-')[-1])
            min_delta = abs(qseq_loc - int(res['send']))
            if min_delta > abs(qseq_loc - int(res['sstart'])):
                min_delta = abs(qseq_loc - int(res['sstart']))
            deltas.append(min_delta)


# CLI stuff to launch other tools 
def main(args):
    ref_len = 29903    
    if os.path.exists(args.samplebam) and os.path.exists(args.blankbam):
        full_match_reads, part_match_reads, primer_reads,  total_reads, any_mapped_reads, sample_depths = get_read_stats(args.samplebam)
        sample_name = os.path.basename(args.samplebam).replace('.sorted.bam', '')

        


        # summary = open( args.output + '.' + sample_name  + '.csv','w')
        # genome_cov = round(full_depth_count / ref_len * 100 , 2)
        # read_count = full_depth_count 
        # genome_depth = round(full_depth_count * 150 / ref_len , 2)
        # summary.write(','.join([str(genome_cov), str(read_count), str(genome_depth)]))
        # log.info(f'Genome coverage: {genome_cov}%')
        # log.info(f'Genome depth: {genome_depth}X')
        # log.info(f'Read count: {read_count}')
        blank_depth = get_full_mapped_reads_depth(args.blankbam)
        plot_mapped_series(sample_depths, blank_depth, args.output, sample_name)
    else: 
        log.error(f'BAM file does not exist: {args.samplebam}')


if __name__ == '__main__':
    start_time = time.time()
    log.setLevel(logging.DEBUG)
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