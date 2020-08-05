import subprocess
import os 
import json 

read_dir = '/home/ubuntu/transfer/incoming/QIB_Sequencing/Covid-19_Seq/200415.coronahit/result.illumina.20200418/ncovIllumina_fastqMergeFourLanes_fastqMergeLanes'
all_bases = 0 
all_reads = 0 
for read in os.listdir(read_dir):
    read_path = os.path.join(read_dir, read)
    all_means = [] 
    total_reads = 0 
    total_phred = 0
    read_name = read.split('_')[0]
    subprocess.call(f'zcat {read_path} | fastp --stdin', shell=True)
    js = json.load(open('fastp.json'))
    all_reads += js['summary']['before_filtering']['total_reads']
    all_bases += js['summary']['before_filtering']['total_bases']
    print(f'{all_reads}\t{all_bases}')

all_quals = [] 
# Get average PHRED score. 
all_total_reads = 0 
all_total_reads_bases = 0
for read in os.listdir(read_dir):
    read_path = os.path.join(read_dir, read)
    all_means = [] 
    total_reads = 0 
    total_phred = 0
    read_name = read.split('_')[0]
# with gzip.open("practicezip.fasta.gz", "r") as handle:
 #   for record in SeqIO.parse(handle, "fasta"):
 
    with gzip.open(read_path, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            #if total_reads > 1000:
                #break     
            total_phred += mean(record.letter_annotations["phred_quality"])
            total_reads += 1 
            all_total_reads += 1
            all_total_reads_bases += len(record)
        print(f'{read_name}\t{round(total_phred / total_reads)}')
        all_quals.append(round(total_phred / total_reads))        
print(mean(all_quals))
print(all_total_reads)
print(all_total_reads_bases)

all_quals = [] 

# for sample in os.listdir(bam_dir):
#     sample_path = os.path.join(bam_dir, sample)
#     if not os.path.isdir(sample_path):
#         continue
#     bams = [os.path.join(sample_path, bam) for bam in os.listdir(sample_path) if bam.endswith('.bam')]
#     for bam_file in bams:
#         samfile = pysam.AlignmentFile(bam_file, "rb")
#         total_phred = 0 
#         total_read = 0 
#         for read in samfile:
#             if total_read > 1000:
#                 break            
#             total_phred += mean(read.query_alignment_qualities)
#             total_read += 1
#         print(f'{sample}\t{round(total_phred / total_read)}')
#         all_quals.append(round(total_phred / total_read))
# print(mean(all_quals))
        
