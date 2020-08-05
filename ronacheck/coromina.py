
# Protoype for cleaning up ilumina seq run for sars2cov2
import pysam
import subprocess
import os 

def check_config(reference_file,  r1, r2, workdir):
    errors = [] 
    if not os.path.exists(workdir):
        os.mkdir(workdir)
    if not os.path.exists(reference_file):
        errors.append('No refererence file at ' + reference_file)
    if not os.path.exists(r1):
        errors.append('No r1 file at ' + r1)
    if not os.path.exists(r2):
        errors.append('No r1 file at ' + r2)
    # Is read mapper available. 
    version_info = subprocess.check_output(['bowtie2', '--version']).decode('utf-8').split('\n')[0]
    if not os.path.basename(version_info.split()[0]) == 'bowtie2-align-s':
        errors.append('bowtie2 is not available. Check your PATH')
    return errors

def index_ref(reference_file, workdir):
    index_dir = os.path.join(workdir, 'index')
    index_file_name  = os.path.join(index_dir, os.path.basename(reference_file))
    index_file = os.path.join(index_dir , os.path.basename(reference_file) + '.1.bt2')
    if not os.path.exists(index_file):
        if not os.path.exists(index_dir):
            os.makedirs(index_dir)
        subprocess.call(['bowtie2-build', reference_file, index_file_name])
    return index_file_name


def map_to_ref(reference_file, r1, r2, workdir):
    # Create index if needed
    ref_index = index_ref(reference_file, workdir)
    r1_name = os.path.basename(r1)
    ref_name = os.path.basename(ref_index)
    bam_file = os.path.join(workdir, f'{ref_name}vs{r1_name}.bam')
    if not os.path.exists(bam_file):
        command = ['bowtie2','-x' , ref_index, '-1', r1, '-2', r2, '|', 'samtools', 'view', '-F', '4','-bS', '-', '>', bam_file]
        subprocess.call(' '.join(command), shell=True)
    return bam_file

def main(reference_file,  r1, r2, workdir, read_len_filter = 110): 

    init_errors = check_config(reference_file,  r1, r2, workdir)
    if len(init_errors) == 0:
        # Map to reference genome 
        bam_file = map_to_ref(reference_file, r1, r2, workdir)
        # Remove based on mapping quality
        sam_file = pysam.AlignmentFile(bam_file, 'rb')
        for read in sam_file:
            print(read)

        # Remove based on mapping length. 

        # Calculate variants on remainder

        # Report regions to mask for real samples.         
    else:
        print('We got errors on startup')
        print('\n\t'.join(init_errors))


ref_file = '/home/ubuntu/scratch/covid/ref/hcov_ref.fasta'
r1 = '/home/ubuntu/scratch/covid/run_2_mismatch_0/COG/Blank_COGPL5_rev_S192_R1_001.fastq.gz'
r2 = '/home/ubuntu/scratch/covid/run_2_mismatch_0/COG/Blank_COGPL5_rev_S192_R2_001.fastq.gz'
workdir = '/home/ubuntu/scratch/covid/coromina_temp'
main(ref_file, r1, r2, workdir)


