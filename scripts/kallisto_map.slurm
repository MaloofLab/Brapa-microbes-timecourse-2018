#!/bin/bash
#SBATCH --partition=dev # partition to submit to
#SBATCH --job-name=Brapa_Kallisto # Job name
#SBATCH --array=0-1 #for this script adjust to match number of fastq files
#SBATCH --nodes=1 # single node, anything more than 1 will not run
#SBATCH --ntasks=01 # equivalent to cpus, stick to around 20 max on gc64, or gc128 nodes
#SBATCH --mem=4000 # in MB, memory pool all cores, default is 2GB per cpu
#SBATCH --time=0-01:00:00  # expected time of completion in hours, minutes, seconds, default 1-day
#SBATCH --output=Kallisto_%A_%a.out # STDOUT
#SBATCH --error=Kallisto_%A_%a.err # STDERR
#SBATCH --mail-user=jnmaloof@ucdavis.edu #
#SBATCH --mail-type=ALL #

# This will be run once for a single process

/bin/hostname

start=`date +%s`

# Load Kallisto

module load kallisto

# Change directory

cd /share/malooflab/Fastqs/Brapa_microbe/timecourse-2018

# Identify each array run

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

# create an array of file names:
filelist=($(ls *gz))

# now pick the file that corresponds to the current array
# note that for this script the number of arrays should equal the number of files
f=${filelist[${SLURM_ARRAY_TASK_ID}]}

# trim off directory info and file extensions:
outdir=$(basename $f .fastq.gz)
echo "file stem: " $outdir

kallisto quant \
    --index /share/malooflab/ref_genomes/B_rapa/V3.0/B_rapa_CDS_V3.0_k31_kallisto_index   \
    --output-dir /share/malooflab/Kazu/Brapa-microbes-timecourse-2018/Kallisto_output/$outdir \
    --plaintext \
    --single \
    -l 250 \
    -s 40 \
    $f

end=`date +%s`
runtime=$((end-start))
echo $runtime seconds to completion
