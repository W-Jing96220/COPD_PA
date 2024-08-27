#!/bin/bash
#SBATCH -J singlecell_exp1127_part6
#SBATCH --account=def-ubclxing # for beluga=rrg
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40 # if ask for entire node (SBATCH --nodes=1) on Graham=32 for Cedar=48 for beluga=40
#SBATCH --mem=191000M # allocates all available memory on the compute node for the job # for graham=128000M for Cedar=192000M for beluga=191000M/771000M
#SBATCH	--mail-user=qun482@usask.ca
#SBATCH	--mail-type=ALL
#SBATCH -t 4:00:00 # Running time of 4day
module load r/4.0.2              # Adjust version and add the gcc module used for installing packages.

# Based on Cluster_Code_ReleaseQuality/originalPA-Java/Biomine_Java_FullStagesTesting.sh

module load python/3.6.3
virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip

pip install --no-index numpy pandas

# Now do my computations here on the local disk using the contents of the extracted archive...
# By default, print in Python is buffered, meaning that it does not write to files or stdout immediately, 
#    and needs to be 'flushed' to force the writing to stdout immediately. -u Force stdin, 
#    stdout and stderr to be totally unbuffered. On systems where it matters, also put stdin, stdout and stderr in binary mode. 
