#!/bin/bash
#SBATCH -J singlecell_exp1106_part2
#SBATCH --account=def-ubclxing # for beluga=rrg
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40 # if ask for entire node (SBATCH --nodes=1) on Graham=32 for Cedar=48 for beluga=40
#SBATCH --mem=191000M # allocates all available memory on the compute node for the job # for graham=128000M for Cedar=192000M for beluga=191000M/771000M
#SBATCH	--mail-user=qun482@usask.ca
#SBATCH	--mail-type=ALL
#SBATCH -t 8:00:00 # Running time of 4day
# Load the required environment modules
module load StdEnv/2020
module load gcc/9.3.0
module load r/4.2.1


module load python/3.10
source ~/ENV/bin/activate myenv


# specify the local installation directory according to currently R module that is loaded and install needed package.
mkdir -p ~/.local/R/$EBVERSIONR/
export R_LIBS=~/.local/R/$EBVERSIONR/
unset LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/cvmfs/soft.computecanada.ca/gentoo/2023/x86-64-v3/usr/lib64:$LD_LIBRARY_PATH


# Run the R scripts
Rscript /lustre06/project/6040537/jingw/singlecell_data_processing_part2.R Macrophage /lustre06/project/6040537/jingw
Rscript /lustre06/project/6040537/jingw/singlecell_data_processing_part2.R Macrophage_Alveolar /lustre06/project/6040537/jingw
Rscript /lustre06/project/6040537/jingw/singlecell_data_processing_part2.R cMonocyte /lustre06/project/6040537/jingw
Rscript /lustre06/project/6040537/jingw/singlecell_data_processing_part2.R ncMonocyte /lustre06/project/6040537/jingw
