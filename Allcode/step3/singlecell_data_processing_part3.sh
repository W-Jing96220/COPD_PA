#!/bin/bash
#SBATCH -J singlecell_exp1106_part3
#SBATCH --account=def-ubclxing
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12 # if ask for entire node (SBATCH --nodes=1) on Graham=32 for Cedar=48 for beluga=40
#SBATCH --mem=64000M # allocates all available memory on the compute node for the job # for graham=128000M for Cedar=192000M for beluga=191000M/771000M
#SBATCH	--mail-user=qun482@usask.ca
#SBATCH	--mail-type=ALL
#SBATCH -t 1:00:00 # Running time of 4day
module load r/4.0.2              # Adjust version and add the gcc module used for installing packages.

# Based on Cluster_Code_ReleaseQuality/single_cell_data_experiment/singlecell_processing_then_convert_to_webgraph_thenPA_onlyControlCOPD.sh

module load python/3.6
source ~/ENV/bin/activate

# specify the local installation directory according to currently R module that is loaded and install needed package.
mkdir -p ~/.local/R/$EBVERSIONR/
export R_LIBS=~/.local/R/$EBVERSIONR/

Rscript /lustre06/project/6040537/jingw/singlecell_data_processing_part3.R Macrophage /lustre06/project/6040537/jingw
Rscript /lustre06/project/6040537/jingw/singlecell_data_processing_part3.R Macrophage_Alveolar /lustre06/project/6040537/jingw
Rscript /lustre06/project/6040537/jingw/singlecell_data_processing_part3.R cMonocyte /lustre06/project/6040537/jingw
Rscript /lustre06/project/6040537/jingw/singlecell_data_processing_part3.R ncMonocyte /lustre06/project/6040537/jingw
