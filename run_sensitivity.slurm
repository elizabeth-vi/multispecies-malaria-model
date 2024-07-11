#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
# CPUS = 48*2 = 96
#SBATCH --cpus-per-task=48
#SBATCH --time=0-24:00:00

# Send yourself an email with job updates:
#SBATCH --mail-user="elizabeth.ivory@unimelb.edu.au"
#SBATCH --mail-type=ALL


# Load required modules
module purge
module load foss/2022a Python/3.10.4 SciPy-bundle/2022.05


# Launch multiple process python code
echo "Starting run"
time srun python3 sensitivity.py