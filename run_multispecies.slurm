#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --time=0-48:00:00

# Send yourself an email when the job:
# aborts abnormally (fails)
# SBATCH --mail-type=FAIL
# begins
# SBATCH --mail-type=BEGIN
# ends successfully
# SBATCH --mail-type=END
# SBATCH --mail-user="ivorye@student.unimelb.edu.au"


# Load required modules
module purge
module load foss/2022a Python/3.10.4 SciPy-bundle/2022.05


# Launch multiple process python code
echo "Starting run of 30 simulations"
time srun -n 30 python3 do_simulations_flexibly.py