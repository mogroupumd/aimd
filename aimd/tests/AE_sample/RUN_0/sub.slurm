#!/bin/bash
#SBATCH --job-name pureLiAlNmd1500K	# job name
#SBATCH -N 2
#SBATCH -n 40				# total number of mpi tasks requested
#SBATCH -t 72:00:00			# run time (hh:mm:ss) 
#SBATCH -A mse-yfmo-hi			# User account 
#SBATCH -o run_aimd.out       		# output and error file name (%j expands to jobID)
#SBATCH --mail-user=account@email.com	# User email
#SBATCH --mail-type=end    		# Email user when job finishes

vasp_exe='/homes/alexre/vasp_5.4.1/vasp_gam'	# path to VASP executable


module load intel/2015.0.3.032
module list
OMP_NUM_THREADS=1
export OMP_NUM_THREADS

cp INCAR.heat INCAR			#For RUN_1 and on, this line should be: 'cp INCAR.MD INCAR'
mpirun -np $SLURM_NTASKS $vasp_exe -x $OMP_NUM_THREADS
