#!/bin/bash
#SBATCH --qos=short_qos
#SBATCH --job-name="25_0.9"
#SBATCH --partition=batch,medium,short
#SBATCH --ntasks=1
#SBATCH --mail-user=kirankri27@gmail.com
#SBATCH --array=1-50
#SBATCH --mem-per-cpu=1240
# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=0-15:00:00

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

# Run the job from the directory where it was launched (default)

# The job command(s):
date=`date +"%H:%M:%S on %d %b %Y"`
echo
echo "============================================================="
echo "Timing: Commenced at $date "
./sens #...${SLURM_ARRAY_TASK_ID}
date=`date +"%H:%M:%S on %d %b %Y"`
echo "Timing: Finished at $date "
echo "============================================================="


