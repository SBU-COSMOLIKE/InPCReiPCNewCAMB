#!/bin/bash

#SBATCH --job-name=MC
#SBATCH --output=MC-%A_%a.out
#SBATCH -e MC-%A_%a.out.err
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --ntasks-per-socket=24
#SBATCH --cpus-per-task=2
#SBATCH --partition=hbm-long-96core
#SBATCH --time=48:00:00
#SBATCH --qos=user-high-100-priority-qos

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo Slurm job NAME is $SLURM_JOB_NAME
echo Slurm job ID is $SLURM_JOBID
echo Number of task is $SLURM_NTASKS
echo Number of cpus per task is $SLURM_CPUS_PER_TASK

cd $SLURM_SUBMIT_DIR
if [ -z "${ROOTDIR}" ]; then
  sleep $(( 5 + SLURM_ARRAY_TASK_ID*25 )) # help avoid different scripts running start_cocoa.sh simultaneously
                                          # clash can still happen when one script runs start_cocoa.sh, the other
                                          # script is reading a file (start_cocoa.sh invokes stop_cocoa.sh)
                                          # which delete links and ROOTDIR so the running script can find file/likelihood
  conda activate cocoa
  source start_cocoa.sh
fi

export OMP_PROC_BIND=close
export OMP_PLACES=cores 
export OMP_DYNAMIC=FALSE
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
if [ -n "$SLURM_CPUS_PER_TASK" ]; then
  export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
else
  export OMP_NUM_THREADS=1
fi

mpirun -n ${SLURM_NTASKS} --oversubscribe --mca pml ^ucx  \
  --mca btl vader,tcp,self --bind-to core:overload-allowed \
  --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS}  \
  python ./projects/InPCReiPC/EXAMPLE_MINIMIZE${SLURM_ARRAY_TASK_ID}.py --nstw 500 \
    --root ./projects/InPCReiPC/ --outroot "EXAMPLE_MINIMIZE${SLURM_ARRAY_TASK_ID}"