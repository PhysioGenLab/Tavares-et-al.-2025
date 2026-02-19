#!/bin/bash
#SBATCH --job-name=TELOMEROxADHD   


#SBATCH --ntasks=1                    
#SBATCH --mem=30gb                     
#SBATCH --time=192:00:00               
#SBATCH --output=/home/LAVA/outputs/LAVA_TELOMERExADHD_%j.log   
pwd; hostname; date

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

echo $SLURM_JOB_ID              
echo $SLURM_SUBMIT_DIR          
echo $SLURM_JOB_NODELIST        
echo $SLURM_NTASKS              

# arg = commandArgs(T); ref.prefix = arg[1]; loc.file = arg[2]; info.file = arg[3]; sample.overlap.    file = arg[4]; phenos = unlist(strsplit(arg[5],";")); out.fname = arg[6]

srun Rscript /home/LAVA/LAVA_Rscript.r "/home/LAVA/g1000_eur" \
"/home/LAVA/LAVA_2500loci.txt" \
"/home/LAVA/input.info.txt" \
"/home/LAVA/sample.overlap_telxadhd.txt" \
