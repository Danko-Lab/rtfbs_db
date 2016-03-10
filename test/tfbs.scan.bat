#!/bin/bash 

#----------------------------------------------------
# Generic SLURM script 
#----------------------------------------------------

#SBATCH -J tfbs.class.test # Job name
#SBATCH -o myjob.%j.out # stdout; %j expands to jobid
#SBATCH -e myjob.%j.err # stderr; skip to combine stdout and stderr
#SBATCH -p normal # queue
#SBATCH -N 1 -n 1 # one node and one task
#SBATCH -t 48:00:00 # max time
#SBATCH --mail-user=zw355@cornell.edu
#SBATCH --mail-type=ALL

module load bedtools

$R --vanilla --no-save < tfbs.scan.R > tfbs.scan.out


