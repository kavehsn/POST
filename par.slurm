#!/bin/csh

# simple template job script for submission of a serial MATLAB program (.m file)
# using a single computational thread (core) job with SLURM directives

# In this script the minimum requirements have been set for SLURM and 
# it only requires editing the last line and replace 

# my_matlab_M_file.m 

# with the name of the actual MATLAB program (.m file) to execute.
# The job can be submitted with the command

# sbatch -p  name_of_serial_queue  name_of_this_slurm_job_script

# If successful SLURM will return a jobID for this job which can be
# used to query its status.

#############################################################################

## All lines that start with #SBATCH will be processed by SLURM.
## Lines in this template script that have white space between # and SBATCH 
## will be ignored. They provide examples of further options that can be
## activated by deleting the white space and replacing any text after the 
## option.

## By default SLURM uses as working directory the directory from where the
## job script is submitted. To change this the standard Linux cd command has
## to be used.

## Name of the job as it will appear when querying jobs with squeue (the
## default is the name of the job script)

# SBATCH  -J  job_name

## By default SLURM combines the standard output and error streams in a single
## file based on the jobID and with extension .out
## These streams can be directed to separate files with these two directives

# SBATCH  -o  out_file_name.o%j
# SBATCH  -e  err_file_name.e%j

## where SLURM will expand %j to the jobID of the job.

## Request email to the user when certain type of events occur to 
## the job

# SBATCH  --mail-type=<type>

## where <type> can be one of BEGIN, END, FAIL, REQUEUE or ALL,
## and send to email address

# SBATCH  --mail-user  email_name@durham.ac.uk

## The default email name is that of the submitting user as known to the system.

## Specify project or account name (currently not required).
##SBATCH -A ITSall

#############################################################################

module purge
module load slurm/current

## Load any other modules required here
module load matlab/current

mkdir -p slurm_task_id_/$SLURM_ARRAY_TASK_ID

## Execute MATLAB program
/bin/rm -f slurm_task_id_$SLURM_ARRAY_TASK_ID
echo $SLURM_ARRAY_TASK_ID > slurm_task_id_$SLURM_ARRAY_TASK_ID
#cat slurm_task_id_$SLURM_ARRAY_TASK_ID

matlab -nodisplay -nodesktop -singleCompThread -r "run ./wrapperPar.m ; quit" \
 < slurm_task_id_$SLURM_ARRAY_TASK_ID

/bin/rm -f slurm_task_id_$SLURM_ARRAY_TASK_ID

rm -rf slurm_task_id_/$SLURM_ARRAY_TASK_ID


