#!/bin/bash
#SBATCH -J create_lensed_map
#SBATCH -N 1
#SBATCH -q regular  # 30min limit for debug. Otherwise, use regular
#SBATCH -L SCRATCH   #Job requires $SCRATCH file system
#SBATCH -C haswell   #Use haswell/knl nodes
#SBATCH -t 8:00:00  # hh:mm:ss
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=nishant.mishra@berkeley.edu
#SBATCH -o /global/cscratch1/sd/nishant/output_messages.out
#SBATCH -e /global/cscratch1/sd/nishant/error_messages.err

count=$1
file=input.txt
if [ -f $file ] ; then
 rm $file
fi
for i in $(seq $1)
do
echo test.py $i >> $file
done

xargs --arg-file=$file --max-lines=1 -P 5 python
#python test.py 2

#python test.py 2 &
#python test.py 4 &

# to submit the job, type:
# sbatch ./submit_cori_nishant.sh

# to check your code is schedule:
# squeue -u nishant

# to kill a job, find the job id number with:
# squeue -u nishant
# then kill it with:
# scancel the_job_id_number_you_found
# comment
