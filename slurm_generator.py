for num in range(0,4):

    file_string = "submit_cori_nishant" + str(num)+ ".sh"

    f = open(file_string, "w+")

    main_string1 = """#!/bin/bash
#SBATCH -J create_lensed_map
#SBATCH -N 3
#SBATCH -q debug  # 30min limit for debug. Otherwise, use regular
#SBATCH -L SCRATCH   #Job requires $SCRATCH file system
#SBATCH -C haswell   #Use haswell/knl nodes
#SBATCH -t 00:30:00  # hh:mm:ss
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=nishant.mishra@berkeley.edu
#SBATCH -o /global/cscratch1/sd/nishant/output_messages.out
#SBATCH -e /global/cscratch1/sd/nishant/error_messages.err


python test.py """

    main_string2 = """

# to submit the job, type:
# sbatch ./submit_cori_nishant.sh

# to check your code is schedule:
# squeue -u nishant

# to kill a job, find the job id number with:
# squeue -u nishant
# then kill it with:
# scancel the_job_id_number_you_found
"""

    new_string = main_string1 + str(num) + main_string2

    f.write(new_string)
    f.close()
