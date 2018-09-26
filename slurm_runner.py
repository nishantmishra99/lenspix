import os

for num in range(0,4):
    string = "sbatch ./submit_cori_nishant" + str(num) +".sh"
    os.system(string)
