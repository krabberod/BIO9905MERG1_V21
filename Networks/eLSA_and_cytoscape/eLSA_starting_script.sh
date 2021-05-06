# This is the commands for running eLSA on abel, should be placed in a slurm jobscript.
# Remember to change <anderkkr> with your own username in the last export PYTHONPATH line. 
# 

module purge
module load python2/2.7.10.gnu
module load R/3.2.3.gnu
export PYTHONPATH=$PYTHONPATH:~/Programs/elsa/lib/python2.7/site-packages/
export PYTHONPATH=$PYTHONPATH:/cluster/home/anderkkr/.local/lib/python2.7/site-packages

lsa_compute <input_file> <output_name> -r 1 -s 24 -d 0 -p perm -x 1000 -f linear -n robustZ;

