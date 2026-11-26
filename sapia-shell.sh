# Source from a project directory before using sapia:
#   activates the prosapia env, loads ./.env, and makes SLURM pass the environment to jobs
source /data/brussel/118/vsc11875/envs/prosapia/bin/activate
export SBATCH_EXPORT=ALL
if [ -f .env ]; then set -a; source .env; set +a; fi
