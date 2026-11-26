# Source from a project directory before using sapia:
#   activates the prosapia env, loads ./.env, makes SLURM pass the environment to jobs,
#   and defines helpers that run the tools repo's scripts in its uv env
source /data/brussel/118/vsc11875/envs/prosapia/bin/activate
export SBATCH_EXPORT=ALL
if [ -f .env ]; then set -a; source .env; set +a; fi

_TOOLS=/data/brussel/118/vsc11875/prosapia-tools
_tools_run() { local s=$1; shift; uv run --frozen --project "$_TOOLS" python "$_TOOLS/scripts/$s" "$@"; }
sapia_annotate()        { _tools_run annotate.py "$@"; }
sapia_build_templates() { _tools_run build_templates.py "$@"; }
