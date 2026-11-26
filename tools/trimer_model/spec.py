from pathlib import Path
from prosapia.core import Tool
from .collect_trimer_model import collect_trimer_model
from .run_trimer_model import build_trimer_model_manifest

TOOL = Tool(
    name="trimer_model",
    action="update",
    description="Build a trimer model by applying the A->B symmetry transform twice to the dimer design.",
    default_sbatch=str(Path(__file__).parent / "trimer_model.sbatch"),
    default_input_column="rfdiffusion3_path",
    build_manifest_fn=build_trimer_model_manifest,
    collect_fn=collect_trimer_model,
)
