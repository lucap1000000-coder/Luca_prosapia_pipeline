from pathlib import Path
from prosapia.core import Tool
from .collect_trimer_sequence import collect_trimer_sequence
from .run_trimer_sequence import add_trimer_sequence_args, build_trimer_sequence_manifest

TOOL = Tool(
    name="trimer_sequence",
    action="update",
    description="Build trimer sequences by appending linker+subunit copies to MPNN/rfd4 sequences.",
    default_sbatch=str(Path(__file__).parent / "trimer_sequence.sbatch"),
    default_input_column="proteinmpnn_sequence",
    add_run_args_fn=add_trimer_sequence_args,
    build_manifest_fn=build_trimer_sequence_manifest,
    collect_fn=collect_trimer_sequence,
)
