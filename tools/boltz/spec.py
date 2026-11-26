from pathlib import Path

from prosapia.core import get_builtin

TOOL = get_builtin("boltz").with_overrides(
    default_sbatch=str(Path(__file__).parent / "boltz.sbatch")
)
