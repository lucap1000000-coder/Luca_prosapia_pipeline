from pathlib import Path

from prosapia.core import get_builtin

TOOL = get_builtin("proteinmpnn").with_overrides(
    default_sbatch=str(Path(__file__).parent / "proteinmpnn.sbatch")
)
