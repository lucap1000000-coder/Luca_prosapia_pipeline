import pandas as pd
from typing import Iterable
from prosapia.core import Collected, CollectCtx, CollectEach, DesignCtx

KEYS = ("trimer_sequence",)

def collect_trimer_sequence(ctx: CollectCtx) -> CollectEach:
    na = {k: pd.NA for k in KEYS}

    def one(d: DesignCtx) -> Iterable[Collected]:
        seq_path = ctx.out_dir / f"{d.name}.txt"
        if seq_path.exists():
            seq = seq_path.read_text().strip()
            yield Collected(status="OK", data={"trimer_sequence": seq})
        else:
            yield Collected(status=f"missing: no .txt for {d.name}", data=na)

    return one
