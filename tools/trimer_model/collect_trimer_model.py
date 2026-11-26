import json
from typing import Iterable

import pandas as pd

from prosapia.core import Collected, CollectCtx, CollectEach, DesignCtx

KEYS = ("trimer_model_path", "trimer_linker_len")


def collect_trimer_model(ctx: CollectCtx) -> CollectEach:
    files = {p.stem: p for p in ctx.out_dir.glob("*.json")}
    na = {k: pd.NA for k in KEYS}

    def one(d: DesignCtx) -> Iterable[Collected]:
        f = files.get(d.name)
        if f is None:
            yield Collected(status=f"missing: no .json for {d.name}", data=na)
            return
        try:
            data = json.loads(f.read_text())
            yield Collected(status="OK", data={
                "trimer_model_path": data.get("trimer_model_path", pd.NA),
                "trimer_linker_len": data.get("linker_len", pd.NA),
            })
        except Exception as e:
            yield Collected(status=f"error: {e}", data=na)

    return one
