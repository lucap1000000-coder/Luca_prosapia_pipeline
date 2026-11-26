#!/usr/bin/env python3
"""Collect RFdiffusion4 outputs into a child table."""

import json
from pathlib import Path
from typing import Iterable

import gemmi
import pandas as pd

from prosapia.core import Collected, CollectCtx, CollectEach, DesignCtx

KEYS = ("path", "sequence", "linker_len")


def collect_rfdiffusion4(ctx: CollectCtx) -> CollectEach:
    out_dir = ctx.out_dir
    samples_dir = out_dir / "inference_samples"

    info_path = ctx.args.run_dir / "protein_info.json"
    subunit_len = None
    if info_path.exists():
        info = json.loads(info_path.read_text())
        subunit_len = info.get("subunit_len")

    cif_map: dict[str, Path] = {}
    if samples_dir.exists():
        for design_dir in samples_dir.iterdir():
            if design_dir.is_dir():
                cifs = list(design_dir.glob("*_generated.cif"))
                if cifs:
                    cif_map[design_dir.name] = cifs[0]

    def one(d: DesignCtx) -> Iterable[Collected]:
        children = {k: v for k, v in cif_map.items() if k.startswith(d.name + "_") or k == d.name}

        if not children:
            return

        for design_name, cif_path in sorted(children.items()):
            try:
                st = gemmi.read_structure(str(cif_path))
                chain_a = st[0]["A"]
                residues = [r for r in chain_a if r.find_atom("CA", "*")]
                seq = gemmi.one_letter_code([r.name for r in residues])
                linker_len = len(residues) - 2 * subunit_len if subunit_len else pd.NA
                yield Collected(
                    name=design_name,
                    parent=d.name,
                    data={
                        "path": str(cif_path),
                        "sequence": seq,
                        "linker_len": linker_len,
                    },
                )
            except Exception as e:
                yield Collected(
                    name=design_name,
                    parent=d.name,
                    status=f"error: {e}",
                    data={k: pd.NA for k in KEYS},
                )

    return one
