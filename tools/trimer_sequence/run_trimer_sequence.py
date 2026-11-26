#!/usr/bin/env python3
"""
Build trimer sequences from MPNN/rfd4 sequences.

Extracts the linker from the input sequence and builds:
  native_seq + linker + native_seq + linker + native_seq

Appends trimer_sequence column to the table.
"""

import json
from argparse import ArgumentParser
from pathlib import Path
from typing import cast

from prosapia.core import CommonArgs, ManifestCtx


class TrimerSequenceArgs(CommonArgs):
    protein_info: str
    sequence_column: str


def add_trimer_sequence_args(parser: ArgumentParser) -> None:
    parser.add_argument(
        "--protein-info",
        type=str,
        default="protein_info.json",
        help="Path to protein_info.json. Default: protein_info.json in run_dir.",
    )
    parser.add_argument(
        "--sequence-column",
        type=str,
        default="proteinmpnn_sequence",
        help="Column with the full fused dimer sequence. Default: proteinmpnn_sequence.",
    )


def build_trimer_sequence_manifest(ctx: ManifestCtx[TrimerSequenceArgs]) -> list[tuple[str, ...]]:
    info_path = Path(ctx.args.protein_info)
    if not info_path.is_absolute():
        info_path = ctx.args.run_dir / info_path
    info = json.loads(info_path.read_text())
    subunit_len = info["subunit_len"]
    native_seq = info["native_sequence"]
    seq_col = ctx.args.sequence_column

    out_dir = ctx.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest_rows = []
    for name in ctx.ready.index:
        name = cast(str, name)
        full_seq = str(ctx.ready.at[name, seq_col])
        linker = full_seq[subunit_len:len(full_seq) - subunit_len]
        if not linker:
            continue
        trimer_seq = native_seq + linker + native_seq + linker + native_seq
        # Write just the sequence to a text file
        seq_path = out_dir / f"{name}.txt"
        seq_path.write_text(trimer_seq)
        manifest_rows.append((name, str(seq_path)))

    return manifest_rows
