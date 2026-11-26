from prosapia.core import ManifestCtx


def build_linker_clash_manifest(ctx: ManifestCtx):
    rows = []
    for name, row in ctx.ready.iterrows():
        subunit_len = int(row["subunit_len"])
        rows.append((name, str(row[ctx.args.input_column]), str(subunit_len)))
    return rows
