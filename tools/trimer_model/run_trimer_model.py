from prosapia.core import ManifestCtx


def build_trimer_model_manifest(ctx: ManifestCtx):
    rows = []
    for name, row in ctx.ready.iterrows():
        subunit_len = int(row["subunit_len"])
        src = ctx.lookup(name, ctx.args.input_column)
        if src is None:
            continue
        rows.append((str(name), str(src), str(subunit_len)))
    return rows
