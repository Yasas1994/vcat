#!/usr/bin/env python3
"""
author: Yasas Wijesekara (yasas.wijesekara@uni-greifswald.de)

This script corrects taxonomy information from the protein-profile search using
cluster LCAs extracted by cluster_lca.py (mmseqs_pprofiles_lca.tsv).
Run after mmseqs_pprofile rule.
"""

import shutil
from pathlib import Path

import click
import polars as pl
import taxopy

DEFAULT_HEADER = (
    "query,target,theader,fident,qlen,tlen,alnlen,mismatch,gapopen,"
    "qstart,qend,tstart,tend,evalue,bits,taxid,taxname,taxlineage"
)


@click.command(context_settings=dict(show_default=True))
@click.option(
    "-d",
    "--database-dir",
    type=click.Path(path_type=Path, exists=True, file_okay=False, dir_okay=True),
    required=True,
    help="Root database directory (contains ictv-taxdump/ and VMR_latest/).",
)
@click.option(
    "-i",
    "--input",
    "prof",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Input mmseqs profile m8 file to correct (will be overwritten in place).",
)
@click.option(
    "--header",
    default=DEFAULT_HEADER,
    show_default=False,
    help="Comma-separated header for the input m8 file.",
)
@click.option(
    "--lca-table",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=None,
    help="Override path to mmseqs_pprofiles_lca.tsv (defaults to <database-dir>/VMR_latest/mmseqs_pprofiles/mmseqs_pprofiles_lca.tsv).",
)
def main(database_dir: Path, prof: Path, header: str, lca_table: Path | None) -> None:
    database_dir = database_dir.resolve()
    prof = prof.resolve()

    taxdump = database_dir / "ictv-taxdump"
    nodes = taxdump / "nodes.dmp"
    names = taxdump / "names.dmp"
    merged = taxdump / "merged.dmp"

    if not nodes.exists() or not names.exists() or not merged.exists():
        raise click.ClickException(
            f"Missing taxdump files under: {taxdump} (need nodes.dmp, names.dmp, merged.dmp)"
        )

    lca_path = (
        lca_table.resolve()
        if lca_table is not None
        else (database_dir / "VMR_latest/mmseqs_pprofiles/mmseqs_pprofiles_lca.tsv").resolve()
    )
    if not lca_path.exists():
        raise click.ClickException(f"LCA table not found: {lca_path}")

    # Load taxonomy DB
    taxdb = taxopy.TaxDb(
        nodes_dmp=str(nodes),
        names_dmp=str(names),
        merged_dmp=str(merged),
    )

    # Load LCA table and annotate with name/lineage
    clustlca = pl.read_csv(lca_path, separator="\t")

    def _taxon_name(x):
        try:
            return taxopy.Taxon(int(x), taxdb).name
        except Exception:
            return "NA"

    def _taxon_lineage(x):
        try:
            return str(taxopy.Taxon(int(x), taxdb))
        except Exception:
            return "NA"

    clustlca = clustlca.with_columns(
        pl.col("taxid").map_elements(_taxon_name, return_dtype=pl.String).alias("taxname"),
        pl.col("taxid").map_elements(_taxon_lineage, return_dtype=pl.String).alias("taxlineage"),
    )

    # Load profile m8 and replace tax columns using the LCA table
    cols = [c.strip() for c in header.split(",") if c.strip()]
    ictv_prof = pl.read_csv(prof, has_header=False, separator="\t", new_columns=cols)

    # Drop existing taxonomy columns if present
    drop_cols = [c for c in ["taxid", "taxname", "taxlineage"] if c in ictv_prof.columns]
    if drop_cols:
        ictv_prof = ictv_prof.drop(drop_cols)

    # Join on target cluster id
    # Keep only needed columns from clustlca; tolerate extra columns in file
    keep_from_lca = [c for c in ["target", "taxid", "taxname", "taxlineage", "root_p1", "taxrank"] if c in clustlca.columns]
    clustlca_small = clustlca.select(keep_from_lca)

    out = ictv_prof.join(
        clustlca_small.drop([c for c in ["root_p1", "taxrank"] if c in clustlca_small.columns]),
        on="target",
        how="inner",
    )

    tmp = prof.with_suffix(prof.suffix + "_tmp")
    out.write_csv(tmp, include_header=False, separator="\t")

    # Replace original atomically-ish
    shutil.move(str(tmp), str(prof))

    click.echo(f"Updated: {prof}")


if __name__ == "__main__":
    main()
