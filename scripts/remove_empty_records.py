#!/usr/bin/env python3
from pathlib import Path

import click
from Bio import SeqIO


@click.command(context_settings=dict(show_default=True))
@click.option(
    "-i",
    "--input",
    "path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Input FASTA file.",
)
@click.option(
    "-o",
    "--output",
    "outpath",
    type=click.Path(path_type=Path, dir_okay=False),
    default=None,
    help="Output FASTA file (default: <input>.tmp).",
)
def main(path: Path, outpath: Path | None) -> None:
    """Remove empty sequences from a FASTA file."""
    seq_records = [rec for rec in SeqIO.parse(path, "fasta") if len(rec.seq) > 0]

    if outpath is None:
        outpath = path.parent / f"{path.name}.tmp"

    with open(outpath, "w") as handle:
        SeqIO.write(seq_records, handle, "fasta")

    click.echo(f"Wrote {len(seq_records)} sequences to {outpath}")


if __name__ == "__main__":
    main()
