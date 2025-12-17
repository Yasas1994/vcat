#!/usr/bin/env python

"""
author: Yasas Wijesekara (yasas.wijesekara@uni-greifswald.de)

This scripts adds lineage information to .m8 files with a taxid column
this will be run after blast_genome rule

"""

#!/usr/bin/env python3
import argparse
import logging
from pathlib import Path
import polars as pl
import taxopy
import re

DEFAULT_IN_HEADER = (
    "query,target,theader,fident,qlen,tlen,alnlen,mismatch,gapopen,"
    "qstart,qend,tstart,tend,evalue,bits,taxid"
)

DEFAULT_OUT_HEADER = (
    "query,target,theader,fident,qlen,tlen,alnlen,mismatch,gapopen,"
    "qstart,qend,tstart,tend,evalue,bits,taxid,taxname,taxlineage"
)

DEFAULT_TAXID = "taxid"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="vcat",
        description="VCAT helper: load taxonomy DB and profile hits with configurable paths.",
    )
    p.add_argument(
        "-d",
        "--database-dir",
        type=Path,
        required=True,
        help="Root database directory (contains ictv-taxdump/ and VMR_latest/).",
    )
    p.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Input m8 file path.",
    )
    p.add_argument(
        "--in_header",
        default=DEFAULT_IN_HEADER,
        help="Comma-separated header for the m8 file (used downstream).",
    )
    p.add_argument(
        "--out_header",
        default=DEFAULT_OUT_HEADER,
        help="Comma-separated header for the m8 output file.",
    )
    p.add_argument(
        "--taxid_col",
        default=DEFAULT_TAXID,
        help="column header with taxid",
    )
    p.add_argument(
        "--taxdump-dir",
        type=Path,
        help="Override path to ictv-taxdump/ (defaults to <database-dir>/ictv-taxdump).",
    )
    p.add_argument(
        "--log-level",
        default="INFO",
        choices=["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"],
        help="Logging verbosity.",
    )
    return p.parse_args()


def setup_logging(level: str) -> logging.Logger:
    logging.basicConfig(
        level=getattr(logging, level),
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler()],
    )
    return logging.getLogger("vcat")


def main() -> None:
    args = parse_args()
    logger = setup_logging(args.log_level)

    # Resolve defaults that depend on database-dir
    database_dir: Path = args.database_dir.resolve()
    taxdump_dir: Path = (args.taxdump_dir or (database_dir / "ictv-taxdump")).resolve()
    table_path: Path = args.input.resolve()
    table_path_out: Path = table_path.with_suffix(table_path.suffix + ".tmp")

    # Basic validations
    for path, desc in [
        (database_dir, "database-dir"),
        (taxdump_dir / "nodes.dmp", "nodes.dmp"),
        (taxdump_dir / "names.dmp", "names.dmp"),
        (taxdump_dir / "merged.dmp", "merged.dmp"),
        (table_path, "hit_table.m8"),
    ]:
        if not path.exists():
            raise FileNotFoundError(f"Required path for {desc} not found: {path}")

    logger.info("Using database dir: %s", database_dir)
    logger.info("Using taxdump dir:   %s", taxdump_dir)
    logger.info("Using Table (m8):     %s", table_path)

    DEFAULT_IN_HEADER = args.in_header  # keep variable name if used downstream

    # Initialize taxonomy database
    taxdb = taxopy.TaxDb(
        nodes_dmp=str(taxdump_dir / "nodes.dmp"),
        names_dmp=str(taxdump_dir / "names.dmp"),
        merged_dmp=str(taxdump_dir / "merged.dmp"),
    )
    logger.info("Loaded TaxDb with %d nodes.", len(taxdb.taxid2parent))

    # Read table.m8
    table = pl.read_csv(
        table_path,
        separator="\t",
        has_header=False,
        new_columns=DEFAULT_IN_HEADER.split(","),
    )
    logger.info(
        "Loaded table.m8: %s rows.",
        f"{table.height:,}",
    )
    taxid2name = {
        i: taxopy.Taxon(i, taxdb=taxdb) for i in table[args.taxid_col].to_list()
    }

    def get_name(x, taxid2name: dict = taxid2name) -> str:
        if c := taxid2name.get(x, None):
            return c.name
        logging.warning("taxid2name.name returned None")
        return "NA"

    def search_pattern(x: str) -> str:
        search = re.search(r"(?<=\|)(\S+)(?=\|)", x)
        if search:
            return search.group(1)

        logging.warning("pattern search returned None")
        return ""

    table = table.with_columns(
        pl.col("fident") / 100,
        pl.col("target").map_elements(
            lambda x: search_pattern(x),
            return_dtype=pl.String,
        ),
    )
    table = table.with_columns(
        taxname=pl.col(args.taxid_col).map_elements(
            lambda x: get_name(x), return_dtype=pl.String
        ),
        taxlineage=pl.col(args.taxid_col).map_elements(
            lambda x: str(taxid2name.get(x)), return_dtype=pl.String
        ),
    )
    table.write_csv(table_path_out, separator="\t", include_header=False)


if __name__ == "__main__":
    main()
