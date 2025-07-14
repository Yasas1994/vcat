#!/usr/bin/env python

"""
author: Yasas Wijesekara (yasas.wijesekara@uni-greifswald.de)

This scripts corrects the taxonomy information from the protein-profile search using
the cluster lcas extracted by cluster_lca.py script  (mmseqs_pprofiles_lca.tsv)

this will be run after mmseqs_pprofile rule

"""

import polars as pl
import taxopy
import shutil
from sys import argv
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,  # Set log level
    format="%(asctime)s - %(levelname)s - %(message)s",  # Log format
    handlers=[
        logging.StreamHandler()  # Print logs to the console
    ],
)

# Create a logger instance
logger = logging.getLogger("[vcat]")

DATABASE_DIR = argv[1]  #'/media/ssd/ICTV-TaxonomyChallenge/vcat/databases'
# input m8 file
PROF = argv[
    2
]  # "/media/ssd/ICTV-TaxonomyChallenge/vcat/databases/VMR_latest/cuttoffs/prof/genomes_fna_prof.m8"
# header of the m8 file
HEADER = "query,target,theader,fident,qlen,tlen,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,taxid,taxname,taxlineage"


def lca(taxa, taxdb, fraction=0.6):
    t = taxa.to_list()
    if len(t) > 1:
        return taxopy.find_majority_vote(t, taxdb=taxdb, fraction=fraction)
    return taxa.iloc[0]


def trim_lineage(row, rank, taxdb):
    taxon = taxopy.Taxon(row["taxid"], taxdb=taxdb)
    # add support for subranks in the future
    # if sub-rank is unavailable round off the to nearest whole rank
    ranks = [
        "realm",
        "subrealm",
        "kingdom",
        "subkingdom",
        "phylum",
        "subphylum",
        "class",
        "subclass",
        "order",
        "suborder",
        "family",
        "subfamily",
        "genus",
        "subgenus",
        "species",
        "subspecies",
    ]
    # case one
    if ranks.index(taxon.rank) <= ranks.index(rank):
        return taxon, taxon.taxid, taxon.rank
    else:
        parent = taxon.parent(taxdb)
        daughter = taxon
        while (
            parent.rank != rank
            and parent.rank != "no rank"
            and ranks.index(parent.rank) >= ranks.index(rank)
        ):
            daughter = parent
            parent = parent.parent(taxdb)
        if parent.rank != "no rank":
            return parent, parent.taxid, parent.rank
        else:
            return daughter, daughter.taxid, daughter.rank


taxdb = taxopy.TaxDb(
    nodes_dmp=f"{DATABASE_DIR}/ictv-taxdump/nodes.dmp",
    names_dmp=f"{DATABASE_DIR}/ictv-taxdump/names.dmp",
    merged_dmp=f"{DATABASE_DIR}/ictv-taxdump/merged.dmp",
)

clustlca = pl.read_csv(
    f"{DATABASE_DIR}/VMR_latest/mmseqs_pprofiles/mmseqs_pprofiles_lca.tsv",
    separator="\t",
)

clustlca = clustlca.with_columns(
    pl.col("taxid")
    .map_elements(lambda x: taxopy.Taxon(x, taxdb).name, return_dtype=str)
    .alias("taxname"),
    pl.col("taxid")
    .map_elements(lambda x: str(taxopy.Taxon(x, taxdb)), return_dtype=str)
    .alias("taxlineage"),
    pl.col("taxid")
    .map_elements(lambda x: taxopy.Taxon(x, taxdb).rank, return_dtype=str)
    .alias("taxrank"),
)


ictv_prof = pl.read_csv(
    PROF, has_header=False, separator="\t", new_columns=HEADER.split(",")
)

ictv_prof = ictv_prof.drop(["taxid", "taxname", "taxlineage"]).join(
    clustlca.drop(["root_p1", "taxrank"]), on="target", how="inner"
)

ictv_prof.write_csv(PROF + "_tmp", include_header=False, separator="\t")

# replace the original with the tmp file
shutil.move(PROF + "_tmp", PROF)
