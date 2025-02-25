#!/usr/bin/env python
"""
author: Yasas Wijesekara (yasas.wijesekara@uni-greifswald.de)

postprocesses genome.m8, protein.m8 and profile.m8 files by calculating
ani, aai, and api and summarizes the summarizes the taxonomy predictions
to a single .tsv file

"""
import polars as pl
import numpy as np
import taxopy
import pyfastx
from sys import argv
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,  # Set log level
    format="%(asctime)s - %(levelname)s - %(message)s",  # Log format
    handlers=[
        logging.StreamHandler()  # Print logs to the console
    ]
)

# Create a logger instance
logger = logging.getLogger("[vcat]")


DATABASE_DIR = argv[1]
NUC = argv[2]
PROT = argv[3]
PROF = argv[4]
OUTFILE = argv[5]
# these should be moved to another place
QCOV_ANI = argv[6]
ANI = argv[7]

# fasta = pyfastx.Fasta(FASTA)
# headers = list(fasta.keys())

# taxdb = taxopy.TaxDb(nodes_dmp=f"{DATABASE_DIR}/ictv-taxdump/nodes.dmp",
#                      names_dmp=f"{DATABASE_DIR}/ictv-taxdump/names.dmp",
#                      merged_dmp=f"{DATABASE_DIR}/ictv-taxdump/merged.dmp")


keys = ['SequenceID', 'Realm (-viria)', 'Realm_score', 'Subrealm (-vira)',
       'Subrealm_score', 'Kingdom (-virae)', 'Kingdom_score',
       'Subkingdom (-virites)', 'Subkingdom_score', 'Phylum (-viricota)',
       'Phylum_score', 'Subphylum (-viricotina)', 'Subphylum_score',
       'Class (-viricetes)', 'Class_score', 'Subclass (-viricetidae)',
       'Subclass_score', 'Order (-virales)', 'Order_score',
       'Suborder (-virineae)', 'Suborder_score', 'Family (-viridae)',
       'Family_score', 'Subfamily (-virinae)', 'Subfamily_score',
       'Genus (-virus)', 'Genus_score', 'Subgenus (-virus)', 'Subgenus_score',
       'Species (binomial)', 'Species_score']

keys_full = ['SequenceID', 'Seqlen', "Score", "Method", 'Realm (-viria)', 'Realm_score', 'Subrealm (-vira)',
       'Subrealm_score', 'Kingdom (-virae)', 'Kingdom_score',
       'Subkingdom (-virites)', 'Subkingdom_score', 'Phylum (-viricota)',
       'Phylum_score', 'Subphylum (-viricotina)', 'Subphylum_score',
       'Class (-viricetes)', 'Class_score', 'Subclass (-viricetidae)',
       'Subclass_score', 'Order (-virales)', 'Order_score',
       'Suborder (-virineae)', 'Suborder_score', 'Family (-viridae)',
       'Family_score', 'Subfamily (-virinae)', 'Subfamily_score',
       'Genus (-virus)', 'Genus_score', 'Subgenus (-virus)', 'Subgenus_score',
       'Species (binomial)', 'Species_score']
n = [    
    pl.col("taxlineage").str.extract(r"_([A-Za-z]+viria);?", 1).alias('Realm (-viria)'),
    pl.lit(None).alias("Realm_score"),
    pl.col("taxlineage").str.extract(r"_([A-Za-z]+vira);", 1).alias('Subrealm (-vira)'),
    pl.lit(None).alias('Subrealm_score'),
    pl.col("taxlineage").str.extract(r"_([A-Za-z]+virae);?", 1).alias('Kingdom (-virae)'),
    pl.lit(None).alias('Kingdom_score'),
    pl.col("taxlineage").str.extract(r"_([A-Za-z]+virites);?", 1).alias('Subkingdom (-virites)'),
    pl.lit(None).alias('Subkingdom_score'),
    pl.col("taxlineage").str.extract(r"_([A-Za-z]+viricota);?", 1).alias('Phylum (-viricota)'),
    pl.lit(None).alias('Phylum_score'),
    pl.col("taxlineage").str.extract(r"_([A-Za-z]+viricotina);?", 1).alias('Subphylum (-viricotina)'),
    pl.lit(None).alias('Subphylum_score'),
    pl.col("taxlineage").str.extract(r"_([A-Za-z]+viricetes);?", 1).alias('Class (-viricetes)'),
    pl.lit(None).alias('Class_score'),
    pl.col("taxlineage").str.extract(r"_([A-Za-z]+viricetidae);?", 1).alias('Subclass (-viricetidae)'),
    pl.lit(None).alias('Subclass_score'),
    pl.col("taxlineage").str.extract(r"_([A-Za-z]+virales);?", 1).alias('Order (-virales)'),
    pl.lit(None).alias('Order_score'),
    pl.col("taxlineage").str.extract(r"_([A-Za-z]+virineae);?", 1).alias('Suborder (-virineae)'),
    pl.lit(None).alias('Suborder_score'),
    pl.col("taxlineage").str.extract(r"_([A-Za-z]+viridae);?", 1).alias('Family (-viridae)'),
    pl.lit(None).alias('Family_score'),
    pl.col("taxlineage").str.extract(r"_([A-Za-z]+virinae);?", 1).alias('Subfamily (-virinae)'),
    pl.lit(None).alias('Subfamily_score'),
    pl.col("taxlineage").str.extract(r"_([A-Za-z]+virus);?", 1).alias('Genus (-virus)'),
    pl.lit(None).alias('Genus_score'),
    pl.col("taxlineage").str.extract(r"_[A-Za-z]+virus;-_([A-Za-z]+virus);?", 1).alias('Subgenus (-virus)'),
    pl.lit(None).alias('Subgenus_score'),
    pl.col("taxlineage").str.extract(r"_([A-Za-z]+\s[A-Za-z0-9]+);?", 1).alias('Species (binomial)'),
    pl.lit(None).alias('Species_score'),]

nuc = pl.read_csv(NUC, separator="\t")
prot = pl.read_csv(PROT, separator="\t")
prof = pl.read_csv(PROF, separator="\t")

nuc = nuc.filter((pl.col('qcov') > 0.7) & (pl.col('ani') > 0.7)).group_by('query').agg(
                        pl.all().top_k_by("tani", 1),
                        ).explode(pl.all().exclude("query")).rename({"query":"seqid", "qlen":"qseqlen"})

nuc = nuc.with_columns(pl.lit("ani").alias("Method")).rename({"tani": "Score", "qseqlen": "Seqlen", "seqid": "SequenceID"})
prot = prot.with_columns(pl.lit("aai").alias("Method")).rename({"taai": "Score", "qseqlen": "Seqlen", "seqid": "SequenceID"})
prof = prof.with_columns(pl.lit("api").alias("Method")).rename({"tapi": "Score", "qseqlen": "Seqlen", "seqid": "SequenceID"})

matched = nuc["SequenceID"].to_list()
prot = prot.filter(~pl.col("SequenceID").is_in(matched))
matched.extend(prot["SequenceID"].to_list())
prof = prof.filter(~pl.col("SequenceID").is_in(matched))
matched.extend(prof["SequenceID"].to_list())

# write a ictv taxonomy challange formatted file - this will be removed later
pl.concat([nuc.with_columns(*n).select(keys),
           prot.with_columns(*n).select(keys),
           prof.with_columns(*n).select(keys)]).write_csv(OUTFILE.rstrip(".tsv")+"_ictv.csv", separator=",")

# write a file with more information
pl.concat([nuc.with_columns(*n).select(keys_full),
           prot.with_columns(*n).select(keys_full),
           prof.with_columns(*n).select(keys_full)]).write_csv(OUTFILE, separator="\t")