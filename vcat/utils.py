import polars as pl
import taxopy
import polars.selectors as cs
import numpy as np
from typing import List, Union
import numpy as np
from tqdm import tqdm
import mmap
import pickle
from collections import defaultdict
from io import StringIO  

# pandas implementation: readable but slow
# def compute_ani(alns):
#     # added pandas support for checkv ani script
#     return round(sum(alns["alnlen"] * alns["fident"]) / sum(alns["alnlen"]), 2)
#     # return round(
#     #     sum(a["len"] * a["pid"] for a in alns) / sum(a["len"] for a in alns), 2
#     # )

# def merge_intervals(intervals):
#     starts = intervals[:,0]
#     ends = np.maximum.accumulate(intervals[:,1])
#     valid = np.zeros(len(intervals) + 1, dtype=bool)
#     valid[0] = True
#     valid[-1] = True
#     valid[1:-1] = starts[1:] >= ends[:-1]
#     return np.vstack((starts[:][valid[:-1]], ends[:][valid[1:]])).T

# def compute_cov(alns):   
#     t = alns[["qstart", "qend"]].to_numpy()
#     # sort both axis of the interval tree and merge overlaps
#     merged = merge_intervals(np.sort(np.sort(t, axis=0), axis=1))
#     # calculate query coverage
#     return np.round(np.diff(merged).sum()/alns["qlen"].iloc[0], 2)


# def ani_summary(infile, outfile, header):
#     mmseqs_nuc = pd.read_table(infile, names=header)
#     output = []
#     groups = mmseqs_nuc.groupby(["query", "target"])
#     for g in tqdm(groups.groups):
#         g_ = groups.get_group(g)
#         output.append(
#             {"query":g[0],
#             "target": g[1],
#             "coverage":compute_cov(g_), 
#             "ani":compute_ani(g_),
#             "lineage": g_["taxlineage"].iloc[0],
#             "qlen": g_["qlen"].iloc[0],
#             "tlen": g_["tlen"].iloc[0]
#             }
#             )
#     ani = pd.DataFrame(output)
#     ani["score"] = ani["ani"] * ani["coverage"]
#     ani.iloc[ani.groupby("query")["score"].idxmax().values].to_csv(outfile, index=None, sep="\t")



# polars implementation: fast!

def merge_intervals(intervals: np.ndarray) -> np.ndarray:
    """Merge overlapping intervals using numpy."""
    starts = intervals[:, 0]
    ends = np.maximum.accumulate(intervals[:, 1])
    valid = np.zeros(len(intervals) + 1, dtype=bool)
    valid[0], valid[-1] = True, True
    valid[1:-1] = starts[1:] >= ends[:-1]
    return np.vstack((starts[valid[:-1]], ends[valid[1:]])).T


def compute_cov(alns: List[pl.Series]) -> float:
    """Compute query coverage using polars and numpy"""

    t = np.stack((alns[0].to_numpy(), alns[1].to_numpy()), axis=1)
    merged = merge_intervals(np.sort(np.sort(t, axis=0), axis=1))
    return np.round(np.sum(np.diff(merged)) / alns[2][0], 3)


# ** Process the grouped data **
def ani_summary(infile: str, outfile: str, header: list) -> Union[int, Exception]:
    try:
        mmseqs_nuc = pl.read_csv(infile,
                            has_header = False,
                            separator="\t",
                            new_columns=header)
        output = (
            mmseqs_nuc
            .group_by(["query", "target"])
            .agg([
                pl.col("qlen").first(),
                pl.col("tlen").first(),
                pl.col("taxlineage").first(),
                # to do: re-implement in native api
                pl.map_groups(exprs=["qstart", "qend", "qlen"],function=lambda df: compute_cov(df), return_dtype=pl.Float64).alias("qcov"),
                # computes ani 
                ((pl.col("alnlen") * pl.col("fident")).sum()/pl.col("alnlen").sum()).round(3).alias("ani"),
                

            ])
        )

        # Compute the final column 'tani' and find the max value per query
        output = output.with_columns(
            (pl.col("ani") * pl.col("qcov")).round(3).alias("tani")
        )

        best_hits = output.filter(
            pl.col("tani") == output.group_by("query").agg(pl.max("tani")).join(output, on="query")["tani"]
        )
        best_hits.write_csv(outfile,  separator="\t")
        return 0
    except Exception as e:
        return e


# aai calculation code

def trim_lineage(taxid: int, taxdb: taxopy.core.TaxDb) -> int:
    # trims lineages to genus level
    if taxdb.taxid2rank[taxid] == "species":
        return taxdb.taxid2parent[taxid]
    return taxid

def get_taxid2taxon_map(df: pl.DataFrame, taxdb: taxopy.core.TaxDb) -> dict[int,taxopy.Taxon]:
    tmap = dict()
    for x in list(df["taxid"].unique()):
        tmap[x] =  taxopy.Taxon(x, taxdb=taxdb)
    return tmap

def aai_summary(input: str, gff: str, dbdir:str, outfile:str, header:list) -> Union[int, Exception]:
    # ps: optimized for speed not for readability
    # to: to run a second round of aai calculations at the next taxonomic rank
    # or to calculate aai at every taxonomic rank until kingdom
    TAAI_CUTOFF = 0.3
    taxdb = taxopy.TaxDb(nodes_dmp=f"{dbdir}/ictv-taxdump/nodes.dmp",
                    names_dmp=f"{dbdir}/ictv-taxdump/names.dmp",
                    merged_dmp=f"{dbdir}/ictv-taxdump/merged.dmp")
    # extract query lengths (genomic seqs)
    leninf = []
    with open(gff) as fh:
        for i in fh:
            if i.startswith("# Sequence Data:"):
                for j in i.split()[-1].split(";"):
                    if j.startswith("seqlen="):
                        s = j.lstrip("seqlen=")
                    elif j.startswith("seqhdr="):
                        h =j.lstrip("seqhdr=").strip('"')
                leninf.append({"seqid": h, "qseqlen" : int(s)})

    seqleninfo = pl.DataFrame(leninf)

    col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    # extract ORF stats
    gff_df = pl.read_csv(gff, 
                        separator="\t", 
                        comment_prefix="#", 
                        has_header=False,
                        new_columns=col_names, 
                        schema_overrides={"start": pl.Int64, "end": pl.Int64, "score": pl.Float64})
    

    gff_df = gff_df.with_columns((pl.col("seqid")+ pl.col("attributes")\
                                  .map_elements(lambda x : f'_{x.split(";")[0].strip("ID=").split("_")[-1]}', 
                                                 return_dtype=str)).alias("query"))
    gff_df = gff_df.join(seqleninfo, on="seqid", how="inner")

    mmseqs_nuc = pl.read_csv(input,
                            has_header = False,
                            separator="\t",
                            new_columns=header)
    mmseqs_nuc = gff_df.select(~cs.by_name(["source", "type", "phase", "attributes"])).join(mmseqs_nuc, on="query", how="full")
    mmseqs_nuc = mmseqs_nuc.with_columns(pl.col("fident").fill_null(0))
    mmseqs_nuc = mmseqs_nuc.with_columns(pl.col("taxid").fill_null(1))
    mmseqs_nuc = mmseqs_nuc.with_columns(pl.col("alnlen").fill_null(0))

    mmseqs_nuc = mmseqs_nuc.with_columns(pl.col("taxid").map_elements(lambda x : trim_lineage(x, taxdb=taxdb),
                                                           return_dtype=int))
    taxonmap = get_taxid2taxon_map(mmseqs_nuc, taxdb=taxdb)

    mmseqs_nuc = mmseqs_nuc.with_columns(pl.col("taxid").map_elements(lambda x: str(taxonmap.get(x, 'no_rank')),
                                                           return_dtype=str).alias("taxlineage"))
    mmseqs_nuc = mmseqs_nuc.with_columns(pl.col("taxid").map_elements(lambda x: taxonmap.get(x, 'no_rank').name,
                                                            return_dtype=str).alias("taxname"))
    
    seqid2qlens = mmseqs_nuc.group_by("seqid")\
                .agg(
                    pl.col("qlen").unique().alias("qlens")
                    )
    mmseqs_nuc = mmseqs_nuc.join(seqid2qlens, on="seqid", how="left")

    try: 
            mmseqs_nuc.filter(pl.col("taxid") > 1).\
                group_by(["seqid","taxid"]).agg(
                                                pl.first("taxlineage"),
                                                pl.first("taxname"),
                                                pl.first("qseqlen"),
                                                ((pl.col("fident") * pl.col("alnlen")).sum()/ pl.col("alnlen").sum()).round(3).alias("aai"),
                                                (pl.col("qlen").unique().len()/ pl.first("qlens").list.len()).round(3).alias("qcov"),
                                                pl.col("qlen").unique().len().alias("mgenes"),
                                                (pl.first("qlens").list.len()).round(3).alias("qgenes"))\
                                            .with_columns(
                                                (pl.col("aai") * pl.col("qcov")).round(3).alias("taai"),
                                                    )\
                                            .group_by("seqid")\
                                        .agg(
                                                pl.all().sort_by('taai').last(),
                                                ) \
                                            .sort('taai').filter(pl.col('taai') > TAAI_CUTOFF)\
                                            .write_csv(outfile,  separator="\t")
            return 0
    except Exception as e:
        return e

def index_m8(input: str) -> defaultdict[list]:
    """
    m8 files with blast results can be insanely large and may require a large amount of
    memory to load.
    to over come this, we index the m8 file and load it chunk wise.
    This approch is slow but effective :)
    """
    index = defaultdict(list)

    with open(input, "r+b") as f:
        mmapped_file = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        total_size = len(mmapped_file)

        pbar = tqdm(unit="lines")
        
        pos = 0
        while pos < total_size:
            end = mmapped_file.find(b"\n", pos)
            if end == -1:  # End of file
                break

            line = mmapped_file[pos:end]  # Process as bytes
            key = line.split(b"\t", 1)[0].rsplit(b"_", 1)[0]  # Extract key as bytes
            index[key.decode()].append(int(pos))  # Store file offset

            pos = end + 1  # Move to the next line
            pbar.update(1)

    pbar.close()

    pbar = tqdm(unit="records")
    index2 = defaultdict(list)
    for k, (_, byte_offsets) in enumerate(index.items()):
        pbar.update(1)
        index2[k] = byte_offsets

    return index2

def load_chunk(input: str,
               index: defaultdict[list], 
               recstart: int,
               recend: int) -> StringIO:
    string_dump = ""
    with open(input, "r+b") as fh:
        mmapped_file = mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ)
        
        records = list(range(recstart, recend))
            # Retrieve lines at given offsets
        for rec in records:
            for pos in index[rec]:
                end = mmapped_file.find(b"\n", pos)  # Find the end of the line
                if end == -1:  # Handle last line case
                    end = len(mmapped_file)
                
                string_dump += mmapped_file[pos:end].decode()+"\n"

    return StringIO(string_dump)
                