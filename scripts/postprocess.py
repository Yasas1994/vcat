#!/usr/bin/env python
import pandas as pd
import numpy as np
import taxopy
import pyfastx
from sys import argv, stderr

DATABASE_DIR = argv[1]
NUC = argv[2]
PROT = argv[3]
PROF = argv[4]
FASTA = argv[5]
OUTFILE = argv[6]

fasta = pyfastx.Fasta(FASTA)
headers = list(fasta.keys())

taxdb = taxopy.TaxDb(nodes_dmp=f"{DATABASE_DIR}/ictv-taxdump/nodes.dmp",
                     names_dmp=f"{DATABASE_DIR}/ictv-taxdump/names.dmp",
                     merged_dmp=f"{DATABASE_DIR}/ictv-taxdump/merged.dmp")

clustlca = pd.read_table(f'{DATABASE_DIR}/VMR_latest/mmseqs_pprofiles/mmseqs_pprofiles_lca.tsv')
clustlca['lca'] = clustlca.apply(lambda x: taxopy.Taxon(x['taxid'], taxdb), axis=1)

def lca(taxa, taxdb=taxdb, fraction=0.6):
    
    t = taxa.to_list()
    if len(t) > 1:
        return taxopy.find_majority_vote(t, taxdb=taxdb, fraction=fraction)
    return taxa.iloc[0]

def trim_lineage(row, rank='class', taxdb=taxdb):
    taxon = taxopy.Taxon(row['taxid'], taxdb=taxdb)
    # add support for subranks in the future
    # if sub-rank is unavailable round off the to nearest whole rank
    ranks = ['realm','subrealm', 'kingdom', 'subkingdom','phylum','subphylum', 'class', 'subclass', 
             'order','suborder', 'family','subfamily', 'genus','subgenus','species', 'subspecies']
    # case one
    if ranks.index(taxon.rank) <= ranks.index(rank):
        return taxon, taxon.taxid, taxon.rank
    else:
        parent = taxon.parent(taxdb)
        daughter = taxon
        while parent.rank != rank and parent.rank != 'no rank' and ranks.index(parent.rank) >= ranks.index(rank):
            daughter = parent
            parent = parent.parent(taxdb)
        if parent.rank != 'no rank':
            return parent, parent.taxid, parent.rank
        else:
            return daughter, daughter.taxid, daughter.rank
    

## mmseqs nucleotide search against ictv genomes
ictv_nuc=pd.read_table(NUC,
                   names='query,target,fident,qlen,tlen,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,taxid,taxname,taxlineage'.split(","))

ictv_nuc['qcov'] = ictv_nuc['alnlen']/ictv_nuc['qlen']

# filter nuc results
idx = ictv_nuc.query('qcov >= 0.85  & fident >= 0.95').groupby("query")['alnlen'].idxmax()
ictv_nuc = ictv_nuc.iloc[idx.values]

# searches against ICTV proteindb using mmseqs taxonomy
ictv_prot=pd.read_table(PROT,
                   names='query,taxid,trank,tname,nfrag,nlfrag,nsfrag,sfract,lineage'.split(","))

ictv_prot = ictv_prot.query('tname != "unclassified"')
ictv_prot = ictv_prot.query('taxid > 1').dropna()


## mmseqs pprofile search against ictv genomes
ictv_prof=pd.read_table(PROF,
                   names='query,target,fident,pident, qlen,tlen,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,taxid,taxname,taxlineage'.split(","))

ictv_prof = pd.merge(left=clustlca, right=ictv_prof, right_on="target", left_on="cluster_rep")
ictv_prof = pd.DataFrame(ictv_prof.groupby("query").apply(lambda x : lca(x['lca']), include_groups=False )).reset_index()
ictv_prof.columns = ["query", "lca"]
ictv_prof['taxid']= ictv_prof.apply(lambda x : x['lca'].taxid, axis=1)
ictv_prof = ictv_prof[ictv_prof["lca"].apply(lambda x : len(x.taxid_lineage)) != 1]

found_in_nuc = set(ictv_nuc["query"])
missing_in_nuc = set(headers) - found_in_nuc
found_in_prot = missing_in_nuc.intersection(set(ictv_prot["query"]))
missing_in_prot = set(missing_in_nuc) - set(found_in_prot)
# filter prot
ictv_prot = ictv_prot[ictv_prot['query'].isin(found_in_prot)]

found_in_prof = missing_in_prot.intersection(set(ictv_prof["query"]))
#filter prof
ictv_prof = ictv_prof[ictv_prof['query'].isin(found_in_prof)]

ictv_nuc[['lineage', 'taxid', 'trank']] = ictv_nuc.apply(lambda x : trim_lineage(x, rank="species"), axis=1, result_type="expand")
ictv_prot[['lineage', 'taxid', 'trank']] = ictv_prot.apply(lambda x : trim_lineage(x, rank="genus"), axis=1, result_type="expand")
ictv_prof[['lineage', 'taxid', 'trank']] = ictv_prof.apply(lambda x : trim_lineage(x, rank="family"), axis=1, result_type="expand")

pd.concat([ictv_nuc[['query','taxid','trank','lineage']],
           ictv_prot[['query','taxid','trank','lineage']],
           ictv_prof[['query','taxid','trank','lineage']]]).to_csv(OUTFILE, index=None)

