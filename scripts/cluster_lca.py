#!/usr/bin/env python
import pandas as pd
import numpy as np
import taxopy
import os
from sys import argv


# You can also use your own set of taxonomy files:
DATABASE_DIR = argv[1]

taxdb = taxopy.TaxDb(nodes_dmp=f"{DATABASE_DIR}/ictv-taxdump/nodes.dmp",
                     names_dmp=f"{DATABASE_DIR}/ictv-taxdump/names.dmp",
                     merged_dmp=f"{DATABASE_DIR}/ictv-taxdump/merged.dmp")

acc2tx = f'{DATABASE_DIR}/VMR_latest/virus_protein.accession2taxid'
clust = f'{DATABASE_DIR}/VMR_latest/mmseqs_pclusters/mmseqs_pclusters.tsv'

acc2tx = pd.read_table(acc2tx)
acc2tx = acc2tx.drop('gi', axis=1)

clusters = pd.read_table(clust, names=["cluster_rep", "cluster_mem"])

def lca(taxa,taxdb=taxdb, fraction=0.6):
    
    t = taxa.to_list()
    if len(t) > 1:
        return taxopy.find_majority_vote(t, taxdb=taxdb, fraction=fraction).taxid
    return taxa.iloc[0].taxid

clusters = pd.merge(left=clusters, right=acc2tx, left_on='cluster_mem', right_on='accession.version')
clusters["lineage"] = clusters.apply(lambda x : taxopy.Taxon(x['taxid'], taxdb=taxdb), axis=1)
clusters['root_p1'] =  clusters.apply(lambda x : x["lineage"].taxid_lineage[-2], axis=1)

clustlca = clusters.groupby(['cluster_rep', 'root_p1']).apply(lambda x : lca(x['lineage']), include_groups=False)
clustlca = pd.DataFrame(clustlca).reset_index()
clustlca.columns = ["cluster_rep","root_p1", "taxid"]

clustlca.to_csv(f'{DATABASE_DIR}/VMR_latest/mmseqs_pprofiles/mmseqs_pprofiles_lca.tsv', index=None, sep="\t")