#!/usr/bin/env python
import pandas as pd
import numpy as np
import taxopy
import pyfastx
from sys import argv, stderr

DATABASE_DIR = argv[1]
INFILE = argv[2]
OUTFILE = argv[3]



taxdb = taxopy.TaxDb(nodes_dmp=f"{DATABASE_DIR}/ictv-taxdump/nodes.dmp",
                     names_dmp=f"{DATABASE_DIR}/ictv-taxdump/names.dmp",
                     merged_dmp=f"{DATABASE_DIR}/ictv-taxdump/merged.dmp")

results = pd.read_csv(INFILE)

rows = []
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
        
for row in results.iterrows():
    t = dict.fromkeys(keys)
    t['SequenceID'] = row[-1]['query']
    i =  taxopy.Taxon(row[-1]['taxid'], taxdb=taxdb).rank_name_dictionary
    for k,v in i.items():
        for c in keys:
            if c.startswith(f"{k.capitalize()} "):
                t[c.capitalize()] = v
    rows.append(t)

pd.DataFrame(rows).to_csv(OUTFILE, index=None)