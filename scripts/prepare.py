#!/usr/bin/env python
"""
author: Yasas Wijesekara (yasas.wijesekara@uni-greifswald.de)

prepares the genome and protein sequence files for mmseqs
databse creation steps
"""
from ete3 import NCBITaxa
import pandas as pd
import numpy as np
from Bio import Entrez
from Bio import GenBank
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import openpyxl
from sys import argv
from tqdm import tqdm
from glob import glob
import taxopy
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


# Specify the path to the custom taxdump database
DATABASE_DIR=argv[1]  # SQLite database path
XLTABLE=argv[2] 
GBDIR=argv[3] 
SEQDIR=argv[4] 


def extract_regions(string):
    r = [list(map(int, i)) for i in re.findall(r"\((\d+)\.{1,3}(\d+)\)", string)]
    return r

def filter_genbank(ids):
    ids_ = []
    for i in ids:
        if '_' in i:
            ids_.append(i)
    if ids_:
        return ids_
    else:
        return ids


def extract_from_hyperlink(hypstr):
    url_match = re.search(r'=HYPERLINK\("([^"]+)"', hypstr)
    url = url_match.group(1) if url_match else None
    text_match = re.search(r',"([^"]+)"\)', hypstr)
    display_text = text_match.group(1) if text_match else None
    return {"url": f'{url}', "display_text": f'{display_text}'}


taxdb = taxopy.TaxDb(nodes_dmp=f"{DATABASE_DIR}/ictv-taxdump/nodes.dmp",
                     names_dmp=f"{DATABASE_DIR}/ictv-taxdump/names.dmp",
                     merged_dmp=f"{DATABASE_DIR}/ictv-taxdump/merged.dmp")

wb = openpyxl.load_workbook(XLTABLE)
sheets = wb.sheetnames
ws = wb[sheets[1]]
headers = [cell.value for cell in ws[1]]
# Convert the worksheet to a Pandas DataFrame
data = ws.values  # Extract cell values
columns = next(data)  # Extract the first row for column headers
ictv = pd.DataFrame(data, columns=columns)

# =HYPERLINK("https://ictv.global/taxonomy/taxondetails?taxnode_id=202308643","Alphalipothrixvirus beppuense")
# =HYPERLINK("https://www.ncbi.nlm.nih.gov/nuccore/MH447526","NCBI Nucleotide")

ictv = ictv[ictv["Accessions Link"].notna()]
#ictv["Species"] = ictv["Species"].apply(lambda x: extract_from_hyperlink(x)["display_text"])
ictv["Taxid"] = ictv["Species"].apply(lambda x : taxopy.taxid_from_name(x, taxdb)[0])
ictv["IDS"] = ictv.apply(lambda x: extract_from_hyperlink(x["Accessions Link"])["url"].split("/")[-1].split(","), axis=1)
ictv["Range"] = ictv.apply(lambda x: extract_regions(x["Virus GENBANK accession"] ), axis=1)




id2newtax = {}
id2range = {}
for row in ictv[['Taxid','Virus name(s)',  'IDS', 'Range']].iterrows(): # ID = genome ID
    for i in row[1]["IDS"]:
        i = i.split("(")[0]
        id2newtax[i] = row[1]['Taxid']
        if row[1]['Range']:
            if len(row[1]['Range']) > 0:
                id2range[i] = row[1]['Range']
logger.info(id2range)
ictv["IDS_NR"] = ictv["IDS"].apply(lambda x : filter_genbank(x)) 
all_nrids =[]
for ii in ictv.iterrows():
    i = ii[1]
    if len(i["IDS_NR"])>0 and i["IDS_NR"][0] != "":
        all_nrids.extend(i["IDS_NR"])
    else:
        pass

all_nrids = set(all_nrids)
all_files = glob(f"{GBDIR}/*.gb")

all_seqs = 0
for file in all_files:    
    all= open(file, 'r')   
    all_=all.read()  
    all.close()
    with open(f"{GBDIR}/tmp.txt",'w') as fh:
        for i  in re.findall(r'LOCUS.+\n', all_):
            fh.write(i)
    all_seqs += len(re.findall(r'LOCUS\s+(\S+)\s+', all_))

proteins = []
for file in tqdm(all_files):
    with open(file) as handle:
        for record in GenBank.parse(handle):
            if record.version.split(".")[0] in all_nrids:
                if record.version.split(".")[0] in id2range:
                    
                    feature_range = id2range[record.version.split(".")[0]][0]
                    #print(record.version.split(".")[0], feature_range)
                else:
                    feature_range = 0, np.inf

                for feature in record.features:    

                    if feature.key == "CDS" :
                        # Get the protein sequence 
                        # filter out host protein sequences and only obtain accessions of phage proteins
                        for qualifier in feature.qualifiers:
                            if qualifier.key =='/protein_id=':
                                protein_id = qualifier.value.strip("\"").rstrip("\"")
                            if qualifier.key =='/product=':
                                product = '['+qualifier.value.strip("\"").rstrip("\"")+']'
                            if qualifier.key =='/translation=':
                                protein_seq = qualifier.value.strip("\"").rstrip("\"")
                        # Create a SeqRecord for the protein
                        #print(feature.location)
                        start, end = list(map(int, re.findall(r"(\d+)\D+(\d+)", feature.location)[0]))
                            # get start and end positions of the feature
                        #print(feature.location, start, end)    
                        # if the feature lies within prophage co-ordinates 
                        if min(feature_range[0], start) == feature_range[0] and max(feature_range[1], end) == feature_range[1]: 

                            protein_record = SeqRecord(
                                Seq(protein_seq),
                                id=protein_id,
                                description=product,
                                name=product
                            )

                            proteins.append(protein_record)

# Write the protein sequences to a file
# YP_011036812.1 |RNA-dependent RNA polymerase [Frijoles virus]
with open(f'{SEQDIR}/proteins.faa', "w") as output_handle:
    SeqIO.write(proteins, output_handle, "fasta")

logger.info(f'{SEQDIR}/proteins.faa created')

# Write accession2taxid file
#accession       accession.version       taxid   gi
with open(f'{SEQDIR}/virus_protein.accession2taxid','w') as fh2:
    fh2.write(f"accession\taccession.version\ttaxid\tgi\n")
    for file in tqdm(all_files):
        with open(file) as handle:
            for record in GenBank.parse(handle):
                if record.version.split(".")[0] in all_nrids:
                    if record.version.split(".")[0] in id2range:
                        feature_range = id2range[record.version.split(".")[0]][0]
                        #print(record.version.split(".")[0], feature_range)
                    else:
                        feature_range = 0, np.inf
                    # for i in record.features[0].qualifiers:
                        # if 'taxon' in i.value:
                        #     taxid = i.value.split("=")[-1].rstrip("\"").strip("\"").split(":")[-1]
                    for feature in record.features:    
                        if feature.key == "CDS":
                            # Get the protein sequence
                            for qualifier in feature.qualifiers:
                                if qualifier.key =='/protein_id=':
                                    protein_id = qualifier.value.strip("\"").rstrip("\"")
                                if qualifier.key =='/product=':
                                    product = '['+qualifier.value.strip("\"").rstrip("\"")+']'
                                if qualifier.key =='/translation=':
                                    protein_seq = qualifier.value.strip("\"").rstrip("\"")
                                # Create a SeqRecord for the protein
                            # Create a SeqRecord for the protein

                            start, end = list(map(int, re.findall(r"(\d+)\D+(\d+)", feature.location)[0]))
                                # get start and end positions of the feature
                            #print(feature.location, start, end)   
                            # if the feature lies within prophage co-ordinates 
                            if min(feature_range[0], start) == feature_range[0] and max(feature_range[1], end) == feature_range[1]: 
                                fh2.write(f"{protein_id.split('.')[0]}\t{protein_id}\t{id2newtax[record.accession[0]]}\t{'-'}\n")
logger.info(f'{SEQDIR}/virus_protein.accession2taxid created')

seq_records = []
for file in tqdm(all_files):
    with open(file) as handle:
        for record in GenBank.parse(handle):
            if record.version.split(".")[0] in all_nrids:

                if record.version.split(".")[0] in id2range:
                    feature_range = id2range[record.version.split(".")[0]][0]
                    logger.info(f'{record.version.split(".")[0]}, {feature_range}')
                else:
                    feature_range = 0, None
                if record.sequence[feature_range[0]:feature_range[1]] != "":
                    seq_records.append(SeqRecord(
                                Seq(record.sequence[feature_range[0]:feature_range[1]]),
                                id=record.accession[0],
                                description=f'[{record.source}] {record.data_file_division}'

                            ))
# Write the genome sequences to a file
with open(f'{SEQDIR}/genomes.fna', "w") as output_handle:
    SeqIO.write(seq_records, output_handle, "fasta")

logger.info(f'{SEQDIR}/genomes.fna created')

# Write accession2taxid file
#accession       accession.version       taxid   gi
with open(f'{SEQDIR}/virus_genome.accession2taxid','w') as fh2:
    fh2.write(f"accession\taccession.version\ttaxid\tgi\n")
    for file in tqdm(all_files):
        with open(file) as handle:
            
            for record in GenBank.parse(handle):
                if record.version.split(".")[0] in all_nrids:
                    fh2.write(f"{record.accession[0].split('.')[0]}\t{record.accession[0].split('.')[0]}\t{id2newtax[record.accession[0]]}\t{'-'}\n")

logger.info(f'{SEQDIR}/virus_genome.accession2taxid created')