#!/usr/bin/env python

"""
author: Yasas Wijesekara (yasas.wijesekara@uni-greifswald.de)

parses the ICTV datasheet and grabs the genomes of viruses directly
from NCBI in gb format

"""

import openpyxl
from pathlib import Path
from sys import argv
import re
import requests
from vcat.color_logger import logger
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

infile = argv[1]
outdir = argv[2]


# Regular expression pattern
pattern = r"https://www\.ncbi\.nlm\.nih\.gov/nuccore/.+"
DBDIR = outdir
directory = Path(outdir)
# Create the directory if it doesn't exist
directory.mkdir(parents=True, exist_ok=True)


# NCBI EFetch API URL for FASTA format
base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

wb = openpyxl.load_workbook(infile)
sheets = wb.sheetnames
ws = wb[sheets[1]]
headers = [cell.value for cell in ws[1]]

dlink = "Accessions Link"
virus_names = "Virus name(s)"
num_links = 0
all_ids = []
with open(f"{DBDIR}/ID.list", "w") as fh:
    # Iterate through the cells in the last column
    for row in range(2, ws.max_row + 1):  # Iterate through all rows
        names = ws.cell(row=row, column=headers.index(virus_names) + 1)
        tmp = "" if None is names.value else names.value
        if "gene transfer agent" not in tmp:
            cell = ws.cell(row=row, column=headers.index(dlink) + 1)
            if cell.value:
                link = cell.value.split('"')[1]
                match = re.search(pattern, link)
                if match:
                    num_links += 1

                    # Extract sequence IDs from the URL
                    ids = link.split("/")[-1]
                    fh.write(ids + "\n")
                    all_ids.extend(ids.split("(")[0].split(","))


for i in range(0, len(all_ids), 200):
    # API parameters
    params = {
        "db": "nuccore",  # Nucleotide database
        "id": ",".join(all_ids[i : i + 200]),  # Sequence IDs
        "rettype": "gb",  # Return type
        "retmode": "text",  # Return mode
    }

    # Send the request
    response = requests.get(base_url, params=params)

    # Check if the request was successful
    if response.status_code == 200:
        # Save the FASTA data to a file
        with open(f"{DBDIR}/sequences_{i}_{i+200}.gb", "w") as file:
            file.write(response.text)
        logger.info(f"FASTA file downloaded successfully as sequences_{i}_{i+200}.gb")
    else:
        logger.info(
            f"Failed to retrieve data. HTTP Status Code: {response.status_code}"
        )


open(f"{DBDIR}/download_complete", "w").close()
