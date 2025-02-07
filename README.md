# vcat

The Virus Contig Annotation Tool (vcat) is a straightforward, homology-based application designed to provide taxonomic annotations for viral contigs.

#### Database Generation Pipeline
The pipeline automates the creation of comprehensive viral nucleotide, protein, and profile databases. The database generation process includes the following steps:

1. <b>ICTV Metadata Parsing and Genome Download:</b>
VCAT parses the ICTV Metadata Resource and downloads all viral genomes and proteins from NCBI GenBank.

2. <b>Custom Taxonomy Creation:</b>
A custom taxdump for the ICTV taxonomy is generated using taxonkit.

3. <b>Nucleotide and Protein Database Construction:</b>
Custom nucleotide and protein databases are created using mmseqs2, incorporating the ICTV-specific taxdump.

4. <b>Protein Clustering and Profile Database Creation:</b>
Viral proteins are clustered, and these clusters are converted into a protein profile database.

#### Contig Annotation Pipeline
`vcat` annotates viral contigs by sequentially comparing query sequences against the nucleotide, protein, and profile databases uisng `mmseqs2`. The annotation process follows these steps:

1. <b>Nucleotide Database Comparison:</b>
Query sequences are first compared to the nucleotide database. If matches with >95% identity and >85% query coverage are found, the lineage of the target sequence is assigned to the query sequence.

2. <b>Protein Database Comparison:</b>
For query sequences with no significant matches at the nucleotide level, open reading frames (ORFs) within the query sequence are extracted and compared to the viral protein database.

3. <b>Profile Database Comparison:</b>
Sequences that remain unannotated at the protein level are then compared to the protein profile database to provide taxonomic annotations.

#### how to install?
```
git clone https://github.com/Yasas1994/vcat.git
cd vcat

# create a conda env
mamba create env -f environment.yml

# install vcat pipeline
pip install .

# test the installation
vcat -h
```

#### how to run?
```
# download and prepare the databases
vcat preparedb -d databases

# run vcat
vcat annotate -i <input>.fasta -o <output>
```


##### Expected runtime ?

It takes ~4hrs to run vcat on the ICTV Taxonomy challenge dataset on a laptop computer.