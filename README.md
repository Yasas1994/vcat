

```

                XXXX       XXXX
              XXXXXXXXXXXXXXXXXXX
             X    XXXXXXXXXXX    X
           XX      XXXXXXXXX      XX
         XX         XXXXXXX          X                           _
        XXXXXXX     XXXXXXX         XXX          __   _____ __ _| |_
        XXXXXXXXXX   XXXXX     XXXXXXXX          \ \ / / __/ _` | __|
         XXXXXXXXXXXXXXXXXXXXXXXXXXXXX            \ V / (_| (_| | |_
         XXXXXXXXXX76636174XXXXXXXXXXX             \_/ \___\__,_|\__|
         XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        XXXXXXXX     XXXXX   XXXXXXXXX  
        XXX         XXXXXXX     XXXXXXX
         X          XXXXXXX         XX
           XX      XXXXXXXXX     X X
             X    XXXXXXXXXXX    X
              XXXXXXXXXXXXXXXXXXX
                XXXX       XXXX

```

The Virus Contig Annotation Tool (vcat) is a straightforward, homology-based application designed to provide taxonomic annotations for viral contigs, mapping reads to virus contigs and much more.

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

---
#### Contig Annotation Pipeline
`vcat` annotates viral contigs by sequentially comparing query sequences against the nucleotide, protein, and profile databases uisng `mmseqs2`. The annotation process follows these steps:

1. <b>Nucleotide Database Comparison:</b>
Query sequences are first compared to the nucleotide database. If matches with >95% identity and >85% query coverage are found, the lineage of the target sequence is assigned to the query sequence.

- Here, ANI is defined as follows. Let, $C$ be a unknown virus contig that needs to be taxonomically annotated. Let, $G_{db}$ be a database of known virus genomes, in out case: ICTV genomes. We search $C$ agaist $G_{db}$ to retrieve all matches from $C$ in $G_{db}$. Let, $M_{db}$ be the list of genomic regions $C$ map to. $M_{db}$ may contain hits from several taxonoically related genomes. We calculate the total length of alignmnets to each taxonomic rank and the similarity of these alignments to $C$ using,

$$\color{white}{X \sim Normal \; (\mu,\sigma^2)}$$

2. <b>Protein Database Comparison:</b>
For query sequences with no significant matches at the nucleotide level, open reading frames (ORFs) within the query sequence are extracted and compared to the viral protein database.

- Here, AAI is defined as follows. extending the definitions from the above section, Let $P_{c}$ be a list of genes on $C$. 


3. <b>Profile Database Comparison:</b>
Sequences that remain unannotated at the protein level are then compared to the protein profile database to provide taxonomic annotations.
---
#### how to install?
```
git clone https://github.com/Yasas1994/vcat.git
cd vcat

# create a conda env
mamba create env -f environment.yml

# install vcat pipeline
pip install .

# test the installation
vcat --help

Usage: vcat [OPTIONS] COMMAND [ARGS]...

  vcat: a command-line tool-kit for adding ICTV taxonomy annotations to virus
  contigs, mapping reads to virus genomes and much more.
  (https://github.com/Yasas1994/vcat)

Options:
  --version   Show the version and exit.
  -h, --help  Show this message and exit.

Commands:
  contigs    run contig annotation workflow
  preparedb  download and build reference databases
  reads      run read annotation workflow
  utils      tool chain for calculating ani, aai and visualizations

```

---
#### download and prepare the databases
---
```
vcat preparedb -d databases
```

---
#### running vcat
---
```

# run vcat contig annotation pipeline
vcat contigs -i <input>.fasta -o outdir

# run vcat read annotation pipeline (comming soon)
vcat reads -i1 <reads1>.fastq [-i2 <reads2.fastq>] -o outdir

# calculate aai of query contigs to ICTV genomes
vcat utils aai  [OPTIONS] -i contigs.m8 -g configs.gff -d [DBDIR]

# calculate ani of query contigs to ICTV genomes
vcat utils ani [OPTIONS] -i contigs.m8

# creating genome comparision plots. i.e query sequence to highly similar ICTV genomes (comming soon)
vcat utils visualize --ani --taxa [taxname] -i contigs.m8 -o outdir

# create phage contig annotation plots (coming soon)
vcat utils visualize --phrogs -i contigs.fasta -o outdir 

# identify provirus (coming soon)
vcat utils provirus -i contigs.fasta -o outdir
```

##### Expected runtime ?

It takes ~4hrs to run vcat on the ICTV Taxonomy challenge dataset on a laptop computer.

---
If you use vcat please cite,
```

[MMSEQS2] MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets.
          Steinegger M & Söding J. 2017. Nat. Biotech. 35, 1026–1028. https://doi.org/10.1038/nbt.3988

[PRODIGAL] Prodigal: prokaryotic gene recognition and translation initiation site identification.
           Hyatt et al. 2010. BMC Bioinformatics 11, 119. https://doi.org/10.1186/1471-2105-11-119.

[TAXONKIT] TaxonKit: A practical and efficient NCBI taxonomy toolkit.
           Shen, W. & Ren, H. J. 2021. Genet. Genomics https://doi:10.1016/j.jgg.2021.03.006.

```