
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

1. **ICTV Metadata Parsing and Genome Download:**
VCAT parses the ICTV Metadata Resource and downloads all viral genomes and proteins from NCBI GenBank.

2. **Custom Taxonomy Creation:**
A custom taxdump for the ICTV taxonomy is generated using taxonkit.

3. **Nucleotide and Protein Database Construction:**
Custom nucleotide and protein databases are created using mmseqs2, incorporating the ICTV-specific taxdump.

4. **Protein Clustering and Profile Database Creation:**
Viral proteins are clustered, and these clusters are converted into a protein profile database.

---

#### Contig Annotation Pipeline

`vcat` annotates viral contigs by sequentially comparing query sequences against the nucleotide, protein, and profile databases uisng `mmseqs2`. The annotation process follows these steps:

1. **Nucleotide Database Comparison:**
Query sequences are first compared to the nucleotide database. If matches with >95% identity and >85% query coverage are found, the lineage of the target sequence is assigned to the query sequence.

- Here, pseudo ANI is defined as follows. Let, $C$ be a unknown virus contig that needs to be taxonomically annotated. Let, $G_{db}$ be a database of known virus genomes, in our case: ICTV genomes. We search $C$ agaist $G_{db}$ to retrieve all matches for $C$ in $G_{db}$. Let, $M_{db}$ be the list of genomic regions $C$ map to. $M_{db}$ may contain hits from several taxonoically related genomes and $M_{{db}_t}$ represents all regions matching to a taxonomic rank. We calculate the total length of alignmnets to each taxonomic rank and the similarity of these alignments to $C$ using,

$$pseudo\text{ }ANI = \frac{\sum M_{{db}_t} \text{ }identidy}{ \sum M_{{db}_t} \text{ }alignment\text{ }length}$$

2. **Protein Database Comparison:**
For query sequences with no significant matches at the nucleotide level, open reading frames (ORFs) within the query sequence are extracted and compared to the viral protein database.

- Here, pseudo AAI is defined as follows. extending the definitions from the above section, Let $P_{c}$ be a list of genes on $C$. $P_{db}$ is the protein database of ICTV viruses. $P_{db}$ is searched for matches for ORFs found on $C$ i.e. $C_{ORF}$. $P_m \in P_c$ represents all ORFs with a match. $I_m$ is the aminnno acid identity of proteins in $P_m$. 

$$ pseudo\text{ }AAI = \frac{1}{n}\sum_{i}^{n}I_{m_{i}}$$


3. **Profile Database Comparison:**
Sequences that remain unannotated at the protein level are then compared to the protein profile database to provide taxonomic annotations.

- Here, pseudo API calculation is simlar to the AAI calculation above. The main difference is the target database is a protein profile database instead of protein database.

---

#### how to install?

```
git clone https://github.com/Yasas1994/vcat.git
cd vcat

# create a conda env
 mamba env create -f environment.yml

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

##### Some additional stuff :)

```
# to view the contig length distribution of your contigs
seqkit fx2tab -lg {input.fasta} | awk -F "\t" '{print $4}' | tail -n +2 | hist -b 100 -s 10

# to view the length distribution of contigs in the ani calculation output (after apply a filter)
cat testing/gut_jaeger/nuc/gut_jaeger_virus_seqs_fasta_genome_ani.tsv | awk -F "\t" '$8 > 0.1' | awk -F "\t" '{print $3}' | tail -n +2 | hist -b 100 -x

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
