


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



The virus contig annotation Tool (vcat) is a straightforward, homology-based application designed to provide taxonomic annotations for viral contigs, mapping reads to virus contigs and much more.

### Changelog
---
<b>0.0.2</b>
* added downloaddb commad to download pre-built databases (VMR_40v1 and VMR_39v4)
* added a definition file to build Apptainer containers
* added subroutines to clean up tmp files generated duting the database build

---

### Database Generation Pipeline

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

### Virus Contig Annotation Pipeline

`vcat` annotates viral contigs by sequentially comparing query sequences against the nucleotide, protein, and profile databases uisng `mmseqs2`. The annotation process follows these steps:

1. **Nucleotide Database Comparison:**
Query sequences are first compared to the nucleotide database. If matches with >95% identity and >85% query coverage are found, the lineage of the target sequence is assigned to the query sequence.

- Here, taxo ANI is defined as follows. Let, $C$ be a unknown virus contig that needs to be taxonomically annotated. Let, $G_{db}$ be a database of known virus genomes, in our case: ICTV genomes. We search $C$ agaist $G_{db}$ to retrieve all matches for $C$ in $G_{db}$. Let, $M_{db}$ be the list of genomic regions $C$ map to. $M_{db}$ may contain hits from several taxonoically related genomes and $M_n$ represents all regions matching to a taxonomic rank. We calculate the total length of alignmnets to each taxonomic rank and the similarity of these alignments to $C$ using,

$$taxo\text{ }ANI = \frac{\sum_{i=1}^{M_n} talnlen_i\times fidentity_i}{\sum_i talnlen_i} \text{ where}$$

$$talnlen_i = \text{length of the alignment } i \text{ between } C \text{ and } G_db$$

$$tfident_i = \text{identity of the alignment } i \text{ between } C \text{ and } G_db
$$

2. **Protein Database Comparison:**
For query sequences with no significant matches at the nucleotide level, open reading frames (ORFs) within the query sequence are extracted and compared to the viral protein database.

- Here, taxo AAI is defined as follows. Extending the definitions from the above section, let $P_c$ be a list of genes (ORFs) on $C$. Let $P_{db}$ be a database of viral proteins from ICTV genomes. $P_{db}$ is searched for matches to the ORFs found in $C$, denoted as $C_{ORF}$. Let $P_m \subseteq P_c$ represent the subset of ORFs from $C$ that have a significant match in $P_{db}$. $I_{m_i}$ is the amino acid identity between the $i$-th matched ORF in $P_m$ and its corresponding hit in $P_{db}$. The taxonomic average amino acid identity is calculated as:

$$ taxo\text{ }AAI = \frac{1}{n}\sum_{i}^{n}I_{m_{i}} \text{ where,}$$

$$n = |P_m|, \text{ the number of ORFs with matches}$$

$$I_{m_i} = \text{amino acid identity of matched ORF } i \text{ between } C \text{ and } P_{db}
$$


3. **Profile Database Comparison:** Sequences that remain unannotated at the protein level are compared to a protein profile (HMM) database to provide taxonomic annotations.

- Here, taxo API is calculated similarly to taxo AAI, except the target database is a protein profile (HMM) database instead of a protein sequence database. Let $H_{db}$ be the profile HMM database derived from ICTV viral protein families. Let $H_m$ represent the subset of ORFs from $C$ that match to profiles in $H_{db}$. Let $S_{m_i}$ be the percent sequence identity (or bit-score-normalized identity, depending on implementation) between the $i$-th ORF and its matched profile.

$$
\text{taxo API} = \frac{1}{n} \sum_{i=1}^{n} S_{m_i} \text{ where,}$$

$$n = |H_m|, \text{ the number of ORFs with profile matches }$$

$$S_{m_i} = \text{similarity score between ORF} i \text{ in } C \text{ and its matched profile in } H_{db}
$$

---

### how to install? (conda)

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
  downloaddb download pre-built reference databases
  preparedb  download and build reference databases
  reads      run read annotation workflow
  utils      tool chain for calculating ani, aai and visualizations

```
### Singularity (Now Apptainer)

```
git clone https://github.com/Yasas1994/vcat.git
cd vcat

# build container
apptainer build vcat.sif Apptainer

# test the container build
apptainer run vcat.sif vcat --help

Usage: vcat [OPTIONS] COMMAND [ARGS]...

  vcat: a command-line tool-kit for adding ICTV taxonomy annotations to virus
  contigs, mapping reads to virus genomes and much more.
  (https://github.com/Yasas1994/vcat)

Options:
  --version   Show the version and exit.
  -h, --help  Show this message and exit.

Commands:
  contigs    run contig annotation workflow
  downloaddb download pre-built reference databases
  preparedb  download sequences and build reference databases
  reads      run read annotation workflow
  utils      tool chain for calculating ani, aai and visualizations

```
---
### downloading pre-built databases
---

```
# pulling pre-built databases from remote server [vmr39v and vmr40v1]
vcat downloaddb --dbversion vmr39v4 -d <path to save the database>
```

---

### downloading and preparing the databases

---

```
vcat preparedb -d <path to save the database>
```

---

### running contig annotation pipeline

---

```
vcat contigs -i <input>.fasta -o outdir
```

results can be found in the results directory within the ouput directory
```
.
├── logs
├── nuc
│   ├── input_genome_ani.tsv
│   └── input_genome.m8
├── prof
│   ├── input_fasta_prof_api.tsv
│   └── input_fasta_prof.m8
├── prot
│   ├── input_fasta.faa
│   ├── input_fasta.gff
│   ├── input_fasta_prot_aai.tsv
│   └── input_fasta_prot.m8
├── results
│   ├── *input_fasta_ictv.csv
│   └── input_fasta.tsv
└── tmp
```

#### Expected runtime ?

It takes ~4hrs to run vcat on the ICTV Taxonomy challenge dataset on a laptop computer.

---
### running other workflows
---
```
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
---
#### Some additional stuff
---
```
# to view the contig length distribution of your contigs
seqkit fx2tab -lg {input.fasta} | awk -F "\t" '{print $4}' | tail -n +2 | hist -b 100 -s 10

# to view the length distribution of contigs in the ani calculation output (after apply a filter)
cat testing/gut_jaeger/nuc/gut_jaeger_virus_seqs_fasta_genome_ani.tsv | awk -F "\t" '$8 > 0.1' | awk -F "\t" '{print $3}' | tail -n +2 | hist -b 100 -x

```


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
