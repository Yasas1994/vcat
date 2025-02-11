import os
# Set the input file and output directory as global variables
configfile: "config.yaml"
# this values are incuded in the snakefile
DBDIR = config["database_dir"]
CWD = os.path.dirname(__file__)
rule all:
    input:
        f"{DBDIR}/VMR_latest/virus_genome.accession2taxid",
        f"{DBDIR}/VMR_latest/virus_protein.accession2taxid",
        f"{DBDIR}/VMR_latest/proteins.faa",
        f"{DBDIR}/VMR_latest/genomes.fna",
        f"{DBDIR}/VMR_latest/bbmap_index/index_done",
        f"{DBDIR}/VMR_latest/mmseqs_proteins/mmseqs_proteins",
        f"{DBDIR}/VMR_latest/mmseqs_genomes/mmseqs_genomes",
        f"{DBDIR}/VMR_latest/mmseqs_pprofiles/mmseqs_pprofiles",
        f"{DBDIR}/VMR_latest/mmseqs_pprofiles/mmseqs_pprofiles_lca.tsv"

# Download the ICTV genomes
rule download_genbank:
    input:
        f"{DBDIR}/ictv.xlsx"
    output:
        directory(f"{DBDIR}/VMR_latest/tmp"),
        f"{DBDIR}/VMR_latest/tmp/download_complete"
    log:
        f"{DBDIR}/logs/download_genbank.log"
    shell:
        """
        download_gb.py {input} {output} 2> {log}
        """

# preprocess and prepare .gb files to build mmseqs datasets
rule prepare:
    input:
        DATADIR=f"{DBDIR}/ictv-taxdump.tar.gz",
        XLTABLE=f"{DBDIR}/ictv.xlsx",
        MAN=f"{DBDIR}/VMR_latest/tmp/download_complete"

    output:
        f"{DBDIR}/VMR_latest/virus_genome.accession2taxid",
        f"{DBDIR}/VMR_latest/virus_protein.accession2taxid",
        f"{DBDIR}/VMR_latest/proteins.faa",
        f"{DBDIR}/VMR_latest/genomes.fna",

    log:
        f"{DBDIR}/logs/prepare.log"
    params:
        DBFILE =f"{DBDIR}/VMR_latest/VMR_latest.sql",
        GBDIR=f"{DBDIR}/VMR_latest/tmp",
        SEQDIR=f"{DBDIR}/VMR_latest",
    shell:
        """
        prepare.py {input.DATADIR} {params.DBFILE} {input.XLTABLE} {params.GBDIR} {params.SEQDIR} 2> {log}
        """
     
rule make_mmseqs_proteindb:
    input:
        seqs=f"{DBDIR}/VMR_latest/proteins.faa",
        map=f"{DBDIR}/VMR_latest/virus_protein.accession2taxid"
    output:
        f"{DBDIR}/VMR_latest/mmseqs_proteins/mmseqs_proteins"
    log:
        f"{DBDIR}/logs/make_mmseqs_proteindb.log"
    params:
        tmp=f"{DBDIR}/tmp",
        taxdump=f"{DBDIR}/ictv-taxdump",
        map=f"{DBDIR}/VMR_latest/mmseqs_proteins/mmseqs_proteins_mapping"
    shell:
        """
        mmseqs createdb --dbtype 1 {input.seqs} {output}  2> {log}
        mmseqs createtaxdb {output} {params.tmp} --ncbi-tax-dump {params.taxdump} --tax-mapping-file {input.map}  >> {log}
        mmseqs nrtotaxmapping {input.map} {output} {params.map}  2> {log}
        """

rule make_mmseqs_genomedb:
    input:
        seqs=f"{DBDIR}/VMR_latest/genomes.fna",
        map=f"{DBDIR}/VMR_latest/virus_genome.accession2taxid"
    output:
        f"{DBDIR}/VMR_latest/mmseqs_genomes/mmseqs_genomes"
    log:
        f"{DBDIR}/logs/make_mmseqs_genomedb.log"
    params:
        tmp=f"{DBDIR}/tmp",
        taxdump=f"{DBDIR}/ictv-taxdump",
        map=f"{DBDIR}/VMR_latest/mmseqs_genomes/mmseqs_genomes_mapping"
    shell:
        """
        mmseqs createdb --dbtype 0 {input.seqs} {output}  2> {log}
        mmseqs createtaxdb {output} {params.tmp} --ncbi-tax-dump {params.taxdump} --tax-mapping-file {input.map}  2> {log}
        mmseqs nrtotaxmapping {input.map} {output} {params.map}  2> {log}
        """

rule make_bbmap_index:
    input:
        seqs=f"{DBDIR}/VMR_latest/genomes.fna",
    output:
        f"{DBDIR}/VMR_latest/bbmap_index/index_done"
    log:
        f"{DBDIR}/logs/make_bbmap_genomedb.log"
    params:
        build = 1
        outdir=f"{DBDIR}/VMR_latest/bbmap_index"
    shell:
        """
        bbmap.sh ref={input} path={params.outdir} build={params.build} &> {log} && touch {output}

        """

rule make_mmseqs_profiledb:
    input:
        proteindb=f"{DBDIR}/VMR_latest/mmseqs_proteins/mmseqs_proteins",
        
    output:
        profiledb=f"{DBDIR}/VMR_latest/mmseqs_pprofiles/mmseqs_pprofiles",
        pclustdir=directory(f"{DBDIR}/VMR_latest/mmseqs_pclusters")

    log:
        f"{DBDIR}/logs/make_mmseqs_profiledb.log"
    params:
        tmp=f"{DBDIR}/tmp",
        pclustdb=f"{DBDIR}/VMR_latest/mmseqs_pclusters/mmseqs_pclusters",
        clusters_tsv=f"{DBDIR}/VMR_latest/mmseqs_pclusters/mmseqs_pclusters.tsv",
        profiledir=f"{DBDIR}/VMR_latest/mmseqs_pprofiles",
        seqrep=f"{DBDIR}/VMR_latest/mmseqs_pclusters/sequenceRepDb",
        seqrep_h=f"{DBDIR}/VMR_latest/mmseqs_pclusters/sequenceRepDb_h",
        protein_h=f"{DBDIR}/VMR_latest/mmseqs_proteins/mmseqs_proteins_h",
        
    shell:
        """
        mkdir -p {output.pclustdir} >> {log}
        mkdir -p {params.profiledir} >> {log}
        mmseqs cluster {input.proteindb} {params.pclustdb} {params.tmp} --cluster-reassign  -s 7 2> {log}
        mmseqs createtsv {input.proteindb} {input.proteindb} {params.pclustdb} {params.clusters_tsv} 2> {log}
        mmseqs createsubdb {params.pclustdb}  {input.proteindb} {params.seqrep}  2> {log}
        mmseqs createsubdb {params.pclustdb}  {params.protein_h} {params.seqrep_h}  2> {log}
        mmseqs result2profile {params.seqrep} {input.proteindb} {params.pclustdb} {output.profiledb}  2> {log}

        """

rule cluster_lca:
    input:
        DATADIR=DBDIR

    output:
        f"{DBDIR}/VMR_latest/mmseqs_pprofiles/mmseqs_pprofiles_lca.tsv",

    log:
        f"{DBDIR}/logs/prepare.log"

    shell:
        """
        cluster_lca.py {input} 2> {log}
        """