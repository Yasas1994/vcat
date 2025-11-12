from pathlib import Path
import multiprocessing
import subprocess
import os

# Extract global variables from the config
DBDIR = config["database_dir"]
INFILES = config["results"]  # This can be a list of files
API_PARAMS = config["api"]
AAI_PARAMS = config["aai"]
ANI_PARAMS = config["ani"]
LEVEOUT_LEVEL = config["level"]
BATCH = config["batch"]

# Define database paths
PROTDB = f"{DBDIR}/VMR_latest/mmseqs_proteins/mmseqs_proteins"
GENOMEDB = f"{DBDIR}/VMR_latest/mmseqs_genomes/mmseqs_genomes"
PROFDB = f"{DBDIR}/VMR_latest/mmseqs_pprofiles/mmseqs_pprofiles"
INDIR = os.path.dirname(INFILES)


# Function to dynamically generate output files
def generate_outputs(wildcards):
    return f"{OUTDIR}/{wildcards.sample}_genome.m8"


# Function to dynamically generate input files
def generate_inputs(wildcards):
    directory_ = os.path.dirname(INFILES)
    file_name, file_extension = os.path.splitext(os.path.basename(INFILES))
    return os.path.join(directory_, f"{wildcards.sample}{file_extension}")


# Rule to define final output
rule all:
    input:
        expand(f"{OUTDIR}/nuc/{{sample}}_{{ext}}_genome_leaveout_{LEVEOUT_LEVEL}_ani.tsv",
        sample=[os.path.splitext(os.path.basename(INFILES))[0]],
        ext=[os.path.basename(INFILES).split(".")[-1]]),
        expand(f"{OUTDIR}/prot/{{sample}}_{{ext}}_prot_leaveout_{LEVEOUT_LEVEL}_aai.tsv",
        sample=[os.path.splitext(os.path.basename(INFILES))[0]],
        ext=[os.path.basename(INFILES).split(".")[-1]]),
        expand(f"{OUTDIR}/prof/{{sample}}_{{ext}}_prof_leaveout_{LEVEOUT_LEVEL}_api.tsv",
        sample=[os.path.splitext(os.path.basename(INFILES))[0]],
        ext=[os.path.basename(INFILES).split(".")[-1]]),        

rule cal_ani:
    input:
        f"{OUTDIR}/nuc/{{sample}}_{{ext}}_genome.m8"
    output:
        f"{OUTDIR}/nuc/{{sample}}_{{ext}}_genome_leaveout_{LEVEOUT_LEVEL}_ani.tsv"
    log:
        f"{OUTDIR}/logs/{{sample}}_{{ext}}_leaveout_ani.log"
    params:
        batch = BATCH,
        ani = ANI_PARAMS,
        db = DBDIR,
        level = LEVEOUT_LEVEL
    shell:
        """
        filter
        vcat utils ani -i {input} \
        --header query,target,theader,fident,qlen,tlen,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,taxid,taxname,taxlineage \
        --level {params.level} --dbdir {params.dbdir} --batch {params.batch} {params.ani} &> {log}
        """

rule cal_aai:
    input:
        m8 = f"{OUTDIR}/prot/{{sample}}_{{ext}}_prot.m8",
        gff = f"{OUTDIR}/prot/{{sample}}_{{ext}}.gff"
    output:
        f"{OUTDIR}/prot/{{sample}}_{{ext}}_leaveout_{LEVEOUT_LEVEL}_prot_aai.tsv"
    log:
        f"{OUTDIR}/logs/{{sample}}_{{ext}}_leaveout_aai.log"
    params:
        batch = BATCH,
        aai = AAI_PARAMS,
        db = DBDIR,
        level = LEVEOUT_LEVEL,
    shell:
        """
        vcat utils aai -i {input.m8} -g {input.gff} -d {params.db} \
        --header query,target,theader,fident,qlen,tlen,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,taxid,taxname,taxlineage \
        --level {params.level} --batch {params.batch} {params.aai} &> {log}
        """

rule cal_api:
    input:
        m8 = f"{OUTDIR}/prof/{{sample}}_{{ext}}_prof.m8",
        gff = f"{OUTDIR}/prot/{{sample}}_{{ext}}.gff"
    output:
        f"{OUTDIR}/prof/{{sample}}_{{ext}}_leaveout_{LEVEOUT_LEVEL}_prof_api.tsv"
    log:
        f"{OUTDIR}/logs/{{sample}}_{{ext}}_leaveout_api.log"
    params:
        batch = BATCH,
        api = API_PARAMS,
        db = DBDIR,
        level = LEVEOUT_LEVEL,
    shell:
        """
        vcat utils api -i {input.m8} -g {input.gff} -d {params.db}\
        --header query,target,theader,fident,qlen,tlen,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,taxid,taxname,taxlineage \
        --level {params.level} --batch {params.batch} {params.api} &> {log}
        """

rule summarize:
    input:
        DATADIR = DBDIR,
        PROF = f"{OUTDIR}/prof/{{sample}}_{{ext}}_leaveout_{LEVEOUT_LEVEL}_prof_api.tsv",
        PROT = f"{OUTDIR}/prot/{{sample}}_{{ext}}_leaveout_{LEVEOUT_LEVEL}_prot_aai.tsv",
        NUC = f"{OUTDIR}/nuc/{{sample}}_{{ext}}_leaveout_{LEVEOUT_LEVEL}_genome_ani.tsv",

    output:
        f"{OUTDIR}/results/{{sample}}_{{ext}}_leaveout_{LEVEOUT_LEVEL}.tsv"
    log:
        f"{OUTDIR}/logs/{{sample}}_{{ext}}_leaveout_summarize.log"

    threads: int(workflow.cores * 0.75)

    params:
        ani=ANI_PARAMS

    shell:
        """ 
        postprocess.py {input.DATADIR} {input.NUC} {input.PROT} {input.PROF} {output} {params.ani} &> {log}
        
        """