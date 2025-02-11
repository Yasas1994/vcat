import click
import subprocess
import multiprocessing
import os
import sys
import yaml
import polars as pl
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm
from vcat.utils import ani_summary, aai_summary, index_m8, load_chunk



def load_configfile(file_path):
    with open(file_path, 'r') as f:
        return yaml.safe_load(f)

def parse_csv(ctx, param, value):
    """Split the comma-separated input into a list."""
    if value:
        return value.split(",")
    return []

from .color_logger import logger

# Define the directory containing the pipeline files
PIPELINE_DIR = os.path.join(os.path.dirname(__file__), "./pipeline")
VERSION = "0.0.1a"
CONFIG = os.path.join(PIPELINE_DIR, "config.yaml")
HEADER = "query,target,theader,fident,qlen,tlen,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,taxid,taxname,taxlineage"

def handle_max_mem(max_mem, profile):
    "Specify maximum virtual memory to use by atlas."
    "For numbers >1 its the memory in GB. "
    "For numbers <1 it's the fraction of available memory."

    if profile is not None:
        if max_mem is not None:
            logger.info(
                "Memory requirements are handled by the profile, I ignore max-mem argument."
            )
        # memory is handled via the profile, user should know what he is doing
        return ""
    else:
        import psutil
        from math import floor

        # calulate max  system meory in GB (float!)
        max_system_memory = psutil.virtual_memory().total / (1024**3)

        if max_mem is None:
            max_mem = 0.95
        if max_mem > 1:
            if max_mem > max_system_memory:
                logger.critical(
                    f"You specified {max_mem} GB as maximum memory, but your system only has {floor(max_system_memory)} GB"
                )
                sys.exit(1)

        else:
            max_mem = max_mem * max_system_memory

        # specify max_mem_string including java mem and max mem

        return f" --resources mem={floor(max_mem)} mem_mb={floor(max_mem*1024)} java_mem={floor(0.85* max_mem)} "


def get_snakefile(file=f"{PIPELINE_DIR }/Snakefile"):
    sf = os.path.join(os.path.dirname(os.path.abspath(__file__)), file)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file; tried %s" % sf)
    return sf

def update_config(config_path, data:dict):

    import yaml

    # Step 1: Read the YAML file
    with open(config_path, "r") as file:
        config = yaml.safe_load(file)  # Load the YAML into a Python dictionary

    # Step 2: Update the configuration
    for k,v in data.items():
        config[k] = v


    # Step 3: Write back to the YAML file
    with open(config_path, "w") as file:
        yaml.safe_dump(config, file, default_flow_style=False)




@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option("0.0.1a")
@click.pass_context
def cli(obj):
    r"""
    vcat: a command-line tool-kit for adding ICTV taxonomy annotations to virus contigs,
    mapping reads to virus genomes and much more.
    (https://github.com/Yasas1994/vcat)"""
    pass

@cli.command(
    context_settings=dict(ignore_unknown_options=True),
    short_help="run contig annotation workflow",
    help="""
    The Virus Contig Annotation Tool (vcat) is a straightforward, homology-based application designed to 
    provide taxonomic annotations for viral contigs. 
    
    The workflow consists of three steps 

    Nucleotide Database Comparison: Query sequences are first compared to the nucleotide database. 
    If matches with >95% ANI and >85% query coverage are found, the lineage of the target 
    sequence is assigned to the query sequence.

    Protein Database Comparison: For query sequences with no significant matches at the nucleotide
    level, open reading frames (ORFs) within the query sequence are extracted and compared to the 
    viral protein database to calculate AAI and the query coverage (number of proteins shared between
    query and the targer)

    Profile Database Comparison: Sequences that remain unannotated at the protein level are then
    compared to the protein profile database to provide taxonomic annotations.

    """,
)
@click.option(
    "-i",
    "--input",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="input file/s to run vcat on",
    required=True,
)
@click.option(
    "-o",
    "--output",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="dir to store vcat results",
    required=True,
)
@click.option(
    "-j",
    "--jobs",
    type=int,
    default=multiprocessing.cpu_count(),
    show_default=True,
    help="use at most this many jobs in parallel (see cluster submission for more details).",
)
@click.option(
    "--profile",
    default=None,
    help="snakemake profile e.g. for cluster execution.",
)
@click.option(
    "-n",
    "--dryrun",
    is_flag=True,
    default=False,
    show_default=True,
    help="Test execution.",
)
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
def contigs(
    input, output, jobs, profile, dryrun, snakemake_args
):
    """
    Runs vcat pipeline on contigs

    Most snakemake arguments can be appended to the command for more info see 'snakemake --help'
    """

    logger.info(f"vcat version: {VERSION}")


    conf = load_configfile(CONFIG)
    db_dir = conf["database_dir"]

    cmd = (
        "snakemake --snakefile {snakefile} "
        "--jobs {jobs} --rerun-incomplete "
        " --configfile {configfile} "
        "--scheduler greedy "
        " --show-failed-logs "
        " --groups group1=1 "
        "--config database_dir='{db_dir}' sample='{input}' output_dir='{output}' {add_args} "
        "{args}"
    ).format(
        snakefile=get_snakefile("./pipeline/Snakefile"),
        jobs=jobs,
        configfile=CONFIG,
        db_dir=db_dir,
        input=input,
        output=output,
        add_args="" if snakemake_args and snakemake_args[0].startswith("-") else "--",
        args=" ".join(snakemake_args),
    )
    logger.debug("Executing: %s" % cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        # removes the traceback
        logger.critical(e)
        exit(1)

@cli.command(
    context_settings=dict(ignore_unknown_options=True),
    short_help="run read annotation workflow",
    help="""
    The Virus Contig Annotation Tool (vcat) is a straightforward, homology-based application designed to 
    provide taxonomy annotations to virus contigs and mapping reads directly to virus genomes.

    usage
    -----
    vcat reads [OPTIONS] -i1 pair1.fastq -i2 pair2.fastq -o mapping_results.tsv

    """,
)
@click.option(
    "-i",
    "--input",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="input read file/s to run vcat on",
    required=True,
)
@click.option(
    "-o",
    "--output",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="dir to store vcat results",
    required=True,
)
@click.option(
    "-j",
    "--jobs",
    type=int,
    default=multiprocessing.cpu_count(),
    show_default=True,
    help="use at most this many jobs in parallel (see cluster submission for more details).",
)
@click.option(
    "--profile",
    default=None,
    help="snakemake profile e.g. for cluster execution.",
)
@click.option(
    "-n",
    "--dryrun",
    is_flag=True,
    default=False,
    show_default=True,
    help="Test execution.",
)
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
def reads(
    input, output, jobs, profile, dryrun, snakemake_args
):
    """
    Runs vcat pipeline on reads

    Most snakemake arguments can be appended to the command for more info see 'snakemake --help'
    """

    logger.info(f"vcat version: {VERSION}")


    conf = load_configfile(CONFIG)
    db_dir = conf["database_dir"]

    cmd = (
        "snakemake --snakefile {snakefile} "
        "--jobs {jobs} --rerun-incomplete "
        " --configfile {configfile} "
        "--scheduler greedy "
        " --show-failed-logs "
        " --groups group1=1 "
        "--config database_dir='{db_dir}' sample='{input}' output_dir='{output}' {add_args} "
        "{args}"
    ).format(
        snakefile=get_snakefile("./pipeline/Snakefile"),
        jobs=jobs,
        configfile=CONFIG,
        db_dir=db_dir,
        input=input,
        output=output,
        add_args="" if snakemake_args and snakemake_args[0].startswith("-") else "--",
        args=" ".join(snakemake_args),
    )
    logger.debug("Executing: %s" % cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        # removes the traceback
        logger.critical(e)
        exit(1)

# Download
@cli.command(
    context_settings=dict(ignore_unknown_options=True),
    short_help="download and build reference databases",
)
@click.option(
    "-d",
    "--db-dir",
    help="location to store databases",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    required=True,
)
@click.option(
    "-j",
    "--jobs",
    default=1,
    type=int,
    show_default=True,
    help="number of simultaneous downloads",
)
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
def preparedb(db_dir, jobs, snakemake_args):
    """Executes a snakemake workflow to downlod and building the databases
    """
    logger.info("Building taxdb")
    cmd = (
        "snakemake --snakefile {snakefile} "
        "--jobs {jobs} --rerun-incomplete "
        " --configfile {configfile} "
        "--scheduler greedy "
        " --show-failed-logs "
        "--config database_dir='{db_dir}' {add_args} "
        "{args}"
    ).format(
        snakefile=get_snakefile("./pipeline/rules/taxdump.smk"),
        jobs=jobs,
        configfile=CONFIG,
        db_dir=db_dir,
        add_args="" if snakemake_args and snakemake_args[0].startswith("-") else "--",
        args=" ".join(snakemake_args),
    )
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        # removes the traceback
        logger.critical(e)
        exit(1)

    logger.info("Building mmseqs databases")

    cmd = (
        "snakemake --snakefile {snakefile} "
        "--jobs {jobs} --rerun-incomplete "
        " --configfile {configfile} "
        "--scheduler greedy "
        " --show-failed-logs "
        "--config database_dir='{db_dir}' {add_args} "
        "{args}"
    ).format(
        snakefile=get_snakefile("./pipeline/rules/createdb.smk"),
        jobs=jobs,
        configfile=CONFIG,
        db_dir=db_dir,
        add_args="" if snakemake_args and snakemake_args[0].startswith("-") else "--",
        args=" ".join(snakemake_args),
    )

    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        # removes the traceback
        logger.critical(e)
        exit(1)

    logger.info(f"Adding {db_dir} to config.yaml")
    update_config(config_path=CONFIG, data={"database_dir" :db_dir})


# utility functions
@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option("0.0.1a")
@click.pass_context
def utils(obj):
    """
    tool chain for calculating ani, aai and visualizations
    """
    pass
@utils.command(
    context_settings=dict(ignore_unknown_options=True),
    help="""calculates ani from mmseqs ICTV genome comparision results
            and writes the results to <input>_ani.tsv

            output includes qcov (i.e synonymous to alignment fraction),
            ani (average nucleotide identity) and tani (ani x qcov)

            usage                                                                      
            -----                                                                      

            vcat utils ani [OPTIONS] -i contigs.fasta

        """,
)
@click.option(
    "-i",
    "--input",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="BLAST tabular (m8) output file",
    required=True
)
@click.option(
    "--header",
    callback=parse_csv,
    help="columnames of the m8 file",
    default = HEADER,
    required=False
)
@click.option(
    "--qcov",
    type=float,
    default=0,
    help="filter results below this qcov cutoff",
    required=False
)
@click.option(
    "--tani",
    type=float,
    default=0,
    help="filter results below this tani cutoff",
    required=False
)
@click.option(
    "--ani",
    type=float,
    default=0,
    help="filter results below this ani cutoff",
    required=False
)
@click.option(
    "--all",
    is_flag=True,
    default=False,
    help="get ani for all hits per query sequence. by default only outputs the besthits",
    required=False
)
@click.option(
    "--batch",
    type=int,
    default=5000,
    help="number of records to process at a time",
    required=False
)
def ani(input, header, ani, tani, qcov, all, batch):
    """
    calculates the average nucleotide identity and coverage of a query sequence 
    to the best 
    """
    
    CHUNK_SIZE = batch
    file_name = os.path.basename(input)
    
    index = index_m8(input)
    tmp_files = []
    with logging_redirect_tqdm():
        for i in tqdm(range(0, len(index), CHUNK_SIZE), ncols=70, ascii=' ='):
            finput = load_chunk(input, index=index, recstart=i, recend=min(i + CHUNK_SIZE, len(index)))
            outfile = os.path.join(os.path.dirname(input), f"{file_name.split('.')[0]}_ani_{min(i + CHUNK_SIZE, len(index))}.tsv")
            status = ani_summary(finput, outfile, all=all, header=header)
            if status == 0:
                logger.info(f"{outfile} updated")
                tmp_files.append(outfile)
                
            else:
                # exit code 1
                logger.error("error occured!")
                logger.exception(status)
                exit(1)

    logger.info("merging temporary files")
    tmp = [pl.read_csv(f, separator="\t") for f in tmp_files]
    df = pl.concat([i for i in tmp if not i.is_empty()])
    outfile = os.path.join(os.path.dirname(input), f"{file_name.split('.')[0]}_ani.tsv")
    df.write_csv(outfile, separator="\t")
    logger.info(f"{outfile} updated")

    # Remove temporary TSV files
    for file in tmp_files:
        os.remove(file)

@utils.command(
    context_settings=dict(ignore_unknown_options=True),
    help="""
            calculates aai from mmseqs ICTV viral protein comparision results
            and writes the results to <input>_aai.tsv

            output includes qcov (i.e number of genes shared between the query
            and target over number of genes on the query),
            aai (average aminoacid identity) and taai (aai x qcov)

            usage                                                                      
            -----                                                                      

            vcat utils aai  [OPTIONS] -i contigs.fasta -g configs.gff -d [DBDIR]

        """,
)
@click.option(
    "-i",
    "--input",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="BLAST tabular (m8) output file",
    required=True
)
@click.option(
    "-g",
    "--gff",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="gff file with gene cordinates of the query sequences",
    required=True
)
@click.option(
    "--header",
    callback=parse_csv,
    help="columnames of the m8 file ",
    default = HEADER,
    required=False
)
@click.option(
    "--qcov",
    type=float,
    default=0,
    help="filter results below this qcov cutoff",
    required=False
)
@click.option(
    "--taai",
    type=float,
    default=0,
    help="filter results below this taai cutoff",
    required=False
)
@click.option(
    "--aai",
    type=float,
    default=0,
    help="filter results below this aai cutoff",
    required=False
)
@click.option(
    "-d",
    "--dbdir",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="path to the ref databases",
    required=True
)
@click.option(
    "--batch",
    type=int,
    default=5000,
    help="number of records to process at a time",
    required=False
)
def aai(input, header, aai, taai, qcov, batch, dbdir, gff):
    """
    calculates the average aminoacid identity and coverage of a query sequence 
    to the best 
    """
    CHUNK_SIZE = batch
    file_name = os.path.basename(input)
    
    index = index_m8(input)
    tmp_files = []
    with logging_redirect_tqdm():
        for i in tqdm(range(0, len(index), CHUNK_SIZE), ncols=70, ascii=' ='):
            finput = load_chunk(input, index=index, recstart=i, recend=min(i + CHUNK_SIZE, len(index)))
            outfile = os.path.join(os.path.dirname(input), f"{file_name.split('.')[0]}_aai_{min(i + CHUNK_SIZE, len(index))}.tsv")
            status = aai_summary(finput, gff, dbdir, outfile, header)
            if status == 0:
                logger.info(f"{outfile} updated")
                tmp_files.append(outfile)
                
            else:
                # exit code 1
                logger.error("error occured!")
                logger.exception(status)
                exit(1)

    logger.info("merging temporary files")
    tmp = [pl.read_csv(f, separator="\t") for f in tmp_files]
    df = pl.concat([i for i in tmp if not i.is_empty()])
    outfile = os.path.join(os.path.dirname(input),
                           f"{file_name.split('.')[0]}{'_all' if all else ''}_aai.tsv")
    df.write_csv(outfile, separator="\t")
    logger.info(f"{outfile} updated")

    # Remove temporary TSV files
    for file in tmp_files:
        os.remove(file)

cli.add_command(utils)

if __name__ == "__main__":
    cli()

