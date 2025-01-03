import click
import subprocess
import multiprocessing
import os
import sys
import yaml

def load_configfile(file_path):
    with open(file_path, 'r') as f:
        return yaml.safe_load(f)


from .color_logger import logger

# Define the directory containing the pipeline files
PIPELINE_DIR = os.path.join(os.path.dirname(__file__), "./pipeline")
VERSION = "0.0.1a"
CONFIG = os.path.join(PIPELINE_DIR, "config.yaml")

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
    """vcat: A command-line tool for taxonomically annotating virus contigs"""
    pass

@cli.command(
    "annotate",
    context_settings=dict(ignore_unknown_options=True),
    short_help="run annotation workflow",
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
def run_workflow(
    input, output, jobs, profile, dryrun, snakemake_args
):
    """Runs vcat pipeline

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
    "preparedb",
    context_settings=dict(ignore_unknown_options=True),
    short_help="download reference files",
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
def run_download(db_dir, jobs, snakemake_args):
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

if __name__ == "__main__":
    cli()

