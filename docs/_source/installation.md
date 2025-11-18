## Installation


currently you can install the development version of vcat with following steps.


#### Installing with conda/mamba
clone the git repositoty and use the environment.yaml file to create a conda environemnt.
alternatively you can also use pixi. Then, install the vcat python package into the conda 
environment you just created.

```bash
git clone https://github.com/Yasas1994/vcat.git
cd vcat
mamba create -f environment.yml
mamba activate vcat
pip install .
vcat --help
```



#### Running vcat with singularity
clone the git reposity and used the Apptainer definition file to build a singularity container.

```bash
git clone https://github.com/Yasas1994/vcat.git
cd vcat
apptainer build vcat.sif Apptainer
apptainer run vcat.sif vcat --help
```


if everything goes soothly, you should see vcat help on the the terminal.

```text
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


#### Downloading pre-built databases

You can download pre-built databases instead of building from scratch. 
Currently, msl39v and masl40v1 are available to download

```bash

vcat downloaddb --dbversion masl40v1 -d <path-to-save-the-database> --cores 1
```

optionally, you can also build it yourself with the following command.

```bash
vcat preparedb -d <path-to-save-the-database>
```
