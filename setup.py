from setuptools import setup

__author__ = "Yasas Wijesekara"
__copyright__ = "Copyright 2024, Yasas Wijesekara"
__email__ = "yasas.Wijesekara@uni-greifswald.de"
__license__ = "BSD-3"

# read the contents of your README file
from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()


setup(
    name="vcat",
    version="0.0.1a",
    url="https://github.com/Yasas94/vcat",
    license=__license__,
    author=__author__,
    author_email=__email__,
    zip_safe=False,
    description="vact - workflows for transferring taxonomic annotations to virus contigs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=["vcat"],
    include_package_data=True,
    package_data={
        "vcat":["pipeline/*", "pipeline/rules/*"],
    },
    data_files=[(".", ["README.md", "LICENSE.txt", "MANIFEST.in"])],
    install_requires=[],
    # install via conda: click, pandas, pyyaml, snakemake
    entry_points={"console_scripts": ["vcat = vcat.cli:cli"]},
    classifiers=["Topic :: Scientific/Engineering :: Bio-Informatics"],
    scripts=["scripts/download_gb.py",
             "scripts/cluster_lca.py",
             "scripts/postprocess.py",
             "scripts/ictv_out.py",
             "scripts/prepare.py"],
)