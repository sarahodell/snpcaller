"""Snakemake wrapper for fastqc"""
### From https://github.com/dib-lab/eelpond/blob/master/rules/fastqc/wrapper.py

from os import path
from tempfile import TemporaryDirectory

from snakemake.shell import shell

log=snakemake.log_fmt_shell(stdout=False,stderr=True)

def basename_without_ext(file_path):
    """ Returns basename of file path without the file extension."""

    base=path.basename(file_path)

    split_ind = 2 if base.endswith(".gz") else 1
    base=".".join(base.split(".")[:-split_ind])

    return base

with TemporaryDirectory() as tempdir:
    shell("fastqc {snakemake.params} --quiet "
          "--outdir {tempdir} snakemake.input[0]}"
          "{log}")

    # Move outputs into proper position.
    output_base = basename_without_ext(snakemake.input[0])
    html_path = path.join(tempdir,output_base + "_fastqc.html")
    zip_path = path.join(tempdir, output_base + "_fastqc.zip")

    if snakemake.output.html != html_path:
        shell("mv {html_path} {snakemake.output.html}")

    if snakemake.output.zip != zip_path:
        shell("mv {zip_path} {snakemake.output.zip}")

