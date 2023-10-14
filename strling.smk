import os
import sys
from warnings import warn
from py.helpers import *
configfile: "config.yaml"


# rule all inputs are ouputs of other rules 
rule_all = [
    # strling output
    expand(
        "{output_dir}/strling/genotypes/{sample}-genotype.tsv",
        output_dir=configfile["output_dir"],
        sample=configfile["samples"]
    ),
    expand(
        "{output_dir}/strling/motifs/{sample}-motif.tsv",
        output_dir=configfile["output_dir"],
        sample=configfile["samples"]
    ),
    # expansionhunter output
    expand(
        "{output_dir}/ehdn/profiles/{sample}.profile.json",
        output_dir=configfile["output_dir"],
        sample=configfile["samples"]
    ),
    # TRGT output
    #?
    # TRASH output
    #?
    # strling pileups
    expand(
        "{output_dir}/strling/{sample}-pileup.tsv",
        output_dir=configfile["output_dir"],
        sample=configfile["samples"]
    ),
