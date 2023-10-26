import os
import sys
from warnings import warn
from py.helpers import *
configfile: "config.yaml"

"""
Download the manifest data, merge bams if necessary, and index.
"""

rule_all = [
 
]


rule download_hifi_bam:
    output: 
        bam = "data/hifi/bams/{sample}.bam",
        bai = "data/hifi/bams/{sample}.bam.bai"
    input: 
        "data/hifi/manifests/{sample}.txt"