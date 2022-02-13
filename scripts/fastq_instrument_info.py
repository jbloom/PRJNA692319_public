#!/usr/bin/env python

"""Information about the machine used to generate FASTQ from header."""


import snakemake
import collections
import re

import pandas as pd

import pysam


# get variables from snakemake
fastq_file = snakemake.input.fastq
output_csv = snakemake.output.csv
fastq_name = snakemake.wildcards.fastq

# For Illumina, match read name in CASAVA 1.8 format, splitting into two
# groups:
#  - the machine / run / flowcell / lane ID
#  - the cluster location (should be unique for all reads)
illumina_regex = re.compile(
    r'(?P<info>[\w\-]+:\d+:[\w\-]+:\d+):(?P<cluster>\d+:\d+:\d+)')
info = collections.defaultdict(int)
clusters = collections.defaultdict(int)

# For BGI-SEQ, get flow-cell serial number as from here:
# https://github.com/powellgenomicslab/BGI_vs_Illumina_Benchmark
bgi_regex = re.compile(r'(?P<info>[A-Z0-9]{10}L\d+)(?P<cluster>C\d+R\d+/[12])')

print(f"Parsing {fastq_file}")

platform = None
platform_regexes = {'illumina': illumina_regex, 'bgi': bgi_regex}
with pysam.FastxFile(snakemake.input.fastq) as f:
    for r in f:
        if platform is None:
            for key, regex in platform_regexes.items():
                if regex.fullmatch(r.name):
                    platform = key
                    break
            else:
                raise ValueError(f"cannot match {r.name=} for either platform")
        m = platform_regexes[platform].fullmatch(r.name)
        if not m:
            raise ValueError(f"cannot match {r.name=} for {platform=}")
        info[m.group('info')] += 1
        clusters[m.group('info') + ':' + m.group('cluster')] += 1

if max(clusters.values()) > 1:
    multi_clusters = pd.Series(clusters).sort_values(False)
    raise ValueError('duplicate clusters\n' + str(multi_clusters))

with open(output_csv, 'w') as f:
    f.write('fastq,info\n')
    for info in info:
        f.write(f"{fastq_name},{info}\n")
