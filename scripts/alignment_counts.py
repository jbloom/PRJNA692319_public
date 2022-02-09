"""Aggregate alignment counts across accessions."""


import subprocess

import pandas as pd


accessions = snakemake.config['accessions']
bams = snakemake.input.bams
high_ident_bams = snakemake.input.high_ident_bams

assert len(accessions) == len(bams) == len(high_ident_bams)

records = []

for accession, bam, high_ident_bam in zip(accessions, bams, high_ident_bams):
    for name, bamfile in [('any', bam), ('high_identity', high_ident_bam)]:
        res = subprocess.run(['samtools', 'view', '-c', bamfile], capture_output=True)
        assert res.returncode == 0
        count = int(res.stdout)
        assert count >= 0
        records.append((accession, name, count))

(pd.DataFrame(records, columns=['accession', 'alignment_type', 'read_count'])
 .to_csv(snakemake.output.csv, index=False)
 )
