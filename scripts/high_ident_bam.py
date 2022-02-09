"""Implements ``snakemake`` rule `high_ident_bam`."""


import pysam


min_length = snakemake.params.min_length
min_identity = snakemake.params.min_identity

print(f"Examining reads in {snakemake.input.bam=}")
print(f"Writing high-identity ones to {snakemake.output.bam=}")
print(f"The high-identity cutoffs are {min_length=} and {min_identity=}")
counts = high_ident_counts = 0
with pysam.AlignmentFile(snakemake.input.bam, 'rb') as bamfile:
    with pysam.AlignmentFile(snakemake.output.bam, 'wb',
                             template=bamfile) as bamfile_out:
        for read in bamfile:
            counts += 1
            length = read.query_alignment_length
            if length >= min_length:
                try:
                    nm = read.get_tag('NM')
                except KeyError:
                    raise KeyError(f"read in {bam=} lacks NM tag")
                assert 0 <= nm, f"{nm=}, {length=}"
                identity = 1 - nm / length
                if identity >= min_identity:
                    high_ident_counts += 1
                    bamfile_out.write(read)
print(f"Retained {high_ident_counts=} of {counts=} reads")
