#!/usr/bin/env python

"""Script for rule ``scripts/trim3_polyA.py``."""

import snakemake
import Bio.SeqIO

seq = Bio.SeqIO.read(snakemake.input.fasta, 'fasta')
seq.seq = seq.seq.rstrip('Aa')
Bio.SeqIO.write([seq], snakemake.output.fasta, 'fasta')
