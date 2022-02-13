#!/usr/bin/env python

"""Analyze variants called by ``ivar`` relative to outgroups."""

import snakemake
import pandas as pd


site_map = pd.read_csv(snakemake.input.site_map)

assert len(snakemake.params.descriptors) == len(snakemake.input.variants)
variants = (pd.concat([pd.read_csv(f, sep='\t').assign(**descriptor)
                       for (f, descriptor) in zip(snakemake.input.variants,
                                                  snakemake.params.descriptors)
                       ],
                      ignore_index=True)
            .rename(columns={'POS': 'site'})
            [['site', 'accession', 'identity', 'REF', 'ALT',
              'REF_DP', 'ALT_DP', 'ALT_FREQ']]
            .merge(site_map, how='left', on='site', validate='many_to_one')
            )
assert (variants['REF'] == variants['reference']).all()
variants = variants.drop(columns='reference')

variants.to_csv(snakemake.output.all_csv, index=False)

outgroups = list(snakemake.config['outgroups'])
(variants
 .query(' | '.join(f"(ALT == `{outgroup}`)" for outgroup in outgroups))
 .to_csv(snakemake.output.to_outgroup_csv, index=False)
 )
