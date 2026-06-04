#!/usr/bin/env python
"""Convert a MetaPhlAn 4 SGB profile into a GTDB-taxonomy profile.

This is a self-contained adaptation of the official MetaPhlAn utility
(metaphlan/utils/sgb_to_gtdb_profile.py, by Aitor Blanco-Miguez). Two
deliberate changes were made so the pipeline can own the mapping file:

  1. The SGB->GTDB assignment table is supplied via -m/--mapping instead of
     being hard-coded relative to the installed MetaPhlAn package. This lets
     the Nextflow pipeline download and cache the table in its storeDir and
     pass that copy in explicitly (see modules/processes/databases.nf).
  2. The util_fun info/error helpers are inlined so the script runs standalone
     from the pipeline bin/ directory without importing the MetaPhlAn package.

The conversion logic (get_gtdb_profile) is otherwise unchanged from upstream:
  - 1:1 mappings are substituted directly,
  - n:1 mappings are binned (summed) into the shared GTDB taxon,
  - higher ranks are aggregated up the lineage from species.

Upstream source:
https://github.com/biobakery/MetaPhlAn/blob/master/metaphlan/utils/sgb_to_gtdb_profile.py
"""

import sys
import time
import argparse as ap


def info(message):
    sys.stdout.write('{}\n'.format(message))
    sys.stdout.flush()


def error(message, exit=False):
    sys.stderr.write('[ERROR] {}\n'.format(message))
    sys.stderr.flush()
    if exit:
        sys.exit(1)


def read_params():
    """Reads and parses the command line arguments of the script."""
    p = ap.ArgumentParser(
        description=__doc__, formatter_class=ap.ArgumentDefaultsHelpFormatter)
    p.add_argument('-i', '--input', type=str, default=None,
                   help="The input MetaPhlAn SGB profile")
    p.add_argument('-o', '--output', type=str, default=None,
                   help="The output GTDB profile")
    p.add_argument('-m', '--mapping', type=str, default=None,
                   help="The SGB-to-GTDB assignment table (TSV)")
    return p.parse_args()


def check_params(args):
    """Checks the mandatory command line arguments of the script."""
    if not args.input:
        error('-i (or --input) must be specified', exit=True)
    if not args.output:
        error('-o (or --output) must be specified', exit=True)
    if not args.mapping:
        error('-m (or --mapping) must be specified', exit=True)


def get_gtdb_profile(mpa_profile, gtdb_profile, mapping_file):
    """Creates the GTDB profile from a MetaPhlAn one."""
    tax_levels = ['d', 'p', 'c', 'o', 'f', 'g', 's']

    def get_parent_taxon(taxon, level):
        parts = taxon.split(';')
        if len(parts) > 1:
            return ';'.join(parts[:-1])
        else:
            return parts[0]

    sgb2gtdb = dict()
    with open(mapping_file, 'r') as read_file:
        for line in read_file:
            line = line.strip().split('\t')
            sgb2gtdb[line[0]] = line[1]
    with open(gtdb_profile, 'w') as wf:
        with open(mpa_profile, 'r') as rf:
            unclassified = 0
            abundances = {x: dict() for x in tax_levels}
            for line in rf:
                if line.startswith('#mpa_'):
                    wf.write(line)
                    wf.write('#clade_name\trelative_abundance\n')
                elif line.startswith('UNCLASSIFIED'):
                    unclassified = float(line.strip().split('\t')[2])
                    wf.write('UNCLASSIFIED\t{}\n'.format(unclassified))
                elif 't__SGB' in line:
                    line = line.strip().split('\t')
                    abundance = float(line[2])
                    sgb_id = line[0].split('|')[-1][3:]
                    gtdb_tax = sgb2gtdb.get(sgb_id, None)
                    if gtdb_tax is None:
                        continue
                    if gtdb_tax not in abundances['s']:
                        abundances['s'][gtdb_tax] = 0
                    abundances['s'][gtdb_tax] += abundance
        tax_levels_rev = list(reversed(tax_levels))
        for i, tax_level in enumerate(tax_levels_rev[:-1]):
            for tax in abundances[tax_level]:
                parent_tax = get_parent_taxon(tax, tax_level)
                new_level = tax_levels_rev[i + 1]
                if parent_tax not in abundances[new_level]:
                    abundances[new_level][parent_tax] = 0
                abundances[new_level][parent_tax] += abundances[tax_level][tax]
        for tax_level in tax_levels:
            for tax in abundances[tax_level]:
                wf.write('{}\t{}\n'.format(tax, abundances[tax_level][tax]))


def main():
    t0 = time.time()
    args = read_params()
    info("Start execution")
    check_params(args)
    get_gtdb_profile(args.input, args.output, args.mapping)
    exec_time = time.time() - t0
    info("Finish execution ({} seconds)".format(round(exec_time, 2)))


if __name__ == "__main__":
    main()
