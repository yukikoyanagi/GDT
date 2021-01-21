#!/usr/bin/env python
#
# File: gdtg_ts.py
#
# Time-stamp: <2016-12-16 12:03:20 au447708>
#
# Description: Computes GDT-graph_TS scores for all decoys in the
# given decoy directory which correspond to the given target.
# seqfile is the file containing residue length for each target,
# and tertiarydir is the directory with tertiary bond files.
# It returns a dict of {decoy name: GDT-G}, unless the optional
# output is specified.
#
# Usage: Needs gdtcdp.py and cdp.py files in the same folder.
#
# Author: Yuki Koyanagi
# History:
#  2016-12-16: Created
#

import os
import argparse
import gdtcdp


# list of cut-off values to use when computing GDT-graph_TS
cutoffs = [25, 50, 100, 200]


def run(targetf, seqfile, decoydir, tertiarydir, output=None):
    ddir = os.path.abspath(decoydir)
    tertdir = os.path.abspath(tertiarydir)

    # load target and decoy graphs in dicts
    targetfile = os.path.abspath(targetf)
    tname, _ = os.path.splitext(os.path.basename(targetfile))
    tertfile = os.path.join(tertiarydir, os.path.basename(targetfile))
    seqf = os.path.abspath(seqfile)

    target = gdtcdp.Protein(name=tname)
    target.from_file(targetfile)
    target.add_tertiary_interactions(tertfile)
    target.addresidueids(seqf)

    decoys = makeprotdict(ddir, tertdir, target)

    # compute gdtg_ts for each decoy and store in a dict
    gdtgs = {}
    for dname, decoy in decoys.iteritems():
        l = [gdtcdp.computegdtg(target, decoy, c) for c in cutoffs]
        gdtg = float(sum(l))/float(len(l))
        gdtgs[dname] = gdtg

    if output:
        with open(os.path.abspath(output), 'w') as outf:
            for dname in gdtgs:
                outf.write('{}\t{}\n'.format(dname, gdtgs[dname]))
        return

    return gdtgs


def makeprotdict(protdir, tertdir, target):
    """
    Returns a dict containing {protein name: graph (gdtcdp.Protein
    object)}. The proteins are decoys corresponding to the given
    target.
    """
    res = {}
    gen = (f for f in os.listdir(protdir)
           if f.startswith(target.name))
    for f in gen:
        k, _ = os.path.splitext(f)
        prot = gdtcdp.Protein(name=k)
        prot.from_file(os.path.join(protdir, f))
        prot.add_tertiary_interactions(os.path.join(tertdir, f))
        prot.residueids = target.residueids
        res[k] = prot
    return res


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('target', help='A target file.')
    parser.add_argument('seqfile', help='A file containing '
                        'seqence length for all targets.')
    parser.add_argument('decoy_dir', help='Directory containing '
                        'decoy files.')
    parser.add_argument('tertiary_dir', help='Directory containing '
                        'tertiary interaction files.')
    parser.add_argument('--output', default=None,
                        help='Outpuf file. Default: None')
    args = parser.parse_args()
    run(args.target, args.seqfile,
        args.decoy_dir, args.tertiary_dir,
        args.output)
