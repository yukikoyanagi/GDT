#!/usr/bin/env python
#
# File: gdtgraph.py
#
# Time-stamp: <2016-11-29 14:53:08 au447708>
#
# Description: Calculate GDT-like score from graphs of decoy and
#  target.
#
# Author: Yuki Koyanagi
# History:
#

import gdtcdp


def compute(targetf, decoyf, tertiary=True,
            targettertf=None, decoytertf=None):
    """
    Compute GDT-Graph score for given target and decoy. targetf and
    decoyf are paths to the target and decoy files.
    """
    target = cdp.Protein.from_file(targetf)
    decoy = cdp.Protein.from_file(decoyf)
    if tertiary:
        assert targettertf and decoytertf
        target.add_tertiary_interactions(targettertf)
        decoy.add_tertiary_interactions(decoytertf)
    pass


def d(t, d):
    """
    Compute distance-measure between the given subgraphs.
    The measure is the number of edges that must be removed/added
    to obtain subgraph t from d.
    """
    pass


def getresidlen(prot, seqf):
    """
    Get residue length for target prot from seqf file.
    """
    pass
