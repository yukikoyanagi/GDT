#!/usr/bin/env python
#
# File: gdtcdp.py
#
# Time-stamp: <2016-12-19 14:37:05 au447708>
#
# Description: Wrapper for cdp.py code for GDT-Graph analysis.
# Class: Protein (inherits cdp.Protein)
#  Attributes:
#   residueids: (1-based) list of residue id's.
#  Methods:
#   subgraph():
#   addresidueids():
#   addedges():
#   bondexists():
#   toresidueids():
#   toatomids():
#   grow():
#   addbond():
#   bbedges():
#   d(): static method
# Computegdtg():
#
# Usage: Needs cdp.py in the same directory. Note cdp.py uses
# internally atom id's, which are 3-based, with the first N atom
# having index 3.
#
# Author: Yuki Koyanagi
# History:
#  2016-12-14: Created
#
import copy
import cdp


class Protein(cdp.Protein):

    def __init__(self, name):
        cdp.Protein.__init__(self, name)
        self.residueids = None

    def subgraph(self, residues):
        """
        Returns a subgraph (i.e. Protein instance) of this graph
        containing only the residues specified.
        """
        subg = Protein(name='{}_sub'.format(self.name))
        subg.residueids = residues
        subg.addedges(self)
        return subg

    def addresidueids(self, seqf):
        if not self.name:
            raise Exception("No name is defined for Protein.")
        name = self.name.split('_')[0]  # works also for decoys
        with open(seqf) as seqfile:
            for line in seqfile:
                if name in line:
                    maxidx = int(line.split()[1])
                    self.residueids = range(1, maxidx+1)

    def addedges(self, supgraph):
        """
        Copy Hbonds and Tbonds from supgraph, if both endpoints are
        in the current graph.
        """
        atomids = self.toatomids(self.residueids)
        for idx in set(atomids) & set(supgraph.vertices.keys()):
            bond = supgraph.vertices[idx]
            if (bond.other_end(idx) in atomids and
                    not self.bondexists(bond)):
                self.addbond(bond)

    def bondexists(self, bond):
        """
        Returns true if a given bond already exists in the graph
        """
        return (bond.donor in self.vertices and
                bond.accptr in self.vertices)

    def toatomids(self, resids):
        """
        Converts the given iterable of residue indices (1-based) to a
        list of atom indices (3(!?)-based)
        """
        atomids = [[i*3, i*3+1, i*3+2] for i in resids]
        return sorted(sum(atomids, []))  # flattens the list(!)

    def toresidueids(self, atomids):
        """
        Converts the given iterable of atom indices (3(!?)-based) to a
        list of residue indices (1-based)
        """
        resids = set([i/3 for i in atomids])  # de-duplicate
        return sorted(list(resids))

    def grow(self, supgraph):
        """
        Grow this graph by one edge inside supgraph
        """
        # Collect outgoing bonds from subgraph
        # Note de-duplication as 'outgoing' bond may terminate
        # inside subgraph.
        vs = list(set(supgraph.vertices) &
                  set(self.toatomids(self.residueids)) -
                  set(self.vertices))
        outbonds = set([supgraph.vertices[i] for i in vs])

        # Construct residueids post-growing
        l = [[i-1, i, i+1] for i in self.residueids]
        x = set(sum(l, []))  # flatten and de-duplicate
        y = x & set(supgraph.residueids)
        newresid = sorted(list(y))

        # Now we do the actual growing
        backup = copy.deepcopy(self)
        try:
            for bond in outbonds:
                self.addbond(bond)
            # Need to add new residues at the endpoints of
            # outgoing bonds
            l = (set(self.vertices.keys()) |
                 set(self.toatomids(newresid)))
            self.residueids = self.toresidueids(l)
        except:
            # Rollback
            self = backup
            raise

    def addbond(self, bond):
        """
        Adds given bond (T or H) to self
        """
        if isinstance(bond, cdp.Hbond):
            self.add_Hbond(bond.linenumber,
                           bond.donor,
                           bond.accptr,
                           bond.length,
                           bond.cluster,
                           bond.flags,
                           bond.residues,
                           bond.so3matrix)
        elif isinstance(bond, cdp.Tbond):
            self.add_tert(bond.linenumber,
                          bond.donor,
                          bond.accptr,
                          bond.d_VDW)

    def bbedges(self):
        """
        Returns a list of backbone-edges (i, i+1) constructed from
        residueids
        """
        return [(i, i+1) for i in self.residueids
                if i+1 in self.residueids]

    @staticmethod
    def d(target, decoy):
        """
        Returns the number of edges that need to be added/removed
        to obtain the target graph from the decoy.
        """
        # bonds are H- or Tbonds
        bonds = (set(target.Hbonds + target.Tbonds) ^
                 set(decoy.Hbonds + decoy.Tbonds))
        # edges are along the backbone
        edges = set(target.bbedges()) ^ set(decoy.bbedges())
        return len(bonds) + len(edges)


def computegdtg(target, decoy, limit):
    """
    Compute GDT-Graph score
    """
    assert isinstance(target, Protein)
    assert isinstance(decoy, Protein)
    l = []
    for i in target.residueids[:-3]:
        tsub = target.subgraph(range(i, i+3))
        dsub = decoy.subgraph(range(i, i+3))
        dist = Protein.d(tsub, dsub)
        if dist > limit:
            l.append(len(dsub.residueids))
            continue
        while True:
            res = len(dsub.residueids)
            tsub.grow(target)
            dsub.grow(decoy)
            if Protein.d(tsub, dsub) > limit:
                l.append(res)
                break
            if tsub.residueids == target.residueids:
                l.append(len(tsub.residueids))
                break
    return float(max(l))/float(len(target.residueids))*100
