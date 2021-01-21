#!/usr/bin/env python3
#
# File: gdtprocessor.py
#
# Time-stamp: <2016-11-08 12:11:46 au447708>
#
# Description: Class for processing GDT 'scores'.
# Methods:
#  __init__: Load data from data_a#_n#.txt files to create
#  np.ndarray object
#  toidx: Converts parameter pair (a,n) to index along the 3rd axis
#  toparam: Converts index along the 3rd axis to parameter pair (a,n)
#  to2d: Converts input vector to column vector with ndim=2
#  extract: Extract data from array for least sq analysis
#  runlstsq: Runs least square analysis
#  evaluate: Evaluate the results from least square analysis
#
# Author: Yuki Koyanagi
# History:
#  2016-10-14 (yk): Created
#  2016-11-08 (yk): Added noarr option to init


import numpy as np
from itertools import groupby
from operator import itemgetter

# Columns to use in data_a#_n#.txt files
cols = (2, 3, 4, 5, 6)
namecols = (0, 1)
# Column in prediction array containing GDT score
gdtcol = 1
# Column in prediction array containing GDT guess
guesscol = -3  # 3rd last


class Gdtprocessor:

    def __init__(self, dir=None,
                 arange=range(4, 13), nrange=range(1, 11),
                 noarr=False):
        if dir is None:
            dir = '.'
        self._arange = arange
        self._nrange = nrange
        if noarr:
            return  # Creation without array
        for idx in [(a, n) for a in arange for n in nrange]:
            fn = '{}/data_a{}_n{}.txt'.format(dir, idx[0], idx[1])
            try:
                self._arr = np.dstack((self._arr,
                                      np.loadtxt(fn, usecols=cols)))
            except AttributeError:
                self._arr = np.loadtxt(fn, usecols=cols)
                self._prot = np.loadtxt(fn, dtype='a8, a15',
                                        usecols=namecols)

    def paramtoidx(self, a, n):
        """
        Each data_a#_n# is stacked along the 3rd axis in the order
        (a1,n1),(a1,n2),(a1,n3),...,(a2,n1),(a2,n2),...,(a_max,n_max)
        This method converts parameter pair a,n to index along the
        3rd axis
        """
        return (len(self._nrange) * self._arange.index(a) +
                self._nrange.index(n))

    def idxtoparam(self, i):
        """
        Converts index along the 3rd axis to a parameter pair (a,n)
        """
        return (self._arange.index(i // len(self._nrange)),
                self._nrange.index(i % len(self._nrange)))

    def to2d(self, A):
        """
        Converts the input vector to a column vector with ndim=2.
        """
        assert isinstance(A, np.ndarray)
        if A.ndim == 1:
            A = np.atleast_2d(A).T
        return A

    def extract(self, avals, nvals):
        """
        Extracts n-vector of dependent variable (Y) and
        nxm matrix of independent variables (Z) from
        self._arr array.
        """
        nrows = self._arr.shape[0]
        ncols = (1 + 1 + len(nvals) + len(avals) +
                 (len(avals) * len(nvals)))
        A = np.ones((nrows, 1))
        B = self._arr[:, 1, [0]]
        C = self._arr[:, 2,
                      [self.paramtoidx(avals[0], n) for n in nvals]]
        D = self._arr[:, 3,
                      [self.paramtoidx(a, nvals[0]) for a in avals]]
        E = self._arr[:, 4,
                      [self.paramtoidx(a, n)
                       for a in avals
                       for n in nvals]]
#        print("\n".join([str(s) for s in
#                         [A.shape, B.shape, C.shape, D.shape, E.shape]]))
        Z = np.hstack((A, B, C, D, E))
        assert Z.shape[1] == ncols
        Y = self._arr[:, 0, 0]
        return Y, Z

    def runlstsq(self, avals, nvals):
        """
        Performs least square computation using dependent variables
        given in parameter combination (avals, nvals).
        Return values:
        result: tuple (solution of lstsq, sum of sqr of residuals,
        rank of solution matrix, singular values of ind. var)
        prediction: np.ndarray in the format of out_all.txt file.
        colmns are;
        1: GDT score
        2-(n+1): dependent variables (in matrix Z)
        -3 (3rd last): predicted GDT value
        -2 (2nd last): GDT - predicted
        -1 (last): (GDT - predicted)^2
        """
        self._avals = avals
        self._nvals = nvals
        Y, Z = self.extract(self._avals, self._nvals)
        try:
            A, res, rank, sing = np.linalg.lstsq(Z, Y)
        except np.LinAlgError:
            print('Least square computation does not converge.\n')
            print('a values:{}, n values:{}\n'.format(avals, nvals))
            return
        result = (A, res, rank, sing)
        pred = np.dot(Z, A)
        resid = Y - pred
        resid2 = resid**2
        p = (self.to2d(A) for A in (Y, Z, pred, resid, resid2))
        prediction = np.hstack(p)
        return result, prediction

    def evaluate(self, prediction):
        """
        Evaluate prediction results. Returns a list of tuples, one
        tuple per protein target. Tuple entries are
        (target, 'winner' decoy, gdt of 'winner' decoy,
        'guess' decoy, gdt of 'guess' decoy,
        winner gdt - guess gdt)
        """
        d = {}
        allpred = [name + tuple(pred)
                   for name, pred
                   in zip(self._prot.tolist(), prediction.tolist())]
        for prot, decoys in groupby(allpred, lambda x: x[0]):
            d[prot] = [decoy[1:] for decoy in decoys]
        for prot in d:
            wdecoy = max(d[prot], key=itemgetter(gdtcol))
            gdecoy = max(d[prot], key=itemgetter(guesscol))
            difprot = ((prot,) +
                       wdecoy[0:2] +
                       gdecoy[0:2] +
                       (wdecoy[1] - gdecoy[1],))
            try:
                gdtdif.append(difprot)
            except NameError:
                gdtdif = [difprot]
        return gdtdif
