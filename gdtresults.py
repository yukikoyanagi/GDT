#!/usr/bin/env python3
#
# File: gdtresults.py
#
# Time-stamp: <2016-10-18 15:44:03 au447708>
#
# Description: Methods for processing text files produced by analysegdt
#  script.
# Methods:
#  loaddata: Load data from specified txt file.
#
# Author: Yuki Koyanagi
# History:
#


from collections import namedtuple


def tofile(lst, f):
    """
    Write list of Record namedtuples to a file f.
    """
    l = ['\t'.join(['a values', 'n values', '2% for all',
                    '10% for all', '2% for >40',
                    '10% for >40'])]
    for rec in lst:
        l.append('\t'.join([str(item) for item in rec]))
    with open(f, 'w') as outf:
        outf.write('\n'.join(l))


def loaddata(f):
    """
    Load data from txt file produced by analysegdt script.
    Returns a list of namedtuples: (avals, nvals, 2%_all, 10%_all,
    2%_>40, 10%_>40)
    """
    Record = namedtuple('Record', 'avals, nvals, pct2_all, '
                        'pct10_all, pct2_40, pct10_40')
    with open(f) as initf:
        for line in initf:
            if line.startswith('Results for parameter'):
                av, nv = parseparams(next(initf))
            elif line.startswith('Guesses within'):
                p2a = parsepct(next(initf))
                p10a = parsepct(next(initf))
            elif line.startswith('for max. GDT'):
                p240 = parsepct(next(initf))
                p1040 = parsepct(next(initf))
                r = Record(av, nv, p2a, p10a, p240, p1040)
                try:
                    recs.append(r)
                except NameError:
                    recs = [r]
    return recs


def parseparams(s):
    """
    Parse param string in the form: a=2,3;n=1,2,3
    Returns 2 lists avals and nvals
    """
    for p in s.strip().split(';'):
        if p.split('=')[0] == 'a':
            av = [int(v) for v in p.split('=')[1].split(',')]
        elif p.split('=')[0] == 'n':
            nv = [int(v) for v in p.split('=')[1].split(',')]
    return av, nv


def parsepct(s):
    """
    Parse string of the form: 2%-pt\t25(28.45%).
    Returns a float of pct value.
    """
    return float(s.split()[1].split('(')[1][:-2])
