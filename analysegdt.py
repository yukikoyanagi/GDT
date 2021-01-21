#!/usr/bin/env python3
#
# File: analysegdt.py
#
# Time-stamp: <2016-10-24 16:36:16 yuki>
#
# Description: Analyse GDT parameter combinations.
# Usage: ./analysegdt.py [parameter file] -d [path to data dir]
# ToDo: Specify output file name.
#
# Author: Yuki Koyanagi
# History:
#

import gdtprocessor
import argparse
from collections import namedtuple
from operator import itemgetter

Param = namedtuple('Param', 'avals, nvals')


def parseparams(p):
    """
    Parse param file and return a list of named tuples
    (avals, nvals). Each item in a tuple is a list of param values.
    """
    with open(p) as pfile:
        for line in pfile:
            for s in line.strip().split(';'):
                if s.split('=')[0] == 'avals':
                    av = [int(v) for v in s.split('=')[1].split(',')]
                elif s.split('=')[0] == 'nvals':
                    nv = [int(v) for v in s.split('=')[1].split(',')]
            param = Param(avals=av, nvals=nv)
            try:
                params.append(param)
            except NameError:
                params = [param]
    return params


def filestoload(params):
    """
    Given list of params, return a namedtuple of parameter
    values to load.
    """
    for param in params:
        try:
            aset = aset.union(set(param.avals))
        except NameError:
            aset = set(param.avals)
        try:
            nset = nset.union(set(param.nvals))
        except NameError:
            nset = set(param.nvals)
    av = sorted(list(aset))
    nv = sorted(list(nset))
    return Param(avals=av, nvals=nv)


def getcount(lst, pt):
    """
    Get no. of items in the lst less than pt.
    Returns count and %
    """
    n = len([l for l in lst if l[-1] < float(pt)])
    return n, float(n)/float(len(lst))*100


def tofile(t, dir):
    transtable = str.maketrans('', '', '[] ')
    for r in t:
        try:
            s.append('#'*40)
        except NameError:
            s = ['#'*40]
        s.append('Results for parameter combination:')
        s.append('a={};n={}'.format(str(r[0]).translate(transtable),
                                    str(r[1]).translate(transtable)))
        s.append('')
        for i, c in enumerate(r[2][0].tolist()):
            s.append('a{}\t{}'.format(i, c))
        s.append('')
        s.append('Sum of squared residuals: '
                 '{}'.format(r[2][1].tolist()))
        s.append('Rank of matrix: {}'.format(r[2][2]))
        s.append('')
        s.append('Deviation from max GDT')
        for row in r[3]:
            s.append('\t'.join(['{:.2f}'.format(item)
                                if isinstance(item, float)
                                else item.decode('utf-8')
                                for item in row]))
        s.append('')
        s.append('Guesses within given %-pt of max. GDT')
        s.append('2%-pt\t{}({:.2f}%)'.format(
            *getcount(r[3], 2)))
        s.append('10%-pt\t{}({:.2f}%)'.format(
            *getcount(r[3], 10)))
        s.append('')
        s.append('for max. GDT > 40 only')
        s.append('2%-pt\t{}({:.2f}%)'.format(
            *getcount([row for row in r[3] if row[2] > 40], 2)))
        s.append('10%-pt\t{}({:.2f}%)'.format(
            *getcount([row for row in r[3] if row[2] > 40], 10)))
        s.append('')
    with open(dir + '/result.txt', 'w') as f:
        f.write('\n'.join(s))


def run(p, d):
    params = parseparams(p)
    uparam = filestoload(params)
    gp = gdtprocessor.Gdtprocessor(dir=d,
                                   arange=uparam.avals,
                                   nrange=uparam.nvals)
    for avals, nvals in params:
        result, prediction = gp.runlstsq(avals, nvals)
        evaluation = sorted(gp.evaluate(prediction),
                            key=itemgetter(0))
        try:
            t.append((avals, nvals, result, evaluation))
        except NameError:
            t = [(avals, nvals, result, evaluation)]
    tofile(t, d)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('paramfile', help='File containing '
                        'list of parameters avals and nvals. '
                        'These must be of the form; '
                        'avals=5,6,7;nvals=1,2,3,4 '
                        'One param combination per line.')
    parser.add_argument('-d', '--data_dir',
                        help='Directory containing '
                        'data_a#_n#.txt files.',
                        default='/home/qgm/QGM/GDT/out/linfit')
    args = parser.parse_args()
    run(args.paramfile, args.data_dir)
