#!/usr/bin/env python
#
# File: limval.py
#
# Time-stamp: <2016-12-22 11:37:42 yuki>
#
# Description: Script to analyse the behaviour of gdt-g score w.r.t.
# cutoff values. Dir's are set for use on spencer.
#
# Author: Yuki Koyanagi
# History:
#

import os
import argparse
import gdtcdp

tdir = os.path.abspath('data/casp10/targets')
#ddir = os.path.abspath('data/casp10/decoytest')
ddir = os.path.abspath('data/casp10/decoysample')
tertdir = os.path.abspath('../tertiary/data/20151102/casp10')
seqf = 'casp10.seqlen.txt'


def addprots(d, tertd, seqf):
    res = []
    for f in os.listdir(d):
        n, _ = os.path.splitext(f)
        prot = gdtcdp.Protein(name=n)
        prot.from_file(os.path.join(d, f))
        prot.add_tertiary_interactions(os.path.join(tertd, f))
        prot.addresidueids(seqf)
        res.append(prot)
    return res


def getgdts(tprots, dprots, lim, out=None):
    res = {}
    for decoy in dprots:
        dname = decoy.name
        for prot in tprots:
            if prot.name == decoy.name.split('_')[0]:
                target = prot
        res[decoy.name] = gdtcdp.computegdtg(target, decoy, lim)
    if out:
        with open(out, 'w') as outf:
            outf.write('\n'.join(['{}\t{}'.format(d, res[d])
                                  for d in res]))
        return
    for d in res:
        print '{}\t{}'.format(d, res[d])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('lim', type=int)
    parser.add_argument('-o', '--output')
    args = parser.parse_args()
    
    tprots = addprots(tdir, tertdir, seqf)
    dprots = addprots(ddir, tertdir, seqf)

    getgdts(tprots, dprots, args.lim, args.output)
