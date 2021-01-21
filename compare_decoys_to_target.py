#!/usr/bin/python
# coding: utf-8

# File: compare_decoys_to_target.py
# Author: Rasmus Villemoes


import sys
import os
import platform
import itertools
import math

import string
import gc
import re

sys.path.append("/home/qgm/QGM/cdp")

from cdp import *

from argparse import ArgumentParser
parser = ArgumentParser()

parser.add_argument("files", nargs='*')
parser.add_argument("--decoy-dir", type=str,
                    help="directory containing decoys corresponding to the given targets")

parser.add_argument("--tert-dir", type=str,
                    help="directory containing tertiary information for targets and decoys")

parser.add_argument("--gdt-file", type=str,
                    help="file with gdt scores for decoys")

parser.add_argument("--const", type=float,
                    help="constant for use in bond_set_score")

parser.add_argument("--nconst", type=int,
                    help="another constant for use in bond_set_score")

parser.add_argument("--Calpha-max", type=float, default=float('inf'),
                    help="cutoff for Calpha distance, default infinity")

parser.add_argument("--Cbeta-max", type=float, default=float('inf'),
                    help="cutoff for Cbeta distance, default infinity")



args = parser.parse_args()
decoy_dir = args.decoy_dir
if decoy_dir is None:
    decoy_dir = "/home/qgm/QGM/GDT/data/casp10/decoys/split"

tert_dir = args.tert_dir
if tert_dir is None:
    tert_dir = "/home/qgm/QGM/tertiary/data/20151102/casp10"


gdt_file = args.gdt_file
if gdt_file is None:
    gdt_file = "/home/qgm/QGM/GDT/data/casp10/summaries/gdt_ts.txt"

GDT = dict()
for l in open(gdt_file):
    d, gdt = l.strip().split()
    GDT[d] = float(gdt)

const = args.const
if const is None:
    const = 0.0

nconst = args.nconst
if nconst is None:
    nconst = 5

Calpha_max = args.Calpha_max
Cbeta_max = args.Cbeta_max


def pdist(a, b):
    return abs(a[0]-b[0]) + abs(a[1]+b[1])


# This is NOT a distance. Larger means better. So it is more like 'score'.
def bond_set_score(b, S):
    def f(x):
        x = abs(x)
        if x > 2*nconst:
            return -1
        return 1 - float(x)/nconst

    if len(S) == 0:
        return -1
    return max([f(b[0]-x[0]) + f(b[1]-x[1]) for x in S])
# Kunne prøve at ændre d(h,D) til følgende

# d(h,D) = max_{h'\in D} f(|i-j|) + f(|k-l|)

# hvor h er fra i til k og h' er fra j til l og f er en linear funktion fra [0, 2n], som er 1 i nul, og -1 i 2n, hvor vi kan varierer n, lad os sige fra 1 til 10 og for argumenter mindre en 2n er den konstant -1. Hvis D er tomt skal den bare være -1.

    # Kludge due to the scoring function not being well-defined when
    # the set of decoy bonds is a subset of the target bonds, or vice
    # versa.
    if len(S) == 0:
        return 0
    dist = min([pdist(b, x) for x in S])
    return math.exp(-const * dist)


def set_scores(target_set, decoy_set):
    common_set = target_set & decoy_set

    if len(target_set) == 0:
        P1 = 100  # kludge, but what else should we set 0/0 to?
    else:
        P1 = float(len(common_set))/float(len(target_set)) * 100

    S2 = 0
    S2 += sum([bond_set_score(b, decoy_set - common_set) for b in target_set - common_set])
    S2 += sum([bond_set_score(b, target_set - common_set) for b in decoy_set - common_set])
    if (len(target_set - common_set) + len(decoy_set - common_set)):
        S2 /= (len(target_set - common_set) + len(decoy_set - common_set))
    else:
        assert(S2 == 0)

    S = P1 + S2*(100-P1)
    return len(target_set), len(decoy_set), P1, S2, S


decoy_files = sorted(os.listdir(decoy_dir))

for f in args.files:
    base, ext = os.path.splitext(os.path.basename(f))
    target = Protein(name=base)
    target.from_file(f)

    if tert_dir:
        tfile = tert_dir + "/" + os.path.basename(f)
        try:
            target.add_tertiary_interactions(tfile, Calpha_max_dist = Calpha_max, Cbeta_max_dist = Cbeta_max)
        except Exception as e:
            # sys.stderr.write("Failed to add tertiary interaction info for %s: %s\n" % (base, e))
            pass

    T_Hbonds = set([(hb.donor, hb.accptr) for hb in target.Hbonds])
    T_Tbonds = set([(tb.donor, tb.accptr) for tb in target.Tbonds])
    # print T_Tbonds

    for d in decoy_files:
        if d.endswith(".txt") and d.startswith(base):
            dbase, dext = os.path.splitext(d)
            decoy = Protein(name=dbase)
            decoy.from_file(decoy_dir + "/" + d)

            if dbase in GDT:
                gdt = GDT[dbase]
            else:
                gdt = float('nan')

            if tert_dir:
                tfile = tert_dir + "/" + os.path.basename(d)
                try:
                    decoy.add_tertiary_interactions(
                        tfile,
                        Calpha_max_dist=Calpha_max,
                        Cbeta_max_dist=Cbeta_max)
                except Exception as e:
                    # sys.stderr.write("Failed to add tertiary interaction info for %s: %s\n" % (dbase, e))
                    pass

            D_Hbonds = set([(hb.donor, hb.accptr) for hb in decoy.Hbonds])
            D_Tbonds = set([(tb.donor, tb.accptr) for tb in decoy.Tbonds])

            T_Hcount, D_Hcount, H_P1, H_S2, H_S = set_scores(T_Hbonds, D_Hbonds)
            T_Tcount, D_Tcount, T_P1, T_S2, T_S = set_scores(T_Tbonds, D_Tbonds)


            columns = [base, dbase, str(gdt)]
            columns += [str(x) for x in [T_Hcount, D_Hcount, H_P1, H_S2, H_S]]
            columns += [str(x) for x in [T_Tcount, D_Tcount, T_P1, T_S2, T_S]]

            print "\t".join(columns)
