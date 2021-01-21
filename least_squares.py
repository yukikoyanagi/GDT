#!/usr/bin/python
#
# File: least_squares.py
# Author: Rasmus Villemoes


# Input: List of rows
#
# y_1 x_11 x_12 x_13 ... x_1n
# y_2 x_21 x_22 x_23 ... x_2n
# y_3 x_31 x_32 x_33 ... x_3n
# ...
# y_m x_m1 x_m2 x_m3 ... x_mn
#
# Every row must contain n+1 columns (determined from the first line
# of input). This is treated as a (normally overdetermined) system of
# linear equations
#
# Y = Z A
#
# where Z is the mx(n+1) matrix consisting of a column of all 1s
# followed by the mxn matrix of x_ij input values, A = [a0 a1 ... an]
# is the n+1 (column) vector of unknowns, and Y is the m-column vector
# of input values.
#
# The output is the, usually unique, vector A such that |ZA-Y|^2 is
# minimal.
#
# TODO: Add option to print the input values, but with a few columns
# appended: Certainly the value which one would get by using the found
# a_j values, but maybe also the quadratic deviation.

import sys
import argparse
import fileinput
import numpy as np

from argparse import ArgumentParser
parser = ArgumentParser()

parser.add_argument("--outfile")

Y = []
Z = []
n = None

args, unk = parser.parse_known_args()

for line in fileinput.input(unk):
    vals = [float(x) for x in line.strip().split()]
    Y.append(vals[0])
    vals[0] = 1.0
    Z.append(vals)

Y = np.array(Y)
Z = np.array(Z)
# print Y
# print Z


A, res, rank, sing = np.linalg.lstsq(Z, Y)

for i,a in enumerate(A):
    print "a%d\t%f" % (i, a)

print "Sum of squared residuals: %f" % res
print "Rank of matrix: %d" % rank
# print sing


if args.outfile:
    f = open(args.outfile, "w")
    print len(Z)
    for r in range(len(Z)):
        row = Z[r,:]
        dotprod = np.dot(row, A)
        resid = Y[r] - dotprod
        f.write("\t".join([str(x) for x in (
            [Y[r]] + list(np.nditer(row)) + [dotprod, resid, resid*resid])]))
        
        # f.write("%f\t" % Y[r])
        # f.write("\t".join([str(x) for x in np.nditer(row)]))
        # f.write("\t%f" % dotprod)
        # f.write()
        f.write("\n")

