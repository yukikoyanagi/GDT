#!/usr/bin/env python
#
# File: selectdecoys.py
#
# Time-stamp: <2017-01-04 14:58:40 au447708>
#
# Author: Yuki Koyanagi
# History:
#
import argparse

def run(sdir, ddir, tsfile):
    decoys = {}
    with open(tsfile) as tsf:
        for line in tsf:

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
