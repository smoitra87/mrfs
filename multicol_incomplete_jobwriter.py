#! /usr/bin/env python


from itertools import product, repeat
import os
import sys
import pickle
from operator import itemgetter
from os.path import split,splitext
import glob

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description="Create Regression jobs")
    parser.add_argument('jobf',type=str)
    parser.add_argument('--startline',type=int,default=1)
    parser.add_argument('--endline',type=int,default=None)
    parser.add_argument('--outf_idx',type=int,default=None, help="Out field idx")
    args = parser.parse_args()

# Read jobfile list
    with open(args.jobf) as fin:
        old_jobs = [line.strip() for line in fin]

    start = args.startline
    end = len(old_jobs) if not args.endline else args.endline

    old_jobs = [j for (idx,j) in enumerate(old_jobs) if idx+1 in range(start,end+1)]


    completed_jobs = []
    for old_job in old_jobs:
        outf = old_job.strip("()").split(",")[args.outf_idx].strip("'")
        if os.path.exists(outf):
            completed_jobs.append(old_job)

    missing_jobs = [j for j in old_jobs if j not in completed_jobs]
    for j in missing_jobs:
        print >>sys.stdout, j


