#!/bin/sh

infile=$1
outfile=$2
gapopen=$3

clustalw2 -infile=${infile} -align -type=protein -ktuple=1 -window=5 -score=percent -topdiags=5 -pairgap=3 -pwmatrix=gonnet -pwdnamatrix=iub -pwgapopen=10 -pwgapext=0.1 -matrix=gonnet -dnamatrix=iub -gapopen=${gapopen} -gapext=0.2 -gapdist=5 -iteration=none -numiter=1 -clustering=NJ -output=fasta -outorder=input -outfile=${outfile}
