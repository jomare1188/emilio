#!/bin/bash

#$ -cwd
#$ -q all.q
#$ -pe smp 1

module load mcl/14-137

in=`ls -1 *_top20CV_mcl.txt`
out=${in/.txt}_I2.0.txt

/usr/bin/time -v mcl $in -I 2 -te 1 --abc -o $out

