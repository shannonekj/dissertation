#!/bin/bash

set -x

infile=$1

sed -n 1p ${infile}  > ${infile}.no0.sort
grep -v $'\t''0'$'\t''0'$'\t''0'$'\t''0'$'\t''0'$'\t''0'$'\t''0'$'\t''0' ${infile} > ${infile}.no0
sort -nrk 9 ${infile}.no0 >> ${infile}.no0.sort
sed -i '$ d' ${infile}.no0.sort
rm ${infile}.no0
lzma ${infile}
