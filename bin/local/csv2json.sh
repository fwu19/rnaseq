#!/usr/bin/env bash

csv=$1; shift
json=$(echo $csv | sed "s/csv$/json/")

(echo "{"; 
cat $csv | awk -F "," -v q='"' -v q2=':' -v q3=',' 'BEGIN {OFS=""} {print "\t"q,$1,q,q2," ",q,$2,q,q3}';
echo "}"
) >$json