#!/usr/bin/env bash

json=$1; shift
csv=$(echo $json | sed "s/json$/csv/")

cat $json | \
egrep -v "{|}" | \
sed "s/\"//g; s/[[:space:]]//g; s/,//; s/:/,/" >$csv

