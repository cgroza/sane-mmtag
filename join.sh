#!/bin/bash

# This script creates the table for input into lift_5mC.py

# gziped files
GRAPH_ALN=$1
BMOD=$1
join -t $'\t' -1 1 -2 1 <(gunzip -c ${GRAPH_ALN} | sort ) <(gunzip -c ${BMOD} | sort ) | gzip

