#!/bin/bash

if [ $# != 1 ]; then
	echo "lack of parameter"
	exit 1
fi

mkdir res_$1
mv *.log res_$1
cp sim.conf res_$1
mv ccdf.csv path.csv pdf.csv peer_pro.csv weight.csv res_$1
