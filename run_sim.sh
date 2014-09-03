#!/bin/bash

if [ $# != 1 ]; then
	echo "input the name for sim results"
	exit 1
fi

for i in $( seq 0 9 )
do
	cd ./sim$i
	./Sim ./graph.csv $i $1 &
	cd ..
done
