#!/bin/bash

for i in $( seq 0 9 )
do
	mkdir ./sim$i
	cp Sim sim.conf mkres.sh ./sim$i
	ln -s ${PWD}/mobility.csv ./sim$i/mobility.csv
	ln -s ${PWD}/graph.csv ./sim$i/graph.csv
done
