#!/bin/bash

for i in $( seq 0 $1 )
do
	cat $PWD/solver_p1 > SOLVER_$i
	cat $PWD/${i}_interest.gams >> SOLVER_$i
	cat $PWD/solver_p2 >> SOLVER_$i
	cat $PWD/${i}_probability.gams >> SOLVER_$i
	cat $PWD/solver_p3 >> SOLVER_$i
done
