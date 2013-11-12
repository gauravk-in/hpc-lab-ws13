#!/bin/sh

for test in serial_pi omp_reduction_pi omp_critical_pi
do
	for size in 1000000 2000000 4000000 8000000
	do
		echo $test" : Problem Size "$size
		for i in {1..4}
		do	
			./$test $size
		done
		echo ""
	done
done	


