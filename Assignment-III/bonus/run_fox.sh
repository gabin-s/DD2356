#!/bin/bash

for i in $(seq 4 16); do
	let 'nproc = i * i'
	let 'M = i * (2048 / i)'
	
	echo nproc=$nproc M=$M

	srun -n $nproc ./fox $M
	[ $? -ne 0 ] && break
done;