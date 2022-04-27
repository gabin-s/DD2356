#!/bin/sh

for task_number in 8 16 32 64 128; do
	echo -n "$task_number "
	srun -n $task_number ./pi_mpi | awk -F "=" "{ print $3}"
	[ $? -ne 0 ] && break
done;
