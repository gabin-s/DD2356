matrix_multiply:
	cc -g -O2 -o matrix_multiply matrix_multiply.c

matrix_multiply_vect:
	cc -O2 -ftree-vectorize -msse2 -mfpmath=sse -fopt-info-vec-missed -o matrix_multiply_vect matrix_multiply.c

stream: stream.c
	cc -O2 -o stream stream.c

spmv: spmv.c
	cc -O2 -o spmv spmv.c

transpose:
	cc -O2 -o transpose transpose.c

# compile benchmark.c
benchmark: benchmark.c
	cc -o benchmark -O2 benchmark.c

# compile benchmark.c to asm
benchmark_asm: benchmark.c
	cc -S -fverbose-asm -O2 benchmark.c

clockgranularity: clockgranularity.c
	cc -o clockgranularity clockgranularity.c

<<<<<<< Updated upstream
transpose: transpose.c
	cc -O2 -o transpose transpose.c

=======
>>>>>>> Stashed changes
