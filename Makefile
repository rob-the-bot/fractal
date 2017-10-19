output: fract.c
	mpicc fract.c -Wall -lm -lgmp -lmpfr -lmpc -fopenmp -o fract -O3

clean:
	rm -r -f *.png fract
