# This is for testing optimized FFT routines....

all: realft realft_sparc_test1 realft_sparc_test2

realft: realft.c printvec.o
	cc -fast -o realft realft.c printvec.o -xlic_lib=sunperf\
           -L/usr/local/recipes/lib -lrecipes_c -lm 

realft_sparc_test1: realft_sparc_test.c printvec.o
	cc -fast -o realft_sparc_test1 realft_sparc_test.c printvec.o\
           -L/usr/local/recipes/lib -lrecipes_c -lm 

realft_sparc_test2: realft_sparc_test.c realft_sparc.o printvec.o
	cc -fast -o realft_sparc_test2 realft_sparc_test.c printvec.o realft_sparc.o\
           -xlic_lib=sunperf -L/usr/local/recipes/lib -lrecipes_c -lm 

realft_sparc.o: realft_sparc.c
	cc -fast  -c realft_sparc.c

printvec.o: printvec.c
	cc -fast  -c printvec.c
