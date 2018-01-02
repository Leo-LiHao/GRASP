# This is for testing optimized FFT routines....

all: realft realft_sgi_test1 realft_sgi_test2

realft: realft.c printvec.o
	cc -mips4 -n32 -o realft realft.c printvec.o -L/home/iris132/mnp/NumericalRecipes/lib -lrecipes_cg -lm -lcomplib.sgimath

realft_sgi_test1: realft_sgi_test.c printvec.o
	cc -mips4 -n32 -o realft_sgi_test1 realft_sgi_test.c printvec.o\
 -L/home/iris132/mnp/NumericalRecipes/lib -lrecipes_cg -lm

realft_sgi_test2: realft_sgi_test.c realft_sgi.o printvec.o
	cc -mips4 -n32 -o realft_sgi_test2 realft_sgi_test.c printvec.o realft_sgi.o\
 -L/home/iris132/mnp/NumericalRecipes/lib -lrecipes_cg -lm -lcomplib.sgimath

realft_sgi.o:
	cc -mips4 -n32 -c realft_sgi.c

printvec.o: printvec.c
	cc -mips4 -n32 -c printvec.c
