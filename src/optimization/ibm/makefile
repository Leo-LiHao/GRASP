#  Makefile for IBM ESSL optimized subroutine library

CC = cc -O3 -qarch=pwrx -qstrict
RECIPES=-L/home/cacr/staff/roy/recipes -lrecipes.RS6K
ESSL=-lesslp2 -lm

all: realft realft_ibm_test1 realft_ibm_test2

# A useful routine for printing out arrays
printvec.o: printvec.c
	$(CC) -c printvec.c

# Compares the Numerical Recipes realft() output format to the ESSL equivalent
realft: realft.c printvec.o
	$(CC) -o realft realft.c printvec.o $(ESSL) $(RECIPES) -lm 

# Replacement for Numerical Recipes realft() that uses ESSL routines
realft_ibm.o: realft_ibm.c
	$(CC) -c realft_ibm.c

# A test program that does 1 forward and 1 reverse FFT using our ESSL based replacement for realft()
realft_ibm_test1: realft_ibm_test.c realft_ibm.o printvec.o
	$(CC) -o realft_ibm_test1 realft_ibm_test.c printvec.o realft_ibm.o $(ESSL) $(RECIPES) -lm

# A test program that does 1 forward and 1 reverse FFT using Numerical Recipes realft()
realft_ibm_test2: realft_ibm_test.c printvec.o
	$(CC) -o realft_ibm_test2 realft_ibm_test.c printvec.o $(ESSL) $(RECIPES) -lm
