# This is for building/testing optimized FFT routines for DEC
# machines using the DEC extended math library.

# The name of the C compiler:
CC = cc

# Flags for the C compiler
# CFLAGS = -std1 -fast -tune host -migrate 
CFLAGS = -std1 -O4 -tune host -migrate 

# The location and name of the numerical recipes library (used for comparison)
# NR_LIB = -L/usr/local/recipes/lib -lrecipes_c
NR_LIB = -L/usr/local/recipes/recipes_c-ansi/lib -lrecipes_c

all: realft realft_dec_test1 realft_dec_test2 realft_dec_time1 realft_dec_time2

# This is a handy helper routine for printout out arrays:
printvec.o: printvec.c
	$(CC) $(CFLAGS) -c printvec.c

# This routine prints out a comparison of NR and DXML FFT conventions:
realft: realft.c printvec.o
	$(CC) $(CFLAGS) -o realft realft.c printvec.o\
		       $(NR_LIB) -lm -ldxml

# This is the optimized replacement for the NR real FFT:
realft_dec.o: realft_dec.c
	$(CC) $(CFLAGS) -c realft_dec.c

# Comment out the previous two lines and uncomment these two if you
# want to use the "one step" FFT routine rather than the three-step one.
# In general the three step should be faster:
#realft_dec.o: realft_dec_1step.c
#	$(CC) $(CFLAGS) -c realft_dec_1step.c

# This prints out the FFT of a given array, using NR:
realft_dec_test1: realft_dec_test.c printvec.o
	$(CC) $(CFLAGS) -o realft_dec_test1 realft_dec_test.c printvec.o\
                       $(NR_LIB) -lm

# This prints out the FFT of a given array, using DXML:
realft_dec_test2: realft_dec_test.c realft_dec.o printvec.o
	$(CC) $(CFLAGS) -o realft_dec_test2 realft_dec_test.c printvec.o\
                        realft_dec.o $(NR_LIB) -lm -ldxml

# This does a timing test, using NR:
realft_dec_time1: realft_dec_time.c printvec.o
	$(CC) $(CFLAGS) -o realft_dec_time1 realft_dec_time.c printvec.o\
                       $(NR_LIB) -lm

# This does a timing test, using DXML:
realft_dec_time2: realft_dec_time.c realft_dec.o printvec.o
	$(CC) $(CFLAGS) -o realft_dec_time2 realft_dec_time.c printvec.o\
                        realft_dec.o $(NR_LIB) -lm -ldxml
