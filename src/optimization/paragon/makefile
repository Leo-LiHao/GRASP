# 	setup paragon first!

all: realft realft_par_test1 realft_par_test2 realft_par_test3
realft: realft.c printvec.o
	icc -o realft realft.c printvec.o -lkmath -L /home/sampson/ballen/recipes -lrecipes.PGON -lm 

realft_par_test1: realft_par_test.c realft_paragon.o printvec.o
	icc -o realft_par_test1 realft_par_test.c printvec.o realft_paragon.o -lkmath\
 -L /home/sampson/ballen/recipes -lrecipes.PGON -lm 

realft_par_test2: realft_par_test.c printvec.o
	icc -o realft_par_test2 realft_par_test.c printvec.o -L /home/sampson/ballen/recipes -lrecipes.PGON -lm 

realft_par_test3: realft_par_test.c realft_paragon_risky.o printvec.o
	icc -o realft_par_test3 realft_par_test.c printvec.o realft_paragon_risky.o -lkmath\
 -L /home/sampson/ballen/recipes -lrecipes.PGON -lm 

realft_paragon.o: realft_paragon.c
	icc -O4 -Mvect -Minfo -c realft_paragon.c

realft_paragon_risky.o: realft_paragon_risky.c
	icc -O4 -Mvect -Minfo -c realft_paragon_risky.c

printvec.o: printvec.c
	icc -O4 -Mvect -Minfo -c printvec.c
