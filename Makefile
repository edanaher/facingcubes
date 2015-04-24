brute: brute.c
	gcc brute.c -obrute -O3 -Wall

brute4: brute.c
	gcc brute.c -obrute4 -O3 -Wall -DDIMENSION=4

brute5: brute.c
	gcc brute.c -obrute5 -O3 -Wall -DDIMENSION=5

brute6: brute.c
	gcc brute.c -obrute5 -O3 -Wall -DDIMENSION=5

debug: brute.c
	gcc brute.c -obrute -g -Wall

profile: brute.c
	gcc brute.c -obrute -O3 -g -lprofiler
