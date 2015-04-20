brute: brute.c
	gcc brute.c -obrute -O3 -Wall

debug: brute.c
	gcc brute.c -obrute -g -Wall

profile: brute.c
	gcc brute.c -obrute -O3 -g -lprofiler
