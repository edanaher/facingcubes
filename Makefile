brute: brute.c
	gcc brute.c -obrute -O2 -Wall

debug: brute.c
	gcc brute.c -obrute -g -Wall

profile: brute.c
	gcc brute.c -obrute -O2 -g -lprofiler
