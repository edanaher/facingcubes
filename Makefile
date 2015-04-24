CC=gcc
CFLAGS=-Wall
OPTFLAGS=-O3

brute: brute.c
	${CC} brute.c -obrute ${OPTFLAGS} ${CFLAGS}

brute4: brute.c
	${CC} brute.c -obrute4 ${OPTFLAGS} ${CFLAGS} -DDIMENSION=4

brute5: brute.c
	${CC} brute.c -obrute5 ${OPTFLAGS} ${CFLAGS} -DDIMENSION=5

brute6: brute.c
	${CC} brute.c -obrute5 ${OPTFLAGS} ${CFLAGS} -DDIMENSION=5

profile: brute.c
	${CC} brute.c -obrute ${OPTFLAGS} -g -lprofiler

debug: brute.c
	${CC} brute.c -obrute -g ${CFLAGS}

