CC=gcc
CFLAGS=-Wall
OPTFLAGS=-O3 ${CFLAGS}

brute: brute.c
	${CC} brute.c -obrute ${OPTFLAGS}

brute4: brute.c
	${CC} brute.c -obrute4 ${OPTFLAGS} -DDIMENSION=4

brute5: brute.c
	${CC} brute.c -obrute5 ${OPTFLAGS} -DDIMENSION=5

brute6: brute.c
	${CC} brute.c -obrute5 ${OPTFLAGS} -DDIMENSION=5

random6: brute.c
	${CC} brute.c -orandom6 ${OPTFLAGS} -DDIMENSION=6 -DRANDOMSKIP=0.3

profile: brute.c
	${CC} brute.c -obrute ${OPTFLAGS} -g -lprofiler

profile5: brute.c
	${CC} brute.c -obrute5 ${OPTFLAGS} -g -lprofiler -DDIMENSION=5

debug: brute.c
	${CC} brute.c -obrute -g ${CFLAGS}

