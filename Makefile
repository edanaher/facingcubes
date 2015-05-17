CC=gcc
CFLAGS=-Wall
OPTFLAGS=-O3 ${CFLAGS}

layout: layout.c
	${CC} layout.c -olayout ${OPTFLAGS}

# Compile flags:
# - DIMENSION: hardcode dimension instead of using ARGV[0]
# - RANDOMCHOOSE: at each step, choose this many options instead of trying all
# - COUNTINGDEPTH: count how many times each depth is hit, up to this depth.
# - ONLYCOUNTDEPTH: return once COUNTINGDEPTH depth is reached.
# - DISPLAYDEPTH: How deep to show progress
# - CACHEDEPTH: how deep to cache.
# - CACHESIZE: the size of the cache
# - CACHEMAPSIZE: the size of the hashmap indexing the cache

brute: brute.c
	${CC} brute.c -obrute ${OPTFLAGS}

brute4: brute.c
	${CC} brute.c -obrute4 ${OPTFLAGS} -DDIMENSION=4

brute5: brute.c
	${CC} brute.c -obrute5 ${OPTFLAGS} -DDIMENSION=5

brute6: brute.c
	${CC} brute.c -obrute6 ${OPTFLAGS} -DDIMENSION=6

random5: brute.c
	${CC} brute.c -orandom5 ${OPTFLAGS} -DDIMENSION=5 -DRANDOMCHOOSE=3

random6: brute.c
	${CC} brute.c -orandom6 ${OPTFLAGS} -DDIMENSION=6 -DRANDOMCHOOSE=2

profile: brute.c
	${CC} brute.c -obrute ${OPTFLAGS} -g -lprofiler

profile5: brute.c
	${CC} brute.c -obrute5 ${OPTFLAGS} -g -lprofiler -DDIMENSION=5

debug: brute.c
	${CC} brute.c -obrute -g ${CFLAGS}

benchmark: brute.c
	${CC} brute.c -obrute6 ${OPTFLAGS} -g -lprofiler -DDIMENSION=6 -DCOUNTINGDEPTH=14 -DCACHEDEPTH=6
