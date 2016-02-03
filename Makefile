CC=gcc
CFLAGS=-Wall
OPTFLAGS=-O3 ${CFLAGS}


# Compile flags:
# - DIMENSION: hardcode dimension instead of using ARGV[0]
# - TIMELIMIT: how many seconds to spend on each case before bailing
# - CACHESIZE: how much space to use for the cache
# - CACHEMAPSIZE: how many elements to put in the cachemap
# - CACHEDEPTH: How many layers of the search tree to cache
# - FASTHASH: Use a less precise but faster hash function
# - DISPLAYDEPTH: How deep to go for display purposes
# - KEEPPROGRESS: Keep this many lines on stderr for later information
# - BIRTHDAYHASH: Store 1 << BIRTHDAYHASH hashes for each value to reduce lookup time by the same value.a
#   						  A straightforward time/space tradeoff.

layout: layout.c buildlayout.c
	${CC} layout.c -olayout ${OPTFLAGS}

layout5: layout.c buildlayout.c
	${CC} layout.c -olayout5 ${OPTFLAGS} -DDIMENSION=5

layout6: layout.c buildlayout.c
	${CC} layout.c -olayout6 ${OPTFLAGS} -DDIMENSION=6 -DCACHEDEPTH=6

layout6stable: layout.c buildlayout.c
	${CC} layout.c -olayout6stable ${OPTFLAGS} -DDIMENSION=6 -DCACHEDEPTH=6 -DKEEPPROGRESS=20 -DDISPLAYDEPTH=5 -DBIRTHDAYHASH=1

layout6timed: layout.c buildlayout.c
	${CC} layout.c -olayout6timed ${OPTFLAGS} -DDIMENSION=6 -DTIMELIMIT=600 -DCACHEDEPTH=4

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
