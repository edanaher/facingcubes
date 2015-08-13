#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>

#ifdef DIMENSION
#define global_dim DIMENSION
#define ncells (1 << DIMENSION)
#define MAXDIMENSION DIMENSION
#else
int global_dim;
int ncells;
#define MAXDIMENSION 6
#endif // DIMENSION

#ifdef TIMELIMIT
unsigned long long global_hist_timeout;
int hist_timeout_counter = 0;
#endif

#ifndef CACHESIZE
#define CACHESIZE (1LL<<32)
#endif

#ifndef CACHEMAPSIZE
// The ten-millionth prime; will use just under 700MB of RAM.
#define CACHEMAPSIZE 179424691
#endif

#ifndef CACHEDEPTH
#define CACHEDEPTH 4
#endif

int histogram[MAXDIMENSION + 1];

// store dimensions by which histogram index they belong in.
// usual trick where first element is length
int dims[MAXDIMENSION + 1][MAXDIMENSION * MAXDIMENSION];

// store the offsets within each dimension; e.g., 6 is 0, 2, 4, 6 => 0b1010101, 1 is 0, 1 => 0b11
long long dimoffsets[1 << MAXDIMENSION];

// masks to check each dimension for unused cubes; bit n of dim d is set if n is a multiple of (1 << d)
// Combined with cellUsed, this allows O(dimensions) check for adjacent unused cells.
long long dimoverlapchecks[MAXDIMENSION];

// Adjacent cubes which will block it in shared dimensions.
long long adjacentCells[1 << MAXDIMENSION];

// A list of cells that have all 0s where dim has 1s.
int baseCellsForDim[1 << MAXDIMENSION][(1 << (MAXDIMENSION-1)) + 1];

// All permutations of up to global_dim elements; used for rotations.
int npermutations;
int permutations[7*720][MAXDIMENSION];

// Actually, let's just store each permutation of each potential value...
int rotates[7*720][1 << MAXDIMENSION];

typedef struct {
  int index;
  int dim;
  int coord;
} placed_cube_t;

typedef struct {
  int len;
  placed_cube_t cubes[1 << (MAXDIMENSION - 1)];
} placed_cubes_t;

// A packed array of placed_cubes;
unsigned char *cache;
unsigned int cachetail;
unsigned int *cachemap;

void initPermutations() {
  int d, i, b, p;
  int nsofar = 1;
  for(i = 0; i < global_dim; i++)
    permutations[0][i] = i;
  for(d = 1; d < global_dim; d++) {
    for(b = 0; b < d; b++)
      for(i = 0; i < nsofar; i++) {
        int cur = i + (b + 1) * nsofar;
        memcpy(permutations[cur], permutations[i], sizeof(permutations[0]));
        permutations[cur][d] = permutations[i][b];
        permutations[cur][b] = d;
      }
    nsofar *= (d + 1);
  }
  npermutations = nsofar;

  for(p = 0; p < npermutations; p++)
    for(i = 0; i < ncells; i++) {
      rotates[p][i] = 0;
      for(b = 0; b < global_dim; b++)
        rotates[p][i] |= ((i >> b) & 1) << permutations[p][b];
    }

  /*for(d = 0; d < npermutations; d++) {
    for(i = 0; i < global_dim; i++)
      printf("%d ", permutations[d][i]);
    printf("\n");
  }*/
}


void initCache() {
  cachetail = 1; // leave zero empty so zero means not in cache.
  cache = malloc(CACHESIZE);
  if(!cache) {
    printf("Failure to allocate cache of size %lld\n", CACHESIZE);
    exit(1);
  }
  cachemap = calloc(CACHEMAPSIZE, sizeof(unsigned int));
  if(!cachemap) {
    printf("Failure to allocate cachemap of size %d * %lu\n", CACHEMAPSIZE, sizeof(unsigned int));
    exit(1);
  }
}

void printPlacedCubes(placed_cubes_t *cubes);

#ifdef FASTHASH
unsigned int hash(placed_cubes_t *cubes) {
  unsigned int h = (cubes->len * 1299721) % CACHEMAPSIZE;
  int i;
  for(i = 0; i < cubes->len; i++) {
    h = (h + (cubes->cubes[i].dim + 1000) * (cubes->cubes[i].coord + 1));
  }
  return h % CACHEMAPSIZE;
}
#else
unsigned int hash(placed_cubes_t *cubes) {
  unsigned int h = (cubes->len * 1299721) % CACHEMAPSIZE;
  int i;
  for(i = 0; i < cubes->len; i++)
    h = (h * 7927 + (cubes->cubes[i].dim * 101) + cubes->cubes[i].coord) % CACHEMAPSIZE;
  return h;
}
#endif

void getSkipDims(placed_cubes_t *cubes, int *skipDims, int *skipDimsMask) {
  int numBits[MAXDIMENSION][2];
  int i, b;
  *skipDims = 0;
  *skipDimsMask = 0;

  for(b = 0; b < global_dim; b++)
    numBits[b][0] = numBits[b][1] = 0;

  for(i = 0; i < cubes->len; i++)
    for(b = 0; b < global_dim; b++)
      if(!(cubes->cubes[i].dim & (1 << b)))
        numBits[b][(cubes->cubes[i].coord >> b) & 1]++;

  for(b = 0; b < global_dim; b++)
    if(numBits[b][0] != numBits[b][1] || (!numBits[b][0] && !numBits[b][1])) {
      *skipDims |= (numBits[b][1] > numBits[b][0]) << b;
      *skipDimsMask |= 1 << b;
    }

  /*printf("numBits:");
  for(b = 0; b < global_dim; b++)
    printf(" (%d %d)", numBits[b][0], numBits[b][1]);
  printf("  %d %d\n", *skipDimsMask, *skipDims);*/
}


#define debugcache(x)
//#define debugcache(x) x

long long int cacheConflicts = 0;
long long int cacheSemiConflicts = 0;
long long int cacheLoad = 0;
void addToCache(placed_cubes_t *placedCubes) {
  int i;
  int skipDims, skipDimsMask;

  getSkipDims(placedCubes, &skipDims, &skipDimsMask);
  debugcache(printf("skipDims is %d/%d\n", skipDims, skipDimsMask));

  placed_cubes_t flippedCubes;
  flippedCubes.len = placedCubes->len;
  for(i = 0; i < placedCubes->len; i++) {
    flippedCubes.cubes[i].index = placedCubes->cubes[i].index;
    flippedCubes.cubes[i].dim = placedCubes->cubes[i].dim;
    flippedCubes.cubes[i].coord = placedCubes->cubes[i].coord;
  }
#ifdef BIRTHDAYHASH
  int flip_dims, netFlipper = 0, lastFlip = 0;

  // If flipping along an axis doesn't change any cubes, don't (re)cache that flip.
  int ignoreFlipDims = 0;
  for(flip_dims = (1 << (global_dim - BIRTHDAYHASH)); flip_dims < (1 << global_dim); flip_dims <<= 1) {
    for(i = 0; i < placedCubes->len; i++) {
      int flipped = flippedCubes.cubes[i].coord ^ (flip_dims & ~flippedCubes.cubes[i].dim);
      if(flipped != flippedCubes.cubes[i].coord)
        break;
    }
    if(i == placedCubes->len)
      ignoreFlipDims |= flip_dims;
  }

  for(flip_dims = 0; flip_dims < (1 << global_dim); flip_dims += (1 << (global_dim - BIRTHDAYHASH))) {
    if(flip_dims & skipDimsMask)
      continue;

    if(flip_dims & ignoreFlipDims)
      continue;
    int realFlipDims = flip_dims | skipDims;
    netFlipper = realFlipDims ^ lastFlip;
    lastFlip = realFlipDims;
    for(i = 0; i < placedCubes->len; i++)
      flippedCubes.cubes[i].coord ^= (netFlipper & ~flippedCubes.cubes[i].dim);
#else
    debugcache(printf("Adding to cache unskipped, "));
    debugcache(printPlacedCubes(&flippedCubes));
    for(i = 0; i < placedCubes->len; i++)
      flippedCubes.cubes[i].coord ^= (skipDims & ~flippedCubes.cubes[i].dim);

#endif

    unsigned int h = hash(&flippedCubes);
    debugcache(printf("Adding to cache: %d, ", h));
    debugcache(printPlacedCubes(&flippedCubes));
    while(cachemap[h]) {
      //printf("Conflict on %d\n", h);
      cacheConflicts++;
      h++;
      if(h > CACHEMAPSIZE)
        h = 0;
    }
    cachemap[h] = cachetail;
    cacheLoad++;

    cache[cachetail++] = placedCubes->len;
    for(i = 0; i < placedCubes->len; i++) {
      cache[cachetail++] = flippedCubes.cubes[i].dim;
      cache[cachetail++] = flippedCubes.cubes[i].coord;
    }
#ifdef BIRTHDAYHASH
  }
#endif
}

void printPlacedCubes(placed_cubes_t *cubes) {
  int i;
  for(i = 0; i < cubes->len; i++)
    printf("{%d %d %d} ", cubes->cubes[i].index, cubes->cubes[i].dim, cubes->cubes[i].coord);
  printf("\n");
}


int checkCacheElement(placed_cubes_t *cubes, unsigned long long e) {
  int j;
  //printf("Checking cache: %d\n", i);
  cacheSemiConflicts++;
  if(cache[e] != cubes->len)
    return 0;
  for(j = 0; j < cubes->len; j++) {
    if(cache[e + 1 + 2*j] != cubes->cubes[j].dim)
      return 0;
    if(cache[e + 2 + 2*j] != cubes->cubes[j].coord)
      return 0;
  }
  //printf("Found in cache: (%d) ", i);
  //printPlacedCubes(cubes);
  cacheSemiConflicts--;
  return 1;
}

int checkCacheRotation(placed_cubes_t *cubes) {
  unsigned int h = hash(cubes);
  while(cachemap[h]) {
    if(checkCacheElement(cubes, cachemap[h]))
      return 1;
    h++;
  }
  return 0;
}

int nrotationsChecked = 0;
int ncacheHits = 0;
int checkCache(placed_cubes_t *placedCubes) {
  int flip_dims = 0, netFlipper = 0, lastFlip = 0;
  int i, perm;
  placed_cubes_t cubes;
  int skipDimsMask, skipDimsMaskBase;
  int skipDims, skipDimsBase;

  getSkipDims(placedCubes, &skipDimsBase, &skipDimsMaskBase);
  cubes.len = placedCubes->len;

  debugcache(printf("checking placement: "));
  debugcache(printPlacedCubes(&placedCubes));
  for(perm = 0; perm < npermutations; perm++) {
    for(i = 0; i < placedCubes->len; i++) {
      cubes.cubes[i].index = placedCubes->cubes[i].index;
      cubes.cubes[i].dim = rotates[perm][placedCubes->cubes[i].dim];
      cubes.cubes[i].coord = rotates[perm][placedCubes->cubes[i].coord];
    }
    skipDims = rotates[perm][skipDimsBase];
    skipDimsMask = rotates[perm][skipDimsMaskBase];

    int firstIndex = cubes.cubes[0].index;
    //printf("%d\n", firstIndex);
    for(i = 0; i < cubes.len && cubes.cubes[i].index == firstIndex; i++)
      if(cubes.cubes[i].dim == (1 << (global_dim - firstIndex)) - 1) {
        //printf("Brekaing: %d == %d\n", cubes.cubes[i].dim, (1 << (global_dim - firstIndex)) - 1);
        break;
      }
    //printf("%d/%d %d=%d\n", i, cubes.len, cubes.cubes[i].index, firstIndex);
    if(i == cubes.len || cubes.cubes[i].index != firstIndex)
      continue;

    int lastIndex = cubes.cubes[0].index;
    int lastDim = cubes.cubes[0].dim;
    int tmp;
    for(i = 1; i < cubes.len; i++) {
      if(cubes.cubes[i].index == lastIndex && lastDim > cubes.cubes[i].dim) {
        cubes.cubes[i-1].dim = cubes.cubes[i].dim;
        cubes.cubes[i].dim = lastDim;
        tmp = cubes.cubes[i-1].coord;
        cubes.cubes[i-1].coord = cubes.cubes[i].coord;
        cubes.cubes[i].coord = tmp;
        if(i > 1)
          i -= 2;
      }
      lastIndex = cubes.cubes[i].index;
      lastDim = cubes.cubes[i].dim;
    }

    nrotationsChecked++;

    debugcache(printf("checking rotation: "));
    debugcache(printPlacedCubes(&cubes));

    lastFlip = 0;
#ifdef BIRTHDAYHASH
    for(flip_dims = 0; flip_dims < (1 << (global_dim - BIRTHDAYHASH)); flip_dims++) {
#else
    for(flip_dims = 0; flip_dims < (1 << global_dim); flip_dims++) {
#endif
      if(flip_dims & skipDimsMask)
        continue;
      int realFlipDims = flip_dims | skipDims;
      netFlipper = realFlipDims ^ lastFlip;
      lastFlip = realFlipDims;
      for(i = 0; i < placedCubes->len; i++)
        cubes.cubes[i].coord ^= (netFlipper & ~cubes.cubes[i].dim);

      debugcache(printf("checking reflection %x: ", realFlipDims));
      debugcache(printPlacedCubes(&cubes));

      if(checkCacheRotation(&cubes)) {
        debugcache(printf("Cache hit for: "));
        debugcache(printPlacedCubes(placedCubes));
        debugcache(printf(" on: "));
        debugcache(printPlacedCubes(&cubes));
        ncacheHits++;
        return 1;
      }
    }
  }
  return 0;
}

unsigned long long global_start_time;
unsigned long long global_current_start_time;
clock_t global_current_clock_time;
unsigned long long currentTime() {
  struct timeval now;
  gettimeofday(&now, NULL);

  return (long long)(now.tv_sec) * 1000000 + now.tv_usec;
}

unsigned long long runningTime() {
  return currentTime() - global_start_time;
}

void arrangeDims() {
  int i, j, b, d;
  for(i = 0; i < ncells; i++) {
    int ones = 0;
    for(b = 0; b <= global_dim; b++)
      if(i & (1 << b))
        ones++;
    dims[ones][++dims[ones][0]] = i;
    for(j = 0; j < ncells; j++)
      if(!(j & ~i))
        dimoffsets[i] |= (1LL << j);
    //printf("%d => %llx\n", i, dimoffsets[i]);
  }

  // These are the "relevant" bits after shifting cellUsed right by (i << dim)
  // and anding with itself.
  // I admit the logic is terrifying, but it is correct if you think about it.
  // Bit fields are strange like this.
  for(i = 0; i < global_dim; i++)
    for(j = 0; j < ncells; j++)
      if(!(j & (1LL << i)))
        dimoverlapchecks[i] |= (1LL << j);

  adjacentCells[i] = 1LL << i;
  for(i = 0; i < ncells; i++)
    for(b = 1; b < ncells; b <<= 1)
      adjacentCells[i] |= (1LL << (i ^ b));

  // When adding a cell along a dimension, only use ones that are 0 when anded with it.
  for(d = 1; d < ncells; d++) {
    baseCellsForDim[d][0] = 0;
    for(i = 0; i < ncells; i++)
      if(!(i & d))
        baseCellsForDim[d][++baseCellsForDim[d][0]] = i;
  }

  /*for(d = 1; d < ncells; d++) {
    printf("%d:", d);
    for(i = 1; i <= baseCellsForDim[d][0]; i++)
      printf(" %d", baseCellsForDim[d][i]);
    printf("\n");
  }*/

  /*for(i = 0; i < ncells; i++)
    printf("%d: %llx\n", i, adjacentCells[i]);*/
  /*for(i = 0; i < global_dim; i++)
    printf("%llx\n", dimoverlapchecks[i]);*/
  /*for(i = 0; i <= global_dim; i++) {
    printf("index %d:", i);
    int j;
    for(j = 0; j <= dims[i][0]; j++)
      printf(" %d", dims[i][j]);
    printf("\n");
  }*/
  /*for(i = 0; i < ncells; i++) {
    printf("%d:", i);
    int j;
    for(j = 0; dimoffsets[i][j]; j++)
      printf(" %d", dimoffsets[i][j]);
    printf("\n");
  }*/
}

// This is quite slow, but it avoids bookkeeping when building, which is far more important.
void printLayout(long long *cellUsedByDim) {
  int i;
  int index, d;
  long long used = 0;
  for(i = 0; i < ncells; i++) {
    for(index = global_dim; index > 0; index--) {
      for(d = 1; d <= dims[index][0]; d++) {
        int dim = dims[index][d];
        if(cellUsedByDim[dim] & (1LL << i)) {
          printf("[%d %d] ", i, dim);
          used |= (dimoffsets[dim] << i);
        }
      }
    }
    if(!(used & (1LL << i)))
      printf("[%d %d] ", i, 0);
  }

  for(index = global_dim; index > 0; index--) {
    for(d = 1; d <= dims[index][0]; d++) {
      int dim = dims[index][d];
      for(i = 0; i < 63; i++)
        if(cellUsedByDim[dim] & (1LL << i))
          printf("{%d %d %d} ", global_dim - index, dim, i);
    }
  }

  //printf("\n");
}

void printHistogram();

void placeCube(int c, int dim, int index, long long *cellUsed, long long *cellUsedByDim, placed_cubes_t *placedCubes) {
  // This cell and all cells in the cube are used.
  *cellUsed |= dimoffsets[dim] << c;

  // And it's used in dim, for checking adjacent cubes in dim.
  cellUsedByDim[dim] |= (1LL << c);

  placedCubes->cubes[placedCubes->len].index = index;
  placedCubes->cubes[placedCubes->len].dim = dim;
  placedCubes->cubes[placedCubes->len].coord = c;
  placedCubes->len++;
}

void placeCubeNoCache(int c, int dim, int index, long long *cellUsed, long long *cellUsedByDim, placed_cubes_t *placedCubes) {
  // This cell and all cells in the cube are used.
  *cellUsed |= dimoffsets[dim] << c;

  // And it's used in dim, for checking adjacent cubes in dim.
  cellUsedByDim[dim] |= (1LL << c);

  placedCubes->len++;
}

void removeCube(int c, int dim, long long *cellUsed, long long *cellUsedByDim, placed_cubes_t *placedCubes) {
  // This cell and all cells in the cube are now unused.
  *cellUsed &= ~(dimoffsets[dim] << c);

  // And it's no longer used in dim, for checking adjacent cubes in dim.
  cellUsedByDim[dim] &= ~(1LL << c);

  placedCubes->len--;
}

long long counts[MAXDIMENSION + 1][1 << (MAXDIMENSION - 2)];

int real;
long long displayTotal;


int buildLayoutNoCache(int index, int count, int d, int c, long long *cellUsed, long long *cellUsedByDim, placed_cubes_t *placedCubes);
int buildLayoutCountNoCache(int index, int count, int d, int c, long long *cellUsed, long long *cellUsedByDim, placed_cubes_t *placedCubes);
#ifdef DISPLAYDEPTH
int buildLayoutNoDisplay(int index, int count, int d, int c, long long *cellUsed, long long *cellUsedByDim, placed_cubes_t *placedCubes);
int buildLayoutNoCacheNoDisplay(int index, int count, int d, int c, long long *cellUsed, long long *cellUsedByDim, placed_cubes_t *placedCubes);
int buildLayoutNoDisplayNoCache(int index, int count, int d, int c, long long *cellUsed, long long *cellUsedByDim, placed_cubes_t *placedCubes);
#endif

#define LAYOUTNAME
#include "buildlayout.c"
#undef LAYOUTNAME

#define LAYOUTNAME Count
#define ONLYCOUNT
#include "buildlayout.c"
#undef ONLYCOUNT
#undef LAYOUTNAME

#define LAYOUTNAME CountNoCache
#define ONLYCOUNT
#define NOCACHE
#include "buildlayout.c"
#undef NOCACHE
#undef ONLYCOUNT
#undef LAYOUTNAME

#define LAYOUTNAME NoCache
#define NOCACHE
#include "buildlayout.c"
#undef NOCACHE
#undef LAYOUTNAME

#ifdef DISPLAYDEPTH
#define NODISPLAY

#define LAYOUTNAME NoDisplay
#include "buildlayout.c"
#undef LAYOUTNAME

// Ew... this will rapidly explode.
#define NOCACHE
#define LAYOUTNAME NoDisplayNoCache
#include "buildlayout.c"
#undef LAYOUTNAME
#define LAYOUTNAME NoCacheNoDisplay
#include "buildlayout.c"
#undef LAYOUTNAME
#undef NOCACHE

#undef NODISPLAY
#endif

void printHistogram() {
  int i;
  for(i = 0; i <= global_dim; i++)
    printf("%d ", histogram[i]);
}

int global_progress = 0;
int total_histograms = 0;
void printProgress(int result) {
  int i, j;
  long long now = runningTime();
  fprintf(stderr, "\n%6lld.%03lld %d/%d ", now / 1000000, (now / 1000) % 1000, ++global_progress, total_histograms);
  for(i = 0; i <= global_dim; i++)
    fprintf(stderr, "%d ", histogram[i]);
  switch(result) {
    case 0: fprintf(stderr, "Failure"); break;
    case 1: fprintf(stderr, "Success"); break;
    case -1: fprintf(stderr, "Timeout"); break;
  }
  for(i = 0; i < global_dim; i++) {
    fprintf(stderr, " ");
    for(j = 0; j < (histogram[i] ? histogram[i] : 1); j++)
      fprintf(stderr, " %lld", counts[i][j]);
  }
  fprintf(stderr, "  %lld", counts[i][0]);
  fprintf(stderr, "  %lld/%dM; %uk/%lldM; +%lld ++%lld =%lld", cacheLoad, CACHEMAPSIZE / 1000000, cachetail, CACHESIZE / 1000000, cacheConflicts, cacheSemiConflicts, (global_current_start_time + (now - global_current_start_time) * total_histograms / global_progress) / 1000000);
  fprintf(stderr, "\n");
}

void buildHistograms(int index) {
  int i, j;
  if(index == global_dim) {
    if(!real) {
      total_histograms++;
      return;
    }
    printHistogram();
    printf(": ");
    fflush(stdout);
#ifdef DISPLAYDEPTH
    int c;
    startBuildLayoutCount();
    fprintf(stderr, "\n");
    cacheLoad = 0;
    nrotationsChecked = 0;
    ncacheHits = 0;
    cacheConflicts = 0;
    cacheSemiConflicts = 0;
    for(i = 0, c = DISPLAYDEPTH; i <= global_dim && c >= histogram[i]; i++)
      c -= histogram[i];
    displayTotal = counts[i][c];
#endif
    int result = startBuildLayout();
    if(!result)
      printf("Failure");
    if(result== -1)
      printf("Timeout");
    for(i = 0; i < global_dim; i++) {
      printf(" ");
      for(j = 0; j < (histogram[i] ? histogram[i] : 1); j++)
        printf(" %lld", counts[i][j]);
    }
    printf("  %lld", counts[i][0]);
    printf("  %lld/%dM; %uk/%lldM; +%lld ++%lld", cacheLoad, CACHEMAPSIZE / 1000000, cachetail / 1000, CACHESIZE / 1000000, cacheConflicts, cacheSemiConflicts);
    printf("\n");
    if(!isatty(STDOUT_FILENO))
      printProgress(result);
    return;
  }

  int total = histogram[index];
  histogram[index + 1] = 0;
  while(histogram[index] >= 0) {
    buildHistograms(index + 1);
    histogram[index]--;
    histogram[index + 1] += 2;
  }
  histogram[index + 1] += 2;
  histogram[index] = total;
}

int main(int argc, char **argv) {
  int i;
  int histargs = 0;

#ifdef DIMENSION
  if(argc != 1 && argc != 2 + global_dim) {
    printf("Usage: %s [histogram]\n", argv[0]);
    printf("Dimension: %d\n", global_dim);
    return 1;
  }
  if(argc == 2 + global_dim)
    histargs = 1;
#else
  if(argc < 2) {
    printf("Usage: %s dimensions [histogram]\n", argv[0]);
    return 1;
  }
  sscanf(argv[1], "%d", &global_dim);
  if(argc != 2 && argc != 3 + global_dim) {
    printf("Usage: %s dimensions [histogram]\n", argv[0]);
    return 1;
  }
  if(argc == 3 + global_dim)
    histargs = 2;
  ncells = 1 << global_dim;
#endif // DIMENSION
  global_start_time = currentTime();

  arrangeDims();
  initCache();
  initPermutations();

  if(histargs) {
    for(i = 0; i <= global_dim; i++)
      histogram[i] = atoi(argv[i + histargs]);
    real = 1;
    buildHistograms(global_dim);
    fprintf(stderr, "Rotations checked: %d\n", nrotationsChecked);
    fprintf(stderr, "cache hits: %d\n", ncacheHits);
    fprintf(stderr, "cache slots used: %lld\n", cacheLoad);
    return 0;
  }
  histogram[0] = 1;
  for(i = 1; i <= global_dim; i++)
    histogram[i] = 0;

  real = 0;
  buildHistograms(0);
  real = 1;
  buildHistograms(0);
  fprintf(stderr, "Rotations checked: %d\n", nrotationsChecked);

  return 0;
}
