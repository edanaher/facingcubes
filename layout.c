#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>

#ifdef DIMENSION
#define global_dim DIMENSION
#define ncells (1 << DIMENSION)
#define MAXDIMENSION DIMENSION
#else
int global_dim;
int ncells;
#define MAXDIMENSION 8
#endif // DIMENSION

#ifdef TIMELIMIT
unsigned long long global_hist_timeout;
int hist_timeout_counter = 0;
#endif

#ifndef CACHESIZE
#define CACHESIZE (1LL<<30)
#endif

#ifndef CACHEMAPSIZE
// The five-millionth prime; will use just under 700MB of RAM.
#define CACHEMAPSIZE 86028121
#endif

#ifndef CACHEDEPTH
#define CACHEDEPTH 4
#endif

int histogram[MAXDIMENSION + 1];

// store dimensions by which histogram index they belong in.
// usual trick where first element is length
int dims[MAXDIMENSION + 1][MAXDIMENSION * MAXDIMENSION];

// store the offsets within each dimension; e.g., 6 is 2, 4, 6
int dimoffsets[1 << MAXDIMENSION][1 + (1 << MAXDIMENSION)];

// All permutations of up to global_dim elements; used for rotations.
int npermutations;
int permutations[7*720][MAXDIMENSION];

typedef struct {
  int index;
  int dim;
  int coord;
} placed_cube_t;

typedef struct {
  int len;
  placed_cube_t cubes[1 << (MAXDIMENSION - 1)];
} placed_cubes_t;

placed_cubes_t placedCubes;

// A packed array of placed_cubes;
unsigned char *cache;
unsigned long long cachetail;
unsigned long long *cachemap;

void initPermutations() {
  int d, i, b;
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
  cachemap = calloc(CACHEMAPSIZE, sizeof(unsigned long long));
  if(!cachemap) {
    printf("Failure to allocate cachemap of size %d * %lu\n", CACHEMAPSIZE, sizeof(unsigned long long));
    exit(1);
  }
}

void printPlacedCubes(placed_cubes_t *cubes);

#ifdef FASTHASH
unsigned int hash(placed_cubes_t *cubes) {
  int h = (cubes->len * 1299721) % CACHEMAPSIZE;
  int i;
  for(i = 0; i < cubes->len; i++) {
    h = (h + (cubes->cubes[i].dim + 1000) * (cubes->cubes[i].coord + 1));
  }
  return h % CACHEMAPSIZE;
}
#else
unsigned int hash(placed_cubes_t *cubes) {
  unsigned int h = (cubes->len * 1299721) % CACHEMAPSIZE;
  int i, j, cur;
  int last = -1;
  for(i = 0; i < cubes->len; i++) {
    for(cur = 0; cur < cubes->len; cur++)
      if(cubes->cubes[cur].coord > last)
        break;

    for(j = cur; j < cubes->len; j++) {
      if(cubes->cubes[j].coord > last && cubes->cubes[j].coord < cubes->cubes[cur].coord)
        cur = j;
    }
    h = (h * 7927 + (cubes->cubes[i].dim * 101) + cubes->cubes[i].coord) % CACHEMAPSIZE;
  }
  return h;
}
#endif

long long int cacheConflicts = 0;
long long int cacheLoad = 0;
void addToCache() {
  int i;

  unsigned int h = hash(&placedCubes);
  //printf("Added to cache %d (%lld): ", h, cacheLoad);
  //printPlacedCubes(&placedCubes);
  while(cachemap[h]) {
    //printf("Conflict on %d\n", h);
    cacheConflicts++;
    h++;
    if(h > CACHEMAPSIZE)
      h = 0;
  }
  cachemap[h] = cachetail;

  cache[cachetail++] = placedCubes.len;
  for(i = 0; i < placedCubes.len; i++) {
    cache[cachetail++] = placedCubes.cubes[i].index;
    cache[cachetail++] = placedCubes.cubes[i].dim;
    cache[cachetail++] = placedCubes.cubes[i].coord;
  }
  cacheLoad++;
}

void printPlacedCubes(placed_cubes_t *cubes) {
  int i;
  for(i = 0; i < cubes->len; i++)
    printf("{%d %d %d} ", cubes->cubes[i].index, cubes->cubes[i].dim, cubes->cubes[i].coord);
  printf("\n");
}


int checkCacheElement(placed_cubes_t *cubes, int e) {
  int j;
  //printf("Checking cache: %d\n", i);
  if(cache[e] != cubes->len)
    return 0;
  for(j = 0; j < cubes->len; j++) {
    if(cache[e + 1 + 3*j] != cubes->cubes[j].index)
      return 0;
    if(cache[e + 2 + 3*j] != cubes->cubes[j].dim)
      return 0;
    if(cache[e + 3 + 3*j] != cubes->cubes[j].coord)
      return 0;
  }
  //printf("Found in cache: (%d) ", i);
  //printPlacedCubes(cubes);
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

int rotate(int n, int perm) {
  int r = 0;
  int b;
  for(b = 0; b < global_dim; b++)
    r |= ((n >> b) & 1) << permutations[perm][b];
  return r;
}

int nrotationsChecked = 0;
int checkCache() {
  int flip_dims = 0, netFlipper = 0;
  int i, perm;
  placed_cubes_t cubes;

  cubes.len = placedCubes.len;
  for(perm = 0; perm < npermutations; perm++) {
    for(i = 0; i < placedCubes.len; i++) {
      cubes.cubes[i].index = placedCubes.cubes[i].index;
      cubes.cubes[i].dim = rotate(placedCubes.cubes[i].dim, perm);
      cubes.cubes[i].coord = rotate(placedCubes.cubes[i].coord, perm);
    }

    //printf("Rotation: ");
    //printPlacedCubes(&cubes);
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

    nrotationsChecked++;

    for(flip_dims = 0; flip_dims < (1 << global_dim); flip_dims++) {
      netFlipper ^= flip_dims;
      for(i = 0; i < placedCubes.len; i++)
        cubes.cubes[i].coord ^= netFlipper & ~cubes.cubes[i].dim;
      //printf("Reflection: ");
      //printPlacedCubes(&cubes);

      if(checkCacheRotation(&cubes))
        return 1;
    }
  }
  return 0;
}

unsigned long long global_start_time;
unsigned long long currentTime() {
  struct timeval now;
  gettimeofday(&now, NULL);

  return (long long)(now.tv_sec) * 1000000 + now.tv_usec;
}

unsigned long long runningTime() {
  return currentTime() - global_start_time;
}

void arrangeDims() {
  int i, j, b;
  for(i = 0; i < ncells; i++) {
    int ones = 0;
    for(b = 0; b <= global_dim; b++)
      if(i & (1 << b))
        ones++;
    dims[ones][++dims[ones][0]] = i;
    int n = 0;
    for(j = 1; j < ncells; j++)
      if(!(j & ~i))
        dimoffsets[i][n++] = j;
    //printf("%d => %d\n", i, ones);
  }
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

// If the cell is in a cube
int cellUsed[1 << MAXDIMENSION];

// cellBlocked[dim][c]: 1 if there is a cube adjacent to c in dim dim (i.e., it would be face-aligned).
int cellBlocked[1 << MAXDIMENSION][1 << MAXDIMENSION];

// This is a bit slow, but it avoids bookkeeping when building, which is far more important.
void printLayout() {
  int i;
  for(i = 0; i < ncells; i++)
    if(!(i & cellUsed[i]))
      printf("[%d %d] ", i, cellUsed[i]);
  for(i = 0; i < placedCubes.len; i++)
    printf("{%d %d %d} ", placedCubes.cubes[i].index, placedCubes.cubes[i].dim, placedCubes.cubes[i].coord);
  //printf("\n");
}

void printHistogram();

void placeCube(int c, int dim, int index) {
  int i, b;
  // This cell and all cells in the cube are used.
  cellUsed[c] = dim;
  for(i = 0; dimoffsets[dim][i]; i++)
    cellUsed[c + dimoffsets[dim][i]] = dim;

  // And its adjacent cubes in this dimension as bad.
  for(b = 1; b < ncells; b <<= 1)
    cellBlocked[dim][c ^ b]++;

  placedCubes.cubes[placedCubes.len].index = index;
  placedCubes.cubes[placedCubes.len].dim = dim;
  placedCubes.cubes[placedCubes.len].coord = c;
  placedCubes.len++;
}

void removeCube(int c, int dim) {
  int i, b;
  // This cell and all cells in the cube are now unused.
  cellUsed[c] = 0;
  for(i = 0; dimoffsets[dim][i]; i++)
    cellUsed[c + dimoffsets[dim][i]] = 0;

  // And its adjacent cubes in this dimension are less bad.
  for(b = 1; b < ncells; b <<= 1)
    cellBlocked[dim][c ^ b]--;
  placedCubes.len--;
}

long long counts[MAXDIMENSION + 1][1 << (MAXDIMENSION - 2)];

// - index: which histogram index
// - count: which element within this histogram index
// - d: current dimension index
int buildLayout(int index, int count, int d, int c, int depth) {
  int dim = dims[global_dim - index][d];
  int b, i, success = 0;

  counts[index][count]++;

  //printf("Entering %d %d %d (%d)\n", index, count, d, dim);
  //printLayout();

#ifdef TIMELIMIT
  // Tradeoff accuracy for less time wasted checking
  if(++hist_timeout_counter > 1000000 * TIMELIMIT) {
    //fprintf(stderr, "bail: %lld/%lld\n", runningTime(), global_hist_timeout);
    if(runningTime() > global_hist_timeout)
      return -1;
    hist_timeout_counter = 0;
  }
#endif

  if(count == histogram[index]) { // Placed all of this index
    //printf("Final check...\n");
    if(index == global_dim - 1) { // Placed all indices except zero-dimensional; check those!
      counts[index+1][0]++;
      for(c = 0; c < ncells; c++)
        if(!cellUsed[c])
          for(b = 1; b < ncells; b <<= 1)
            if(!cellUsed[c ^ b]) { // There's an adjacent (face-aligned) 0-cell; this fails.
              //printf("Final check failed: %d %d\n", c, c^b);
              return 0;
            }
      printf("Success: ");
      printLayout();
      return 1;
    } else { // Not finished; start next index
      return buildLayout(index + 1, 0, 1, 0);
    }
  } else if(d > dims[global_dim - index][0]) { // Out of dimensions for this index; give up.
    return 0;
  }

  for(; c < ncells; c++) {
    //printf("Checking [%d %d]\n", c, dim);
    // Only consider cubes which are 0 in the dimension we're checking.
    if(c & dim)
      continue;
    //printf("Checked [%d %d]: %d\n", c, dim, cellBlocked[dim][c]);

    // Shortcut if the cell itself is used or it would be face-aligned in this dimension.
    if(cellUsed[c] || cellBlocked[dim][c])
      continue;

    // For each cell, if it differs from c only by dimensions in dim, it's in
    // the cube.  If it's occupied, then this position won't work.
    for(i = 0; dimoffsets[dim][i]; i++)
      if(cellUsed[c + dimoffsets[dim][i]])
        break;
    if(dimoffsets[dim][i])
      continue;

    // Now we know this cube is safe.  Let's go!
    placeCube(c, dim, index);
    //printf("Placed cube %d %d %d\n", index, dim, c);
    if(placedCubes.len <= CACHEDEPTH) {
      if(checkCache()) { // Already saw it, must have failed.
        //printf("cached\n");
        removeCube(c, dim);
        continue;
      }
      addToCache();
    }
    success = buildLayout(index, count + 1, d, c + 1);
    removeCube(c, dim);

    //printf("Leaving\n");
    //printLayout();
    if(success)
      return success;
  }

  // Finally, let's consider this dimension finished and try the next one:
  return buildLayout(index, count, d + 1, 0);
}

int startBuildLayout() {
  int i, j;
  int index;
#ifdef TIMELIMIT
  global_hist_timeout = runningTime() + TIMELIMIT * 1000000;
#endif
  for(i = 0; i < global_dim; i++)
    for(j = 0; j < histogram[i]; j++)
      counts[i][j] = 0;
  counts[global_dim][0] = 0;
  cachetail = 1;
  // This is sparse, so it's faster to free and re-alloc than to zero out.
  free(cachemap);
  cachemap = calloc(CACHEMAPSIZE, sizeof(unsigned long long));
  for(index = 0; !histogram[index]; index++);
  if(index == global_dim) // Don't add a simple cube; it makes things sad.
    return buildLayout(index - 1, 0, 1, 0);
  placeCube(0, dims[global_dim - index][1], index);
  int success = buildLayout(index, 1, 1, 1);
  removeCube(0, dims[global_dim - index][1]);
  return success;
}

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
  fprintf(stderr, "%6lld.%03lld %d/%d  ", now / 1000000, (now / 1000) % 1000, ++global_progress, total_histograms);
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
  fprintf(stderr, "  %lld/%dM; %lld/%lldM; +%lld", cacheLoad, CACHEMAPSIZE / 1000000, cachetail, CACHESIZE / 1000000, cacheConflicts);
  fprintf(stderr, "\n");
}

void buildHistograms(int index, int real) {
  int i, j;
  if(index == global_dim) {
    if(!real) {
      total_histograms++;
      return;
    }
    printHistogram();
    printf(": ");
    fflush(stdout);
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
    printf("  %lld/%dM; %lld/%lldM; +%lld", cacheLoad, CACHEMAPSIZE / 1000000, cachetail, CACHESIZE / 1000000, cacheConflicts);
    printf("\n");
    if(!isatty(STDOUT_FILENO))
      printProgress(result);
    return;
  }

  int total = histogram[index];
  histogram[index + 1] = 0;
  while(histogram[index] >= 0) {
    buildHistograms(index + 1, real);
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
    buildHistograms(global_dim, 1);
    return 0;
  }
  histogram[0] = 1;
  for(i = 1; i <= global_dim; i++)
    histogram[i] = 0;

  buildHistograms(0, 0);
  buildHistograms(0, 1);
  fprintf(stderr, "Rotations checked: %d\n", nrotationsChecked);

  return 0;
}
