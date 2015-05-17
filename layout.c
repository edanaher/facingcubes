#include <stdio.h>
#include <unistd.h>
#include <sys/time.h>

#define MAXDIMENSION 8

int global_dim;
int histogram[MAXDIMENSION + 1];

// store dimensions by which histogram index they belong in.
// usual trick where first element is length
int dims[MAXDIMENSION + 1][MAXDIMENSION * MAXDIMENSION];

// store the offsets within each dimension; e.g., 6 is 2, 4, 6
int dimoffsets[1 << MAXDIMENSION][1 + (1 << MAXDIMENSION)];


long long global_start_time;
long long currentTime() {
  struct timeval now;
  gettimeofday(&now, NULL);

  return (long long)(now.tv_sec) * 1000000 + now.tv_usec;
}

long long runningTime() {
  return currentTime() - global_start_time;
}

void arrangeDims() {
  int i, j, b;
  for(i = 0; i < (1 << global_dim); i++) {
    int ones = 0;
    for(b = 0; b <= global_dim; b++)
      if(i & (1 << b))
        ones++;
    dims[ones][++dims[ones][0]] = i;
    int n = 0;
    for(j = 1; j < (1 << global_dim); j++)
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
  /*for(i = 0; i < (1 << global_dim); i++) {
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
  for(i = 0; i < (1 << global_dim); i++)
    if(!(i & cellUsed[i]))
      printf("[%d %d] ", i, cellUsed[i]);
  printf("\n");
}

void printHistogram();

void placeCube(int c, int dim) {
  int i, b;
  // This cell and all cells in the cube are used.
  cellUsed[c] = dim;
  for(i = 0; dimoffsets[dim][i]; i++)
    cellUsed[c + dimoffsets[dim][i]] = dim;

  // And its adjacent cubes in this dimension as bad.
  for(b = 1; b < (1 << global_dim); b <<= 1)
    cellBlocked[dim][c ^ b]++;
}

void removeCube(int c, int dim) {
  int i, b;
  // This cell and all cells in the cube are now unused.
  cellUsed[c] = 0;
  for(i = 0; dimoffsets[dim][i]; i++)
    cellUsed[c + dimoffsets[dim][i]] = 0;

  // And its adjacent cubes in this dimension are less bad.
  for(b = 1; b < (1 << global_dim); b <<= 1)
    cellBlocked[dim][c ^ b]--;
}

// - index: which histogram index
// - count: which element within this histogram index
// - d: current dimension index
int buildLayout(int index, int count, int d, int c) {
  int dim = dims[global_dim - index][d];
  int b, i, success = 0;

  //printf("Entering %d %d %d (%d)\n", index, count, d, dim);
  //printLayout();

  if(count == histogram[index]) { // Placed all of this index
    //printf("Final check...\n");
    if(index == global_dim - 1) { // Placed all indices except zero-dimensional; check those!
      for(c = 0; c < (1 << global_dim); c++)
        if(!cellUsed[c])
          for(b = 1; b < (1 << global_dim); b <<= 1)
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

  for(; c < (1 << global_dim); c++) {
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
    placeCube(c, dim);
    success = buildLayout(index, count + 1, d, c + 1);
    removeCube(c, dim);

    //printf("Leaving\n");
    //printLayout();
    if(success)
      return 1;
  }

  // Finally, let's consider this dimension finished and try the next one:
  return buildLayout(index, count, d + 1, 0);
}

int startBuildLayout() {
  int index;
  for(index = 0; !histogram[index]; index++);
  placeCube(0, dims[global_dim - index][1]);
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
void printProgress() {
  int i;
  long long now = runningTime();
  fprintf(stderr, "%6lld.%03lld %d/%d  ", now / 1000000, (now / 1000) % 1000, ++global_progress, total_histograms);
  for(i = 0; i <= global_dim; i++)
    fprintf(stderr, "%d ", histogram[i]);
  fprintf(stderr, "\n");
}

void buildHistograms(int index, int real) {
  if(index == global_dim) {
    if(!real) {
      total_histograms++;
      return;
    }
    if(!isatty(STDOUT_FILENO))
      printProgress();
    printHistogram();
    printf(": ");
    fflush(stdout);
    if(!startBuildLayout())
      printf("Failure\n");
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
  int dim;
  int i;

  if(argc != 2) {
    printf("Usage: %s [dimensions]\n", argv[0]);
    return 1;
  }
  sscanf(argv[1], "%d", &dim);
  global_dim = dim;
  global_start_time = currentTime();

  arrangeDims();
  histogram[0] = 1;
  for(i = 1; i <= global_dim; i++)
    histogram[i] = 0;
  buildHistograms(0, 0);
  buildHistograms(0, 1);

  return 0;
}
