#include <stdio.h>

#define MAXDIMENSION 8

int global_dim;
int histogram[MAXDIMENSION + 1];

// store dimensions by which histogram index they belong in.
// usual trick where first element is length
int dims[MAXDIMENSION + 1][MAXDIMENSION * MAXDIMENSION];

void arrangeDims() {
  int i, b;
  for(i = 0; i < (1 << global_dim); i++) {
    int ones = 0;
    for(b = 0; b <= global_dim; b++)
      if(i & (1 << b))
        ones++;
    dims[ones][++dims[ones][0]] = i;
    //printf("%d => %d\n", i, ones);
  }
  /*for(i = 0; i <= global_dim; i++) {
    printf("index %d:", i);
    int j;
    for(j = 0; j <= dims[i][0]; j++)
      printf(" %d", dims[i][j]);
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

// - index: which histogram index
// - count: which element within this histogram index
// - d: current dimension index
int buildLayout(int index, int count, int d, int c) {
  int dim = dims[global_dim - index][d];
  int b, c2, success = 0;

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

    // This is a stupid way to do this, but it's easy: For each cell, if it
    // differs from c only by dimensions in dim, it's in the cube.  If it's
    // occupied, then this position won't work.
    for(c2 = 0; c2 < (1 << global_dim); c2++) {
      if(!((c ^ c2) & (~dim)) && cellUsed[c2])
        break;
    }
    if(c2 != (1 << global_dim))
      continue;

    // Now we know this cube is safe.  Let's place it...
    // Again, this is a stupid way to do this.
    for(c2 = 0; c2 < (1 << global_dim); c2++)
      if(!((c ^ c2) & (~dim)))
        cellUsed[c2] = dim;

    // And mark its adjacent cubes in this dimension as bad.
    for(b = 1; b < (1 << global_dim); b <<= 1) {
      //printf("Blocking [%d %d]\n", c ^ b, dim);
      cellBlocked[dim][c ^ b]++;
    }

    //printf("Adding [%d %d]\n", c, dim);

    // And recurse
    success = buildLayout(index, count + 1, d, c + 1);

    //printf("Removing [%d %d]\n", c, dim);

    // And remove the cube:
    // Unmark its adjacent cubes in this dimension as bad.
    for(b = 1; b < (1 << global_dim); b <<= 1)
      cellBlocked[dim][c ^ b]--;

    // And mark this cube as unused.
    // Again, this is a stupid way to do this.
    for(c2 = 0; c2 < (1 << global_dim); c2++)
      if(!((c ^ c2) & (~dim)))
        cellUsed[c2] = 0;
    //printf("Leaving\n");
    //printLayout();
    if(success)
      return 1;
  }

  // Finally, let's consider this dimension finished and try the next one:
  return buildLayout(index, count, d + 1, 0);
}

void printHistogram() {
  int i;
  for(i = 0; i <= global_dim; i++)
    printf("%d ", histogram[i]);
}

void buildHistograms(int index) {
  if(index == global_dim) {
    printHistogram();
    printf(": ");
    if(!buildLayout(0, 0, 1, 0))
      printf("Failure\n");
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
  int dim;
  int i;

  if(argc != 2) {
    printf("Usage: %s [dimensions]\n", argv[0]);
    return 1;
  }
  sscanf(argv[1], "%d", &dim);
  global_dim = dim;

  arrangeDims();
  histogram[0] = 1;
  for(i = 1; i <= global_dim; i++)
    histogram[i] = 0;
  buildHistograms(0);

  return 0;
}
