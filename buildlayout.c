#define concat(x,y) x ## y
#define buildLayoutName(x) concat(buildLayout, x)
#define startBuildLayoutName(x) concat(startBuildLayout, x)

// - index: which histogram index
// - count: which element within this histogram index
// - d: current dimension index
int buildLayoutName(LAYOUTNAME)(int index, int count, int d, int c) {
  int dim = dims[global_dim - index][d];
  int b, i, j, success = 0;

  counts[index][count]++;

  //printf("Entering %d %d %d (%d)\n", index, count, d, dim);
  //printLayout();

#ifdef DISPLAYDEPTH
#ifndef ONLYCOUNT
  if(placedCubes.len == DISPLAYDEPTH) {
    long long now = runningTime();
    fprintf(stderr, "%6lld.%03lld %lld/%lld (%lld +%lld)\033[1A\n", now / 1000000, (now / 1000) % 1000, counts[index][count], displayTotal, cacheLoad, cacheConflicts);
  }
#else
  if(placedCubes.len == DISPLAYDEPTH) {
    if(!(counts[index][count] % 10000)) {
      long long now = runningTime();
      fprintf(stderr, "%6lld.%03lld %lld (%lld +%lld)\033[1A\n", now / 1000000, (now / 1000) % 1000, counts[index][count], cacheLoad, cacheConflicts);
    }
  }
#endif
#endif

#ifdef TIMELIMIT
  // Tradeoff accuracy for less time wasted checking
  if(++hist_timeout_counter > 1000000 * TIMELIMIT) {
    //fprintf(stderr, "bail: %lld/%lld\n", runningTime(), global_hist_timeout);
    if(runningTime() > global_hist_timeout)
      return -1;
    hist_timeout_counter = 0;
  }
#endif

  if(checkTilingAt[index][count]) {
    for(c = 0; c < ncells; c++)
      if(!cellIsUsed(c)) {
        for(b = 1; b < ncells; b <<= 1)
          if(!cellIsUsed(c ^ b)) { // There's an adjacent (face-aligned) 0-cell; this fails.
            break;
          }
        if(b != ncells)
          break;
      }
    if(c == ncells) {
      printf("\nindex/count: %d/%d\n", index, count);
      printf("Success: ");
      printLayout(index, count);
      printf("\n");
      checkTilingAt[index][count] = 0;
      int keepGoing = 0;
      for(i = 0; i < global_dim; i++)
        for(j = 0; j < histogram[i]; j++) {
          if(checkTilingAt[i][j])
            keepGoing = 1;
        }
      if(!keepGoing && !checkTilingAt[global_dim - 1][histogram[global_dim - 1]])
        return 1;
    }
  }

  if(count == histogram[index]) { // Placed all of this index
    if(index == global_dim - 1) { // Placed all indices except zero-dimensional; check those!
    } else { // Not finished; start next index
      return buildLayoutName(LAYOUTNAME)(index + 1, 0, 1, 0);
    }
    return 0;
  } else if(d > dims[global_dim - index][0]) { // Out of dimensions for this index; give up.
    return 0;
  }

#ifdef DISPLAYDEPTH
#ifdef ONLYCOUNT
  if(placedCubes.len != DISPLAYDEPTH)
#endif
#endif
  for(; c < ncells; c++) {
    //printf("Checking [%d %d]\n", c, dim);
    // Only consider cubes which are 0 in the dimension we're checking.
    if(c & dim)
      continue;
    //printf("Checked [%d %d]: %d\n", c, dim, cellBlocked[dim][c]);

    // Shortcut if the cell itself is used or it would be face-aligned in this dimension.
    if(cellIsUsed(c) || cellBlocked[dim][c])
      continue;

    // For each cell, if it differs from c only by dimensions in dim, it's in
    // the cube.  If it's occupied, then this position won't work.
    for(i = 0; dimoffsets[dim][i]; i++)
      if(cellIsUsed(c + dimoffsets[dim][i]))
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
    success = buildLayoutName(LAYOUTNAME)(index, count + 1, d, c + 1);
    removeCube(c, dim);

    //printf("Leaving\n");
    //printLayout();
    if(success)
      return success;
  }

  // Finally, let's consider this dimension finished and try the next one:
  return buildLayoutName(LAYOUTNAME)(index, count, d + 1, 0);
}

int startBuildLayoutName(LAYOUTNAME)() {
  int i, j;
  int index;
#ifdef TIMELIMIT
  global_hist_timeout = runningTime() + TIMELIMIT * 1000000;
#endif
  for(i = 0; i < global_dim; i++)
    for(j = 0; j < histogram[i]; j++) {
      counts[i][j] = 0;
      checkTilingAt[i][j] = 1;
    }
  checkTilingAt[global_dim - 1][histogram[global_dim - 1]] = 1;
  counts[global_dim][0] = 0;
  cachetail = 1;
  // This is sparse, so it's faster to free and re-alloc than to zero out.
  free(cachemap);
  cachemap = calloc(CACHEMAPSIZE, sizeof(unsigned long long));
  for(index = 0; !histogram[index]; index++);
  if(index == global_dim) // Don't add a simple cube; it makes things sad.
    return buildLayoutName(LAYOUTNAME)(index - 1, 0, 1, 0);
  placeCube(0, dims[global_dim - index][1], index);
  int success = buildLayoutName(LAYOUTNAME)(index, 1, 1, 1);
  removeCube(0, dims[global_dim - index][1]);
  return success;
}

