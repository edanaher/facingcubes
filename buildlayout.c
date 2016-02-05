#define concat(x,y) x ## y
#define concat3(x,y,z) x ## y ## z
#define buildLayoutName(x) concat(buildLayout, x)
#define buildLayoutNameSuffix(x,y) concat3(buildLayout, x, y)
#define startBuildLayoutName(x) concat(startBuildLayout, x)

// - index: which histogram index
// - count: which element within this histogram index
// - d: current dimension index
int buildLayoutName(LAYOUTNAME)(int index, int count, int d, int c, long long *cellUsed, long long *cellUsedByDim, placed_cubes_t *placedCubes) {
  int dim = dims[global_dim - index][d];
  int b, success = 0;

#if !defined(NOCACHE)
  if(placedCubes->len >= CACHEDEPTH)
    return buildLayoutNameSuffix(LAYOUTNAME,NoCache)(index, count, d, c, cellUsed, cellUsedByDim, placedCubes);
#endif

  counts[index][count]++;

  //printf("Entering %d %d %d (%d)\n", index, count, d, dim);
  //printLayout(cellUsedByDim);

#if defined(DISPLAYDEPTH) && !defined(NODISPLAY)
#ifndef ONLYCOUNT
  if(placedCubes->len == DISPLAYDEPTH && count < histogram[index]) {
    long long now = runningTime();
    clock_t clocknow = clock();
    fprintf(stderr, "%6ld.%03ld/%6lld.%03lld %lld/%lld (%lld +%lld ++%lld; %uM) =%lld/%lld  \033[1A\n", clocknow / CLOCKS_PER_SEC, clocknow / (CLOCKS_PER_SEC / 1000) % 1000, now / 1000000, (now / 1000) % 1000, counts[index][count], displayTotal, cacheLoad, cacheConflicts, cacheSemiConflicts, cachetail/(1<<20), (global_current_clock_time + (clocknow - global_current_clock_time) * displayTotal / counts[index][count]) / CLOCKS_PER_SEC, (global_current_start_time + (now - global_current_start_time) * displayTotal / counts[index][count]) / 1000000);
#ifdef KEEPPROGRESS
    if(counts[index][count] % (displayTotal / KEEPPROGRESS + 1) == 0)
      fprintf(stderr, "\n");
#endif
  }

  if(placedCubes->len > DISPLAYDEPTH) {
    counts[index][count]--;
    return buildLayoutNameSuffix(LAYOUTNAME,NoDisplay)(index, count, d, c, cellUsed, cellUsedByDim, placedCubes);
  }
#else
  if(placedCubes->len == DISPLAYDEPTH) {
    if(!(counts[index][count] % 10000)) {
      long long now = runningTime();
      clock_t clocknow = clock();
      fprintf(stderr, "%6ld.%03ld/%6lld.%03lld %lld (%lld +%lld ++%lld)\033[1A\n", clocknow / CLOCKS_PER_SEC, clocknow / (CLOCKS_PER_SEC / 1000) % 1000, now / 1000000, (now / 1000) % 1000, counts[index][count], cacheLoad, cacheConflicts, cacheSemiConflicts);
    }
  }
#endif
#endif

#ifdef TIMELIMIT
  // Tradeoff accuracy for less time wasted checking
  if(++hist_timeout_counter > 100000 * TIMELIMIT) {
    clock_t clocknow = clock();
    if(clocknow * 1000000 / CLOCKS_PER_SEC > global_hist_timeout) {
#ifdef TIMELIMITEXPECTED
      int i;
      for(i = 0, c = DISPLAYDEPTH; i <= global_dim && c >= histogram[i]; i++)
        c -= histogram[i];
      if((global_current_clock_time + (clocknow - global_current_clock_time) * displayTotal / counts[i][c]) / CLOCKS_PER_SEC > TIMELIMITEXPECTED)
#endif /* TIMELIMITEXPECTED */
        return -1;
    }
    hist_timeout_counter = 0;
  }
#endif /* TIMELIMIT */

  if(count == histogram[index]) { // Placed all of this index
    //printf("Final check...\n");
    if(index == global_dim - 1) { // Placed all indices except zero-dimensional; check those!
      counts[index+1][0]++;
      for(b = 0; b < global_dim; b++)
        if(((~*cellUsed) >> (1 << b)) & (~*cellUsed) & dimoverlapchecks[b])
          return 0;
      printf("Success: ");
      printLayout(cellUsedByDim);
      return 1;
    } else { // Not finished; start next index
      int matching = maxMatching(cellUsed);
      if(matching < histogram[global_dim-1] || matching / 2 > histogram[global_dim-1])
        return 0;

      return buildLayoutName(LAYOUTNAME)(index + 1, 0, 1, 1, cellUsed, cellUsedByDim, placedCubes);
    }
  } else if(d > dims[global_dim - index][0]) { // Out of dimensions for this index; give up.
    return 0;
  }

#if defined(DISPLAYDEPTH) && !defined(NODISPLAY)
#ifdef ONLYCOUNT
  if(placedCubes->len != DISPLAYDEPTH)
#endif
#endif
  for(; c <= baseCellsForDim[dim][0]; c++) {
    int cube = baseCellsForDim[dim][c];
    //printf("Checking [%d %d]\n", cube, dim);

    // For each cell, if it differs from cube only by dimensions in dim, it's in
    // the cube.  If it's occupied, then this position won't work.
    if(*cellUsed & (dimoffsets[dim] << cube))
      continue;

    // Shortcut if the cell would be face-aligned in this dimension (including with itself).
    if(cellUsedByDim[dim] & adjacentCells[cube])
      continue;

    // Now we know this cube is safe.  Let's go!
#ifdef NOCACHE
    placeCubeNoCache(cube, dim, index, cellUsed, cellUsedByDim, placedCubes);
#else
    placeCube(cube, dim, index, cellUsed, cellUsedByDim, placedCubes);
#endif
    //printf("Placed cube %d %d %d\n", index, dim, cube);
#if !defined(NOCACHE)
      if(checkCache(placedCubes)) { // Already saw it, must have failed.
        //printf("cached\n");
        removeCube(cube, dim, cellUsed, cellUsedByDim, placedCubes);
        continue;
      }
      addToCache(placedCubes);
#endif
    success = buildLayoutName(LAYOUTNAME)(index, count + 1, d, c + 1, cellUsed, cellUsedByDim, placedCubes);
    removeCube(cube, dim, cellUsed, cellUsedByDim, placedCubes);

    //printf("Leaving\n");
    //printLayout(cellUsedByDim);
    if(success)
      return success;
  }

  // Finally, let's consider this dimension finished and try the next one:
  return buildLayoutName(LAYOUTNAME)(index, count, d + 1, 1, cellUsed, cellUsedByDim, placedCubes);
}

int startBuildLayoutName(LAYOUTNAME)() {
  int i, j;
  int index;

  // If the cell is in a cube
  long long cellUsed = 0;

  // cellUsedByDim[dim]: Track cells used in each dimension to avoid faceing cubes.
  long long cellUsedByDim[1 << MAXDIMENSION];
  for(i = 0; i < ncells; i++)
    cellUsedByDim[i] = 0;

  placed_cubes_t placedCubes;
  placedCubes.len = 0;

  global_current_start_time = runningTime();
  global_current_clock_time = clock();
#ifdef TIMELIMIT
  global_hist_timeout = runningTime() + TIMELIMIT * 1000000;
#endif
  for(i = 0; i < global_dim; i++)
    for(j = 0; j < histogram[i]; j++)
      counts[i][j] = 0;
  counts[global_dim][0] = 0;
  matchingPruned = 0;
  cachetail = 1;
  // This is sparse, so it's faster to free and re-alloc than to zero out.
  free(cachemap);
  cachemap = calloc(CACHEMAPSIZE, sizeof(unsigned long long));
  for(index = 0; !histogram[index]; index++);
  if(index == global_dim) // Don't add a simple cube; it makes things sad.
    return buildLayoutName(LAYOUTNAME)(index - 1, 0, 1, 1, &cellUsed, cellUsedByDim, &placedCubes);
  placeCube(0, dims[global_dim - index][1], index, &cellUsed, cellUsedByDim, &placedCubes);
  int success = buildLayoutName(LAYOUTNAME)(index, 1, 1, 2, &cellUsed, cellUsedByDim, &placedCubes);
  removeCube(0, dims[global_dim - index][1], &cellUsed, cellUsedByDim, &placedCubes);
  return success;
}

