#define concat(x,y) x ## y
#define concat3(x,y,z) x ## y ## z
#define buildLayoutName(x) concat(buildLayout, x)
#define buildLayoutNameSuffix(x,y) concat3(buildLayout, x, y)
#define startBuildLayoutName(x) concat(startBuildLayout, x)

// - index: which histogram index
// - count: which element within this histogram index
// - d: current dimension index
int buildLayoutName(LAYOUTNAME)(int index, int count, int d, int c) {
  int dim = dims[global_dim - index][d];
  int b, i, success = 0;

#if !defined(NOCACHE) && !defined(ONLYCOUNT)
  if(placedCubes.len >= CACHEDEPTH)
    return buildLayoutNameSuffix(LAYOUTNAME,NoCache)(index, count, d, c);
#endif

  counts[index][count]++;

  //printf("Entering %d %d %d (%d)\n", index, count, d, dim);
  //printLayout();

#if defined(DISPLAYDEPTH) && !defined(NODISPLAY)
#ifndef ONLYCOUNT
  if(placedCubes.len == DISPLAYDEPTH) {
    long long now = runningTime();
    fprintf(stderr, "%6lld.%03lld %lld/%lld (%lld +%lld)\033[1A\n", now / 1000000, (now / 1000) % 1000, counts[index][count], displayTotal, cacheLoad, cacheConflicts);
  }
  if(placedCubes.len > DISPLAYDEPTH) {
    counts[index][count]--;
    return buildLayoutNameSuffix(LAYOUTNAME,NoDisplay)(index, count, d, c);
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

  if(count == histogram[index]) { // Placed all of this index
    //printf("Final check...\n");
    if(index == global_dim - 1) { // Placed all indices except zero-dimensional; check those!
      counts[index+1][0]++;
      for(b = 0; b < global_dim; b++)
        if(((~cellUsed) >> (1 << b)) & (~cellUsed) & dimoverlapchecks[b])
          break;
      if(b != global_dim)
        return 0;
      printf("Success: ");
      printLayout();
      return 1;
    } else { // Not finished; start next index
      return buildLayoutName(LAYOUTNAME)(index + 1, 0, 1, 0);
    }
  } else if(d > dims[global_dim - index][0]) { // Out of dimensions for this index; give up.
    return 0;
  }

#if defined(DISPLAYDEPTH) && !defined(NODISPLAY)
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
    if(cellUsed & (dimoffsets[dim] << c))
      continue;

    // Now we know this cube is safe.  Let's go!
    placeCube(c, dim, index);
    //printf("Placed cube %d %d %d\n", index, dim, c);
#if !defined(NOCACHE)
      if(checkCache()) { // Already saw it, must have failed.
        //printf("cached\n");
        removeCube(c, dim);
        continue;
      }
      addToCache();
#endif
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
    for(j = 0; j < histogram[i]; j++)
      counts[i][j] = 0;
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

