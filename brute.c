#include <stdio.h>
#include <signal.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

typedef struct  {
  int dims;
  int len;
  int *pairs;
  int *matched;
} pairing;

typedef struct {
  int len;
  pairing *pairings;
} dimpairing;

typedef struct {
  int *cur;
  int len;
  long long all[10000][2];
} dist_t;

int fatal(char *msg) {
  printf("Fatal error: %s\n", msg);
  exit(1);
}


#ifdef DIMENSION
#define global_dim DIMENSION
#else
int global_dim;
#endif

void addPair(int dims, int pair, dimpairing *pairings, int d) {
  int p;
  for(p = 0; p < pairings[d].len; p++)
    if(dims == pairings[d].pairings[p].dims)
      break;
  if(p == pairings[d].len) {
    pairings[d].pairings[p].dims = dims;
    pairings[d].pairings[p].len = 0;
    pairings[d].len++;
  }
  int len = pairings[d].pairings[p].len;
  pairings[d].pairings[p].pairs[len] = pair & ~dims;
  pairings[d].pairings[p].matched[len] = -1;
  pairings[d].pairings[p].len++;
}

void removePair(int dims, int pair, dimpairing *pairings, int d) {
  int p;
  for(p = pairings[d].len - 1; p >= 0; p--)
    if(dims == pairings[d].pairings[p].dims)
      break;
  if(p == -1)
    fatal("didn't find dimension while removing pair");
  int len = pairings[d].pairings[p].len;
  if(pairings[d].pairings[p].pairs[len - 1] != (pair & ~dims))
    fatal("removing pair not at end of pairs list");
  pairings[d].pairings[p].len--;
  if(!pairings[d].pairings[p].len) {
    if(p != pairings[d].len - 1)
      fatal("removing paired dimension not at end of dimensions list");
    pairings[d].len--;
  }
}

void printPairings(dimpairing *pairings) {
  int d, i, j;
  for(d = 0; d < global_dim; d++) {
    for(i = 0; i < pairings[d].len; i++) {
      printf("  %d dim %d:", d, pairings[d].pairings[i].dims);
      for(j = 0; j < pairings[d].pairings[i].len; j++) {
        printf(" %s%d[0m", pairings[d].pairings[i].matched[j] == -1 ? "" : "[30;1m", pairings[d].pairings[i].pairs[j]);
      }
      printf("\n");
    }
  }
}

void printMatches(int *matching, dimpairing *pairings, int ncubes) {
  int matched = 0;
  int i;
  for(i = 0; i < ncubes; i++)
    if(matching[i] > i)
      matched++;
  printf("Matched: %d; ", matched);
  for(i = 0; i < ncubes; i++)
    if(matching[i] > i)
      printf("%d-%d  ", i, matching[i]);
  printf("\n");
  printPairings(pairings);
}

inline int adjacent(int a, int b) {
  int diff = a ^ b;
  return !(diff & (diff - 1));
}

dist_t distribution;

void printDistribution(dimpairing *pairings) {
  int i;
  printf("Distribution: ");
  for(i = 0; i <= global_dim; i++)
    printf("%d ", distribution.cur[i]);
  printf("\n");
}

void updateDistributions(dimpairing *pairings) {
  int i, d = 0;
  long long dist = 0, mult = 1;

  for(i = 0; i <= global_dim; i++) {
    dist += distribution.cur[i] * mult;
    mult *= (1 << global_dim);
  }

  int min = 0, max = distribution.len;
  d = (min + max) / 2;
  while(min < max) {
    if(distribution.all[d][0] == dist)
      break;
    if(distribution.all[d][0] < dist)
      min = d + 1;
    else
      max = d;
    d = (min + max) / 2;
  }

  if(distribution.all[d][0] != dist) {
    memmove(distribution.all[d+1], distribution.all[d],
            sizeof(distribution.all[0]) * distribution.len - d + 1);
    distribution.all[d][0] = dist;
    distribution.all[d][1] = 0;
    distribution.len++;
  }

  distribution.all[d][1]++;
}

void printDistributions(char *filename, int dim) {
  int i, d;
  long long t;
  int *dist = malloc(sizeof(int) * (dim + 1));
  FILE *fout = fopen(filename, "w");

  for(d = 0; d < distribution.len; d++) {
    t = distribution.all[d][0];
    for(i = 0; i <= dim; i++) {
      dist[i] = t % (1 << dim);
      t /= (1 << dim);
    }

    fprintf(fout, "%-12lld ", distribution.all[d][1]);
    for(i = 0; i <= dim; i++)
      fprintf(fout, " %d", dist[i]);
    fprintf(fout, "\n");
  }

  fclose(fout);
}

void signalHandler(int sig) {
  char filename[100];
  sprintf(filename, "results/brute.%d.wip", global_dim);
  printDistributions(filename, global_dim);
}

int mergeMatches(int curdim, int curp, dimpairing *pairings, int cur, pairing *curPairings) {
  //pairing curPairings = pairings[curdim].pairings[curp];
  for(; cur < curPairings->len && curPairings->matched[cur] != -1; cur++);
  if(cur == curPairings->len) {  // No more unmatched pairs
    if(curp < pairings[curdim].len)
      return mergeMatches(curdim, curp + 1, pairings, 0, &pairings[curdim].pairings[curp + 1]);
    else if(curdim + 1 < global_dim)
      return mergeMatches(curdim + 1, 0, pairings, 0, pairings[curdim + 1].pairings);
    else {
      //printf("===== Complete =========\n");
      //printPairings(pairings);
      //printDistribution(pairings);
      updateDistributions(pairings);
      return 1;
    }
  }

  int coords = curPairings->pairs[cur];
  int other;

  int mustMatch = 0;
  for(other = 0; other < cur; other++)
    if(curPairings->matched[other] == -1 && adjacent(curPairings->pairs[other], coords))
      mustMatch = 1;

  int nmatches = 0;
  for(other = cur + 1; other < curPairings->len; other++) {
    int otherCoords = curPairings->pairs[other];
    if(adjacent(otherCoords, coords) && curPairings->matched[other] == -1) {
      int newDim = coords ^ otherCoords;
      // only merge "upwards"; if we have dimensions 1 and 2, merge the two 1's
      // into 2's, but don't merge the two 2's into 1's.  They have the same
      // result, so it's a waste of time to do both.
      // And this direction works nicely with assuming 0-1 pairing, since that
      // will always merge if possible.
      if(((newDim - 1) & curPairings->dims) != curPairings->dims)
        continue;

      curPairings->matched[cur] = otherCoords;
      curPairings->matched[other] = coords;

      distribution.cur[curdim + 1] -= 2;
      distribution.cur[curdim + 2]++;
      int dims = curPairings->dims | newDim;
      addPair(dims, coords, pairings, curdim + 1);

      nmatches += mergeMatches(curdim, curp, pairings, cur + 1, curPairings);

      removePair(dims, coords, pairings, curdim + 1);
      distribution.cur[curdim + 1] += 2;
      distribution.cur[curdim + 2]--;

      curPairings->matched[cur] = -1;
      curPairings->matched[other] = -1;
    }
  }

  if(!mustMatch)
    nmatches += mergeMatches(curdim, curp, pairings, cur + 1, curPairings);
  return nmatches;
}

int buildMatches(int *matching, dimpairing *pairings, int cur) {
  int ncubes = 1 << global_dim;
  int b, other;
  int nmatches = 0;

  for(; cur < ncubes && matching[cur] != -1; cur++);
  if(cur == ncubes) { // No more unmatched cubes
    //printMatches(matching, pairings, ncubes);
    mergeMatches(0, 0, pairings, 0, pairings[0].pairings);
    return 1;
  }

  // We always match lower numbers to higher numbers to avoid duplication.
  // So if there's a potential lower match, we *must* match this one to a
  // higher one in order for the matching to be maximal.
  int mustMatch = 0;
  for(b = 1; b < ncubes; b <<= 1) {
    other = cur ^ b;

    // If a "wrong way" match is available, this cube must be matched.
    if(other < cur && matching[other] == -1)
      mustMatch = 1;

    // If a "right way" match is available, try it.
    if(other > cur && matching[other] == -1) {
      matching[cur] = other;
      matching[other] = cur;
      distribution.cur[0] -= 2;
      distribution.cur[1]++;
      addPair(b, cur, pairings, 0);
      nmatches += buildMatches(matching, pairings, cur + 1);
      removePair(b, cur, pairings, 0);
      distribution.cur[1]--;
      distribution.cur[0] += 2;
      matching[other] = -1;
    }
  }

  matching[cur] = -1;
  // And if we don't have to match it, try leaving it unmatched.
  if(!mustMatch) {
    nmatches += buildMatches(matching, pairings, cur + 1);
  }

  return nmatches;
}

int main(int argc, char **argv) {
#ifdef DIMENSION
  if(argc != 1) {
    printf("Usage: %s\n", argv[0]);
    printf("Dimension: %d\n", global_dim);
    return 1;
  }
#else
  if(argc != 2) {
    printf("Usage: %s [dimensions]\n", argv[0]);
    return 1;
  }
#endif // DIMENSION

  printf("pid is %d\n", getpid());

#ifndef DIMENSION
  int dim;
  sscanf(argv[1], "%d", &dim);
  global_dim = dim;
#endif
  signal(SIGUSR2, signalHandler);
  int ncubes = 1 << global_dim;
  int i, d;

  // Map of each cube to its match in the top-level matching.
  int *matching = malloc(sizeof(int) * ncubes);
  for(i = 0; i < ncubes; i++)
    matching[i] = -1;

  // pairings[d][p] is the list of all pairings with d+1 shared dimensions that
  // share pairings[d][p].dims: - dims: a // bitmask with 1s on the shared
  // dimensions)
  // - pairs: a list of all merged pairs merged along that dimension
  //          (normalized with shared dimensions set to 0
  // - len: the length of pairs.
  int maxPairings = ncubes;
  dimpairing *pairings = malloc(sizeof(dimpairing) * global_dim);
  for(d = 0; d < global_dim; d++) {
    pairings[d].len = 0;
    pairings[d].pairings = malloc(sizeof(pairing *) * maxPairings);
    for(i = 0; i < maxPairings; i++) {
      pairings[d].pairings[i].pairs = malloc(sizeof(int) * maxPairings);
      pairings[d].pairings[i].matched = malloc(sizeof(int) * maxPairings);
    }
  }
  distribution.len = 0;
  distribution.cur = malloc(sizeof(distribution.cur[0]) * global_dim + 1);
  for(d = 0; d <= global_dim; d++)
    distribution.cur[d] = 0;

  // Assume that 0 connects to 1.  This is a trivial symmetry that gives a
  // ~d-fold speedup for d dimensions.
  matching[0] = 1;
  matching[1] = 0;
  pairings[0].len = 1;
  pairings[0].pairings[0].len = 1;
  pairings[0].pairings[0].dims = 1;
  pairings[0].pairings[0].pairs[0] = 0;
  pairings[0].pairings[0].matched[0] = -1;
  distribution.cur[0] = (1 << global_dim) - 2;
  distribution.cur[1] = 1;
  int matches = buildMatches(matching, pairings, 0);
  printf("Top-level checks: %d\n", matches);
  char filename[100];
  sprintf(filename, "results/brute.%d.out", global_dim);
  printDistributions(filename, global_dim);

  return 0;
}
