#include <stdio.h>
#include <signal.h>
#include <stdlib.h>
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


int fatal(char *msg) {
  printf("Fatal error: %s\n", msg);
  exit(1);
}

int global_dim;

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

void printPairings(int dim, dimpairing *pairings) {
  int d, i, j;
  for(d = 0; d < dim; d++) {
    for(i = 0; i < pairings[d].len; i++) {
      printf("  %d dim %d:", d, pairings[d].pairings[i].dims);
      for(j = 0; j < pairings[d].pairings[i].len; j++) {
        printf(" %s%d[0m", pairings[d].pairings[i].matched[j] == -1 ? "" : "[30;1m", pairings[d].pairings[i].pairs[j]);
      }
      printf("\n");
    }
  }
}

void printMatches(int dim, int *matching, dimpairing *pairings, int ncubes) {
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
  printPairings(dim, pairings);
}

inline int adjacent(int a, int b) {
  int diff = a ^ b;
  return !(diff & (diff - 1));
}

int distribution[100];

void printDistribution(int dim, dimpairing *pairings) {
  int i;
  printf("Distribution: ");
  for(i = 0; i <= dim; i++)
    printf("%d ", distribution[i]);
  printf("\n");
}

int ndistributions = 0;
long long distributions[10000][2];

void updateDistributions(int dim, dimpairing *pairings) {
  int i, d;
  long long dist = 0, mult = 1;

  for(i = 0; i <= dim; i++) {
    dist += distribution[i] * mult;
    mult *= (1 << dim);
  }

  for(d = 0; d < ndistributions; d++)
    if(distributions[d][0] == dist)
      break;
  if(d == ndistributions) {
    ndistributions++;
    distributions[d][0] = dist;
    distributions[d][1] = 0;
  }

  distributions[d][1]++;
}

void printDistributions(char *filename, int dim) {
  int i, d;
  long long t;
  int *dist = malloc(sizeof(int) * (dim + 1));
  FILE *fout = fopen(filename, "w");

  for(d = 0; d < ndistributions; d++) {
    t = distributions[d][0];
    for(i = 0; i <= dim; i++) {
      dist[i] = t % (1 << dim);
      t /= (1 << dim);
    }

    fprintf(fout, "%-12lld ", distributions[d][1]);
    for(i = 0; i <= dim; i++)
      fprintf(fout, " %d", dist[i]);
    fprintf(fout, "\n");
  }

  fclose(fout);
}

void signalHandler(int sig) {
  printDistributions("brute.wip", global_dim);
}

int mergeMatches(int dim, int curdim, int curp, dimpairing *pairings, int cur) {
  pairing curPairings = pairings[curdim].pairings[curp];
  for(; cur < curPairings.len && curPairings.matched[cur] != -1; cur++);
  if(cur == curPairings.len) {  // No more unmatched pairs
    if(curp < pairings[curdim].len)
      return mergeMatches(dim, curdim, curp + 1, pairings, 0);
    else if(curdim + 1 < dim)
      return mergeMatches(dim, curdim + 1, 0, pairings, 0);
    else {
      //printf("===== Complete =========\n");
      //printPairings(dim, pairings);
      //printDistribution(dim, pairings);
      updateDistributions(dim, pairings);
      return 1;
    }
  }

  int coords = curPairings.pairs[cur];
  int other;

  int mustMatch = 0;
  for(other = 0; other < cur; other++)
    if(curPairings.matched[other] == -1 && adjacent(curPairings.pairs[other], coords))
      mustMatch = 1;

  int nmatches = 0;
  for(other = cur + 1; other < curPairings.len; other++) {
    int otherCoords = curPairings.pairs[other];
    if(adjacent(otherCoords, coords) && curPairings.matched[other] == -1) {
      curPairings.matched[cur] = otherCoords;
      curPairings.matched[other] = coords;

      distribution[curdim + 1] -= 2;
      distribution[curdim + 2]++;
      int dims = curPairings.dims | (coords ^ otherCoords);
      addPair(dims, coords, pairings, curdim + 1);

      nmatches += mergeMatches(dim, curdim, curp, pairings, cur + 1);

      removePair(dims, coords, pairings, curdim + 1);
      distribution[curdim + 1] += 2;
      distribution[curdim + 2]--;

      curPairings.matched[cur] = -1;
      curPairings.matched[other] = -1;
    }
  }

  if(!mustMatch)
    nmatches += mergeMatches(dim, curdim, curp, pairings, cur + 1);
  return nmatches;
}

void mirror(dimpairing *pairings, int mirror) {
  int p, i;
  for(p = 0; p < pairings[0].len; p++) {
    int flipper = mirror & ~pairings[0].pairings[p].dims;
    for(i = 0; i < pairings[0].pairings[p].len; i++) {
      pairings[0].pairings[p].pairs[i] ^= flipper;
    }
  }
}

int hash(dimpairing *pairings, int mirror) {
  int p, i;
  int h = 0;
  for(p = 0; p < pairings[0].len; p++) {
    int flipper = mirror & ~pairings[0].pairings[p].dims;
    for(i = 0; i < pairings[0].pairings[p].len; i++)
      h += ((pairings[0].pairings[p].pairs[i] ^ flipper) + 1) * pairings[p].len;
    h *= pairings[0].pairings[p].dims;
  }
  return h;
}

int checkCanonical(int dim, dimpairing *pairings) {
  int myHash = hash(pairings, 0);
  int i;
  for(i = 1; i < (1 << dim); i = i << 1) {
    int h = hash(pairings, i);
    if(h < myHash)
      return 0;
  }
  return 1;
}

int buildMatches(int dim, int *matching, dimpairing *pairings, int cur) {
  int ncubes = 1 << dim;
  int b, other;
  int nmatches = 0;

  for(; cur < ncubes && matching[cur] != -1; cur++);
  if(cur == ncubes) { // No more unmatched cubes
    //printMatches(dim, matching, pairings, ncubes);
    if(!checkCanonical(dim, pairings))
      return 0;
    mergeMatches(dim, 0, 0, pairings, 0);
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
      distribution[0] -= 2;
      distribution[1]++;
      addPair(b, cur, pairings, 0);
      nmatches += buildMatches(dim, matching, pairings, cur + 1);
      removePair(b, cur, pairings, 0);
      distribution[1]--;
      distribution[0] += 2;
      matching[other] = -1;
    }
  }

  matching[cur] = -1;
  // And if we don't have to match it, try leaving it unmatched.
  if(!mustMatch) {
    nmatches += buildMatches(dim, matching, pairings, cur + 1);
  }

  return nmatches;
}

int main(int argc, char **argv) {
  if(argc != 2) {
    printf("Usage: %s [dimensions]\n", argv[0]);
    return 1;
  }

  printf("pid is %d\n", getpid());

  int dim;
  sscanf(argv[1], "%d", &dim);
  global_dim = dim;
  signal(SIGUSR2, signalHandler);
  int ncubes = 1 << dim;
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
  dimpairing *pairings = malloc(sizeof(dimpairing) * dim);
  for(d = 0; d < dim; d++) {
    pairings[d].len = 0;
    pairings[d].pairings = malloc(sizeof(pairing *) * maxPairings);
    for(i = 0; i < maxPairings; i++) {
      pairings[d].pairings[i].pairs = malloc(sizeof(int) * maxPairings);
      pairings[d].pairings[i].matched = malloc(sizeof(int) * maxPairings);
    }
  }
  for(d = 0; d <= dim; d++)
    distribution[d] = 0;

  // Assume that 0 connects to 1.  This is a trivial symmetry that gives a
  // ~d-fold speedup for d dimensions.
  matching[0] = 1;
  matching[1] = 0;
  pairings[0].len = 1;
  pairings[0].pairings[0].len = 1;
  pairings[0].pairings[0].dims = 1;
  pairings[0].pairings[0].pairs[0] = 0;
  pairings[0].pairings[0].matched[0] = -1;
  distribution[0] = (1 << dim) - 2;
  distribution[1] = 1;
  int nmatches = buildMatches(dim, matching, pairings, 0);
  printf("%d top-level matchings checked\n", nmatches);
  printDistributions("brute.out", dim);

  return 0;
}
