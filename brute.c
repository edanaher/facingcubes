#include <stdio.h>
#include <stdlib.h>

typedef struct  {
  int dims;
  int len;
  int *pairs;
} pairing;

typedef struct {
  int len;
  pairing *pairings;
} dimpairing;


int fatal(char *msg) {
  printf("Fatal error: %s\n", msg);
  exit(1);
}

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
  pairings[d].pairings[p].len++;
}

void removePair(int dims, int pair, dimpairing *pairings, int d) {
  int p;
  for(p = 0; p < pairings[d].len; p++)
    if(dims == pairings[d].pairings[p].dims)
      break;
  if(p == pairings[d].len)
    fatal("didn't find dimension while removing pair");
  int len = pairings[d].pairings[p].len;
  if(pairings[d].pairings[p].pairs[len - 1] != pair & ~dims)
    fatal("removing pair not at end of pairs list");
  pairings[d].pairings[p].len--;
  if(!pairings[d].pairings[p].len) {
    if(p != pairings[d].len - 1)
      fatal("removing paired dimension not at end of dimensions list");
    pairings[d].len--;
  }
}

void printMatches(int dim, int *matching, dimpairing *pairings, int ncubes) {
  int matched = 0;
  int i, j, d;
  for(i = 0; i < ncubes; i++)
    if(matching[i] > i)
      matched++;
  printf("Matched: %d; ", matched);
  for(i = 0; i < ncubes; i++)
    if(matching[i] > i)
      printf("%d-%d  ", i, matching[i]);
  printf("\n");
  for(d = 0; d < dim; d++) {
    for(i = 0; i < pairings[d].len; i++) {
      printf("  %d dim %d:", d, pairings[d].pairings[i].dims);
      for(j = 0; j < pairings[d].pairings[i].len; j++)
        printf(" %d", pairings[d].pairings[i].pairs[j]);
      printf("\n");
    }
  }
}

int buildMatches(int dim, int *matching, int npairdims, pairing *pairings, int cur) {
  int ncubes = 1 << dim;
  int b, other;
  int nmatches = 0;

  for(; cur < ncubes && matching[cur] != -1; cur++);
  if(cur == ncubes) { // No more unmatched cubes
    printMatches(dim, matching, pairings, ncubes);
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
      addPair(b, cur, pairings, 0);
      nmatches += buildMatches(dim, matching, npairdims, pairings, cur + 1);
      removePair(b, cur, pairings, 0);
      matching[other] = -1;
    }
  }

  matching[cur] = -1;
  // And if we don't have to match it, try leaving it unmatched.
  if(!mustMatch) {
    nmatches += buildMatches(dim, matching, npairdims, pairings, cur + 1);
  }

  return nmatches;
}

int main(int argc, char **argv) {
  if(argc != 2) {
    printf("Usage: %s [dimensions]\n", argv[0]);
    return 1;
  }

  int dim;
  sscanf(argv[1], "%d", &dim);
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
  int maxPairings = ncubes / 2;
  dimpairing *pairings = malloc(sizeof(dimpairing) * dim);
  for(d = 0; d < dim; d++) {
    pairings[d].len = 0;
    pairings[d].pairings = malloc(sizeof(pairing *) * maxPairings);
    for(i = 0; i < maxPairings; i++)
      pairings[d].pairings[i].pairs = malloc(sizeof(int) * maxPairings);
  }

  int nmatches = buildMatches(dim, matching, 0, pairings, 0);
  printf("Found %d matches\n", nmatches);
}
