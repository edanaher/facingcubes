#include <stdio.h>
#include <stdlib.h>

typedef struct  {
  int dims;
  int len;
  int *pairs;
} pairing;

int fatal(char *msg) {
  printf("Fatal error: %s\n", msg);
  exit(1);
}

int addPair(int dims, int pair, int npairdims, pairing *pairings) {
  int p;
  for(p = 0; p < npairdims; p++)
    if(dims == pairings[p].dims)
      break;
  if(p == npairdims) {
    pairings[p].dims = dims;
    pairings[p].len = 0;
    npairdims++;
  }
  int len = pairings[p].len;
  pairings[p].pairs[len] = pair & ~dims;
  pairings[p].len++;
  return npairdims;
}

int removePair(int dims, int pair, int npairdims, pairing *pairings) {
  int p;
  for(p = 0; p < npairdims; p++)
    if(dims == pairings[p].dims)
      break;
  if(p == npairdims)
    fatal("didn't find dimension while removing pair");
  int len = pairings[p].len;
  if(pairings[p].pairs[len - 1] != pair & ~dims)
    fatal("removing pair not at end of pairs list");
  pairings[p].len--;
  if(!pairings[p].len) {
    if(p != npairdims - 1)
      fatal("removing paired dimension not at end of dimensions list");
    npairdims--;
  }
  return npairdims;
}

void printMatches(int *matching, int npairdims, pairing *pairings, int ncubes) {
  int matched = 0;
  int i, j;
  for(i = 0; i < ncubes; i++)
    if(matching[i] > i)
      matched++;
  printf("Matched: %d; ", matched);
  for(i = 0; i < ncubes; i++)
    if(matching[i] > i)
      printf("%d-%d  ", i, matching[i]);
  printf("\n");
  for(i = 0; i < npairdims; i++) {
    printf("  dim %d:", pairings[i].dims);
    for(j = 0; j < pairings[i].len; j++)
      printf(" %d", pairings[i].pairs[j]);
    printf("\n");
  }
}

int buildMatches(int dim, int *matching, int npairdims, pairing *pairings, int cur) {
  int ncubes = 1 << dim;
  int b, other;
  int nmatches = 0;

  for(; cur < ncubes && matching[cur] != -1; cur++);
  if(cur == ncubes) { // No more unmatched cubes
    printMatches(matching, npairdims, pairings, ncubes);
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
      npairdims = addPair(b, cur, npairdims, pairings);
      nmatches += buildMatches(dim, matching, npairdims, pairings, cur + 1);
      npairdims = removePair(b, cur, npairdims, pairings);
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
  int i;

  // Map of each cube to its match in the top-level matching.
  int *matching = malloc(sizeof(int) * ncubes);
  for(i = 0; i < ncubes; i++)
    matching[i] = -1;

  // pairings[d] is the list of all pairings that share pairings[d].dims:
  // - dims: a // bitmask with 1s on the shared dimensions)
  // - pairs: a list of all merged pairs merged along that dimension
  //          (normalized with shared dimensions set to 0
  // - len: the length of pairs.
  int maxPairings = ncubes / 2;
  pairing *pairings = malloc(sizeof(int *) * maxPairings);
  for(i = 0; i < maxPairings; i++)
    pairings[i].pairs = malloc(sizeof(int) * maxPairings);

  int nmatches = buildMatches(dim, matching, 0, pairings, 0);
  printf("Found %d matches\n", nmatches);
}
