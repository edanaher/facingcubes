#include <stdio.h>
#include <stdlib.h>

int buildMatches(int dim, int *matching, int matched, int cur) {
  int ncubes = 1 << dim;
  int b, other;
  int nmatches = 0;

  for(; cur < ncubes && matching[cur] != -1; cur++);
  if(cur == ncubes) { // No more unmatched cubes
    /*printf("Matched: %d; ", matched);
    int i;
    for(i = 0; i < ncubes; i++)
      if(matching[i] > i)
        printf("%d-%d  ", i, matching[i]);
    printf("\n");*/
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
      nmatches += buildMatches(dim, matching, matched + 1, cur + 1);
      matching[other] = -1;
    }
  }

  matching[cur] = -1;
  // And if we don't have to match it, try leaving it unmatched.
  if(!mustMatch) {
    nmatches += buildMatches(dim, matching, matched, cur + 1);
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

  int *matching = malloc(sizeof(int) * ncubes);
  for(i = 0; i < ncubes; i++)
    matching[i] = -1;

  int nmatches = buildMatches(dim, matching, 0, 0);
  printf("Found %d matches\n", nmatches);
}
