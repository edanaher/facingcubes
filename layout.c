#include <stdio.h>

#define MAXDIMENSION 8

int global_dim;
int histogram[MAXDIMENSION + 1];

void printHistogram() {
  int i;
  for(i = 0; i <= global_dim; i++)
    printf("%d ", histogram[i]);
  printf("\n");
}

void buildHistograms(int index) {
  if(index == global_dim) {
    printHistogram();
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

  histogram[0] = 1;
  for(i = 1; i <= global_dim; i++)
    histogram[i] = 0;
  buildHistograms(0);

  return 0;
}
