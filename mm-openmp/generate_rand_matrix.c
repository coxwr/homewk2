#include <stdio.h>
#include <stdlib.h> // srand and rand
#include <time.h>
#include <stdlib.h>

int main(int argc, char **argv) {
  time_t t;
  /* Intializes random number generator */
  srand((unsigned)time(&t));

  int i, r, j, n = strtoull(argv[1], NULL, 10);
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      r = rand() % 200;
      printf("%d ", r);
    }
    printf("\n");
  }

  return 0;
}
