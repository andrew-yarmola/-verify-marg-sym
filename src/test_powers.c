#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "arithmetic/roundoff.h"

char *double_to_hex(double x)
{
  static char buf[100];
  union {
    double d;
    long l[2];
  } u;
  u.d = x;
  snprintf(buf, sizeof(buf)/sizeof(char), "'%0lx', '%0lx' (%.18f)", u.l[0], u.l[1], u.d);
  return buf;
}

int main() {
  initialize_roundoff();
  printf("Hex of powers used in box dimensions\n");
  char expected[] =
    "pow(2, 0/4) = '3ff0000000000000', '0' (1.000000000000000000)\n"
    "pow(2, 1/4) = '3ff306fe0a31b715', '3fd0000000000000' (1.189207115002721027)\n"
    "pow(2, 2/4) = '3ff6a09e667f3bcd', '3fe0000000000000' (1.414213562373095145)\n"
    "pow(2, 3/4) = '3ffae89f995ad3ad', '3fe8000000000000' (1.681792830507429004)\n";
  static char buf[1000];
  int offset = 0;
  for (int i = 0; i < 4; i++) {
    offset += snprintf(buf + offset, sizeof(buf)/sizeof(char),
        "pow(2, %d/4) = %s\n", i, double_to_hex(pow(2, i/4.0)));
  }
  if (!roundoff_ok()) {
    printf("Error --  unexpected roundoff error!\n");
    exit(-2);
  }
  if (strncmp(expected, buf, strlen(expected)) == 0) {
    printf("%s", buf);
    printf("Powers above are correct.\n");
    exit(0);
  } else {
    printf("%s", buf);
    printf("Error --  powers above are not correct!\n");
  }
  exit(-1);
}
