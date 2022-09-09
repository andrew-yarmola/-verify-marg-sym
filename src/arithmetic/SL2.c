#include "SL2.h"

void print_type(const SL2AJCC &M) {
  fprintf(stderr, "a = ");
  print_type(M.a);
  fprintf(stderr, "b = ");
  print_type(M.b);
  fprintf(stderr, "c = ");
  print_type(M.c);
  fprintf(stderr, "d = ");
  print_type(M.d);
}

const SL2AJCC operator*(const SL2AJCC &M, const SL2AJCC &N) {
  return SL2AJCC(
      M.a*N.a+M.b*N.c, M.a*N.b+M.b*N.d,
      M.c*N.a+M.d*N.c, M.c*N.b+M.d*N.d);
};

const SL2AJCC inverse(const SL2AJCC &M) {
  return SL2AJCC(M.d,-M.b,-M.c,M.a);
};

const AJCC dist(const SL2AJCC &M1, const SL2AJCC &M2) {
  return abs(M1.a - M2.a) + abs(M1.b - M2.b) + abs(M1.c - M2.c) + abs(M1.d - M2.d); 
};

const SL2AJCC pow(const SL2AJCC &M, int n) {
  SL2AJCC A; // identity
  if (n == 0) { return A; }
  SL2AJCC B;
  if (n < 0) { 
    B = inverse(M);
    n = -n;
  } else {
    B = M;
  }
  while (n > 1) {
    if (n & 1) { // n odd
      A = A*B;
    }
    B = B*B;
    n /= 2; // int division
  } 
  return A*B;
};
