#ifndef _SL2AJCC_
#define _SL2AJCC_

#include "types.hh"

struct SLAJCC {
  SLAJCC() : a(1), b(0), c(0), d(1) {}
  SLAJCC(const AJCC &aa, const AJCC &bb, const AJCC &cc, const AJCC &dd) :
    a(aa), b(bb), c(cc), d(dd) {}
  AJCC a, b, c, d;
};

const SLAJCC operator*(const SLAJCC &M, const SLAJCC &N) {
  return SLAJCC(
      M.a*N.a+M.b*N.c, M.a*N.b+M.b*N.d,
      M.c*N.a+M.d*N.c, M.c*N.b+M.d*N.d);
};

const SLAJCC inverse(const SLAJCC &M) {
  return SLAJCC(M.d,-M.b,-M.c,M.a);
};

const AJCC dist(const SLAJCC &M1, const SLAJCC &M2) {
  return abs(M1.a - M2.a) + abs(M1.b - M2.b) + abs(M1.c - M2.c) + abs(M1.d - M2.d); 
};

inline const SLAJCC pow(const SLAJCC &M, int n) {
  SLAJCC A; // identity
  if (n == 0) { return A; }
  SLAJCC B;
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

#endif // _SL2AJCC_
