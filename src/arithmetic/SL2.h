#ifndef _SL2AJCC_
#define _SL2AJCC_

#include "AJCC.h"

struct SL2AJCC {
  AJCC a, b, c, d;
  SL2AJCC(const AJCC &aa, const AJCC &bb, const AJCC &cc, const AJCC &dd) :
    a(aa), b(bb), c(cc), d(dd) {}
  SL2AJCC() : a(1), b(0), c(0), d(1) {}
};

void print_type(const SL2AJCC &M);

void print_type(const char* desc, const SL2AJCC &M);

const SL2AJCC operator*(const SL2AJCC &M, const SL2AJCC &N);

const SL2AJCC inverse(const SL2AJCC &M);

const AJCC dist(const SL2AJCC &M1, const SL2AJCC &M2);

const SL2AJCC pow(const SL2AJCC &M, int n);

#endif // _SL2AJCC_
