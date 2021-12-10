#ifndef _box_h_
#define _box_h_
#include "AJCC.h"

// Our root box has center (0,0,0,0) and sides of length
// 2 * (2, 2^(3/4), 2^(2/4), 2^(1/4)). The last is > 2.37
#define DIM 4
#define SIZ 2

typedef struct {
  XComplex sinhL2;
  XComplex sinhD2;
} ParamsX;

typedef struct {
  T sinhL2;
  T sinhD2;
  // derived parameters
  T sinhsqL2;
  T sinhsqD2;
  T coshsqL2;
  T coshsqD2;
  T coshL2;
  T coshD2;
  T expD2;
  T expmD2;
  T twocoshreD2;
  T twosinhreD2;
  T coshreD;
  T coshmu;
  T coshreL;
  T sinhreL;
  T cosimL;
} ParamsAJCC;

typedef struct {
  double center_digits[DIM];
  double size_digits[DIM];
  double center[DIM];
  double size[DIM];
  ACJParams cover;
  XParams nearer; // all values closer to 0 than in box or 0 if box overlaps
} Box;

Box build_box(char* where);

#endif // _box_h_
