#include <stdlib.h>
#include "box.h"


double shape[DIM];
static bool shape_initialized = false; 

void compute_center_and_size(Box& box)
{
  for (int i = 0; i < DIM; ++i) {
    // GMT paper page 419 of Annals
    // size guarantees that :
    // center - size <= true_center - true_size
    // center + size >= true_center + true_size
    // where box operations are floating point. 
    box.center[i] = shape[i] * box.center_digits[i];
    box.size[i]= (1 + 2 * EPS) * (box.size_digits[i] * shape[i]
        + HALFEPS * fabs(box.center_digits[i]));
  }
}

void fill_derived(AJCCParams& p) {
  AJCC one = AJCC(1);
  p.sinhsqL2 = p.sinhL2 * p.sinhL2;
  p.coshsqL2 = p.sinhsqL2 + one;
  p.coshL2 = sqrt(p.coshsqL2); // TODO check if correct branch
  p.sinhsqD2 = p.sinhD2 * p.sinhD2;
  p.coshsqD2 = p.sinhsqD2 + one;
  p.coshD2 = sqrt(p.coshsqD2); // TODO check if correct branch
  
  p.expD2 = p.coshD2 + p.sinhD2; 
  p.expmD2 = p.coshD2 - p.sinhD2;

  p.twocoshreD2 = abs(p.expD2) + abs(p.expmD2);
  p.twosinhreD2 = abs(p.expD2) - abs(p.expmD2);
  p.coshreD = abs(p.sinhsqD2) + abs(p.coshsqD2);
  p.coshreL = abs(p.sinhsqL2) + abs(p.coshsqL2);
  p.sinhreL = sqrt(p.coshreL * p.coshreL - one);
  p.cosimL = abs(p.coshsqL2) - abs(p.sinhsqL2);
}

void compute_cover(Box& box)
{
  box.cover.sinhL2 = AJCC(XComplex(box.center[1], box.center[3]), 
                     XComplex(box.size[1], box.size[3]), 0, 
                     0, 0,
                     0);

  box.cover.sinhD2 = AJCC(XComplex(box.center[0], box.center[2]), 
                     0, XComplex(box.size[0], box.size[2]),
                     0, 0,
                     0);

  fill_derived(box.cover);
}

void compute_nearer(Box& box)
{
  double m[DIM];
  for (int i = 0; i < DIM; ++i) {
    m[i] = 0; // inconclusive cases
    if (box.center_digits[i] > 0 && // center is positive 
        box.center_digits[i] > box.size_digits[i] &&  // true diff is positive
        box.center[i]        > box.size[i]) { // machine diff is >= 0
      // Want lower bound on true_center - true_size.  Assume no overflow or underflow 
      // Note, sign(center_digits) == sign(center), unless center == 0.
      // Also, size is always >= 0. 
      // GMT paper page 419 of Annals gives with true arithmetic
      //      center - size <= true_center - true_size
      // Now, in machine arthimetric, by IEEE, if 
      //      center > size then center (-) size >= 0.
      // Lemma 7.0 gives,
      //      (1-EPS)(*)( center (-) size ) <= center - size <= true_center - size. 
      m[i] = (1 - EPS) * (box.center[i] - box.size[i]);
    } else if (box.center_digits[i] < 0 && // center is negative
        box.center_digits[i] < -box.size_digits[i] && // true sum is negative
        box.center[i]        < -box.size[i]) {  // machine sum is negative
      // Want upper bound on true_center - true_size.  Assume no overflow or underflow
      // Note, sign(center_digits) == sign(center), unless center == 0.
      // Also, size is always >= 0. 
      // GMT paper page 419 of Annals gives with true arithmetic
      //      true_center + true_size <= center + size
      // Now, in machine arthimetric, by IEEE, if 
      //      -center > size then (-center) (-) size >= 0.
      // Lemma 7.0 gives,
      //      (1-EPS)(*)( (-center) (-) size ) <= -center - size <= -true_center - true_size.
      // So,
      //      -((1-EPS)(*)( (-center) (-) size )) >= true_center + true_size.
      // Note, negation is exact for machine numbers
      m[i] = -(( 1 - EPS) * ((-box.center[i]) - box.size[i]));
    }
  }

  box.nearer.sinhL2 = XComplex(m[1], m[3]);
  box.nearer.sinhD2 = XComplex(m[0], m[2]);
}

Box build_box(const char* where) {
  Box box;
  // Global scaling of boxes - runs once
  if (!shape_initialized) {
    for (int i = 0; i < DIM; ++i) {
      shape[i] = pow(2, -i / float(DIM));
    }
    shape_initialized = true;
  } 
  for (int i = 0; i < DIM; ++i) {
    box.center_digits[i] = 0;
    box.size_digits[i] = SCL;
  }
  size_t pos = 0;
  size_t idx = 0;
  int dir;
  while (where[idx] != '\0') {
    if (where[idx] == '0') {
      dir = 0;
    } else if (where[idx] == '1') {
      dir = 1;
    } else {
      fprintf(stderr, "Fatal: boxcode is invalid %s\n", where);
      exit(DIM);
    }
    box.size_digits[pos] *= 0.5;
    box.center_digits[pos] += (2 * dir - 1) * box.size_digits[pos];
    ++pos;
    if (pos == DIM) {
      pos = 0;
    }
    ++idx;
  }
  compute_center_and_size(box);
  compute_cover(box);
  compute_nearer(box);
  return box;    
}

void print_box(const Box& box) {
  fprintf(stderr, "sinh(L/2) = %f + i %f with size %f, absLB %f, and absUB %f\n",
                                  box.cover.sinhL2.f.re, box.cover.sinhL2.f.im, box.cover.sinhL2.size, absLB(box.cover.sinhL2), absUB(box.cover.sinhL2));
  fprintf(stderr, "cosh(L/2) = %f + i %f with size %f, absLB %f, and absUB %f\n",
                                  box.cover.coshL2.f.re, box.cover.coshL2.f.im, box.cover.coshL2.size, absLB(box.cover.coshL2), absUB(box.cover.coshL2));
  fprintf(stderr, "sinh(D/2) = %f + i %f with size %f, absLB %f, and absUB %f\n",
                                  box.cover.sinhD2.f.re, box.cover.sinhD2.f.im, box.cover.sinhD2.size, absLB(box.cover.sinhD2), absUB(box.cover.sinhD2));
  fprintf(stderr, "cosh(D/2) = %f + i %f with size %f, absLB %f, and absUB %f\n",
                                  box.cover.coshD2.f.re, box.cover.coshD2.f.im, box.cover.coshD2.size, absLB(box.cover.coshD2), absUB(box.cover.coshD2));
}
