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
    box.size[i]= (1 + 2 * EPS) * (box.size_digits[i] * shape[i] + HALFEPS * fabs(box.center_digits[i]));
  }
}

void compute_cover(Box& box)
{
  box.cover.lattice = ACJ(
      XComplex(box.center[3], box.center[0]),
      XComplex(box.size[3], box.size[0]),
      0.,
      0.
      );
  box.cover.loxodromic_sqrt = ACJ(
      XComplex(box.center[4], box.center[1]),
      0.,
      XComplex(box.size[4], box.size[1]),
      0.
      );
  box.cover.parabolic = ACJ(
      XComplex(box.center[5], box.center[2]),
      0.,
      0.,
      XComplex(box.size[5], box.size[2])
      );
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
      // Note, sign(center_digits) == sign(center), unless center == 0. Also, size is always >= 0. 
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
      // Note, sign(center_digits) == sign(center), unless center == 0. Also, size is always >= 0. 
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

  box.nearer.lattice = XComplex(m[3], m[0]);
  box.nearer.loxodromic_sqrt = XComplex(m[4], m[1]);
  box.nearer.parabolic = XComplex(m[5], m[2]);
}

void compute_further(Box& box)
{
  double m[DIM];
  for (int i = 0; i < DIM; ++i) {
    m[i] = 0; // inconclusive cases
    if (box.center_digits[i] > -box.size_digits[i]) { // true sum is positive 
      // Want upper bound of true_center + true_size. Assume no overflow or underflow
      // Note, sign(center_digits) == sign(center), unless center == 0. Also, size is always >= 0. 
      // GMT paper page 419 of Annals gives with true arithmetic
      //      true_center + true_size <= center + size
      // By IEEE (+) and (-) resepct <= and >=, so center (+) size >=0 and
      // Lemma 7.0 for floating point arithmetic gives and upper bound
      //      (1+EPS)(*)(center (+) size) >= center + size >= true_center + true_size
      m[i] = (1 + EPS) * (box.center[i] + box.size[i]);
    } else { // true sum is <= 0
      // Want lower bound of true_center - true_size. Assume no overflow or underflow
      // Note, sign(center_digits) == sign(center), unless center == 0 
      // GMT paper page 419 of Annals gives with true arithmetic
      //      center - size <= true_center - true_size
      // By IEEE, (+) and (-) respects <= and >=, and negation is exact.
      // Thus, (-center) (+) size >=0 and Lemma 7.0 for floating point arithmetic gives
      //        (1+EPS)(*)( (-center) (+) size) ) >= (-center) + size
      // So,
      //      -((1+EPS)(*)( (-center) (+) size) ))<= center - size <= true_center - true_size
      m[i] = -((1 + EPS) * ((-box.center[i]) + box.size[i]));
    }
  }

  box.further.lattice = XComplex(m[3], m[0]);
  box.further.loxodromic_sqrt = XComplex(m[4], m[1]);
  box.further.parabolic = XComplex(m[5], m[2]);
}

void compute_greater(Box& box)
{
  double m[DIM];
  for (int i = 0; i < DIM; ++i) {
    m[i] = 0; // inconclusive cases
    if (box.center_digits[i] > -box.size_digits[i]) { // true sum is positive
      // Want upper bound of true_center + true_size. Assume no overflow or underflow
      // Note, sign(center_digits) == sign(center), unless center == 0. Also, size is always >= 0. 
      // GMT paper page 419 of Annals gives with true arithmetic
      //      true_center + true_size <= center + size.
      // Notice that center + size >= true_center + true_size > 0.
      // By IEEE, center (+) size >=0, as it's guanrateed to evaluate to nearest representable.
      // Lemma 7.0 for floating point arithmetic gives and upper bound
      //      (1+EPS)(*)(center (+) size) >= center + size >= true_center + true_size
      m[i] = (1 + EPS) * (box.center[i] + box.size[i]);
    } else if (box.center_digits[i] < -box.size_digits[i] && // true sum is negative
        box.center[i]        < -box.size[i]) { // machine sum is <= 0
      // Want upper bound of true_center + true_size. Assume no overflow or underflow
      // Note, sign(center_digits) == sign(center), unless center == 0. Also, size is always >= 0. 
      // GMT paper page 419 of Annals gives with true arithmetic
      //      true_center + true_size <= center + size.
      // Notice that center + size < 0.
      // By IEEE, center (+) size <= 0, as it's guanrateed to evaluate to nearest representable.
      // Lemma 7.0 for floating point arithmetic gives a bound
      //      (1-EPS)(*)| center (+) size | < | center + size |
      // So,
      //      -((1-EPS)(*)(-(center (+) size))) >= center + size >= true_center + true_size
      m[i] = -((1 - EPS) * (-(box.center[i] + box.size[i])));
    }
  }

  box.greater.lattice = XComplex(m[3], m[0]);
  box.greater.loxodromic_sqrt = XComplex(m[4], m[1]);
  box.greater.parabolic = XComplex(m[5], m[2]);
}

Box build_box(char* where) {
  Box box;
  // Global scaling of boxes - runs once
  if (!shape_initialized) {
    for (int i = 0; i < DIM; ++i) {
      shape[i] = pow(2, -i / DIM.0);
    }
    shape_initialized = true;
  } 
  for (int i = 0; i < DIM; ++i) {
    box.center_digits[i] = 0;
    box.size_digits[i] = 8;
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
  compute_further(box);
  compute_greater(box);
  return box;    
}


void fill_derived(ParamsAJCC& p) {
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

double shape[DIM];
static bool shape_initialized = false; 
Box::Box() {
  if (!shape_initialized) {
    shape_initialized = true;
    for (int i = 0; i < DIM; ++i) {
      shape[i] = pow(2, -i / float(DIM));
    }
  }
  for (int i = 0; i < DIM; ++i) {
    center_digits[i] = 0;
    size_digits[i] = SCL;
  }
  pos = 0;
  compute_center_and_size();
  compute_nearer();
  compute_cover();
}

Box Box::child(int dir) const
{
  Box child(*this);
  child.size_digits[pos] *= 0.5;
  child.center_digits[pos] += (2*dir-1)*child.size_digits[pos];
  ++child.pos;
  if (child.pos == DIM) { child.pos = 0; }

  child.name = name;
  child.name.append(1, '0'+dir);

  child.qr = qr;
  child.short_words_cache.clear();

  child.compute_center_and_size();
  child.compute_nearer();
  child.compute_cover();
  return child;
}

Box get_box(std::string code) {
  Box box;
  for (char dir : code) {
    if (dir == '0') {
      box = box.child(0);
    } else if (dir == '1') {
      box = box.child(1);
    }
  }
  return box;
}

// This is now special to the box mapping
std::string Box::desc() {
  AJCC sinhL2 = _cover.sinhL2;
  AJCC sinhD2 = _cover.sinhD2;
  AJCC coshL2 = _cover.coshL2;
  AJCC coshD2 = _cover.coshD2;
  Complex c_sinhL2 = _center.sinhL2;
  Complex c_sinhD2 = _center.sinhD2;

  char _desc[10000];
  sprintf(_desc, "%s\n", name.c_str());
  sprintf(_desc + strlen(_desc), "sinh(L/2) = %f + i %f with size %f, absLB %f, and absUB %f\n",
                                  sinhL2.f.re, sinhL2.f.im, sinhL2.size, absLB(sinhL2), absUB(sinhL2));
  sprintf(_desc + strlen(_desc), "cosh(L/2) = %f + i %f with size %f, absLB %f, and absUB %f\n",
                                  coshL2.f.re, coshL2.f.im, coshL2.size, absLB(coshL2), absUB(coshL2));
  sprintf(_desc + strlen(_desc), "sinh(D/2) = %f + i %f with size %f, absLB %f, and absUB %f\n",
                                  sinhD2.f.re, sinhD2.f.im, sinhD2.size, absLB(sinhD2), absUB(sinhD2));
  sprintf(_desc + strlen(_desc), "cosh(D/2) = %f + i %f with size %f, absLB %f, and absUB %f\n",
                                  coshD2.f.re, coshD2.f.im, coshD2.size, absLB(coshD2), absUB(coshD2));

  sprintf(_desc + strlen(_desc),
      "Center\n    sinh(L/2) %f + i %f | sinh(D/2) %f + i %f\n",
      c_sinhL2.real(), c_sinhL2.imag(),  c_sinhD2.real(), c_sinhD2.imag());

  std::string s(_desc);
  return s;
}

void Box::compute_center_and_size()
{
  for (int i = 0; i < DIM; ++i) {
    // GMT paper page 419 of Annals
    // box_size guarantees that :
    // box_center - box_size <= true_center - true_size
    // box_center + box_size >= true_center + true_size
    // where box operations are floating point. 
    box_center[i] = shape[i]*center_digits[i];
    box_size[i]= (1+2*EPS)*(size_digits[i]*shape[i]+HALFEPS*fabs(center_digits[i]));
  }
  _center.sinhL2 = Complex(box_center[1], box_center[3]);
  _center.sinhD2 = Complex(box_center[0], box_center[2]);

  fill_derived(_center);

  _x_center = construct_x(_center); 
  _y_center = construct_y(_center); 
  _center.coshmu = cosh_move_j(_x_center); 
}

void Box::compute_cover()
{
  // Let A = { (z0,z1) \in C^2 | |zi| <= 1 }
  // Our parameters are functions on A with the following defintions
  // sinh(L/2) = (s[1] + i s[3]) z0 + (c[1] + i c[3]) 
  // sinh(D/2) = (s[0] + i s[2]) z1 + (c[0] + i c[2])

  _cover.sinhL2 = AJCC(XComplex(box_center[1], box_center[3]), 
                     XComplex(box_size[1], box_size[3]), 0, 0);

  _cover.sinhD2 = AJCC(XComplex(box_center[0], box_center[2]), 
                     0, XComplex(box_size[0], box_size[2]), 0);

  fill_derived(_cover);

  _x_cover = construct_x(_cover); 
  _y_cover = construct_y(_cover); 
  _cover.coshmu = cosh_move_j(_x_cover); 
}

void Box::compute_nearer()
{
	double m[DIM];
	for (int i = 0; i < DIM; ++i) {
        m[i] = 0; // inconclusive cases
        if (center_digits[i] > 0 && // center is positive 
            center_digits[i] > size_digits[i] &&  // true diff is positive
            box_center[i]    > box_size[i]) { // machine diff is >= 0
            // Want lower bound on true_center - true_size.  Assume no overflow or underflow 
            // Note, sign(center_digits) == sign(box_center), unless box_center == 0. Also, box_size is always >= 0. 
            // GMT paper page 419 of Annals gives with true arithmetic
            //      box_center - box_size <= true_center - true_size
            // Now, in machine arthimetric, by IEEE, if 
            //      box_center > box_size then box_center (-) box_size >= 0.
            // Lemma 7 gives,
            //      (1-EPS)(*)( box_center (-) box_size ) <= box_center - box_size <= true_center - box_size. 
            m[i] = (1-EPS)*(box_center[i] - box_size[i]);
        } else if (center_digits[i] < 0 && // center is negative
                   center_digits[i] < -size_digits[i] && // true sum is negative
                   box_center[i]    < -box_size[i]) {  // machine sum is negative
            // Want upper bound on true_center - true_size.  Assume no overflow or underflow
            // Note, sign(center_digits) == sign(box_center), unless box_center == 0. Also, box_size is always >= 0. 
            // GMT paper page 419 of Annals gives with true arithmetic
            //      true_center + true_size <= box_center + box_size
            // Now, in machine arthimetric, by IEEE, if 
            //      -box_center > box_size then (-box_center) (-) box_size >= 0.
            // Lemma 7 gives,
            //      (1-EPS)(*)( (-box_center) (-) box_size ) <= -box_center - box_size <= -true_center - true_size.
            // So,
            //      -((1-EPS)(*)( (-box_center) (-) box_size )) >= true_center + true_size.
            // Note, negation is exact for machine numbers
            m[i] = -((1-EPS)*((-box_center[i]) - box_size[i]));
        }
	}
	
  _nearer.sinhL2 = Complex(m[1], m[3]);
  _nearer.sinhD2 = Complex(m[0], m[2]);

  fill_derived(_nearer);
}
