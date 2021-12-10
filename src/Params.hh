#ifndef __Params_h
#define __Params_h
#include <math.h>
#include <string>
#include "SL2.hh"
#include "Generators.hh"
#include "roundoff.h"
#include "types.hh"
#include "assert.h"

// Exp of distance from ortho endpts to the special point on axis
template<typename T>
double exp_dist_to_ortho_x(const T& zm, const T& zp, const Params<T>& p) {
  T one(1.0); 
  T z = ((p.expD2 * zm + one) * (p.expD2 * zp + one))/
    ((p.expD2 * zm - one) * (p.expD2 * zp - one));
  return absUB(sqrt(z));
}

// Exp of distance from ortho endpts to the special point on axis
template<typename T>
double exp_dist_to_ortho_y(const T& zm, const T& zp, const Params<T>& p) {
  T z = ((p.expD2 + zm) * (p.expD2 + zp))/
    ((p.expD2 - zm) * (p.expD2 - zp));
  return absUB(sqrt(z));
}

// Eliminate bad boxes that can't generate non-elementay groups
template<typename T>
const T jorgensen_xy(const Params<T>& p) {
  T z = p.sinhL2 * p.sinhperp; 
  return (abs_sqrd(z) + 1) * abs_sqrd(p.sinhL2) * 4;
}

// Eliminate bad boxes that can't generate non-elementay groups
template<typename T>
const T jorgensen_yx(const Params<T>& p) {
  T z = p.sinhL2 * p.sinhperp; 
  return (abs_sqrd(z) + 1) * abs_sqrd(p.sinhL2) * 4;
}

template<typename T>
const T jorgensen_xw(const SL2<T>& w, const Params<T>& p) {
  T shLx2 = p.sinhL2;
  T td = w.a - w.d;
  T z = w.c * p.expmD2 - w.b * p.expD2;
  return (abs(td * td - z * z) + 4) * abs_sqrd(shLx2);
}

template<typename T>
const T jorgensen_wx(const SL2<T>& w, const Params<T>& p) {
  T shLx2 = p.sinhL2;
  T tr = w.a + w.d;
  T td = w.a - w.d;
  T z = w.c * p.expmD2 - w.b * p.expD2;
  return abs(tr * tr - 4) + abs(td * td - z * z) * abs_sqrd(shLx2);
}

template<typename T>
const T jorgensen_yw(const SL2<T>& w, const Params<T>& p) {
  T shLy2 = p.sinhL2;
  T td = w.a - w.d;
  T z = w.c * p.expD2 - w.b * p.expmD2;
  return (abs(td * td - z * z) + 4) * abs_sqrd(shLy2);
}

template<typename T>
const T jorgensen_wy(const SL2<T>& w, const Params<T>& p) {
  T shLy2 = p.sinhL2;
  T tr = w.a + w.d;
  T td = w.a - w.d;
  T z = w.c * p.expD2 - w.b * p.expmD2;
  return abs(tr * tr - 4) + abs(td * td - z * z) * abs_sqrd(shLy2);
}

// Complex distance between axis(x) and w(axis(x)) 
template<typename T>
const T four_sinh_perp2_sq_ax_wax(const SL2<T>& w, const Params<T>& p) {
  T td = w.a - w.d;
  T zm = w.c * p.expmD2 - w.b * p.expD2;
  // formula by using crossratios
  T four_sinh_sq_perp2 = td * td - zm * zm;  
  return four_sinh_sq_perp2; 
}

// Distance between axis(x) and w(axis(x)) 
template<typename T>
const T four_cosh_dist_ax_wax(const SL2<T>& w, const Params<T>& p) {
  T four_sinh_sq_perp2 = four_sinh_perp2_sq_ax_wax(w, p);
  return abs(four_sinh_sq_perp2 + 4) + abs(four_sinh_sq_perp2);
}

// Complex distance between axis(y) and w(axis(y)) 
template<typename T>
const T four_sinh_perp2_sq_ay_way(const SL2<T>& w, const Params<T>& p) {
  T td = w.a - w.d;
  T zm = w.c * p.expD2 - w.b * p.expmD2;
  // formula by using crossratios
  T four_sinh_sq_perp2 = td * td - zm * zm;  
  return four_sinh_sq_perp2; 
}

// Distance between axis(y) and w(axis(y)) 
template<typename T>
const T four_cosh_dist_ay_way(const SL2<T>& w, const Params<T>& p) {
  T four_sinh_sq_perp2 = four_sinh_perp2_sq_ay_way(w, p); 
  return abs(four_sinh_sq_perp2 + 4) + abs(four_sinh_sq_perp2);
}

// Complex distance between axis(x) and w(axis(y)) 
template<typename T>
const T four_sinh_perp2_sq_ax_way(const SL2<T>& w, const Params<T>& p) {
  T z = ((w.a * w.a) * p.expD2  - (w.b * w.b) * p.expmD2) * p.expD2 +
    ((w.d * w.d) * p.expmD2 - (w.c * w.c) * p.expD2 ) * p.expmD2;
  return z - 2;
}

// Distance between axis(x) and w(axis(y)) 
template<typename T>
const T four_cosh_dist_ax_way(const SL2<T>& w, const Params<T>& p) {
  T z = ((w.a * w.a) * p.expD2  - (w.b * w.b) * p.expmD2) * p.expD2 +
    ((w.d * w.d) * p.expmD2 - (w.c * w.c) * p.expD2 ) * p.expmD2;
  return  abs(z - 2) + abs(z + 2);
}

// Complex distance between axis(y) and w(axis(x)) 
template<typename T>
const T four_sinh_perp2_sq_ay_wax(const SL2<T>& w, const Params<T>& p) {
  T z = ((w.a * w.a) * p.expmD2 - (w.b * w.b) * p.expD2 ) * p.expmD2 +
    ((w.d * w.d) * p.expD2  - (w.c * w.c) * p.expmD2) * p.expD2;
  return  z - 2;
}

// Distance between axis(y) and w(axis(x)) 
template<typename T>
const T four_cosh_dist_ay_wax(const SL2<T>& w, const Params<T>& p) {
  T z = ((w.a * w.a) * p.expmD2 - (w.b * w.b) * p.expD2 ) * p.expmD2 +
    ((w.d * w.d) * p.expD2  - (w.c * w.c) * p.expmD2) * p.expD2;
  return  abs(z - 2) + abs(z + 2);
}

// Complex distance between w(axis(x)) and {0, infty} 
template<typename T>
const T two_sinh_perp2_sq_wax_zero_inf(const SL2<T>& w, const Params<T>& p) {
  T one = T(1);
  return (w.b * w.d) * p.expD2 - (w.a * w.c) * p.expmD2 - one; 
}

// Compley distance between w(ayis(y)) and {0, infty} 
template<typename T>
const T two_sinh_perp2_sq_way_zero_inf(const SL2<T>& w, const Params<T>& p) {
  T one = T(1);
  return (w.b * w.d) * p.expmD2 - (w.a * w.c) * p.expD2 - one; 
}

// Complex distance between w(axis(x)) and {-1, 1} 
template<typename T>
const T four_sinh_perp2_sq_wax_mp_one(const SL2<T>& w, const Params<T>& p) {
  T two = T(2);
  return (w.d * w.d - w.b * w.b) * p.expD2 + (w.a * w.a - w.c * w.c) * p.expmD2 - two; 
}

// Compley distance between w(ayis(y)) and {-1, 1} 
template<typename T>
const T four_sinh_perp2_sq_way_mp_one(const SL2<T>& w, const Params<T>& p) {
  T two = T(2);
  return (w.d * w.d - w.b * w.b) * p.expmD2 + (w.a * w.a - w.c * w.c) * p.expD2 - two; 
}

// Complex distance between w(axis(x)) and {-I, I} 
template<typename T>
const T four_sinh_perp2_sq_wax_mp_iye(const SL2<T>& w, const Params<T>& p) {
  T two = T(2);
  T iye = eye(T(1));
  return ((w.d * w.d + w.b * w.b) * p.expD2 - (w.a * w.a + w.c * w.c) * p.expmD2) * iye - two; 
}

// Complex distance between w(axis(y)) and {-I, I} 
template<typename T>
const T four_sinh_perp2_sq_way_mp_iye(const SL2<T>& w, const Params<T>& p) {
  T two = T(2);
  T iye = eye(T(1));
  return ((w.d * w.d + w.b * w.b) * p.expmD2 - (w.a * w.a + w.c * w.c) * p.expD2) * iye - two; 
}

#endif // __Params_h
