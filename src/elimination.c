#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "elimination.h"

extern double g_max_area;

char impossible[134][32] = {
    "xy",
    "Xy",
    "xyXY",
    "xYXy",
    "xYXY",
    "XyXy",
    "XYxY",
    "XYXY",
    "yxyx",
    "yXYx",
    "yXYX",
    "YxYx",
    "YXyX",
    "xyXyx",
    "XyxyX",
    "yXYxy",
    "xxyxxy",
    "xxYxxY",
    "XyXXXy",
    "xyxxyx",
    "xYxxYx",
    "xyxyxy",
    "xYXyXY",
    "XyxyXY",
    "XyxYxy",
    "XyXyXy",
    "XyXYxY",
    "XYxYXy",
    "XYXyxy",
    "XYXYXY",
    "XyyXyy",
    "YXXYXX",
    "YxyXyx",
    "YxYXyX",
    "YXyxyX",
    "yyxyxY",
    "yyXyXY",
    "xxyXyxx",
    "YXXYYXX",
    "yxyXXyx",
    "yxYXXYx",
    "YxyXXyx",
    "YxYXXYx",
    "YxyXyxY",
    "YxYXYxY",
    "YXyxyXY",
    "YXYxYXY",
    "YYxyxYY",
    "XXXYYYXX",
    "xYXXyXXY",
    "XYXXyXXY",
    "xyXyxyXy",
    "xYxyxyxY",
    "xYxyxYxy",
    "xYxyXyxY",
    "xYxYXYxY",
    "xYXYxYXY",
    "XyxyXyxy",
    "XYxYXYxY",
    "XYXyxyXY",
    "XYXyXyXY",
    "XYXYxYXY",
    "xYYxyxYY",
    "XyyyXyxy",
    "yXXYXXYY",
    "yxxyXyxx",
    "yxxYxYxx",
    "yXXYXYXX",
    "YxxyyxxY",
    "YxxYYxxY",
    "YXXyyXXY",
    "YXXYYXXY",
    "yxyxxxyx",
    "yxyXyxyX",
    "yXyxyXyx",
    "YxYXYxYX",
    "yXyyxyyX",
    "YYxYxxYx",
    "YYXYXXYX",
    "yyyXyxyX",
    "yyyxyxYY",
    "xxyxxyxxy",
    "xxYXXyXXY",
    "XXYXXYXXY",
    "xxYxyXyxY",
    "XXYxYXYxY",
    "XXYXyxyXY",
    "xxYYxyxYY",
    "xyxxyxxyx",
    "xYXXYYXXY",
    "XYxxYYxxY",
    "xYYxYYxYY",
    "XyyXyyXyy",
    "yxxyxxyxx",
    "YXXYXXYXX",
    "YxYYxYYxY",
    "yyxxyXyxx",
    "yyXyyXyyX",
    "XXXYxYXYxY",
    "xxYxxYYxxY",
    "xxYXXYYXXY",
    "XXyxxyyxxy",
    "XXYxxYYxxY",
    "XXYXXYYXXY",
    "xxYxYxxYxY",
    "xxYXYxxYXY",
    "XXyXyXXyXy",
    "XXYXYXXYXY",
    "xYxxYxYxxY",
    "XYXXYXYXXY",
    "xYXYxxYXYx",
    "XyXyXXyXyX",
    "xYxYxYxYxY",
    "XyXyyxyyXy",
    "XYXYYxYYXY",
    "XYYXYxYXYY",
    "YxxxYYxxxY",
    "YXXXyyXXXY",
    "YXXXYYXXXY",
    "yXXyXYXyXX",
    "YxYXyyXYxY",
    "YXyXYYXyXY",
    "YXYxyyxYXY",
    "YYXyXYYXyX",
    "XXyyXYxYXyy",
    "XXYYXyxyXYY",
    "xxxYxYxxxYxY",
    "xxxYXYxxxYXY",
    "XXXYxYXXXYxY",
    "XXXYXYXXXYXY",
    "XyXYXyXyXYXy",
    "XYXyXYXYXyXY",
    "xxxxyyxyyyxxxxyyxyyy",
    "xxYYYYxxxYxxYYYYxxxY"
};

char symmetric[14][32] = {
    "xxYxyxyxY",
    "XXYXyXyXY",
    "xxyyxxYxY",
    "XXyyXXYXY",
    "XXYYXXyXy",
    "XXYYxYxYY",
    "xYxyxxyxY",
    "XyXyXYXXY",
    "XyxyyxyXy",
    "XyyxxyyXy",
    "YxYXYXYxY",
    "YXYxYYxYX",
    "yxyyxyXyX",
    "XXXYYXXYXXYY"
};

// Helper functions

// If inequalities is false, crash the program
void check(bool inequalities, const char* where)
{
  if (!inequalities) {
    fprintf(stderr, "Fatal: verify error at %s\n", where);
    exit(3);
  }
}

// Replaces contents of code with the parsed word or pair
void parse_word(char* code)
{
  char buf[MAX_CODE_LEN];
  strncpy(buf, code, MAX_CODE_LEN);
  char * start = strchr(buf,'(');
  char * end = strchr(buf,')');
  size_t len = end - start - 1;
  strncpy(code, start+1, len);
  code[len] = '\0'; 
}

// Returns a word pair encoded in code
word_pair get_word_pair(const char* code)
{
  word_pair pair;
  char * start = strchr(code,'(');
  char * mid = strchr(code,',');
  char * end = strchr(code,')');
  size_t first_len = mid - start - 1;
  size_t second_len = end - mid - 1;
  strncpy(pair.first, start+1, first_len);
  strncpy(pair.second, mid+1, second_len);
  pair.first[first_len] = '\0';
  pair.second[second_len] = '\0';
  return pair;
}

// Generators and words

SL2AJCC construct_x(const AJCCParams& params) {
  return SL2AJCC(
      params.coshL2, params.expmD2 * params.sinhL2,
      params.expD2 * params.sinhL2, params.coshL2);
};

SL2AJCC construct_y(const AJCCParams& params) {
  return SL2AJCC(
      params.coshL2, params.expD2 * params.sinhL2,
      params.expmD2 * params.sinhL2, params.coshL2);
};

// SL2AJCC matrix constrution from parameters and a word.
// Note, floating point arithmetic is not commutative and so neither is AJCC.
// We have chosen a right to left order of multiplication for validation to pass.
// Different orders of multiplication produce slightly different rounding
// errors. Since our data is greedy, this order is necessary for all words.
SL2AJCC construct_word(const AJCCParams& params, const char* word) {
  SL2AJCC w; // identity
  SL2AJCC x = construct_x(params);
  SL2AJCC y = construct_y(params);

  char h;
  int x_pow = 0;
  int y_pow = 0;	
  size_t pos;
  for (pos = strlen(word); pos > 0; --pos) {
    char h = word[pos-1];
    switch(h) {
      case 'x': ++x_pow; break;
      case 'X': --x_pow; break;
      case 'y': ++y_pow; break;
      case 'Y': --y_pow; break;
    }
    if (y_pow != 0 && x_pow != 0) {
      if (h == 'y' || h == 'Y' ) {
        w = pow(x, x_pow) * w;
        x_pow = 0;
      } else {
        w = pow(y, y_pow) * w;
        y_pow = 0;
      }
    }
  }
  // Only one of these should be true
  if (x_pow != 0) { 
    w = pow(x, x_pow) * w;
  }
  if (y_pow != 0) { 
    w = pow(y, y_pow) * w;
  }
  return w;
}

// Elimination Tools

const AJCC jorgensen_xw(const SL2AJCC& w, const AJCCParams& p) {
  AJCC shLx2 = p.sinhL2;
  AJCC td = w.a - w.d;
  AJCC z = w.c * p.expmD2 - w.b * p.expD2;
  return (abs(td * td - z * z) + 4) * abs_sqrd(shLx2);
}

const AJCC jorgensen_wx(const SL2AJCC& w, const AJCCParams& p) {
  AJCC shLx2 = p.sinhL2;
  AJCC tr = w.a + w.d;
  AJCC td = w.a - w.d;
  AJCC z = w.c * p.expmD2 - w.b * p.expD2;
  return abs(tr * tr - 4) + abs(td * td - z * z) * abs_sqrd(shLx2);
}

const AJCC jorgensen_yw(const SL2AJCC& w, const AJCCParams& p) {
  AJCC shLy2 = p.sinhL2;
  AJCC td = w.a - w.d;
  AJCC z = w.c * p.expD2 - w.b * p.expmD2;
  return (abs(td * td - z * z) + 4) * abs_sqrd(shLy2);
}

const AJCC jorgensen_wy(const SL2AJCC& w, const AJCCParams& p) {
  AJCC shLy2 = p.sinhL2;
  AJCC tr = w.a + w.d;
  AJCC td = w.a - w.d;
  AJCC z = w.c * p.expD2 - w.b * p.expmD2;
  return abs(tr * tr - 4) + abs(td * td - z * z) * abs_sqrd(shLy2);
}

// 4 sinh^2(half complex distance between axis(x) and w(axis(x))) 
const AJCC four_sinh_perp2_sq_ax_wax(const SL2AJCC& w, const AJCCParams& p) {
  AJCC td = w.a - w.d;
  AJCC zm = w.c * p.expmD2 - w.b * p.expD2;
  // formula by using crossratios
  AJCC four_sinh_sq_perp2 = td * td - zm * zm;  
  return four_sinh_sq_perp2; 
}

// 4 cosh(real distance between axis(x) and w(axis(x))) 
const AJCC four_cosh_dist_ax_wax(const SL2AJCC& w, const AJCCParams& p) {
  AJCC four_sinh_sq_perp2 = four_sinh_perp2_sq_ax_wax(w, p);
  return abs(four_sinh_sq_perp2 + 4) + abs(four_sinh_sq_perp2);
}

// 4 sinh^2(half complex distance between axis(y) and w(axis(y))) 
const AJCC four_sinh_perp2_sq_ay_way(const SL2AJCC& w, const AJCCParams& p) {
  AJCC td = w.a - w.d;
  AJCC zm = w.c * p.expD2 - w.b * p.expmD2;
  // formula by using crossratios
  AJCC four_sinh_sq_perp2 = td * td - zm * zm;  
  return four_sinh_sq_perp2; 
}

// 4 cosh(real distance between axis(y) and w(axis(y))) 
const AJCC four_cosh_dist_ay_way(const SL2AJCC& w, const AJCCParams& p) {
  AJCC four_sinh_sq_perp2 = four_sinh_perp2_sq_ay_way(w, p); 
  return abs(four_sinh_sq_perp2 + 4) + abs(four_sinh_sq_perp2);
}

// 4 sinh^2(half complex distance between axis(x) and w(axis(y))) 
const AJCC four_sinh_perp2_sq_ax_way(const SL2AJCC& w, const AJCCParams& p) {
  AJCC z = ((w.a * w.a) * p.expD2  - (w.b * w.b) * p.expmD2) * p.expD2 +
    ((w.d * w.d) * p.expmD2 - (w.c * w.c) * p.expD2 ) * p.expmD2;
  return z - 2;
}

// 4 cosh(real distance between axis(x) and w(axis(y))) 
const AJCC four_cosh_dist_ax_way(const SL2AJCC& w, const AJCCParams& p) {
  AJCC z = ((w.a * w.a) * p.expD2  - (w.b * w.b) * p.expmD2) * p.expD2 +
    ((w.d * w.d) * p.expmD2 - (w.c * w.c) * p.expD2 ) * p.expmD2;
  return  abs(z - 2) + abs(z + 2);
}

// 4 sinh^2(half complex distance between axis(y) and w(axis(x))) 
const AJCC four_sinh_perp2_sq_ay_wax(const SL2AJCC& w, const AJCCParams& p) {
  return four_sinh_perp2_sq_ax_way(inverse(w), p);
}

// 4 cosh(real distance between axis(x) and w(axis(y))) 
const AJCC four_cosh_dist_ay_wax(const SL2AJCC& w, const AJCCParams& p) {
  return four_cosh_dist_ax_way(inverse(w), p);
}

// 2 sinh^2(half complex distance between w(axis(x)) and (0, inf)) 
inline const AJCC two_sinh_perp2_sq_wax_zero_inf(const SL2AJCC& w, const AJCCParams& p) {
  AJCC one = T(1);
  return (w.b * w.d) * p.expD2 - (w.a * w.c) * p.expmD2 - one; 
}

// 2 sinh^2(half complex distance between w(axis(y)) and (0, inf)) 
inline const AJCC two_sinh_perp2_sq_way_zero_inf(const SL2AJCC& w, const AJCCParams& p) {
  AJCC one = T(1);
  return (w.b * w.d) * p.expmD2 - (w.a * w.c) * p.expD2 - one; 
}

// 2 sinh^2(half complex distance between w(axis(x)) and (-1, 1)) 
inline const AJCC four_sinh_perp2_sq_wax_mp_one(const SL2AJCC& w, const AJCCParams& p) {
  AJCC two = T(2);
  return (w.d * w.d - w.b * w.b) * p.expD2 + (w.a * w.a - w.c * w.c) * p.expmD2 - two; 
}

// 2 sinh^2(half complex distance between w(axis(y)) and (-1, 1)) 
inline const AJCC four_sinh_perp2_sq_way_mp_one(const SL2AJCC& w, const AJCCParams& p) {
  AJCC two = T(2);
  return (w.d * w.d - w.b * w.b) * p.expmD2 + (w.a * w.a - w.c * w.c) * p.expD2 - two; 
}

// 2 sinh^2(half complex distance between w(axis(x)) and (-i, i)) 
inline const AJCC four_sinh_perp2_sq_wax_mp_iye(const SL2AJCC& w, const AJCCParams& p) {
  AJCC two = T(2);
  AJCC iye = eye(T(1));
  return ((w.d * w.d + w.b * w.b) * p.expD2 - (w.a * w.a + w.c * w.c) * p.expmD2) * iye - two; 
}

// 2 sinh^2(half complex distance between w(axis(y)) and (-i, i)) 
inline const AJCC four_sinh_perp2_sq_way_mp_iye(const SL2AJCC& w, const AJCCParams& p) {
  AJCC two = T(2);
  AJCC iye = eye(T(1));
  return ((w.d * w.d + w.b * w.b) * p.expmD2 - (w.a * w.a + w.c * w.c) * p.expD2) * iye - two; 
}

inline const bool moves_y_axis_too_close_to_x(const SL2AJCC& w, const AJCCParams& p) {
  AJCC diff = p.coshreD * 4 - four_cosh_dist_ax_way(w, p);
  // We know that diff is away from zero and the diff should be conj symmetrix, so
  // we only test if the real part is to one side of the bound
  return strictly_pos(diff);
}

inline const bool moves_x_axis_too_close_to_y(const SL2AJCC& w, const AJCCParams& p) {
  AJCC diff = p.coshreD * 4 - four_cosh_dist_ay_wax(w, p);
  // We know that diff is away from zero and the diff should be conj symmetrix, so
  // we only test if the real part is to one side of the bound
  return strictly_pos(diff);
}

inline const bool moved_y_axis_not_x_axis(const SL2AJCC& w, const AJCCParams& p) {
  AJCC fsp2sq = four_sinh_perp2_sq_ax_way(w, p);
  return absLB(fsp2sq) > 0; 
}

inline const bool moved_x_axis_not_y_axis(const SL2AJCC& w, const AJCCParams& p) {
  AJCC fsp2sq = four_sinh_perp2_sq_ay_wax(w, p);
  return absLB(fsp2sq) > 0; 
}

const AJCC cosh_move_j(const SL2AJCC& w) {
  AJCC q = abs_sqrd(w.c) + abs_sqrd(w.d);
  AJCC z = w.a * conj(w.c) + w.b * conj(w.d);
  return (abs_sqrd(z) + (q - 1) * (q - 1))/(q * 2) + 1; 
}

inline const bool not_identity(const SL2AJCC& w) {
  return absLB(w.b) > 0 ||  absLB(w.c) > 0 ||
    ((absLB(w.a-1) > 0 || absLB(w.d-1) > 0) && (absLB(w.a+1) > 0 || absLB(w.d+1) > 0));
}

inline const bool tube_hits_axis_two(const T& two_sinh_p2sq, const T& two_cosh_re_tube) {
  AJCC tcd = two_cosh_dist(two_sinh_p2sq);
  // sinh(I Pi/4)^2 = -1/2 which means axes meet othrogonally
  return strictly_pos(two_cosh_re_tube - tcd) && absLB(two_sinh_p2sq + 1) > 0;
}

inline const bool tube_hits_axis_four(const T& four_sinh_p2sq, const T& two_cosh_re_tube) {
  AJCC fcd = four_cosh_dist(four_sinh_p2sq);
  // sinh(I Pi/4)^2 = -1/2 which means axes meet othrogonally
  return strictly_pos(two_cosh_re_tube * 2 - fcd) && absLB(four_sinh_p2sq + 2) > 0;
}

inline const bool wx_hits_sym_axis(const SL2AJCC& w, const AJCCParams& p) {
  AJCC tsp2sq_inf = two_sinh_perp2_sq_wax_zero_inf(w, p);
  AJCC fsp2sq_one = four_sinh_perp2_sq_wax_mp_one(w, p);
  AJCC fsp2sq_iye = four_sinh_perp2_sq_wax_mp_iye(w, p);
  return (tube_hits_axis_two(tsp2sq_inf, p.twocoshreD2) ||
      tube_hits_axis_four(fsp2sq_one, p.twocoshreD2) || 
      tube_hits_axis_four(fsp2sq_iye, p.twocoshreD2));
}

inline const bool wy_hits_sym_axis(const SL2AJCC& w, const AJCCParams& p) {
  AJCC tsp2sq_inf = two_sinh_perp2_sq_way_zero_inf(w, p);
  AJCC fsp2sq_one = four_sinh_perp2_sq_way_mp_one(w, p);
  AJCC fsp2sq_iye = four_sinh_perp2_sq_way_mp_iye(w, p);
  return (tube_hits_axis_two(tsp2sq_inf, p.twocoshreD2) ||
      tube_hits_axis_four(fsp2sq_one, p.twocoshreD2) || 
      tube_hits_axis_four(fsp2sq_iye, p.twocoshreD2));
}

inline const bool does_not_fix_sym_axis(const SL2AJCC& w) {
  AJCC one(1);
  AJCC zero(0);
  AJCC iye(0,1);
  return (absLB(w.b) > 0 && absLB(w.d) > 0) 
    || (absLB(w.a) > 0 && absLB(w.c) > 0) ||
    (absLB(mobius(w, one) - one) > 0 && absLB(mobius(w, one) + one) > 0) || 
    (absLB(mobius(w, iye) - iye) > 0 && absLB(mobius(w, iye) + iye) > 0) || 
    (absLB(mobius(w, -one) - one) > 0 && absLB(mobius(w, -one) + one) > 0) || 
    (absLB(mobius(w, -iye) - iye) > 0 && absLB(mobius(w, -iye) + iye) > 0); 
}

inline const bool cant_fix_x_axis(const SL2AJCC& w, const AJCCParams& p) {
  AJCC fsp2sq = four_sinh_perp2_sq_ax_wax(w, p);
  // return absLB(fsp2sq) > 0 && absLB(fsp2sq + 4) > 0; 
  return absLB(fsp2sq) > 0; 
}

inline const bool must_fix_x_axis(const SL2AJCC& w, const AJCCParams& p) {
  AJCC diff = p.coshreD * 4 - four_cosh_dist_ax_wax(w, p);
  // We know that diff is away from zero and the diff should be conj symmetrix, so
  // we only test if the real part is to one side of the bound
  return strictly_pos(diff);
}

inline const bool cant_fix_y_axis(const SL2AJCC& w, const AJCCParams& p) {
  AJCC fsp2sq = four_sinh_perp2_sq_ay_way(w, p);
  // return absLB(fsp2sq) > 0 && absLB(fsp2sq + 4) > 0; 
  return absLB(fsp2sq) > 0; 
}

inline const bool must_fix_y_axis(const SL2AJCC& w, const AJCCParams& p) {
  AJCC diff = p.coshreD * 4 - four_cosh_dist_ay_way(w, p);
  // We know that diff is away from zero and the diff should be conj symmetrix, so
  // we only test if the real part is to one side of the bound
  return strictly_pos(diff);
}

inline const bool inside_var_nbd_x(const SL2AJCC& w, const AJCCParams& params) {
  return absUB(jorgensen_wx(w, params)) < 1 || absUB(jorgensen_xw(w, params)) < 1 || must_fix_x_axis(w, params);
}

inline const bool inside_var_nbd_y(const SL2AJCC& w, const AJCCParams& params) {
  return absUB(jorgensen_wy(w, params)) < 1 || absUB(jorgensen_yw(w, params)) < 1 || must_fix_y_axis(w, params);
}

inline bool move_less_than_marg(const SL2AJCC& w, const AJCCParams& p) {
  AJCC diff = p.coshmu - cosh_move_j(w);
  return strictly_pos(diff); 
}

inline bool non_cyclic_power(const SL2AJCC& w, const SL2AJCC& x_or_y) {
  // Assume word fixes the same axis as x or y, so it must live in a cyclic group with x or y.
  // Here we check that this is impossible in this box. Must use margulis
  // number to check cut off for roots of x or y
  SL2AJCC commutator = x_or_y * w * inverse(w * x_or_y);
  // TODO Test powers when coshmu > coshsdx + sinhsdx
  return not_identity(commutator); 
}

#define MAX_ROOTS 8
AJCC worst_primitive_cosh_re_len(const T& ch_o, const T& cs_o, const T& four_cosh_tube_diam_UB) {
  // Assumed ch and cs are real valued jets for cosh(Re(L)) and cos(Im(L))
  AJCC ch_prev = ch_o;
  AJCC cs_prev = cs_o;
  for (int i = 0; i < MAX_ROOTS; ++i) {
    AJCC ch = sqrt((ch_prev + 1) / 2);
    AJCC cs = sqrt((cs_prev + 1) / 2); // note, - pi <= Im(L) <= pi, so sign is +
    if (meyerhoff_k_test(ch, cs, four_cosh_tube_diam_UB)) {
      return ch_prev;
    }
    ch_prev = ch;
    cs_prev = cs;
  }
  // no luck
  AJCC zero(0);
  return zero; 
}

// Our compact parameter space has the following bounds:
// TODO
void verify_out_of_bounds(const char* where, char bounds_code)
{
  Box box = build_box(where);
  switch(bounds_code) {
    case '0':	{ // 1.0052 < cosh(0.104) <= cosh(mu) <= 0.846
                check(
                    absUB(box.cover.coshmu) < g_cosh_marg_lower_bound ||
                    absLB(box.cover.coshmu) > g_cosh_marg_upper_bound,
                    where);
              }
    case '1': {
                check(
                    absLB(box.cover.twocoshreD2) > g_cosh_r_bound * 2 ||
                    strictly_pos(re(-box.nearer.sinhL2)) ||
                    strictly_pos(re(-box.nearer.sinhD2)) ||
                    strictly_pos(one - box.cover.coshreD) ||
                    strictly_pos(one - box.cover.coshreL),
                    where);
              }
    default: {
               check(false, where);
             }
  }
}

// Meyerhoff k-test
#define MAX_MEYER 8
bool meyerhoff_k_test(const T& ch_o, const T& cs_o, const T& four_cosh_tube_diam_UB) {
  // Assumed ch and cs are real valued jets for cosh(Re(L)) and cos(Im(L))
  AJCC ch_prev = T(1);
  AJCC cs_prev = T(1);
  AJCC ch = ch_o;
  AJCC cs = cs_o;
  AJCC temp, four_cosh_tube_diam_LB;
  AJCC meyer_k = T(1024); // arbitray large enough number
  int count = 0;
  while (absUB(ch * ch) < 2 && count < MAX_MEYER) {
    temp = ch - cs; 
    if (strictly_pos(meyer_k - temp) && absUB((temp + 1) * (temp + 1)) < 2) {
      meyer_k = temp;
      // See Meyerhoff paper on volume lowerbounds for hyperbolic 3-manifolds
      four_cosh_tube_diam_LB = sqrt(-(meyer_k * 32) + 16) / meyer_k;
      if (strictly_pos(four_cosh_tube_diam_LB - four_cosh_tube_diam_UB)) {
        return true; // box can be killed
      }
    } 
    // Use Chebyshev recurrance relation
    temp = (ch_o * 2) * ch - ch_prev;
    ch_prev = ch;
    ch = temp;  
    AJCC temp = (cs_o * 2) * cs - cs_prev;
    cs_prev = cs;
    cs = temp;
    count +=1;
  }
  return false; // inconclusive
}

// Meyerhoff tube bound.
// Checks if embeded tube about axis(x) is more than (dist(axis(x), y axis(x))
// and similarly for axis(y)
void verify_meyerhoff(const char* where) {
  SL2AJCC x = construct_x(box.cover);
  SL2AJCC y = construct_y(bpx.cover);
  AJCC four_cosh_x_tube_UB = four_cosh_dist_ax_wax(y, box.cover);
  AJCC four_cosh_y_tube_UB = four_cosh_dist_ay_way(x, box.cover);
  check(
      meyerhoff_k_test(cover.coshreL, cover.cosimL, four_cosh_x_tube_UB) ||
      meyerhoff_k_test(cover.coshreL, cover.cosimL, four_cosh_y_tube_UB),
      where);
}

void verify_x_hits_y(const char* where, const char* word) {
  Box box = build_box(where);
  SL2AJCC w = construct_word(box.cover, word);
  check(moves_x_axis_too_close_to_y(w, box.cover) &&
        moved_x_axis_not_y_axis(w, box.cover),
        where);
}

void verify_y_hits_x(const char* where, const char* word) {
  Box box = build_box(where);
  SL2AJCC w = construct_word(box.cover, word);
  check(moves_y_axis_too_close_to_x(w, box.cover) &&
        moved_y_axis_not_x_axis(w, box.cover),
        where);
}

void verify_x_hits_x(const char* where, const char* word) {
  Box box = build_box(where);
  SL2AJCC w = construct_word(box.cover, word);
  check(inside_var_nbd_x(w, box.cover) &&
        cant_fix_x_axis(w, box.cover),
        where);
}

void verify_y_hits_y(const char* where, const char* word) {
  Box box = build_box(where);
  SL2AJCC w = construct_word(box.cover, word);
  check(inside_var_nbd_y(w, box.cover) &&
        cant_fix_y_axis(w, box.cover),
        where);
}

void verify_x_not_cyclic(const char* where, const char* word) {
  Box box = build_box(where);
  SL2AJCC x = construct_x(box.cover);
  SL2AJCC w = construct_word(box.cover, word);
  check(inside_var_nbd_x(w, box.cover) &&
        non_cyclic_power(w, x),
        where);
}

void verify_y_not_cyclic(const char* where, const char* word) {
  Box box = build_box(where);
  SL2AJCC x = construct_x(box.cover);
  SL2AJCC w = construct_word(box.cover, word);
  check(inside_var_nbd_x(w, box.cover) &&
        non_cyclic_power(w, x),
        where);
}

void verify_move(const char* where, const char* word) {
  Box box = build_box(where);
  SL2AJCC w = construct_word(box.cover, word);
  check(not_identity(w) &&
        move_less_than_marg(w, p) &&
        ((x_power(word) == 0 || y_power(word) == 0) || does_not_fix_sym_axis(w)),
        where);
}

void verify_w_ax_hits_sym_axis(const char* where, const char* word) {
  Box box = build_box(where);
  SL2AJCC w = construct_word(box.cover, word);
  check(wx_hits_sym_axis(w, box.cover),
        where);
} 
 
void verify_w_ay_hits_sym_axis(const char* where, const char* word) {
  Box box = build_box(where);
  SL2AJCC w = construct_word(box.cover, word);
  check(wy_hits_sym_axis(w, box.cover),
        where);
} 

inline bool word_in_array(const char* word, const char** array) {
  int len = sizeof(array)/sizeof(array[0]);
  for (int i = 0; i < len; ++i) {
    if (strncmp(word, array[i], 32) == 0) {
      return true;
    }
  }
  return false;
}

bool is_proven_elementary(const char* word, const AJCCParams& p) {
  SL2AJCC x = construct_x(p);
  SL2AJCC y = construct_y(p);
  SL2AJCC w = construct_word(word, p);
  if (y_power(word) > 0 && inside_var_nbd_x(w, p)) {
    AJCC four_cosh_x_tube_UB = four_cosh_dist_ax_wax(y, p);
    AJCC cosh_prim_re_len = worst_primitive_cosh_re_len(p.coshreL, p.cosimL, four_cosh_x_tube_UB); 
    AJCC diff = cosh_prim_re_len * 4 - four_cosh_re_length(w);
      if (strictly_pos(diff)) {
        return true;
      }      
    }
  if (x_power(word) > 0 && inside_var_nbd_y(w, p)) {
    AJCC four_cosh_y_tube_UB = four_cosh_dist_ay_way(x, p);
    AJCC cosh_prim_re_len = worst_primitive_cosh_re_len(p.coshreL, p.cosimL, four_cosh_y_tube_UB); 
    AJCC diff = cosh_prim_re_len * 4 - four_cosh_re_length(w);
      if (strictly_pos(diff)) {
        return true;
      }      
    }
  }
  return false;
}

void verify_impossible_relator(const char* where, const char* word) {
  Box box = build_box(where);
  check(word_in_array(word, impossible) &&
        is_proven_elementary(word, box.cover),
        where); 
}

void verify_symmetric_relator(const char* where, const char* word) {
  Box box = build_box(where);
  check(word_in_array(word, symmetric) &&
        is_proven_elementary(word, box.cover),
        where); 
}

void verify_vol3(const char* where, const char* first, const char* second) {
  Box box = build_box(where);
  check(strncmp(first, "XYXYxYXYXy") == 0 &&
        strncmp(second, "yXYxYxyxYxYY") == 0 &&
        is_proven_elementary(first, box.cover),
        is_proven_elementary(second, box.cover),
        where); 
}
