#ifndef _TestCollection_ 
#define _TestCollection_
#include <unordered_map>
#include <string>
#include <vector>
#include "types.hh"
#include "Box.h"
#include "SL2.hh"
#include "IsomH3.hh"
#include "Params.hh"
#include "RelatorTest.hh"

extern bool g_debug;

struct RelatorTest;

struct TestCollection {
  int size();
  box_state evaluate_center(int index, Box& box);
  TestResult evaluate_box(int index, Box& box);
  TestResult evaluate_AJCC(word_pair& pair, Box& box);
  TestResult evaluate_qrs(Box& box);
  const std::string get_name(int index);
  word_pair get_pair(int index);
  int add(word_pair pair);
  int add(std::string pair);
  void load(const char* file_path);
  void load_relator_test(const char* impos_path,
                         const char* bad_rel_path);
  RelatorTest *relator_test;
private:
  word_pair parse_word_pair(std::string buf);
  std::map<word_pair, int> pair_index;
  std::vector<word_pair> pair_vector;
  box_state evaluate_approx(word_pair pair, const Box& params);
  bool ready_for_elliptics_test(SL2<AJCC>& w);
  bool only_elliptics(SL2<AJCC>& w, Params<AJCC>& params);
};

template<typename T>
inline const bool inside_var_nbd(const SL2<T>& w1, const SL2<T>& w2) {
  // Show that w1 and w2 are commuting non-parabolics when discrete
  // So either we have a relator or they are elliptic
  return (absUB(jorgensen(w1, w2)) < 1 || absUB(jorgensen(w2, w1)) < 1) && (not_parabolic(w1) || not_parabolic(w2));
}

template<typename T>
inline const bool inside_var_nbd_ne(const SL2<T>& w1, const SL2<T>& w2) {
  // Show that w1 and w2 are commuting loxodromics when discrete, so we have a relator
  // Note,we must test both as elliptic can commute with loxodromic
  return (absUB(jorgensen(w1, w2)) < 1 || absUB(jorgensen(w2, w1)) < 1) && (not_elliptic_or_parabolic(w1) && not_elliptic_or_parabolic(w2));
}

template<typename T>
inline const bool not_parabolic(const SL2<T>& w) {
  T tr = w.a + w.d;
  return absLB(im(tr)) > 0 || (absLB(re(tr) - 2) > 0 && absLB(re(tr) + 2) > 0);
}

template<typename T>
inline const bool not_elliptic_or_parabolic(const SL2<T>& w) {
  T tr = w.a + w.d;
  return absLB(im(tr)) > 0 || absLB(re(tr)) > 2; 
}

template<typename T>
inline const bool not_identity(const SL2<T>& w) {
  return absLB(w.b) > 0 ||  absLB(w.c) > 0 ||
    ((absLB(w.a-1) > 0 || absLB(w.d-1) > 0) && (absLB(w.a+1) > 0 || absLB(w.d+1) > 0));
}

template<typename T>
inline const bool tube_hits_axis_two(const T& two_sinh_p2sq, const T& two_cosh_re_tube) {
  T tcd = two_cosh_dist(two_sinh_p2sq);
  // sinh(I Pi/4)^2 = -1/2 which means axes meet othrogonally
  return strictly_pos(two_cosh_re_tube - tcd) && absLB(two_sinh_p2sq + 1) > 0;
}

template<typename T>
inline const bool tube_hits_axis_four(const T& four_sinh_p2sq, const T& two_cosh_re_tube) {
  T fcd = four_cosh_dist(four_sinh_p2sq);
  // sinh(I Pi/4)^2 = -1/2 which means axes meet othrogonally
  return strictly_pos(two_cosh_re_tube * 2 - fcd) && absLB(four_sinh_p2sq + 2) > 0;
}

template<typename T>
inline const bool wx_hits_sym_axis(const SL2<T>& w, const Params<T>& p) {
  T tsp2sq_inf = two_sinh_perp2_sq_wax_zero_inf(w, p);
  T fsp2sq_one = four_sinh_perp2_sq_wax_mp_one(w, p);
  T fsp2sq_iye = four_sinh_perp2_sq_wax_mp_iye(w, p);
  return (tube_hits_axis_two(tsp2sq_inf, p.twocoshreD2) ||
          tube_hits_axis_four(fsp2sq_one, p.twocoshreD2) || 
          tube_hits_axis_four(fsp2sq_iye, p.twocoshreD2));
}

template<typename T>
inline const bool wy_hits_sym_axis(const SL2<T>& w, const Params<T>& p) {
  T tsp2sq_inf = two_sinh_perp2_sq_way_zero_inf(w, p);
  T fsp2sq_one = four_sinh_perp2_sq_way_mp_one(w, p);
  T fsp2sq_iye = four_sinh_perp2_sq_way_mp_iye(w, p);
  if (g_debug) {
    if (tube_hits_axis_two(tsp2sq_inf, p.twocoshreD2)) {
      fprintf(stderr, "Hit zero inf\n");
      print_type(tsp2sq_inf);
    }
    if (tube_hits_axis_four(fsp2sq_one, p.twocoshreD2)) {
      fprintf(stderr, "Hit +- one\n");
      print_type(fsp2sq_one);
    }
    if (tube_hits_axis_four(fsp2sq_iye, p.twocoshreD2)){
      fprintf(stderr, "Hit +- iye\n");
      print_type(fsp2sq_iye);
    }
  }
  return (tube_hits_axis_two(tsp2sq_inf, p.twocoshreD2) ||
          tube_hits_axis_four(fsp2sq_one, p.twocoshreD2) || 
          tube_hits_axis_four(fsp2sq_iye, p.twocoshreD2));
}

template<typename T>
inline const bool really_cant_fix_x_axis(const SL2<T>& w, const Params<T>& p) {
  T fsp2sq = four_sinh_perp2_sq_ax_wax(w, p);
  if (g_debug && absLB(fsp2sq) > 0 && absLB(fsp2sq + 4) > 0) {
    fprintf(stderr, "Realy can't fix x_axis");
    print_type(fsp2sq);
    fprintf(stderr, "LB values %f and %f\n", absLB(fsp2sq), absLB(fsp2sq + 4));
    fprintf(stderr, "LB away from %d and %d\n", absLB(fsp2sq) > 0, absLB(fsp2sq + 4) > 0);
  }
  return absLB(fsp2sq) > 0 && absLB(fsp2sq + 4) > 0; 
}

template<typename T>
inline const bool cant_fix_x_axis(const SL2<T>& w, const Params<T>& p) {
  T fsp2sq = four_sinh_perp2_sq_ax_wax(w, p);
  if (g_debug && 
      std::is_same<T, AJCC>::value && absLB(fsp2sq) > 0 && absLB(fsp2sq + 4) > 0) {
    fprintf(stderr, "********** CANNOT FIX X AXIS ***********\n");
    T fsp2sq = four_sinh_perp2_sq_ax_wax(w, p);
    print_type(fsp2sq);
    fprintf(stderr, "Can't fix x axis LB values %f and %f\n",
      absLB(fsp2sq), absLB(fsp2sq + 4));
    fprintf(stderr,"Can't fix x axis LB away from %d and %d\n",
      absLB(fsp2sq) > 0, absLB(fsp2sq + 4) > 0);
    fprintf(stderr, "*******************************\n");
  }
  return absLB(fsp2sq) > 0 && absLB(fsp2sq + 4) > 0; 
}

template<typename T>
inline const bool must_fix_x_axis(const SL2<T>& w, const Params<T>& p) {
  // The "must" part is only valid for AJCC tests
  T diff = p.coshreD * 4 - four_cosh_dist_ax_wax(w, p);
  if (g_debug && 
      std::is_same<T, AJCC>::value && strictly_pos(diff)) {
    fprintf(stderr, "********** MUST FIX X AXIS ***********\n");
    print_SL2(w);
    print_type("4 cosh 2 dx:", p.coshreD * 4);
    print_type("4 cosh dist ax wax:", four_cosh_dist_ax_wax(w, p));
    print_type("diff:", diff);
    fprintf(stderr, "*******************************\n");
  }
  // We know that diff is away from zero and the diff should be conj symmetrix, so
  // we only test if the real part is to one side of the bound
  return strictly_pos(diff);
}

#define LERR 0.00000000001

template<typename T>
inline const bool really_cant_fix_y_axis(const SL2<T>& w, const Params<T>& p) {
  T fsp2sq = four_sinh_perp2_sq_ay_way(w, p);
  if (g_debug && absLB(fsp2sq) > LERR && absLB(fsp2sq + 4) > LERR) {
      fprintf(stderr, "Realy can't fix y_axis");
    print_type(fsp2sq);
    fprintf(stderr, "LB values %f and %f\n", absLB(fsp2sq), absLB(fsp2sq + 4));
    fprintf(stderr, "LB away from %d and %d\n", absLB(fsp2sq) > 0, absLB(fsp2sq + 4) > 0);
  }
  return absLB(fsp2sq) > LERR && absLB(fsp2sq + 4) > LERR; 
}

template<typename T>
inline const bool cant_fix_y_axis(const SL2<T>& w, const Params<T>& p) {
  T fsp2sq = four_sinh_perp2_sq_ay_way(w, p);
  if (g_debug && 
      std::is_same<T, AJCC>::value && absLB(fsp2sq) > 0 && absLB(fsp2sq + 4) > 0) {
   fprintf(stderr, "********** CANNOT FIX Y AXIS ***********\n");
   print_type(fsp2sq);
   fprintf(stderr, "LB values %f and %f\n", absLB(fsp2sq), absLB(fsp2sq + 4));
   fprintf(stderr, "LB away from %d and %d\n", absLB(fsp2sq) > 0, absLB(fsp2sq + 4) > 0);
  }
  return absLB(fsp2sq) > 0 && absLB(fsp2sq + 4) > 0; 
}

template<typename T>
inline const bool must_fix_y_axis(const SL2<T>& w, const Params<T>& p) {
  // The "must" part is only valid for AJCC tests
  T diff = p.coshreD * 4 - four_cosh_dist_ay_way(w, p);
  if (g_debug && 
      std::is_same<T, AJCC>::value && strictly_pos(diff)) {
    fprintf(stderr, "********** MUST FIX Y AXIS ***********\n");
    print_SL2(w);
    print_type("4 cosh 2 dy:", p.coshreD * 4);
    print_type("4 cosh dist ay way:", four_cosh_dist_ay_way(w, p));
    print_type("diff:", diff);
    fprintf(stderr, "*******************************\n");
  }
  // We know that diff is away from zero and the diff should be conj symmetrix, so
  // we only test if the real part is to one side of the bound
  return strictly_pos(diff);
}

template<typename T>
inline const bool inside_var_nbd_x(const SL2<T>& w, const Params<T>& params) {
  // The second test may only work when x has trace close to +/- 2
  if (g_debug && (absUB(jorgensen_wx(w, params)) < 1 ||
      absUB(jorgensen_xw(w, params)) < 1 || must_fix_x_axis(w, params))) {
      fprintf(stderr, "UB Jwx %f, UB Jxw %f, must_fix %d\n", absUB(jorgensen_wx(w, params)),
        absUB(jorgensen_xw(w, params)), must_fix_x_axis(w, params));
  }
  return absUB(jorgensen_wx(w, params)) < 1 || absUB(jorgensen_xw(w, params)) < 1 || must_fix_x_axis(w, params);
}

template<typename T>
inline const bool inside_var_nbd_y(const SL2<T>& w, const Params<T>& params) {
  // The second test may only work when y has trace close to +/- 2
  if (g_debug && (absUB(jorgensen_wy(w, params)) < 1 ||
    absUB(jorgensen_yw(w, params)) < 1 || must_fix_y_axis(w, params))) {
    fprintf(stderr, "UB Jwy %f, UB Jyw %f, must_fix %d\n", absUB(jorgensen_wy(w, params)),
      absUB(jorgensen_yw(w, params)), must_fix_y_axis(w, params));
  }
  return absUB(jorgensen_wy(w, params)) < 1 || absUB(jorgensen_yw(w, params)) < 1 || must_fix_y_axis(w, params);
}


template<typename T>
inline const bool moves_y_axis_too_close_to_x(const SL2<T>& w, const Params<T>& p) {
  T diff = p.coshreD * 4 - four_cosh_dist_ax_way(w, p);
  // We know that diff is away from zero and the diff should be conj symmetrix, so
  // we only test if the real part is to one side of the bound
  if (g_debug && std::is_same<T, AJCC>::value && strictly_pos(diff)) {
    fprintf(stderr, "****************************************\n");
    fprintf(stderr, "MOVES Y TOO CLOSE TO X\n");
    print_SL2(w);
    print_type("4cosh(dx+dy):", p.coshreD * 4);
    print_type("4coshd(dist(x-axis, w(y-axis))):", four_cosh_dist_ax_way(w, p)); 
    T z = ((w.a * w.a) * p.expD2  - (w.b * w.b) * p.expmD2) * p.expD2 +
      ((w.d * w.d) * p.expmD2 - (w.c * w.c) * p.expD2 ) * p.expmD2;
    print_type("4 sinh^2(dist/2) + 2:", z);
    print_type("|4 sinh^2(dist/2)|:", abs(z - 2));
    print_type("|4 cosh^2(dist/2)|:", abs(z + 2));
    print_type("4 cosh(dist):",  abs(z - 2) + abs(z + 2));
    print_type("diff:", diff);
    fprintf(stderr, "diff is positive: %d\n", strictly_pos(diff));
    fprintf(stderr, "****************************************\n");
  }
  return strictly_pos(diff);
}

template<typename T>
inline const bool moves_x_axis_too_close_to_y(const SL2<T>& w, const Params<T>& p) {
  T diff = p.coshreD * 4 - four_cosh_dist_ay_wax(w, p);
  // We know that diff is away from zero and the diff should be conj symmetrix, so
  // we only test if the real part is to one side of the bound
  if (g_debug && std::is_same<T, AJCC>::value && strictly_pos(diff)) {
    fprintf(stderr, "****************************************\n");
    fprintf(stderr, "MOVES X TOO CLOSE TO Y\n");
    print_SL2(w);
    print_type("4cosh(dx+dy):", p.coshreD * 4);
    print_type("4coshd(dist(y-axis, w(x-axis))):", four_cosh_dist_ay_wax(w, p)); 
    T z = ((w.a * w.a) * p.expmD2 - (w.b * w.b) * p.expD2 ) * p.expmD2 +
      ((w.d * w.d) * p.expD2  - (w.c * w.c) * p.expmD2) * p.expD2;
    print_type("4 sinh^2(dist/2) + 2:", z);
    print_type("|4 sinh^2(dist/2)|:", abs(z - 2));
    print_type("|4 cosh^2(dist/2)|:", abs(z + 2));
    print_type("4 cosh(dist):",  abs(z - 2) + abs(z + 2));
    print_type("diff:", diff);
    fprintf(stderr, "diff is positive: %d\n", strictly_pos(diff));
    fprintf(stderr, "****************************************\n");
  }
  return strictly_pos(diff);
}

template<typename T>
inline const bool moved_y_axis_not_x_axis(const SL2<T>& w, const Params<T>& p) {
  T fsp2sq = four_sinh_perp2_sq_ax_way(w, p);
  return absLB(fsp2sq) > 0 && absLB(fsp2sq + 4) > 0; 
}

template<typename T>
inline const bool moved_x_axis_not_y_axis(const SL2<T>& w, const Params<T>& p) {
  T fsp2sq = four_sinh_perp2_sq_ay_wax(w, p);
  return absLB(fsp2sq) > 0 && absLB(fsp2sq + 4) > 0; 
}

template<typename T>
inline bool margulis_smaller_than_xy(const SL2<T>& w1, const SL2<T>& w2, const Params<T>& p) {
  T diff = p.coshmu * 4 - four_cosh_margulis_simple(w1, w2).first;
  // We know that diff is away from zero and the diff should be conj symmetrix, so
  // we only test if the real part is to one side of the bound
  return strictly_pos(diff);
}

template<typename T>
inline bool move_less_than_marg(const SL2<T>& w, const Params<T>& p) {
  T diff = p.coshmu - cosh_move_j(w);
  return strictly_pos(diff); 
}

template<typename T>
inline bool non_cylic_power(const SL2<T>& w, const SL2<T>& x_or_y) {
  // Assume word fixes the same axis as x or y, so it must live in a cyclic group with x or y.
  // Here we check that this is impossible in this box. Must use margulis
  // number to check cut off for roots of x or y
  SL2<T> commutator = x_or_y * w * inverse(w * x_or_y);
  if (g_debug && std::is_same<T, AJCC>::value && not_identity(commutator)) {
    fprintf(stderr, "****************************************\n");
    fprintf(stderr, "NOT CYCLIC POWER\n");
    fprintf(stderr, "x or y\n");
    print_SL2(x_or_y);
    fprintf(stderr, "(x or y)^2\n");
    print_SL2(x_or_y * x_or_y);
    fprintf(stderr, "w\n");
    print_SL2(w);
    fprintf(stderr, "commutator\n");
    print_SL2(commutator);
    fprintf(stderr, "|b| == 0: %d, |c| == 0: %d, |a-1| == 0: %d, |d-1| == 0: %d, |a+1| == 0: %d, |d+1| == 0: %d\n", absLB(commutator.b) == 0, absLB(commutator.c) == 0, absLB(commutator.a-1) == 0,absLB(commutator.d-1) == 0, absLB(commutator.a+1) == 0, absLB(commutator.d+1) == 0);
    fprintf(stderr, "****************************************\n");
  }
  // TODO Test powers when coshmu > coshsdx + sinhsdx
  return not_identity(commutator); 
}

// Meyerhoff K Test
// We stop computing if we fail the test
#define MAX_MEYER 8
template<typename T>
bool meyerhoff_k_test(const T& ch_o, const T& cs_o, const T& four_cosh_tube_diam_UB) {
  // Assumed ch and cs are real valued jets for cosh(Re(L)) and cos(Im(L))
  T ch_prev = T(1);
  T cs_prev = T(1);
  T ch = ch_o;
  T cs = cs_o;
  T temp, four_cosh_tube_diam_LB;
  T meyer_k = T(1024); // arbitray large enough number
  int count = 0;
  while (absUB(ch * ch) < 2 && count < MAX_MEYER) {
    temp = ch - cs; 
    if (strictly_pos(meyer_k - temp) && absUB((temp + 1) * (temp + 1)) < 2) {
      meyer_k = temp;
      // See Meyerhoff paper on volume lowerbounds for hyperbolic 3-manifolds
      four_cosh_tube_diam_LB = sqrt(-(meyer_k * 32) + 16) / meyer_k;
      if (strictly_pos(four_cosh_tube_diam_LB - four_cosh_tube_diam_UB)) {
        if (g_debug) {
          fprintf(stderr, "Meyer k %f with 4 cosh tube diam LB %f and UB %f\n",
              absLB(meyer_k), absUB(four_cosh_tube_diam_LB), absLB(four_cosh_tube_diam_UB));
        }
        return true; // box can be killed
      }
    } 
    // Use Chebyshev recurrance relation
    temp = (ch_o * 2) * ch - ch_prev;
    ch_prev = ch;
    ch = temp;  
    T temp = (cs_o * 2) * cs - cs_prev;
    cs_prev = cs;
    cs = temp;
    count +=1;
  }
  return false; // inconclusive
}

#define MAX_ROOTS 8
template<typename T>
T worst_primitive_cosh_re_len(const T& ch_o, const T& cs_o, const T& four_cosh_tube_diam_UB) {
  // Assumed ch and cs are real valued jets for cosh(Re(L)) and cos(Im(L))
  T ch_prev = ch_o;
  T cs_prev = cs_o;
  for (int i = 0; i < MAX_ROOTS; ++i) {
    T ch = sqrt((ch_prev + 1) / 2);
    T cs = sqrt((cs_prev + 1) / 2); // note, - pi <= Im(L) <= pi, so sign is +
    if (meyerhoff_k_test(ch, cs, four_cosh_tube_diam_UB)) {
      return ch_prev;
    }
    ch_prev = ch;
    cs_prev = cs;
  }
  // no luck
  T zero(0);
  return zero; 
}

template<typename T>
T cosh_marg_lower_bound(const T& two_sinh_r) {
  T s = two_sinh_r;
  T a8 = powT(s, 8) * (-0.0000014461700558); 
  T a7 = powT(s, 7) *   0.0000365880448817; 
  T a6 = powT(s, 6) * (-0.0003163830157272);
  T a5 = powT(s, 5) *   0.0005316504647188;
  T a4 = powT(s, 4) *   0.0086912125268823;
  T a3 = powT(s, 3) * (-0.061949675652791);
  T a2 = powT(s, 2) *   0.151649220047696;
  T a1 = s          * (-0.01513801009421);
  double a0 = 0.9999999;
  return ((a8 + (a1 + a0)) + (a4 + a5)) + ((a7 + a2) + (a6 + a3)); 
}

bool proven_is_good(const std::string& proven, const Box& box);

#define MAX_ID_SHIFT 5
template<typename T>
std::string proven_identity(std::string word, const Params<T>& p) {
  SL2<T> x = construct_x(p);
  SL2<T> y = construct_y(p);
  if (g_debug) {
    fprintf(stderr, "Testing proven identity for word: %s .\n", word.c_str());
  }
  SL2<T> w = construct_word(word, p);
  std::string new_word;
  if (y_power(word) > 0 && inside_var_nbd_x(w, p)) {
    T four_cosh_x_tube_UB = four_cosh_dist_ax_wax(y, p);
    T cosh_prim_re_len = worst_primitive_cosh_re_len(p.coshreL, p.cosimL, four_cosh_x_tube_UB); 
    for (auto s : {"x", "X"}) {
      new_word = x_strip(word);
      for (int i = 0; i < MAX_ID_SHIFT; ++i) {
        SL2<T> new_w = construct_word(new_word, p); // order matters
        T diff = cosh_prim_re_len * 4 - four_cosh_re_length(new_w);
        if (strictly_pos(diff)) {
          if (g_debug) {
            fprintf(stderr, "Found proven identity: %s .\n", new_word.c_str());
          }
          return new_word;
        }      
        new_word = s + new_word;
      }
    }
  }
  if (x_power(word) > 0 && inside_var_nbd_y(w, p)) {
    T four_cosh_y_tube_UB = four_cosh_dist_ay_way(x, p);
    T cosh_prim_re_len = worst_primitive_cosh_re_len(p.coshreL, p.cosimL, four_cosh_y_tube_UB); 
    for (auto s : {"y", "Y"}) {
      new_word = y_strip(word);
      for (int i = 0; i < MAX_ID_SHIFT; ++i) {
        SL2<T> new_w = construct_word(new_word, p); // order matters
        T diff = cosh_prim_re_len * 4 - four_cosh_re_length(new_w);
        if (strictly_pos(diff)) {
          if (g_debug) {
            fprintf(stderr, "Found proven identity: %s .\n", new_word.c_str());
          }
          return new_word;
        }      
        new_word = s + new_word;
      }
    }
  }
  return "";
}

#endif //_TestCollection_
