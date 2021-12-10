#ifndef __types_h
#define __types_h
#include <complex>
#include <vector>
#include <string.h>
#include <algorithm>

typedef std::complex<double> Complex;
typedef std::pair<double, double> float_pair;
typedef std::pair<std::string, std::string> word_pair;

inline const Complex eye(const Complex&x) { return Complex(0,1); };
inline const Complex operator+(const Complex&x,double y) { return x + Complex(y,0); };
inline const Complex operator-(const Complex&x,double y) { return x - Complex(y,0); };
inline const Complex operator*(const Complex&x,double y) { return x * Complex(y,0); };
inline const Complex operator/(const Complex&x,double y) { return x / Complex(y,0); };
inline const Complex abs(const Complex& x) { return Complex(hypot(x.real(), x.imag()), 0); };
inline const Complex re(const Complex& x) { return Complex(x.real(), 0); };
inline const Complex im(const Complex& x) { return Complex(x.imag(), 0); };
inline const Complex abs_sqrd(const Complex& x) { return Complex(x.real()*x.real()+x.imag()*x.imag(), 0); };

typedef enum _box_state
{
  killed_bounds = 1,
  killed_impossible_relator = 2, // var_nbd of w^k but w not id
  killed_x_hits_y = 3, // w(axis(x)) closer to axis(y) than dx + dy but > 0
  killed_y_hits_x = 4, // w(axis(y)) closer to axis(x) than dx + dy but > 0
  killed_x_hits_x = 5, // w(axis(x)) closer to axis(x) than 2dx but further than 0 (i.e. w not x^k)
  killed_y_hits_y = 6, // w(axis(y)) closer to axis(y) than 2dy but further than 0 (i.e. w not y^k)
  killed_x_not_cyclic = 7, // w(axis(x)) closer to axis(x) than 2dx and provably w not in same cyclic group
  killed_y_not_cyclic = 8, // w(axis(y)) closer to axis(y) than 2dy and provably w not in same cyclic group
  killed_move = 9, // w(1j) moved less than marg, note power not produces by word search
  killed_marg = 10, // w1 and w2 have (simple) margulis less than mu TODO should we do powers?
  killed_nbd_x = 11, // w and x fail Jorgensen and don't commute 
  killed_nbd_y = 12, // w and y fail Jorgensen and don't commute
  killed_nbd = 13, // w1 and w2 fail Jorgensen and one is not parabolic
  killed_w_ax_hits_sym_axis = 14, // only used in symmetric search
  killed_w_ay_hits_sym_axis = 15, // only used in symmetric search
  killed_via_sym = 16,
  proven_relator = 17,
  out_of_bounds_center = 50,
  maybe_killed_center = 51,
  var_nbd_x = 52, // w and x fail Jorgensen 
  var_nbd_y = 53, // w and y fail Jorgensen
  var_nbd = 54, // w1 and w2 fail Jorgensen
  open_with_qr = 99,
  open = -1,
} 
box_state;

struct TestResult {
  int index;
  box_state state;
  word_pair words;
};

Complex parse_complex(const std::string &complex_str);
Complex shift_imag_around_zero(const Complex &a);
Complex shift_imag_around_pi(const Complex &a);

double absUB(const Complex& x);
double absLB(const Complex& x);

std::string repeat(std::string s, int n);
void split_string(const std::string &str, const std::string &delims, std::vector<std::string> &out);

int x_power(std::string w);
int y_power(std::string w);
int syllables(std::string w);
std::string cyclic_strip(std::string w);
std::string x_strip(std::string w);
std::string y_strip(std::string w);
std::string x_rstrip(std::string w);
std::string y_rstrip(std::string w);

bool x_power_sort(std::string a, std::string b);
bool y_power_sort(std::string a, std::string b);

template<typename T>
inline T max_local(const T& x, const T& y) {
  return ((x+y)+abs(x-y))/2;
} 

template<typename T>
inline T min_local(const T& x, const T& y) {
  return ((x+y)-abs(x-y))/2;
} 

template<typename T> struct Params {
  T sinhL2;
  T sinhD2;
  T sinhsqL2; // derived parameter
  T sinhsqD2; // derived parameter
  T coshsqL2; // derived parameter
  T coshsqD2; // derived parameter
  T coshL2; // derived parameter
  T coshD2; // derived parameter
  T expD2; // derived parameter
  T expmD2; // derived parameter
  T twocoshreD2; // derived parameter
  T twosinhreD2; // derived parameter
  T coshreD; // derived parameter
  T coshmu; //derived parameter MUST BE SET LATER
  T coshreL; // derived
  T sinhreL; // derived
  T cosimL; //derived
};

template<typename T>
void fill_derived(Params<T>& p) {
  T one = T(1);
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

inline double re_center(const Complex& x) {
  return x.real();
}

template<typename T>
void print_type(T& x); 

template<typename T>
void print_center(const T& x); 

template<typename T>
bool sort_comp(const T& a, const T& b); 

template<typename T>
inline bool strictly_pos(const T& diff) {
  return absLB(re(diff)) > 0 && re_center(diff) > 0;
}

template<typename T>
void print_type(const char desc[], const T& x) {
  fprintf(stderr, "%s\n", desc);
  print_type(x);
} 

template<typename T>
void print_center(const char desc[], const T& x) {
  fprintf(stderr, "%s\n", desc);
  print_center(x);
} 

template<typename T>
inline const T powT(const T& x, int n) {
  T one(1);
  T a = one;
  if (n == 0) { return a; }
  T b;
  if (n < 0) { 
    b = one / x;
    n = -n;
  } else {
    b = x;
  }
  while (n > 1) {
    if (n & 1) { // n odd
      a = a*b;
    }
    b = b*b;
    n /= 2; // int division
  } 
  return a*b;
};

std::string double_to_hex(double x);

#endif // __types_h
