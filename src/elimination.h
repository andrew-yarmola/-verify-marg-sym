#ifndef _elimination_h_
#define _elimination_h_
#include "box.h"

#define MAX_DEPTH 256
#define MAX_CODE_LEN 512

// Helper functions

void check(bool inequalities, const char* where);

void parse_word(char* code);

typedef struct {
    char first[MAX_CODE_LEN];
    char second[MAX_CODE_LEN];
} word_pair;

word_pair get_word_pair(const char* code);

// Elimination functions

SL2ACJ construct_x(const ACJParams& params);

ACJ construct_T(const ACJParams& params, int M, int N);

SL2ACJ construct_word(const ACJParams& params, const char* word);

void verify_out_of_bounds(const char* where, char bounds_code);

void verify_meyerhoff(const char* where);

void verify_marg_lower_bound(const char* where);

void verify_x_hits_y(const char* where, const char* word);

void verify_y_hits_x(const char* where, const char* word);

void verify_x_hits_x(const char* where, const char* word);

void verify_y_hits_y(const char* where, const char* word);

void verify_x_not_cyclic(const char* where, const char* word);

void verify_y_not_cyclic(const char* where, const char* word);

void verify_move(const char* where, const char* word);

void verify_w_ax_hits_sym_axis(const char* where, const char* word);

void verify_w_ay_hits_sym_axis(const char* where, const char* word);

void verify_impossible_relator(const char* where, const char* word);

void verify_symmetric_relator(const char* where, const char* word);

void verify_vol3(const char* where, const char* first, const char* second);

#endif // _elimination_h_
