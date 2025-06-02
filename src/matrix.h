/*------------------------------------------------------------------------*/
/*! \file matrix.h

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#ifndef TALISMAN_SRC_MAT_
#define TALISMAN_SRC_MAT_
/*------------------------------------------------------------------------*/
#include <set>
#include <sstream>
#include <string.h>
#include <vector>

#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>
#include <flint/fmpz.h>

#include "gmp.h"
#include "gmpxx.h"

/*------------------------------------------------------------------------*/

inline bool
normalize_row(fmpq_mat_t mat, int i) {
  fmpz_t den;
  bool res = false;
  fmpz_init_set_si(den,1);
  for(long j = 0; j < fmpq_mat_ncols(mat); j++)
    fmpz_lcm(den, den, fmpq_mat_entry_den(mat, i, j));

  if(fmpz_popcnt(den) == 1) {
    for(long j = 0; j < fmpq_mat_ncols(mat); j++)
      fmpq_mul_fmpz(fmpq_mat_entry(mat,i,j), fmpq_mat_entry(mat,i,j), den);
  }
  res = true;
  fmpz_clear(den);
  return res;
}  

inline bool
is_denom_free(fmpq_mat_t mat, int i) {
  bool res = true;
  for(long j = 0; j < fmpq_mat_ncols(mat); j++) {
    if(fmpz_cmp_ui(fmpq_mat_entry_den(mat, i, j), 1) != 0) {
      res = false;
      break;
    }
  }
  return res;
}

inline bool
row_is_zero(fmpq_mat_t mat, int i) {
  for(long j = 0; j < fmpq_mat_ncols(mat); j++)
    if(!fmpq_is_zero(fmpq_mat_entry(mat,i,j)))
      return false;
  return true;
}
/*------------------------------------------------------------------------*/
inline std::vector<size_t>
rref(fmpq_mat_t mat) {
  
  fmpq_mat_rref(mat, mat);

  std::vector<size_t> pivots;
  for(auto i = 0; i < fmpq_mat_nrows(mat); i++) {
    bool piv = false;
    for(auto j = 0; j < fmpq_mat_ncols(mat); j++) {
      if(!piv && !fmpq_is_zero(fmpq_mat_entry(mat, i, j))) {
        piv = true;
        pivots.push_back(static_cast<size_t>(j));
      }
    }
  }

  return pivots;
}

/*------------------------------------------------------------------------*/
inline void
kernel(fmpq_mat_t M, fmpq_mat_t K) {

  // compute rref
  std::vector<size_t> pivots = rref(M);

  int n = fmpq_mat_ncols(M);
  fmpq_mat_t M_extended;
  fmpq_mat_init(M_extended, n, n);
  // insert pivot rows
  for(size_t i = 0; i < pivots.size(); i++) {
    for(int j = 0; j < n; j++)
      fmpq_set(fmpq_mat_entry(M_extended, pivots[i], j),
               fmpq_mat_entry(M, i, j));
  }

  // insert other rows
  for(int i = 0; i < n; i++) {
    if(fmpq_is_zero(fmpq_mat_entry(M_extended, i, i)))
      fmpq_set_si(fmpq_mat_entry(M_extended, i, i), -1, 1);
  }

  // extract kernel
  fmpq_mat_init(K, n - pivots.size(), n);
  size_t r = 0;
  for(int i = 0; i < n; i++) {
    if(fmpq_cmp_si(fmpq_mat_entry(M_extended, i, i), -1) == 0) {
      for(int j = 0; j < n; j++)
        if(!fmpq_is_zero(fmpq_mat_entry(M_extended, j, i)))
          fmpq_set(fmpq_mat_entry(K, r, j), fmpq_mat_entry(M_extended, j, i));
      r++;
    }
  }

  fmpq_mat_clear(M_extended);

  // bring kernel into rref
  rref(K);
  fmpq_mat_neg(K, K);
}

#endif// TALISMAN_SRC_MAT_
