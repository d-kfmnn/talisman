/*------------------------------------------------------------------------*/
/*! \file reductionmethods.h
    \brief contains all reduction methods of polynomials

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#ifndef TALISMAN_SRC_REDUCTIONMETHODS_H_ 
#define TALISMAN_SRC_REDUCTIONMETHODS_H_
/*------------------------------------------------------------------------*/
#include <string.h>
#include <vector>
#include <set>

#include "pac.h"

extern Polynomial * x_spec;
extern std::map<Term *, Polynomial*> van_poly;
extern std::map<Term *, Polynomial*>  dual_van_poly;
/*------------------------------------------------------------------------*/
Polynomial *flip_var_in_poly(Polynomial *p1, Var *v, bool rem_van);

//flip all dual variables
Polynomial *unflip_poly(Polynomial *p);
Polynomial* remove_vanishing_monomials(Polynomial* p, std::vector<Polynomial*> *used_van_poly = nullptr);

Polynomial* unflip_poly_and_remove_van_mon(Polynomial* p);
Polynomial *reduce_by_one_poly(Polynomial *p1, Polynomial *p2, bool non_lin_rewriting = false);
Polynomial *substitute_linear_poly(Polynomial *p1, Polynomial *p2);
Polynomial *mod_poly(Polynomial *p1, int exp);

bool reduce_to_zero(Polynomial * p, std::set<Polynomial*> G);
/*---------------------------------------------------------------------------*/

#endif // TALISMAN_SRC_REDUCTIONMETHODS_H_
