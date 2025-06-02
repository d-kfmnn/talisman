/*------------------------------------------------------------------------*/
/*! \file reduction.cpp
    \brief contains functions used in the polynomial solver

  This file contains all functions used for preprocessing the Gröbner basis
  and for reducing the specification by the slices.

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#include "reduction.h"

#include <algorithm>
#include <list>

#include "propagate.h"
#include "variable.h"
/*------------------------------------------------------------------------*/
// Global variables
int proof = 0;
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
Polynomial *non_linear_reduction(Polynomial *rem) {

  Gate *g = gate(rem->get_lt()->get_var_num());
  while (!g->get_input()) {
    Polynomial *gc = g->get_gate_constraint();

    Polynomial *gc_unflip = unflip_poly(gc);
    gc = gc_unflip;

    if (verbose >= 2 && gc) {
      msg_nl("non-linear reducing by ");
      gc->print(stdout);
    }
 

    non_linear_count++;
    Polynomial *tmp = reduce_by_one_poly(rem, gc, 1);
    delete (rem);
    delete (gc);

  

    g->set_elim();

    rem = tmp;

    if (!rem) {
      msg("remainder is 0");
      return 0;
    }

    if (verbose >= 3) {
      msg_nl("remainder is ");
      rem->print(stdout);
      msg(" ");
    }
    if(rem->len() == 1 && !rem->get_lt()) break; // we need to abort when rem is a constant
    g = gate(rem->get_lt()->get_var_num());
  }
  
  if (rem) {
    Polynomial *mod_tmp = mod_poly(rem, NN);
    delete (rem);
    rem = mod_tmp;
  }
  return rem;
}

/*------------------------------------------------------------------------*/
Polynomial *reduce(Polynomial *spec) {
  print_hline();
  msg("starting reduction");
  assert(spec->degree() == 1);

  Polynomial *rem = spec;
  if (verbose > 1) {
    msg_nl("spec is: ");
    rem->print(stdout);
  }

  Gate *g = gate(rem->get_lt()->get_var_num());
  while (!g->get_input()) {
    if (g->get_gate_constraint()->degree() > 1) {
      Polynomial *p = remove_vanishing_monomials(g->get_gate_constraint());
      g->update_gate_poly(p);
    }

    if (g->get_gate_constraint()->degree() > 1) {
      Polynomial *p = unflip_poly_and_remove_van_mon(g->get_gate_constraint());
      g->update_gate_poly(p);
    }

    if (g->get_gate_constraint()->degree() > 1) {
      linearize_via_fglm_or_gap(g);
      if (!g->get_gate_constraint())
        die(2, "g lost gate constraint");
      if (verbose >= 3)
        g->get_gate_constraint()->print(stdout);
    }
    Polynomial *gc = g->get_gate_constraint();


    if (gc->degree() > 1) {
      msg_nl("failed to linearize gate poly: ");
      gc->print(stdout);

      msg("switching to non-linear rewriting");
      return non_linear_reduction(rem);
    }

    if (verbose >= 2 && gc) {
      msg_nl("linear reducing by ");
      gc->print(stdout);
    }

    Polynomial *tmp = substitute_linear_poly(rem, gc);
    linear_count++;
    g->set_elim();
    delete (rem);
    rem = tmp;
   
    if (rem) {
      Polynomial *mod_tmp = mod_poly(rem, NN);
      delete (rem);
      rem = mod_tmp;
    }

    if (!rem) {
      msg("remainder is 0");
      return 0;
    }

    if (verbose > 2) {
      msg_nl("remainder is ");
      rem->print(stdout);
      msg(" ");
    }
    g = gate(rem->get_lt()->get_var_num());
  }

  return rem;
}

/*------------------------------------------------------------------------*/
