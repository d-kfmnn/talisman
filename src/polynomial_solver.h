/*------------------------------------------------------------------------*/
/*! \file polynomial_solver.h
    \brief contains the polynomial solving routine

  Part of TalisMan
  Copyright(C) 2025 TalisMan-Developers
*/
/*------------------------------------------------------------------------*/
#ifndef TALISMAN_SRC_POLYNOMIAL_SOLVER_H_
#define TALISMAN_SRC_POLYNOMIAL_SOLVER_H_
/*------------------------------------------------------------------------*/
#include "preprocessing.h"
#include "vanishing_constraints.h"
#include "reduction.h"
#include "witness.h"
#include "extensions.h"
#include "substitution.h"
/*------------------------------------------------------------------------*/
// / If final remainder is not equal to zero a counter example is generated and
// / printed to file <input_name>.wit, default is true, can be turned of
//  using command line input
extern bool gen_witness;
/*------------------------------------------------------------------------*/

/**
    Calls the preprocessing, slicing, reduction routines
    If certify is true, files need to be provided to store the proof.

    @param inp_f name of input file
    @param out_f1 name of first output file
    @param out_f2 name of second output file
    @param out_f3 name of third output file
    @param certify true when modus -certify is used

    @returns boolean whether circuit is correct (1) or incorrect (0)
*/
bool verify(const char *inp_f = 0, Polynomial *spec = 0, const char *out_f1 = 0, const char *out_f2 = 0,
            const char *out_f3 = 0);

#endif // TALISMAN_SRC_POLYNOMIAL_SOLVER_H_
